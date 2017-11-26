#!/usr/bin/env amm

import $ivy.`com.github.tototoshi::scala-csv:1.3.4`
import $ivy.`org.apache.commons:commons-math3:3.6.1`

import ammonite.ops._
import com.github.tototoshi.csv._

// http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics

import scala.collection.mutable.ListBuffer
import scala.collection.mutable.{ Map => MMap }
import scala.io.Source

@main
def main(input: Path, output: Path, mapping: Path, mappingFuture: Path, year: Int = 2015, threshold: Double = 0, bufSize: Int = 4194304) = {

  // --- grid id to category mapping

  val categories = {
    // open mapping
    val file = Source.fromFile(mapping.toIO, bufSize)
    val reader = CSVReader.open(file)

    // grid_id -> category
    val categories: Map[Int, String] = {
      for {
        (row, index) <- {
          val it = reader.iterator.zipWithIndex
          it.next // drop header
          it
        }
        category = row(0) if category != "NA"
      } yield index -> category
    }.toMap

    reader.close()

    categories
  }

  val categoriesFuture = {
    // open mapping
    val file = Source.fromFile(mappingFuture.toIO, bufSize)
    val reader = CSVReader.open(file)

    // grid_id -> category
    val categories: Map[Int, String] = {
      for {
        (row, index) <- {
          val it = reader.iterator.zipWithIndex
          it.next // drop header
          it
        }
        category = row(0) if category != "NA"
      } yield index -> category
    }.toMap

    reader.close()

    categories
  }

  // --- input

  val readers: Array[CSVReader] = {
    val files = for {
      file <- ls! input |? (_.ext == "csv")
    } yield {
      Source.fromFile(file.toIO, bufSize)
    }

    files.toArray.map(file => CSVReader.open(file))
  }

  // final input data structure

  // file   row    col    cell content
  // list   list   list   string

  val iterators: Array[Iterator[Seq[String]]] = readers.map(_.iterator)

  // --- drop headers

  println("dropping headers ...")

  iterators.foreach(_.next)

  // --- stats per year per category

  //                     cat     bstrap sum
  val year_c_b_sum: MMap[String, Array[Double]] =
    MMap[String, Array[Double]]()

  val future_c_b_sum: MMap[String, Array[Double]] =
    MMap[String, Array[Double]]()

  // --- read data and fill stats

  var grid_id = 1

  while (iterators.forall(_.hasNext)) {

    categories get grid_id match {
      case Some(category) =>
        println(s"""processing line $grid_id with $category ... """)

        // 1000 x 16 x Double
        val bootstraps: Array[Seq[Double]] = iterators.map(_.next.map(_.toDouble))

        for {
          bstrap_id <- 0 until 1000
          bootstrap: Seq[Double] = bootstraps(bstrap_id)
        } {
          val value: Double = bootstrap(year - 1999)

          if (value > threshold) {
            { // year
              val c_b_sum: MMap[String, Array[Double]] =
                year_c_b_sum

              val bstraps = c_b_sum get category match {
                case Some(bstraps) =>
                  bstraps

                case None =>
                  val bstraps = Array.fill(1000)(0.0)
                  c_b_sum.update(category, bstraps)
                  bstraps
              }

              bstraps(bstrap_id) = bstraps(bstrap_id) + value
            }

            { // future
              val c_b_sum: MMap[String, Array[Double]] =
                future_c_b_sum

              val bstraps = c_b_sum get category match {
                case Some(bstraps) =>
                  bstraps

                case None =>
                  val bstraps = Array.fill(1000)(0.0)
                  c_b_sum.update(category, bstraps)
                  bstraps
              }

              if (categoriesFuture.isDefinedAt(grid_id)) {
                bstraps(bstrap_id) = bstraps(bstrap_id) + value
              }
            }
          }
        }

      case None =>
        println(s"""skipping $grid_id - that cell doesn't have that category ... """)
        iterators.map(_.next)
    }

    grid_id += 1
  }

  readers.foreach(_.close())

  // --- calculate stats and write out

  val early_c_b_sum = year_c_b_sum
  val later_c_b_sum = future_c_b_sum

  val stats = for (category <- early_c_b_sum.keys) yield {
    val early_b_sum: Array[Double] = early_c_b_sum(category)
    val later_b_sum: Array[Double] = later_c_b_sum(category)

    val diff_b_sum: Array[Double] = for {
      (early, later) <- early_b_sum zip later_b_sum
    } yield (early - later)

    val ds = diff_b_sum.foldLeft(new DescriptiveStatistics) {
      (ds, value) => {
        ds addValue value
        ds
      }
    }

    category -> ds
  }

  val writer = CSVWriter.open(output.toIO)

  val output_header = Seq (
    "category",
    "percentile_2.5",
    "percentile_25",
    "median",
    "percentile_75",
    "percentile_97.5"
  )

  writer.writeRow(output_header)

  //            category stats
  // stats: Map[String,  DescriptiveStatistics]

  for ((category, ds) <- stats) {
    println(s"""writing $category ...""")

    val values = Seq (
      category,
      ds.getPercentile(2.5),
      ds.getPercentile(25),
      ds.getPercentile(50),
      ds.getPercentile(75),
      ds.getPercentile(97.5)
    )

    writer.writeRow(values)
  }

  writer.close()
}
