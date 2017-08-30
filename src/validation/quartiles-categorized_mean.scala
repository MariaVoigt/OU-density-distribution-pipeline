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
def main(input: Path, output: Path, mapping: Path, years: Seq[Int] = 1999 to 2015, threshold: Double = 0.1, bufSize: Int = 4194304) = {

  // --- grid id to category mapping

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

  //                 year      cat     bstrap sum
  val y_c_b_sum: Map[Int, MMap[String, Array[(Int, Double)]]] =
    years.map(_ -> MMap[String, Array[(Int, Double)]]()).toMap

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
          year <- years
        } {
          val value: Double = bootstrap(year - 1999)

          if (value > threshold) {
            val c_b_sum: MMap[String, Array[(Int, Double)]] =
              y_c_b_sum(year)

            val bstraps = c_b_sum get category match {
              case Some(bstraps) =>
                bstraps

              case None =>
                val bstraps = Array.fill(1000)((0, 0.0))
                c_b_sum.update(category, bstraps)
                bstraps
            }

            val (prev_count: Int, prev_value: Double) = bstraps(bstrap_id)

            bstraps(bstrap_id) = (prev_count + 1, prev_value + value)
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

  val stats = {
    y_c_b_sum mapValues { c_b_sum =>
      c_b_sum mapValues { (b_sum: Array[(Int, Double)]) =>
        val b_means: Array[Double] = b_sum map {
          case (count, sum) => sum / count
        }

        b_means.foldLeft(new DescriptiveStatistics) {
          (ds, value) => {
            ds addValue value
            ds
          }
        }
      }
    }
  }

  val writer = CSVWriter.open(output.toIO)

  val output_header = Seq (
    "year",
    "category",
    "percentile_2.5",
    "percentile_25",
    "median",
    "percentile_75",
    "percentile_97.5"
  )

  writer.writeRow(output_header)

  //            year     category stats
  // stats: Map[Int, Map[String,  DescriptiveStatistics]]

  for {
    (year, catstats) <- stats
    (category, ds) <- catstats
  } {
    println(s"""writing $year $category ...""")

    val values = Seq (
      year,
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
