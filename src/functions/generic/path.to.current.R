# Function to find the most current file,
# given that the files were generated with Sys.Date()
# takes a path, a pattern(unique to file), and a file-type (e.g. csv or rds)
library(stringr)
path.to.current <- function(input, file_pattern, file_type){
  files <- dir(input, pattern = file_pattern, full.names=TRUE)
  files <- files[grepl(file_type, files)]
  # file-ending + .
  length_end <- nchar(file_type) + 1
  # 13 because lenght of date
  dates <- c(str_sub(files , start = unique(nchar(files)) - 13,
                     end = unique(nchar(files) - length_end)))
  path <- files[which.max(as.Date(dates))]
  return(path)
}
