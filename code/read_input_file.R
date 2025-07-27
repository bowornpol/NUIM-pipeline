read_input_file <- function(file_path, file_type = c("csv", "tsv"), ...) {
  file_type <- match.arg(file_type)

  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }

  if (file_type == "csv") {
    data <- read.csv(file_path, ...)
  } else if (file_type == "tsv") {
    data <- read.delim(file_path, ...)
  } else {
    stop("Invalid file_type. Must be 'csv' or 'tsv'.")
  }
  return(data)
}
