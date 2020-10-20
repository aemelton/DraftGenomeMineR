# AEM, 2020
# IDK WHY NOTHING WORKS FOR READING THE FILE INTO MY DELL :(

processFile = function(filepath) {
  res <- NULL
  con <- file(filepath, "r")
  while ( TRUE ) {
    line <- readLines(con, n = 1, encoding = "UTF-8")
    if ( length(line) == 0 ) {
      break
    }
    res <- c(res, line)
  }
  return(res)
  close(con)
  suppressWarnings()
}
