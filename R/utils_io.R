#' Resolve an input path
#'
#' If `path` exists as given, return it; otherwise join it with `base_dir`.
#'
#' @param path A file path.
#' @param base_dir A base directory.
#'
#' @return A resolved file path.
#' @noRd
resolve_path <- function(path, base_dir = NULL) {
  if (!is.null(path) && file.exists(path)) {
    return(path)
  }
  if (!is.null(base_dir)) {
    candidate <- file.path(base_dir, path)
    if (file.exists(candidate)) {
      return(candidate)
    }
  }
  stop("File not found: ", path, call. = FALSE)
}

#' Read a TSV file as data.table
#'
#' @param path Path to a TSV file.
#'
#' @return A data.table.
#' @importFrom data.table fread
#' @noRd
read_tsv_dt <- function(path) {
  data.table::fread(path, sep = "\t")
}