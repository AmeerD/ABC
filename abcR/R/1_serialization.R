# Reader for prepped survey files (qs2 format). A file saved with the former
# qs package is reported rather than misread, as the two formats differ.
read_prepped_survey <- function(path) {
  if (!file.exists(path)) {
    legacy_path <- sub("\\.qs2$", ".qs", path)
    if (file.exists(legacy_path)) {
      stop(
        "Legacy qs file found: ", legacy_path,
        ". abcR now reads qs2 files: convert it with qs::qread() followed by ",
        "qs2::qs_save() (renaming the file is not enough). See the README ",
        "section 'Prepped survey files'.",
        call. = FALSE
      )
    }
    stop("Preprocessed survey file not found: ", path, call. = FALSE)
  }
  qs2::qs_read(path, validate_checksum = TRUE)
}
