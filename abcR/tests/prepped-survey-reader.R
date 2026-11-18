# Tests for the pre-processed survey reader: qs2 files, legacy qs files and
# missing files. No Stan sampling is needed.
local({
  reader <- getFromNamespace("read_prepped_survey", "abcR")
  dir <- tempfile("abcR-prepped-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  survey <- data.frame(
    age = c(10L, 25L, 45L),
    level = factor(c("prim", "prim", "usec")),
    value = c(0, 1, NA_real_)
  )
  error_message <- function(expr) {
    tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  }

  # no file at all
  stopifnot(grepl("not found", error_message(reader("XXX", dir))))

  # a qs2 file is read back unchanged
  qs2::qs_save(survey, file.path(dir, "AAA_prepped.qs2"))
  stopifnot(identical(reader("AAA", dir), survey))

  # a legacy qs file: read when qs is installed, otherwise a conversion hint
  legacy <- file.path(dir, "BBB_prepped.qs")
  if (requireNamespace("qs", quietly = TRUE)) {
    qs::qsave(survey, legacy)
    stopifnot(identical(reader("BBB", dir), survey))
  } else {
    file.create(legacy)
    stopifnot(grepl("qs2::qs_save", error_message(reader("BBB", dir)), fixed = TRUE))
  }

  # when both exist the qs2 file is used
  file.create(file.path(dir, "AAA_prepped.qs"))
  stopifnot(identical(reader("AAA", dir), survey))
})
