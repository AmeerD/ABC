# Regression tests for the qs2 migration. No Stan sampling is required.
local({
  reader <- getFromNamespace("read_prepped_survey", "abcR")
  directory <- tempfile("abcr-qs2-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  path <- file.path(directory, "survey_prepped.qs2")
  survey <- data.frame(
    age = c(10L, 25L, 45L),
    level = factor(c("primary", "primary", "upper secondary")),
    value = c(0, 1, NA_real_)
  )
  attr(survey, "source") <- "synthetic regression fixture"
  qs2::qs_save(survey, path)
  stopifnot(identical(reader(path), survey))

  expect_error <- function(expr, text) {
    message <- tryCatch({ force(expr); NA_character_ },
                        error = function(e) conditionMessage(e))
    stopifnot(!is.na(message), grepl(text, message, fixed = TRUE))
  }
  unlink(path)
  expect_error(reader(path), "Preprocessed survey file not found")
  legacy <- sub("\\.qs2$", ".qs", path)
  file.create(legacy)
  expect_error(reader(path), "renaming the file is not enough")
  # Prefer a valid qs2 file when a legacy file also exists.
  qs2::qs_save(survey, path)
  stopifnot(identical(reader(path), survey))
  # Reject invalid bytes rather than returning an object.
  writeBin(charToRaw("invalid qs2 fixture"), path)
  expect_error(reader(path), "")
})
