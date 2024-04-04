library(testthat)
test_that("check_this_var function works correctly", {

  library(stcompr)
  set_verb_level(0)

  # Test NULL input
  expect_error(check_this_var(NULL), " should not be NULL.")

  # length
  expect_error(check_this_var(x=c("bla", "foo"), type="int", null_accepted = TRUE), "should be of length 1")

  # Test NA input
  expect_error(check_this_var(NA), "Error: |-- STOP : NA  should not be NA")

  # Test NaN input
  x <- NaN
  expect_error(check_this_var(x), "Error: |-- STOP : x  should not be nan")

  # Test Infinite input
  x <- Inf
  expect_error(check_this_var(x), "Error: |-- STOP : x should not be infinite")

  # Test empty string input
  x <- ""
  expect_error(check_this_var(x), "Error: |-- STOP : x should not be an empty string.")

  # Test character input
  x <- "hello"
  expect_silent(check_this_var(x, type = "char"))

  # Test numeric input
  expect_silent(check_this_var(42, type = "num"))

  # Test logical input
  expect_silent(check_this_var(TRUE, type = "bool"))

  # Test integer input
  expect_silent(check_this_var(5L, type = "int"))

  # Test unknown format input
  expect_error(check_this_var("unknown", type = "unknown"), "'arg' should be one of ")

  # Test character input with NaN
  x <- NaN
  expect_error(check_this_var(x, type = "char"), "x  should not be nan")

  # Test numeric input with NA
  x <- c(1, NA, 3)
  expect_error(check_this_var(x, type = "num"), "x  should be of length 1")

  # Test logical input with NULL
  x <- TRUE
  expect_silent(check_this_var(x, type = "bool"))

  })


# Test cases for check_this_file function

test_that("check_this_file function works correctly", {

  library(stcompr)
  library(testthat)
  set_verb_level(0)

  # Test invalid file_path (NULL)
  expect_error(check_this_file(NULL), "Error: |-- STOP : file_path  should not be NULL")

  # Test writing to an existing file without force
  temp_existing_file_write <- tempfile(fileext = ".txt")
  out <- file.create(temp_existing_file_write)
  expect_error(check_this_file(temp_existing_file_write, mode = "write"), "already exists. Use force = TRUE to overwrite.")

  # Test writing to an existing file with force
  expect_silent(check_this_file(temp_existing_file_write, mode = "write", force = TRUE))

  # Test writing to a non-existing file, create directory
  temp_new_file_write <- tempfile(fileext = ".txt")
  set_verb_level(0)
  expect_no_error(check_this_file(temp_new_file_write, mode = "write", force = TRUE))
  set_verb_level(1)

  # Test reading from an existing file
  temp_existing_file_read <- tempfile(fileext = ".txt")
  tf <- file.create(temp_existing_file_read)
  expect_silent(check_this_file(temp_existing_file_read, mode = "read"))

  # Test reading from a non-existing file
  temp_non_existing_file_read <- tempfile(fileext = ".txt")
  expect_error(check_this_file(temp_non_existing_file_read, mode = "read"), "does not exist.")
})


