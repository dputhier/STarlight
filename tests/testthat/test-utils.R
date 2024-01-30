
library(testthat)


test_that("check_var function works correctly", {
  # Test NULL input
  expect_error(check_var(NULL), " should not be NULL.")

  # Test NA input
  expect_error(check_var(NA), " should not be NA.")

  # Test NaN input
  expect_error(check_var(NaN), " should not be nan.")

  # Test Infinite input
  expect_error(check_var(Inf), " should not be infinite.")

  # Test empty string input
  expect_error(check_var(""), "should not be an empty string.")

  # Test character input
  expect_silent(check_var("hello", type = "char"))

  # Test numeric input
  expect_silent(check_var(42, type = "num"))

  # Test logical input
  expect_silent(check_var(TRUE, type = "bool"))

  # Test integer input
  expect_silent(check_var(5L, type = "int"))

  # Test unknown format input
  expect_error(check_var("unknown", type = "unknown"), " has unknown format...")

  # Test character input with NaN
  expect_error(check_var("NaN", type = "char"), " should not be nan.")

  # Test numeric input with NA
  expect_error(check_var(c(1, NA, 3), type = "num"), " should not be NA.")

  # Test logical input with NULL
  expect_error(check_var(c(TRUE, NULL, FALSE), type = "bool"), " should not be NULL.")

  # Test integer input with empty string
  expect_error(check_var(5L, type = "int", x_name = ""), " should not be an empty string.")

  # Test character input with Infinite
  expect_error(check_var("hello", type = "char", x_name = "Inf"), "should not be infinite.")

  # Test numeric input with unknown format
  expect_error(check_var(c(1, 2, 3), type = "unknown"), " has unknown format...")

  # Test logical input with empty string
  expect_error(check_var(TRUE, type = "bool", x_name = ""), " should not be an empty string.")

  # Test integer input with NaN
  expect_error(check_var(5L, type = "int", x_name = "NaN"), " should not be nan.")
})


# Test cases for check_file function

test_that("check_file function works correctly", {
  # Test invalid file_path (NULL)
  expect_error(check_file(NULL), "Error: 'file_path' parameter cannot be NULL.")

  # Test writing to an existing file without force
  temp_existing_file_write <- tempfile(fileext = ".txt")
  file.create(temp_existing_file_write)
  expect_error(check_file(temp_existing_file_write, mode = "write"), "already exists. Use force = TRUE to overwrite.")

  # Test writing to an existing file with force
  expect_silent(check_file(temp_existing_file_write, mode = "write", force = TRUE))

  # Test writing to a non-existing file, create directory
  temp_new_file_write <- tempfile(fileext = ".txt")
  expect_silent(check_file(temp_new_file_write, mode = "write"))

  # Test reading from an existing file
  temp_existing_file_read <- tempfile(fileext = ".txt")
  file.create(temp_existing_file_read)
  expect_silent(check_file(temp_existing_file_read, mode = "read"))

  # Test reading from a non-existing file
  temp_non_existing_file_read <- tempfile(fileext = ".txt")
  expect_error(check_file(temp_non_existing_file_read, mode = "read"), "does not exist.")
})


