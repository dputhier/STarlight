#################################################################
##    set_verb_level
#################################################################
#' Set the verbosity level for the STarlight package
#'
#' This function sets the verbosity level for the STarlight package,
#' which controls the amount of information that is printed to the console by
#' the \code{\link{print_this_msg}} function. The verbosity level can be set to
#'  any non-negative integer, with higher values indicating more detailed output.
#'  By default, the verbosity level is set to 1.
#'
#' @param verbosity_value A non-negative integer indicating the verbosity level to be set.
#'
#' @return NULL
#'
#'
#' @examples
#' # Set verbosity level to 2
#' set_verb_level(2)
#'
#' # Set verbosity level to 0
#' set_verb_level(0)

# 0 : No message
# 1 : Display only INFO type message
# 2 : Display both INFO and DEBUG type message
#' @export
set_verb_level <- function(verbosity_value) {
  if (!is.null(verbosity_value) &
      verbosity_value >= 0 & is.numeric(verbosity_value)) {
    options(STarlight_verbosity = verbosity_value)
  }
}

#################################################################
##    get_verb_level()
#################################################################
#' Get the current verbosity level.
#'
#' This function get the verbosity level of the STarlight package which
#' controls the amount of information that is printed to the console by
#' the \code{\link{print_this_msg}} function.
#'
#'
#' @return A vector
#'
#'
#' @examples
#' get_verb_level()
#' @export
get_verb_level <- function() {
  if (is.null(unlist(options()["STarlight_verbosity"]))) {
    options(STarlight_verbosity = 1)
  }
  return(options()$STarlight_verbosity)
}

#################################################################
##    print_this_msg
#################################################################
#' Print a message based on the level of verbosity
#'
#' @param ... The variables that will be pasted in the message.
#' @param msg_type The type of message, one of "INFO", "DEBUG",
#' "WARNING" or "STOP" (may also accept "ERROR" which is
#' translated as "STOP").
#'
#' @return None
#'
#' @examples
#' opt_warn <- options()$warn
#' set_verb_level(1)
#' print_this_msg("Hello world!", msg_type = "INFO")
#' set_verb_level(3)
#' print_this_msg("Debugging message", msg_type = "DEBUG")
#' set_verb_level(0)
#' print_this_msg("Hello world!", msg_type = "INFO")
#' print_this_msg("Debugging message", msg_type = "DEBUG")
#' options(warn=-1)
#' print_this_msg("A warning message not displayed", msg_type = "WARNING")
#' options(warn=opt_warn)
#' @keywords internal
#' @export
print_this_msg <-
  function(...,
           msg_type = c("INFO", "DEBUG", "WARNING", "STOP", "ERROR")) {

    msg_type <- match.arg(msg_type)

    if(msg_type == "ERROR")
      msg_type <- "STOP"

    msg <- list(...)
    msg <- paste0(msg, collapse = " ")

    if (is.null(unlist(options()["STarlight_verbosity"]))) {
      options(STarlight_verbosity = 1)
    }

    if (msg_type == "DEBUG"){
      if (unname(unlist(options()["STarlight_verbosity"]) > 1))
        cat(paste("|-- DEBUG : ", msg, "\n"))

    }else if (msg_type == "WARNING"){
      warning("|-- WARNING : ", msg, call. = FALSE)

    }else if (msg_type == "STOP"){
      stop(paste0("|-- STOP : ", msg), call. = FALSE)

    }else{
      if (unname(unlist(options()["STarlight_verbosity"]) > 0))
        cat(paste("|-- INFO : ", msg, "\n"))
    }
  }

#################################################################
##    print_this_stat
#################################################################
#' @title Debugging statistics (vector, matrix or dataframe)
#' @description
#' Mostly a debugging function that will print some summary
#' statistics about a numeric vector, matrix or dataframe.
#'
#' @param msg The message to users.
#' @param data  The vector (numeric) for which the stats are to be produced
#' (a vector, )
#' @param msg_type The type of message, one of "INFO", "DEBUG", or "WARNING"
#' @param round_val Round the values in its first argument to the specified number
#'  of decimal. Set argument to -1 for no rounding
#' @return None
#'
#' @examples
#' opt_warn <- options()$warn
#' print_this_stat("My data", 1:10, msg_type="INFO")
#' set_verb_level(3)
#' print_this_stat("My data", matrix(rnorm(10), nc=2), msg_type="DEBUG")
#' set_verb_level(0)
#' print_this_stat("My data", matrix(rnorm(10), nc=2), msg_type="DEBUG")
#' @keywords internal
#' @export
print_this_stat <-
  function(msg,
           data,
           round_val = 2,
           msg_type = c("INFO", "DEBUG", "WARNING")) {

    msg_type <- match.arg(arg = msg_type, c("DEBUG", "WARNING", "INFO"))

    if (inherits(data, "data.frame")) {
      data <- as.matrix(data)
    }

    data <- as.vector(data)

    if (!is.numeric(data)) {

      print_this_msg("Can't print stats from non numeric object", msg_type = "WARNING")
      stats="No Statistics"
    }else{
      stats <- summary(data)
      names(stats) <- c("Min", "Q1", "Med", "Mean", "Q3", "Max")

      if (round_val > 0 & is.numeric(round_val)) {
        stats <- round(stats, round_val)
      }

      stats <- paste(names(stats), stats, sep = ":", collapse = " ")
    }

    print_this_msg(msg, ": ", stats, msg_type = msg_type)

  }

#################################################################
##    A simple function to create a random string
#################################################################
#' @title Generate a random string of letters and numbers
#' @description
#' This function generates a random string of 10 characters, consisting of 3 uppercase letters,
#' 4 digits, and 3 lowercase letters.
#'
#' @returns A character string of length 10, consisting of random letters and numbers
#'
#' @export
#'
#' @examples
#' rand_string()
#' @keywords internal
rand_string <- function() {
  v <- c(
    sample(LETTERS, 3, replace = TRUE),
    sample(0:9, 4, replace = TRUE),
    sample(letters, 3, replace = TRUE)
  )
  return(paste0(sample(v), collapse = ""))
}


# -------------------------------------------------------------------------
# Check file existence and optionally create directory   ------------------
# -------------------------------------------------------------------------
#' Check file existence and optionally create directory
#'
#' This function checks whether a file exists and, depending on the mode, performs additional checks.
#' If the mode is set to "write", it checks whether the file exists. If it exists and force is FALSE,
#' an error is raised. If force is TRUE, the existing file is overwritten or created if it doesn't exist.
#' If the mode is set to "read", it checks whether the file exists. If it doesn't, an error is raised.
#'
#' @param file_path The path to the file to be checked. If NULL, an error is raised.
#' @param mode The mode of operation. Choose from "write" or "read".
#' @param force Logical indicating whether to force the operation (default is FALSE).
#' @return NULL
#'
#' @examples
#' \dontrun{
#' # Example with a temporary file for writing
#' temp_file_write <- tempfile(fileext = ".txt")
#' check_this_file(temp_file_write, mode = "write")  # No error if file exists, force = TRUE
#'
#' # Example creating a file with file.create() and then checking it for writing
#' temp_file_create <- tempfile(fileext = ".txt")
#' file.create(temp_file_create)
#' check_this_file(temp_file_create, mode = "write", force = TRUE)  # No error if file exists, force = TRUE
#'
#' # Example checking a non-existing file for reading
#' temp_file_read <- tempfile(fileext = ".txt")
#' check_this_file(temp_file_read, mode = "read")  # Raises an error if file does not exist
#' }
#'
#' @export
check_this_file <- function(file_path, mode = c("write", "read"), force = FALSE) {

  check_this_var(x = file_path, type = "char")
  check_this_var(x = force, type = "bool")

  mode <- match.arg(mode)

  if (mode == "write") {
    if (file.exists(file_path)) {
      if (!force) {
        stop(paste("Error: File", file_path, "already exists. Use force = TRUE to overwrite.", sep = " "))
      } else {
        unlink(file_path)
        out <- file.create(file_path)
      }
    } else {
      dir_path <- dirname(file_path)
      if (!dir.exists(dir_path)) {
        print_this_msg("Creating directory:", dir_path)
        out <-  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      }
      print_this_msg("Creating file:", file_path)
      out <-  file.create(file_path)
    }
  } else {
    if (!file.exists(file_path)) {
      print_this_msg("File", file_path, "does not exist.", msg_type = "STOP")
    }
  }
}

# -------------------------------------------------------------------------
# check_this_var is intended to check variable type           ------------------
# -------------------------------------------------------------------------
#' Check variable type.
#'
#' This function checks the type of a variable and raises an error if it does not match the expected type.
#'
#' @param x The variable to be checked.
#' @param type The expected type. Choose from "char", "num", "bool", or "int".
#' @param x_name The name of the variable as a string.
#' @param null_accepted Whether NULL are accepted.
#' @param calling_function The name of the calling function (automatically obtained).
#' @examples
#' \dontrun{
#' # Example checking a character variable
#' check_this_var(x = "hello", type = "char")
#'
#' # Example checking a numeric variable
#' check_this_var(x = 42, type = "num")
#' }
#'
#' @keywords internal
#' @importFrom envnames get_fun_name
#' @export
check_this_var <- function(x,
                      type = c("char", "num", "bool", "int"),
                      x_name = deparse(substitute(x)),
                      null_accepted=FALSE,
                      calling_function = envnames::get_fun_name()) {



  fun_info <- paste0("(", calling_function, ")")

  type <- match.arg(type)

  if(length(x) > 1)
    print_this_msg(x_name, " should be of length 1.", msg_type="STOP")

  if(is.null(x)){
    if(!null_accepted)
      print_this_msg(x_name, " should not be NULL.", msg_type="STOP")
  }else if(is.nan(x)){
    print_this_msg(x_name, " should not be nan", msg_type="STOP")
  }else if(is.infinite(x)){
    print_this_msg(x_name, "should not be infinite", msg_type="STOP")
  }else if(is.na(x)){
    print_this_msg(x_name, " should not be NA", msg_type="STOP")
  }else if(x == ""){
    print_this_msg(x_name, "should not be an empty string.", msg_type="STOP")
  }

  if(!is.null(x)){
    if (type == "char") {
      if (!is.character(x))
        print_this_msg(x_name, " should be a character", fun_info, msg_type = "STOP")
    } else if (type == "int") {
      if(is.numeric(x)){
        if(!x %% 1 == 0){
          print_this_msg(x_name, " should be an integer", fun_info, msg_type = "STOP")
        }
      }else{
        print_this_msg(x_name, " should be an integer", fun_info, msg_type = "STOP")
      }

    } else if (type == "num") {
      if (!is.numeric(x))
        print_this_msg(x_name, " should be a numeric", fun_info, msg_type = "STOP")
    } else if (type == "bool") {
      if (!is.logical(x))
        print_this_msg(x_name, " should be a logical", fun_info, msg_type = "STOP")
    }else {
      print_this_msg(x_name, " has unknown format...", fun_info, msg_type = "STOP")
    }
  }
}

# -------------------------------------------------------------------------
# Show available methods --------------------------------------------------
# -------------------------------------------------------------------------
#' @title List the method for the STGrid object.
#' @description
#' List the method for the STGrid object.
#' @param class The name of the class. Default to "STGrid".
#' @param where The package name. Default to "STGrid".
#' @returns A vector with the available methods.
#' @importFrom methods showMethods
#' @importFrom utils capture.output
#' @examples
#' show_st_methods()
#' @export show_st_methods
show_st_methods <- function(class="STGrid",
                         where="STarlight"){
  class_method <- utils::capture.output(methods::showMethods(class=class,
                                                             where = paste0("package:", where)))

  class_method <- class_method[grep("Function", class_method, perl=T)]
  class_method <- gsub("Function: ([^ ]+) \\(.*", "\\1", class_method)

  return(class_method)
}

# -------------------------------------------------------------------------
# Reload the package     --------------------------------------------------
# -------------------------------------------------------------------------
#' @title Reload STarlight (used for development).
#' @description
#' Reload STarlight (used for development).
#' @returns NULL
#' @examples
#' reload_pac()
#' @export reload_pac
reload_pac <- function(){
  tryCatch(detach("package:STarlight", unload = TRUE))
 library(STarlight)
}

# -------------------------------------------------------------------------
# example_dataset()   ------------------------------------------------
# -------------------------------------------------------------------------
#' @title Load/download a Seurat or ClusterSet example dataset.
#' @description
#' Load/download a Seurat or ClusterSet example dataset.
#' @param dataset The name of the dataset.
#' @param timeout Set the timout for download (options(timeout=timeout))
#' @param force Reload the dataset even if already in globalEnv.
#' @returns Load the dataset.
#' @examples
#' example_dataset()
#' example_dataset("11284296/files/cmp_xen")
#' @export
example_dataset <- function(dataset=c("11284233/files/Xenium_Mouse_Brain_Coronal_7g",
                                      "11284296/files/cmp_xen"),
                                 timeout=NULL,
                                 force=FALSE){

  dataset <- match.arg(dataset)

  if(!is.null(timeout))
    options(timeout=timeout)

  file_data <- gsub(".*\\/", "", dataset)
  dir_path <- file.path(path.expand('~'), ".STarlight", "datasets")

  if(!dir.exists(dir_path)){
    print_this_msg(paste0("Creating a path for dataset installation: ",
                     dir_path),
              msg_type = "INFO")
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)

  }

  old_path <- getwd()
  setwd(dir_path)

  if(!file.exists(paste0(file_data, ".rda"))){
    utils::download.file(url=paste0("https://zenodo.org/record/",
                             dataset, ".rda"),
                  destfile = paste0(file_data, ".rda"))
  }

  dataset_short <- gsub(".*/", "", dataset)

  if(!dataset_short %in% ls(envir = globalenv()) | force){
    load(file.path(dir_path, paste0(file_data, ".rda")), envir = .GlobalEnv)
    print_this_msg(paste0("Dataset ", dataset, " has been loaded."),
              msg_type = "INFO")

  }else{
    print_this_msg(paste0("Dataset ", dataset, " was already loaded."),
              msg_type = "INFO")
  }

  setwd(old_path)

}


# -------------------------------------------------------------------------
# control_list()   ------------------------------------------------
# -------------------------------------------------------------------------

#' @title Generate a Control Gene List with Similar Distribution
#' @description
#' This function generates a control gene list with a distribution of expression values similar to the provided gene list.
#'
#' @param expression_value A named vector containing all gene expression values and their associated gene names.
#' @param gene_list A character vector of gene names for which a control gene list made of different genes with similar distribution is to be found.
#' @return A character vector of control gene names with similar expression distribution to the input \code{gene_list}.
#' @details The function works by calculating the absolute differences between the expression values of the genes in \code{gene_list}
#' and those of all other genes in the \code{expression_value} vector. For each gene in the \code{gene_list}, the gene with the
#' smallest difference in expression value from the remaining genes is selected. This process ensures that the selected control genes
#' have expression values that are closely matched to the target genes, thereby maintaining a similar distribution of expression values
#' in the control list. Ensuring the distribution are the same (e.g histogram, boxplot...) is important as the function may fail depending on the way the
#' gene list is sampled from the global distribution.
#' @examples
#' # Continuous values
#' set.seed(1)
#' a <- sort(c(rnorm(100, mean=0), rnorm(100, mean=4)))
#' names(a) <- paste0("gene_", 1:length(a))
#' b <- names(a[sample(1:100, size = 10, replace = FALSE)])
#' control <- control_list(a, b)
#' boxplot(list(a, a[b], a[control]))
#' all(!b %in% control)
#'
#' # Discrete values
#' set.seed(1)
#' a <- c(rpois(100, 10), rpois(100, 30))
#' names(a) <- paste0("gene_", 1:length(a))
#' b <- names(a[sample(50:130, size = 10, replace = FALSE)])
#' control <- control_list(a, b)
#' boxplot(list(a, a[b], a[control]))
#' all(!b %in% control)
#' @export
control_list <- function(expression_value=NULL,
                         gene_list=NULL){

  if(is.null(expression_value) | is.null(gene_list))
    print_this_msg("Please provide a set of expression values and a gene list",
                   msg_type = "STOP")

  if(length(unique(names(expression_value))) != length(names(expression_value)))
    print_this_msg("The provided vector of expression values contains duplicate names.",
                   msg_type = "STOP")

  if(length(unique(gene_list)) != length(gene_list))
    print_this_msg("The provided gene list contains duplicates.",
                   msg_type = "STOP")

  if(!all(sapply(gene_list, "nchar") > 0))
    print_this_msg("Some gene names are empty.",
                   msg_type = "STOP")

  if(!all(sapply(names(expression_value), "nchar") > 0))
    print_this_msg("Some gene names are empty.",
                   msg_type = "STOP")

  if(!all(gene_list %in% names(expression_value)))
    print_this_msg("Some genes from gene_list are not part of expression_value...",
                   msg_type = "STOP")

  other_genes <- names(expression_value[!names(expression_value) %in% gene_list])

  if(length(other_genes) < length(gene_list))
    print_this_msg("The number of remaining genes is not sufficient to find a matched distribution...",
                   msg_type = "STOP")

  diff <- abs(outer(expression_value[gene_list],
                    expression_value[other_genes], "-"))

  f <- vector()
  g <- vector()

  for(i in 1:nrow(diff)){
    if(length(f) > 0){
      if(length(diff[i, -f]) == 0)
        print_this_msg("Not enough remaining expression values.",
                       msg_type = "STOP")

      tmp <- which(diff[i, -f] == min(diff[i, -f]))

      if(length(tmp) > 1)
        tmp <- tmp[sample(length(tmp))][1]

      g[i] <- colnames(diff[i, -f, drop=FALSE])[tmp]
      f[i] <- which(colnames(diff) == g[i])

    }else{

      f[i] <- which(diff[i, ] == min(diff[i, ]))[1]
      g[i] <- colnames(diff)[f[i]]

    }
  }


  ks <- ks.test(expression_value[gene_list], expression_value[g])$p.value

  print_this_msg("Kolmogorov-Smirnov test p-value :",  ks)

  return(g)
}



# -------------------------------------------------------------------------
# st_gg_theming.           ------------------------------------------------
# -------------------------------------------------------------------------

#' Custom ggplot2 Theming
#'
#' This function defines a custom theme for ggplot2 plots with specific adjustments
#' to the legend key dimensions, legend text, axis titles, and strip text sizes.
#'
#' @return A ggplot2 theme object with predefined settings for legend, axis, and strip text.
#'
#' @details
#' The `st_gg_theming` function customizes the appearance of ggplot2 plots.
#'
#' @keywords internal
#' @importFrom ggplot2 theme
#' @export
st_gg_theming <- function(){

  gg_theme <- ggplot2::theme(
    legend.key.width = unit(0.10, "in"),
    legend.key.height = unit(0.10, "in"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 7))

  return(gg_theme)

}
