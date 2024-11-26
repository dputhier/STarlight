#' @title Create a report from an STGrid object.
#' @description
#' Create a report from an STGrid object.
#' @param st_grid_list A list of STGrid objects.
#' @param st_grid_list A list of feat_list to restrict the analysis.
#' @param title A title for the report.
#' @param subtitle A subtitle for the report.
#' @param subtitle A subtitle for the report.
#' @param out_dir A directory where to store the bookdown output.
#' @param author A character string corresponding to one or severe al authors.
#' @param date A date.
#' @param experimenters A data.frame giving informations about experimenters (any number of row/columns accepted).
#' @param sample_info A data.frame giving informations about samples (any number of row/columns accepted).
#' @param image_height A multiplier for image size.
#' @param rmd_dir A path where to find the templates for creating the book.
#' @param corrplot_params Parameters for corrplot() function.
#' @param hc_tree_params Parameters for hc_tree() function.
#' @param plot_rip_k_params Parameters for plot_rip_k() function.
#' @param spatial_image_params Parameters for spatial_image() function.
#' @param cmp_counts_st_params Parameters for cmp_counts_st() function.
#' @param rm_tmpdir Whether to delete temporary directory.
#' @param quiet Whether to run bookdown::render_book() quietly.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_2 <- xen[x_bins, y_bins]
#' st_list <- list(xen_1=xen_1, xen_2=xen_2)
#' st_report(st_list)
#'
#' @export st_report
st_report <- function(st_grid_list = NULL,
                      feat_list = NULL,
                      title = "Spatial transcriptomics report",
                      subtitle = "An example experiment",
                      out_dir = file.path(fs::path_home(), "st_book"),
                      author = "Unknown",
                      date = format(Sys.time(), '%d %B %Y'),
                      experimenters = data.frame(
                        Firstname = NA,
                        Lastname = NA,
                        Laboratory = NA,
                        Affiliation = NA
                      ),
                      sample_info = data.frame(Species = NA,
                                               Age = NA,
                                               row.names = row.names(st_grid_list)),
                      image_height = 3,
                      rmd_dir = file.path(system.file(package = "STarlight"), "rmarkdown"),
                      corrplot_params = list(type="upper",
                                             order="original",
                                             col=RColorBrewer::brewer.pal(n=8, name="RdYlBu"),
                                             tl.cex=0.4, diag = FALSE),
                      hc_tree_params = list(class_nb = 20, offset = 8),
                      plot_rip_k_params = list(ncol = 4),
                      spatial_image_params = list(ncol = 5, features =
                                                    NULL),
                      cmp_counts_st_params = list(fill_color = "#7845FF", transform =
                                                    "log10"),
                      rm_tmpdir = TRUE,
                      quiet=FALSE) {

  verb_level <- get_verb_level()

  if (inherits(st_grid_list, "STGrid"))
    st_grid_list <- list(st_grid_list)

  if (is.null(names(st_grid_list)))
    names(st_grid_list) <- 1:length(st_grid_list)

  check_st_list(st_grid_list, feat_list = spatial_image_params$features)

  print_this_msg("Retrieving full list of features.")
  all_feat <- table(unlist(lapply(st_grid_list, feat_names)))
  all_feat <- names(all_feat)[all_feat == length(st_grid_list)]

  if(hc_tree_params$class_nb > length(all_feat))
    hc_tree_params$class_nb <- length(all_feat)

  if (length(all_feat) == 0)
    print_this_msg("There is no shared features between objets", msg_type = "WARNING")

  if (is.null(spatial_image_params$features)) {
    spatial_image_params$features <- all_feat
  } else{
    if (any(!spatial_image_params$features %in% all_feat)) {
      print_this_msg(
        "Some features (see spatial_image_params arg) were not found all the STGrid objects."
      )
    }
  }

  tmp_dir <- tempdir(check = FALSE)
  tmp_dir <- paste0(tmp_dir, format(Sys.time(), "%a-%b-%e-%H-%M-%S-%Y"))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  print_this_msg("Created temporary directory:", msg_type = "DEBUG")
  print_this_msg(tmp_dir)

  print_this_msg("Copying file in temporary directory.", msg_type = "DEBUG")

  out <- file.copy(rmd_dir, tmp_dir, recursive = TRUE)

  tmp_dir <- file.path(tmp_dir, "rmarkdown")

  # name of the parameterised report rmarkdown
  sample_rmd <- file.path(tmp_dir, "sample.rmd")

  print_this_msg("Preparing parameters for the report.")

  print_this_msg("Looping through parameterised Reports.")

  n <- 1

  all_knited_files <- mapply(
    function(x, n, sample_name) {

      print_this_msg("Preparing rmd files for objects", n, msg_type = "DEBUG")
      fig_path <- paste0("figure_", n, "_")

      cur_rmd <- gsub(".rmd$" , paste0("_", sprintf("%04d", n), ".rmd"), sample_rmd)

      print_this_msg("file : ", cur_rmd, msg_type = "DEBUG")

      code_rmd  <- readLines(sample_rmd)
      code_rmd  <- gsub("ST_GRID", paste0("st_grid_list[[", n, "]]"), x = code_rmd)
      code_rmd  <- gsub("SAMPLE_NAME", sample_name, x = code_rmd)
      code_rmd  <- gsub("FIG_PATH", fig_path, x = code_rmd)
      code_rmd  <- gsub("SMP_NUMBER", n, x = code_rmd)
      writeLines(code_rmd, con = cur_rmd)

      print_this_msg("Preparation of rmd files for objects", n, "finished.", msg_type = "DEBUG")

      n <- n + 1

      # return our new file name
      return(basename(cur_rmd))

    },
    x = st_grid_list,
    n = 1:length(st_grid_list),
    sample_name = names(st_grid_list)
  )

  print_this_msg("Deleting sample.Rmd.")

  unlink(sample_rmd)

  print_this_msg("Preparing index.rmd file.")

  code_rmd  <- readLines(file.path(tmp_dir, "index.Rmd"))
  code_rmd  <- gsub("REPORT_TITLE", title, x = code_rmd)
  code_rmd  <- gsub("REPORT_AUT", title, x = code_rmd)
  code_rmd  <- gsub("REPORT_SUBTITLE", subtitle, x = code_rmd)
  code_rmd  <- gsub("REPORT_DATE", date, x = code_rmd)
  code_rmd  <- gsub("FIG_PATH", "all_", x = code_rmd)
  writeLines(code_rmd, con = file.path(tmp_dir, "index.Rmd"))


  print_this_msg("preparing _bookdown.yml")

  all_knited_files_merged <- paste(sapply(all_knited_files, shQuote), collapse = ", ")

  code_yml  <- readLines(file.path(tmp_dir, "_bookdown.yml"), warn = FALSE)

  code_yml  <- gsub(
    "rmd_files:",
    paste0(
      "rmd_files: [ 'index.Rmd', ",
      all_knited_files_merged,
      "]"
    ),
    x = code_yml
  )

  writeLines(code_yml, con = file.path(tmp_dir, "_bookdown.yml"))
  code_yml  <- readLines(file.path(tmp_dir, "_bookdown.yml"), warn = FALSE)

  options(knitr.duplicate.label = "allow")

  bookdown::render_book(tmp_dir, quiet=quiet)
  set_verb_level(verb_level)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(tmp_dir, "_book"), out_dir, recursive=TRUE)

  print_this_msg("Results have been copied to:", out_dir)

  if(rm_tmpdir)
    unlink(tmp_dir, recursive=TRUE)

}
