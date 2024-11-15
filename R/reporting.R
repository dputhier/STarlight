#' @title Create a report from an STGrid object.
#' @description
#' Create a report from an STGrid object.
#' @param corrplot_text_size Numeric, for the size of text label (variable names) for the coorplot diagram.
#' @examples
#' st_report(exp_1, spatial_image_params=list(ncol=4, features=feat_names(exp_1)[1:10]))
#' @export st_report
st_report <- function(st_grid_list = NULL,
                      feat_list = NULL,
                      title = "Spatial transcriptomics report",
                      subtitle = "Sample 1",
                      out_file = "/Users/puthier/Documents/basic_informations.html",
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
                      ncol = 4,
                      image_height = 0.8,
                      #rmd_dir = paste0(system.file(package = "mypackage"), "/rmd/template.Rmd"),
                      rmd_dir = file.path("/Users/puthier/Documents/git/project_dev/STarlight/inst/rmd/"),

                      corrplot_text_size = 0.4,
                      hc_tree_params = list(class_nb = 5, offset = 8),
                      plot_rip_k_params = list(ncol = 4),
                      spatial_image_params = list(ncol = 3, features =
                                                    NULL),
                      cmp_counts_st_params = list(fill_color = "#7845FF", transform =
                                                    "log10"),
                      rm_tmpdir = TRUE,
                      as_job = FALSE) {

  verb_level <- get_verb_level()

  if (inherits(st_grid_list, "STGrid"))
    st_grid_list <- list(st_grid_list)

  if (is.null(names(st_grid_list)))
    names(st_grid_list) <- 1:length(st_grid_list)

  check_st_list(st_grid_list, feat_list = spatial_image_params$features)

  print_this_msg("Retrieving full list of features.")
  all_feat <- table(unlist(lapply(st_grid_list, feat_names)))
  all_feat <- names(all_feat)[all_feat == length(st_grid_list)]

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

  tmp_dir <- file.path(tmp_dir, "rmd")


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
      #code_rmd  <- gsub("RDS_PATH", params_as_rds, x = code_rmd)
      #code_rmd  <- gsub("FIG_PATH", fig_path, x = code_rmd)

      writeLines(code_rmd, con = cur_rmd)

      print_this_msg("Preparation of rmd files for objects",
                     n,
                     "finished.",
                     msg_type = "DEBUG")

      n <- n + 1

      knited_file <- gsub(".rmd$" , ".knit.md", cur_rmd)
      knited_file_rmd <- gsub(".knit.md$" , ".knit.Rmd", knited_file)

      rmarkdown::render(
        input = cur_rmd,
        output_file = knited_file_rmd,
        params = list(
          all = list(
            st_grid = x,
            sample_name = sample_name,
            date = date,
            sample_info = sample_info,
            image_height = image_height,
            corrplot_text_size = corrplot_text_size,
            hc_tree_params = hc_tree_params,
            spatial_image_params = spatial_image_params,
            cmp_counts_st_params = cmp_counts_st_params,
            plot_rip_k_params = plot_rip_k_params,
            fig_path = fig_path
          )
        ),
        encoding     = 'UTF-8',
        run_pandoc = FALSE,
        quiet = TRUE
      )

      set_verb_level(verb_level)

      unlink(cur_rmd)




      file.rename(knited_file, knited_file_rmd)

      # return our new file name
      return(basename(knited_file_rmd))

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

  rmarkdown::render(
    input = file.path(tmp_dir, "index.Rmd"),
    output_file = file.path(tmp_dir, "index.knit.md"),
    params = list(
      all = list(
        st_grid_list = st_grid_list,
        experimenters = experimenters,
        sample_info = sample_info
      )
    ),
    encoding = 'UTF-8',
    run_pandoc = FALSE,
    quiet = TRUE,
    knit_root_dir=tmp_dir
  )

  set_verb_level(verb_level)

  unlink(file.path(tmp_dir, "index.Rmd"))

  file.rename(file.path(tmp_dir, "index.knit.md"), file.path(tmp_dir, "index.Rmd"))

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

  panel_set_code <- paste(readLines(file.path(tmp_dir,
                                              "panel_set_code.txt"),
                                    warn = FALSE),
                          collapse = "\n")

  print_this_msg("preparing panelsets")

  code_rmd <- readLines(file.path(tmp_dir, "index.Rmd"))
  code_rmd  <- gsub("PANELSET_CODE",  panel_set_code, x = code_rmd)
  writeLines(code_rmd, con = file.path(tmp_dir, "index.Rmd"))


  options(knitr.duplicate.label = "allow")

  bookdown::render_book(tmp_dir, quiet=TRUE)

  set_verb_level(verb_level)
}
