#' @title Create a report from an STGrid object.
#' @description
#' Create a report from an STGrid object.
#' @param corrplot_text_size Numeric, for the size of text label (variable names) for the coorplot diagram.
#' @examples
#' st_report(exp_1, spatial_image_params=list(ncol=4, features=feat_names(exp_1)[1:10]))
#' @export st_report
st_report <- function(st_grid_list=NULL,
                        feat_list=NULL,
                        title="Spatial transcriptomics report",
                        subtitle="Sample 1",
                        out_file="/Users/puthier/Documents/basic_informations.html",
                        author="Unknown",
                        date=format(Sys.time(), '%d %B %Y'),
                        experimenters=data.frame(Firstname=NA,
                                                 Lastname=NA,
                                                 Laboratory=NA,
                                                 Affiliation=NA),
                        sample_info=data.frame(Species=NA,
                                               Age=NA,
                                               row.names = row.names(st_grid_list)),
                        ncol=4,
                        image_height=0.8,
                        #qmd_dir = paste0(system.file(package = "mypackage"), "/rmd/template.Rmd"),
                        qmd_dir=file.path("/Users/puthier/Documents/git/project_dev/STarlight/inst/rmd/"),

                        corrplot_text_size=0.4,
                        hc_tree_params=list(class_nb = 5,
                                            offset=8),
                        plot_rip_k_params=list(ncol=4),
                        spatial_image_params=list(ncol=3,
                                                  features=NULL),
                        cmp_counts_st_params=list(fill_color="#7845FF",
                                                  transform="log10"),
                        rm_tmpdir=TRUE,
                        as_job=FALSE){


  if(inherits(st_grid_list, "STGrid"))
    st_grid_list <- list(st_grid_list)

  if(is.null(names(st_grid_list)))
    names(st_grid_list) <- 1:length(st_grid_list)

  check_st_list(st_grid_list, feat_list = spatial_image_params$features)

  print_this_msg("Retrieving full list of features.")
  all_feat <- table(unlist(lapply(st_grid_list, feat_names)))
  all_feat <- names(all_feat)[all_feat == length(st_grid_list)]

  if(length(all_feat) == 0)
    print_this_msg("There is no shared features between objets",
                 msg_type = "WARNING")

  if(is.null(spatial_image_params$features)){
    spatial_image_params$features <- all_feat
  }else{

    if(any(!spatial_image_params$features %in% all_feat)){
      print_this_msg("Some features (see spatial_image_params arg) were not found all the STGrid objects.")
    }
  }

  tmp_dir <- tempdir(check = FALSE)
  tmp_dir <- paste0(tmp_dir, format(Sys.time(), "%a-%b-%e-%H-%M-%S-%Y"))
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  print_this_msg("Created temporary directory:", msg_type = "DEBUG")
  print_this_msg(tmp_dir)

  print_this_msg("Copying file in temporary directory.", msg_type = "DEBUG")

  out <- file.copy(qmd_dir, tmp_dir, recursive = TRUE)

  tmp_dir <- file.path(tmp_dir, "rmd")


  # name of the parameterised report rmarkdown
  sample_qmd <- file.path(tmp_dir, "sample.qmd")

  print_this_msg("Preparing parameters for the report.")

  print_this_msg("Looping through parameterised Reports.")

  n <- 1

  quarto_files <- mapply(
    function(x, n, sample_name){

      print_this_msg("Preparing qmd files for objects", n, msg_type = "DEBUG")

      fig_path <- paste0("figure_", n, "_")

      params_as_list <- list(st_grid=x,
                             title = title,
                             sample_name=sample_name,
                             subtitle=subtitle,
                             date = date,
                             experimenters=experimenters,
                             sample_info=sample_info,
                             image_height=image_height,
                             corrplot_text_size=corrplot_text_size,
                             hc_tree_params=hc_tree_params,
                             spatial_image_params=spatial_image_params,
                             cmp_counts_st_params=cmp_counts_st_params,
                             plot_rip_k_params=plot_rip_k_params,
                             fig_path=fig_path)

      dir.create(file.path(tmp_dir, "rds"), recursive = TRUE, showWarnings = FALSE)
      params_as_rds <- file.path(tmp_dir, "rds", paste0("sample_", n, ".rds"))

      print_this_msg("Saving object", n, "as rds.", msg_type = "DEBUG")
      saveRDS(params_as_list, file = params_as_rds)

      cur_qmd <- gsub(".qmd$" , paste0("_", sprintf("%04d", n), ".qmd"), sample_qmd)

      print_this_msg("file : ", cur_qmd, msg_type = "DEBUG")

      code_qmd  <- readLines(sample_qmd)
      code_qmd  <- gsub("RDS_PATH", params_as_rds, x = code_qmd)
      code_qmd  <- gsub("FIG_PATH", fig_path, x = code_qmd)

      writeLines(code_qmd, con=cur_qmd)

      print_this_msg("Preparation of qmd files for objects", n, "finished.", msg_type = "DEBUG")

      n <- n + 1

      # return our new file name
      return(basename(cur_qmd))

    },
    x = st_grid_list,
    n = 1:length(st_grid_list),
    sample_name=names(st_grid_list)
  )

  print_this_msg("Preparing index.qmd file.")
  params_as_rds <- file.path(tmp_dir, "rds", "all_sample.rds")
  saveRDS(list(st_grid_list=st_grid_list,
               sample_info=sample_info,
               image_height=image_height,
               experimenters=experimenters,
               fig_path="figure_all"), file = params_as_rds)
  code_qmd  <- readLines(file.path(tmp_dir, "index.qmd"))
  code_qmd  <- gsub("RDS_PATH", params_as_rds, x = code_qmd)
  code_qmd  <- gsub("FIG_PATH", "all_", x = code_qmd)

  writeLines(code_qmd, con=file.path(tmp_dir, "index.qmd"))

  quarto_files <- append(quarto_files, "index.qmd", after=0)
  quarto_yaml <- file.path(tmp_dir, "params.yml")
  quarto_yaml  <- readLines(quarto_yaml)

  quarto_yaml  <- gsub("CHAPTERS_VAR",
               paste('   - text: "',
                     c("Introduction", names(st_grid_list)),'"',
                     "\n     file:",
                     quarto_files, "\n",
                     collapse = ""),
               quarto_yaml)
  quarto_yaml  <- gsub("AUTHOR",
                       author,
                       quarto_yaml)

  quarto_yaml  <- gsub("DATE",
                       Sys.Date(),
                       quarto_yaml)


  writeLines(quarto_yaml, con=file.path(tmp_dir, "_quarto.yml"))

  unlink(sample_qmd)
  unlink(file.path(tmp_dir, "sample.Rmd"))
  unlink(file.path(tmp_dir, "params.yml"))

  quarto::quarto_render(tmp_dir, as_job=as_job)

}
