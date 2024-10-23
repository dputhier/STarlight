#' @title Create a report from an STGrid object.
#' @description
#' Create a report from an STGrid object.
#' @param corrplot_text_size Numeric, for the size of text label (variable names) for the coorplot diagram.
#' @examples
#' st_report(exp_1, spatial_image_params=list(ncol=4, features=feat_names(exp_1)[1:10]))
#'
st_report <- function(st_grid=NULL,
                      feat_list=NULL,
                      title="Spatial transcriptomics report",
                      subtitle="Sample 1",
                      output_format=c("html_document",
                                      "pdf_document",
                                      "beamer_presentation"),
                      out_file="/Users/puthier/Documents/basic_informations.html",
                      date=format(Sys.time(), '%d %B %Y'),
                      authors=list(c(Firstname="NA",
                                    Lastname="NA",
                                    Laboratory="NA",
                                    Affiliation="NA"
                                    )),
                      sample_info=c(Species=NA,
                                    Age=NA),
                      ncol=4,
                      image_height=3,
                      #rmd_path = paste0(system.file(package = "mypackage"), "/rmd/template.Rmd"),
                      rmd_path=file.path("/Users/puthier/Documents/git/project_dev/STarlight",
                                         "inst/rmd/sample.Rmd"),

                      corrplot_text_size=0.5,
                      hc_tree_params=list(class_nb = 5,
                                          offset=5),
                      plot_rip_k_params=list(ncol=4),
                      spatial_image_params=list(ncol=4,
                                                features=NULL),
                      cmp_counts_st_params=list(fill_color="#7845FF",
                                                transform="log10")){

  output_format <- match.arg(output_format)

  check_st_list(list(st_grid), feat_list = spatial_image_params$features)

  if(is.null(spatial_image_params$features)){
    spatial_image_params$features <- feat_names(st_grid)
  }

  rmarkdown::render(
    input = rmd_path,
    output_file = out_file,
    params = list(st_grid=st_grid,
                  title = title,
                  subtitle=subtitle,
                  date = date,
                  authors=authors,
                  sample_info=sample_info,
                  image_height=image_height,
                  corrplot_text_size=corrplot_text_size,
                  hc_tree_params=hc_tree_params,
                  spatial_image_params=spatial_image_params,
                  cmp_counts_st_params=cmp_counts_st_params,
                  plot_rip_k_params=plot_rip_k_params),
    encoding     = 'UTF-8'
  )

  #ext_var <- "test.html"
  #bookdown::render_book("index.Rmd", output_file = ext_var, clean_envir = FALSE)

}
