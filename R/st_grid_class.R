# -------------------------------------------------------------------------
## Spatial transcriptomic analysis class (STGrid) definition
# -------------------------------------------------------------------------

#' Spatial transcriptomic analysis class (STGrid).
#'
#' A class for representing spatial transcriptomic analysis results as x/y coordinates.
#' This class stores transcript molecules or cell coordinates as a grid.
#'
#' @slot coord A data frames representing spatial coordinates (empty by default).
#' @slot bin_mat A data frames representing the matrix (empty by default).
#' @slot y_max Numeric value representing the maximum y-coordinate (default is 0).
#' @slot y_min Numeric value representing the minimum y-coordinate (default is 0).
#' @slot x_min Numeric value representing the minimum x-coordinate (default is 0).
#' @slot x_max Numeric value representing the maximum x-coordinate (default is 0).
#' @slot path Character value representing the file path to the original data.
#' @slot method The technology used.
#' @slot meta A list containing meta-information about bins.
#' @slot bin_size Numeric value representing the bin size.
#' @slot bin_x Character vector containing the names of bins/windows along the x-axis.
#' @slot bin_y Character vector containing the names of bins/windows along the y-axis.
#' @slot ripley_k_function The result of Ripley's K function computation.
#' @slot control A regular expression indicating the string motif related to controls.
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g
#' show_st_methods(class = "STGrid")
#' @export
setClass(
  "STGrid",
  slots = c(
    coord = "data.frame",
    bin_mat = "data.frame",
    y_max = "numeric",
    y_min = "numeric",
    x_min = "numeric",
    x_max = "numeric",
    path = "character",
    method = "character",
    meta = "data.frame",
    bin_size = "numeric",
    bin_x = "character",
    bin_y = "character",
    ripley_k_function = "data.frame",
    moran_index = "data.frame",
    control = "character"
  ),
  prototype = list(
    coord = data.frame(),
    bin_mat = data.frame(),
    y_max = 0,
    y_min = 0,
    x_min = 0,
    x_max = 0,
    path = character(),
    method = character(),
    meta = data.frame(),
    bin_size = 0,
    bin_x = character(),
    bin_y = character(),
    ripley_k_function = data.frame(),
    moran_index = data.frame(),
    control = character()
  )
)

# -------------------------------------------------------------------------
#      NCOL/NROW/DIM METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------
#' @title
#' ncol
#' @description
#' The number of columns of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
setMethod(
  f = "ncol",
  signature = "STGrid",
  definition = function(x)
    length(bin_y(x))
)

#' @title
#' nrow
#' @description
#' The number of rows of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
setMethod(
  f = "nrow",
  signature = "STGrid",
  definition = function(x)
    length(bin_x(x))
)

#' @title Dimension of a STGrid object.
#' dim
#' @description
#' The dimension of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' dim(Xenium_Mouse_Brain_Coronal_7g)
#'
#' @keywords internal
setMethod("dim", signature(x = "STGrid"),
          function(x)
            c(length(bin_x(x)), length(bin_y(x)))
)

#' @title Column names of a STGrid object.
#' @description
#' The column names of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' col_names(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setGeneric(
  name = "col_names",
  def = function(x)
    standardGeneric("col_names")
)

#' @title Column names of a STGrid object.
#' @description
#' The column names of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' col_names(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setMethod(
  f = "col_names",
  signature = "STGrid",
  definition = function(x)
    x@bin_x
)

#' @title Check all specified features exist in  an STGrid object.
#' @description
#' Check all specified features exist in  an STGrid object.
#' @param object The STGrid object
#' @param feat_list The list of features to be tested.
#' @param no_meta Whether meta features should not be considered. Default is FALSE.
#' @examples
#' example_dataset()
#' check_features_exist(Xenium_Mouse_Brain_Coronal_7g, "Ano1")
#' check_features_exist(Xenium_Mouse_Brain_Coronal_7g, "foo")
#' @keywords internal
#' @export
setGeneric(
  name = "check_features_exist",
  def = function(object, feat_list, no_meta=FALSE)
    standardGeneric("check_features_exist")
)

#' @title Check all specified features exist in  an STGrid object.
#' @description
#' Check all specified features exist in  an STGrid object.
#' @param object The STGrid object
#' @param feat_list The list of features to be tested.
#' @param no_meta Whether meta features should not be considered. Default is FALSE.
#' @examples
#' example_dataset()
#' check_features_exist(Xenium_Mouse_Brain_Coronal_7g, "Ano1")
#' check_features_exist(Xenium_Mouse_Brain_Coronal_7g, "foo")
#' @export
setMethod(
  f = "check_features_exist",
  signature = "STGrid",
  definition = function(object, feat_list, no_meta=FALSE)
    if(no_meta){
      return(all(feat_list %in% c(feat_names(object))))
    }else{
      return(all(feat_list %in% c(feat_names(object), meta_names(object))))
    }

)


#' @title Row names of a STGrid object.
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' row_names(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setGeneric("row_names",
           function(x)
             standardGeneric("row_names"))

#' @title Row names of a STGrid object.
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' example_dataset()
#' row_names(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("row_names", "STGrid",
          function(x)
            x@bin_y
)

#' @title X bins of a STGrid object.
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' bin_x(Xenium_Mouse_Brain_Coronal_7g)
#' @export bin_x
setGeneric("bin_x",
           function(x)
             standardGeneric("bin_x"))

#' @title X bins of a STGrid object.
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' bin_x(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("bin_x", "STGrid",
          function(x)
            x@bin_x)

#' @title Y bins of a STGrid object.
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' bin_y(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setGeneric("bin_y",
           function(x)
             standardGeneric("bin_y"))

#' @title Y bins of a STGrid object.
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' bin_y(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("bin_y", "STGrid",
          function(x)
            x@bin_y)

#' @title Number of bins (x axis) of a STGrid object.
#' @description
#' Number of bins (x axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' nbin_x(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setGeneric("nbin_x",
           function(x)
             standardGeneric("nbin_x"))

#' @title Number of bins (x axis) of a STGrid object.
#' @description
#' Number of bins (x axis) of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' nbin_x(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setMethod("nbin_x", "STGrid",
          function(x)
            length(bin_x(x)))

#' @title Number of bins (y axis) of a STGrid object.
#' @description
#' Number of bins (y axis) of a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' nbin_y(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setGeneric("nbin_y",
           function(x)
             standardGeneric("nbin_y"))

#' @title Number of bins (y axis) of a STGrid object.
#' @description
#'  Number of bins (y axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' nbin_y(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("nbin_y", "STGrid",
          function(x)
            length(bin_y(x)))
#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @examples
#' example_dataset()
#' head(coord(Xenium_Mouse_Brain_Coronal_7g))
#' @keywords internal
#' @export
setGeneric("coord",
           function(object,
                    as_factor = FALSE)
             standardGeneric("coord"))

#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @examples
#' example_dataset()
#' head(coord(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setMethod("coord", "STGrid",
          function(object,
                   as_factor = FALSE) {
            if (!as_factor) {
              return(object@coord)
            } else{
              tmp <- object@coord
              tmp$bin_x <- factor(tmp$bin_x,
                                  levels = bin_x(object),
                                  ordered = TRUE)
              tmp$bin_y <- factor(tmp$bin_y,
                                  levels = bin_y(object),
                                  ordered = TRUE)
              return(tmp)
            }
          })

#' @title The binned matrix stored in a STGrid object
#' @description
#' The binned matrix stored in a STGrid object
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @param feat_list Whether to subset to some features.
#' @param melt_tab Whether to melt.
#' @param del_bin Whether to delete the bin_x, bin_y columns.
#' @param pseudo_count Whether to add a pseudo_count.
#' @param transform Whether the count should be transformed.
#' @param axis_as_number If melt_tab is TRUE, whether bin_x and bin_y should be transformed to numerics.
#' @keywords internal
#' @examples
#' example_dataset()
#' head(bin_mat(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setGeneric("bin_mat",
           function(object,
                    as_factor = FALSE,
                    feat_list = character(),
                    melt_tab = FALSE,
                    del_bin = FALSE,
                    pseudo_count=0,
                    transform = c("None", "log2", "log10", "log"),
                    axis_as_number=FALSE)
             standardGeneric("bin_mat"))

#' @title The binned matrix stored in a STGrid object.
#' @description
#' The binned matrix stored in a STGrid object.
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @param feat_list Whether to subset to some features.
#' @param melt_tab Whether to melt.
#' @param del_bin Whether to delete the bin_x, bin_y columns.
#' @param pseudo_count Whether to add a pseudo_count.
#' @param transform Whether the count should be transformed.
#' @param axis_as_number If melt_tab is TRUE, whether bin_x and bin_y should be transformed to numerics.
#' @examples
#' example_dataset()
#' head(bin_mat(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setMethod("bin_mat", "STGrid",
          function(object,
                   as_factor = FALSE,
                   feat_list = character(),
                   melt_tab = FALSE,
                   del_bin = FALSE,
                   pseudo_count=0,
                   transform = c("None", "log2", "log10", "log"),
                   axis_as_number=FALSE) {

            transform <- match.arg(transform)

            if (length(feat_list) == 0) {
              feat_list <- feat_names(object)
            }

            if(any(colnames(object@meta) %in% feat_list)){
              print_this_msg("Using a feature from meta slot.", msg_type = "DEBUG")
              this_bin_mat <- object@bin_mat
              for(i in feat_list[feat_list %in% colnames(object@meta)])
                this_bin_mat[[i]] <- object@meta[[i]]
            }else{
              this_bin_mat <- object@bin_mat
            }

            if (del_bin) {
              this_bin_mat <- this_bin_mat[, feat_list, drop=FALSE]
            } else{
              this_bin_mat <- this_bin_mat[, c("bin_x", "bin_y", feat_list), drop=FALSE]
            }

            if (as_factor) {
              this_bin_mat$bin_x <- factor(this_bin_mat$bin_x,
                                           levels = bin_x(object),
                                           ordered = TRUE)

              this_bin_mat$bin_y <- factor(this_bin_mat$bin_y,
                                           levels = bin_y(object),
                                           ordered = TRUE)
            }

            if(transform == "log2"){
              this_bin_mat[, feat_list] <- log2(this_bin_mat[, feat_list] + pseudo_count)
            }else if(transform == "log10"){
              this_bin_mat[, feat_list] <- log10(this_bin_mat[, feat_list] + pseudo_count)
            }else if(transform == "log"){
              this_bin_mat[, feat_list] <- log(this_bin_mat[, feat_list] + pseudo_count)
            }

            if (melt_tab) {
              print_this_msg("Melting data.frame...", msg_type = "DEBUG")
              if (del_bin) {
                this_bin_mat <- reshape2::melt(this_bin_mat)
                colnames(this_bin_mat) <- c("feature",
                                            "value")
              } else{
                this_bin_mat <- reshape2::melt(this_bin_mat,
                                               id.vars = c("bin_x", "bin_y"))
                colnames(this_bin_mat) <- c("bin_x",
                                            "bin_y",
                                            "feature",
                                            "value")
              }

            }

          if(axis_as_number){
            if(!del_bin){
              if(melt_tab){

                bin_x_2_num <- setNames(1:length(bin_x(object)),
                                        bin_x(object))
                bin_y_2_num <- setNames(1:length(bin_y(object)),
                                        bin_y(object))

                this_bin_mat$bin_x <- bin_x_2_num[this_bin_mat$bin_x]
                this_bin_mat$bin_y <- bin_y_2_num[this_bin_mat$bin_y]

                if(as_factor){
                  this_bin_mat$bin_x <- factor(this_bin_mat$bin_x,
                                                  levels=1:length(bin_x(object)),
                                                  ordered=TRUE)
                  this_bin_mat$bin_y <- factor(this_bin_mat$bin_y,
                                                  levels=1:length(bin_y(object)),
                                                  ordered=TRUE)
                }

              }
            }
          }

            return(this_bin_mat)
          })

#' @title The number of x/y items (e.g. molecules or cells) stored in a STGrid object
#' @description
#' Returns the number of x/y items (e.g. molecules or cells) stored in a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' nb_items(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setGeneric("nb_items",
           function(x)
             standardGeneric("nb_items"))

#' @title The number of x/y items (e.g. molecules or cells) stored in a STGrid object
#' @description
#' Returns the number of x/y items (e.g. molecules or cells) stored in a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' nb_items(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("nb_items", "STGrid",
          function(x)
            nrow(x@coord))

#' @title The number of items (molecules/cells) per features (gene/cell-type).
#' @description
#' Returns The number of items (molecules/cells) per features (gene/cell-type).
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' tab(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setGeneric("tab",
           function(x)
             standardGeneric("tab"))

#' @title The number of items (molecules/cells) per features (gene/cell-type).
#' @description
#' Returns The number of items (molecules/cells) per features (gene/cell-type).
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' tab(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("tab", "STGrid",
          function(x)
            sort(table(x@coord$feature)))


#' @title The number of features stored in a STGrid object
#' @description
#' The number of features stored in a STGrid object
#' @param object The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' nb_feat(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setGeneric("nb_feat",
           function(object) {
             standardGeneric("nb_feat")
           })


#' @title The number of features stored in a STGrid object
#' @description
#' The number of features stored in a STGrid object
#' @param object The STGrid object
#' @examples
#' example_dataset()
#' nb_feat(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("nb_feat", signature(object = "STGrid"),
          function(object) {
            length(feat_names(object))
          })

#' @title Write object coordinates.
#' @description
#' Write object coordinates.
#' @param object The STGrid object
#' @param file_path The path where to store the file.
#' @keywords internal
#' @examples
#' example_dataset()
#' tmp_f <- tempfile()
#' write_coord(Xenium_Mouse_Brain_Coronal_7g, file_path=tmp_f)
#' @export
setGeneric("write_coord",
           function(object=NULL, file_path=NULL) {
             standardGeneric("write_coord")
           })


#' @title Write object coordinates.
#' @description
#' Write object coordinates.
#' @param object The STGrid object
#' @param file_path The path where to store the file.
#' @examples
#' example_dataset()
#' tmp_f <- tempfile()
#' write_coord(Xenium_Mouse_Brain_Coronal_7g, file_path=tmp_f)
#' @export
setMethod("write_coord", signature(object = "STGrid"),
          function(object=NULL,
                   file_path=NULL) {
            if(is.null(object))
              print_this_msg("Please provide an STGrid object.")

            if(is.null(file_path)){
              print_this_msg("Please provide a file path.")
            }else{
              check_this_file(file_path, mode = "write", force = TRUE)
            }
            print_this_msg("Writing to ", file_path)
            data.table::fwrite(object@coord,
                               file = file_path,
                               quote = FALSE, sep = "\t")
})

#' @title The features stored in a STGrid object
#' @description
#' The features stored in a STGrid object
#' @param object The STGrid object
#' @param del_control Whether to delete controls.
#' @keywords internal
#' @examples
#' example_dataset()
#' feat_names(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setGeneric("feat_names",
           function(object,
                    del_control = FALSE)
             standardGeneric("feat_names"))

#' @title The features stored in an STGrid object
#' @description
#' The features stored in a STGrid object
#' @param object The STGrid object.
#' @param del_control Whether to delete controls.
#' @examples
#' example_dataset()
#' feat_names(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("feat_names", signature(object = "STGrid"),
          function(object,
                   del_control = FALSE) {

            fn <- grep(
                "^bin_[xy]$",
                colnames(object@bin_mat),
                invert = TRUE,
                perl = TRUE,
                value = TRUE
            )


            if (del_control) {
              fn <-
                fn[!fn %in% grep(object@control, fn, value = TRUE)]
            }

            fn
          })

#' @title The names of meta informations stored in a STGrid object
#' @description
#' The names of meta informations stored in a STGrid object
#' @param object The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' meta_names(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setGeneric("meta_names",
           function(object)
             standardGeneric("meta_names"))

#' @title The names of meta informations stored in a STGrid object
#' @description
#' The names of meta informations stored in a STGrid object
#' @param object The STGrid object
#' @examples
#' example_dataset()
#' meta_names(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("meta_names",
          signature(object = "STGrid"),
          function(object) {
            colnames(object@meta)
          })

#' @title The size of the bins stored in an STGrid object
#' @description
#' The size of the bins stored in an STGrid object
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' bin_size(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export bin_size
setGeneric("bin_size",
           function(x)
             standardGeneric("bin_size"))

#' @title The size of the bins stored in an STGrid object
#' @description
#' The size of the bins stored in an STGrid object
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' bin_size(Xenium_Mouse_Brain_Coronal_7g)
#' @export bin_size
setMethod("bin_size", "STGrid",
          function(x)
            x@bin_size)

#' @title Extract meta slot of an STGrid object
#' @description
#' Extract meta slot of an STGrid object
#' @param object The STGrid object
#' @examples
#' example_dataset()
#' head(meta(Xenium_Mouse_Brain_Coronal_7g))
#' @keywords internal
#' @export meta
setGeneric("meta",
           function(object)
             standardGeneric("meta"))

#' @title Extract meta slot of an STGrid object
#' @description
#' Extract meta slot of an STGrid object
#' @param object The STGrid object
#' @examples
#' example_dataset()
#' head(meta(Xenium_Mouse_Brain_Coronal_7g))
#' @export meta
setMethod("meta", "STGrid",
          function(object)
            object@meta)
#' @title Remove control (Blank-*) features from a STGrid object.
#' @description
#'  Remove control (Blank-*) features.
#' @param object The STGrid object
#' @param regexp A regexp to search for internal controls.
#' @examples
#' # Here the example_dataset does not contain any
#' # control...
#' example_dataset()
#' rm_controls(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setGeneric("rm_controls",
           function(object,
                    regexp = NULL)
             standardGeneric("rm_controls"))

#' @title Remove control (Blank-*) features from a STGrid object.
#' @description
#'  Remove control (Blank-*) features.
#' @param object The STGrid object
#' @param regexp A regexp to search for internal controls.
#' @examples
#' # Here the example_dataset does not contain any
#' # control...
#' example_dataset()
#' rm_controls(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("rm_controls", "STGrid",
          function(object,
                   regexp = NULL) {
            if (is.null(regexp)) {
              regexp <- object@control
            }
            pos <- -grep(regexp, feat_names(object), perl = TRUE)
            if(length(pos) > 0){
              return(object[pos, ])
            }else{
              return(object)
            }

          })

#' @title The number of items (e.g molecules) per feature of an STGrid object.
#' @description
#' Returns the number of items (e.g molecules) per feature of an STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' table_st(Xenium_Mouse_Brain_Coronal_7g)
#' @keywords internal
#' @export
setGeneric("table_st",
           function(x)
             standardGeneric("table_st"))

#' @title The number of items (e.g molecules) per feature of an STGrid object.
#' @description
#' Returns the number of items (e.g molecules) per feature of an STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' table_st(Xenium_Mouse_Brain_Coronal_7g)
#' @export
setMethod("table_st", signature("STGrid"),
          function(x){
            table(coord(x)$feature)
          }
)


#' @title  Get Ripley's K function slot from a STGrid object.
#' @description
#'   Get Ripley's K function slot from a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g <- compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g, verbose=FALSE)
#' head(ripley_k_function(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setGeneric("ripley_k_function",
           function(object)
             standardGeneric("ripley_k_function"))

#' @title Get Ripley's K function slot from a STGrid object.
#' @description
#'  Get Ripley's K function slot from a STGrid object.
#' @param x The STGrid object
#' @export bin_size
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g <- compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g, verbose=FALSE)
#' head(ripley_k_function(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setMethod("ripley_k_function", "STGrid",
          function(object)
            object@ripley_k_function)


#' @title  Get Moran's index values from corresponding slot from a STGrid object.
#' @description
#'   Get Moran's index values from corresponding slot from a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g <- compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g, verbose=FALSE)
#' head(moran_index(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setGeneric("moran_index",
           function(object)
             standardGeneric("moran_index"))

#' @title  Get Moran's index values from corresponding slot from a STGrid object.
#' @description
#'   Get Moran's index values from corresponding slot from a STGrid object.
#' @param x The STGrid object
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g <- compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g, verbose=FALSE)
#' head(moran_index(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setMethod("moran_index", "STGrid",
          function(object)
            object@moran_index)

# -------------------------------------------------------------------------
##    REDEFINE SHOW() METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------

#' @title
#' The show method of a STGrid object
#' @description
#' The show method of a STGrid object
#' @param object A STGrid object.
#' @keywords internal
#' @examples
#' example_dataset()
#' show(Xenium_Mouse_Brain_Coronal_7g)
#' @import methods
#' @export
setMethod("show", signature("STGrid"),
          function(object) {
            print_this_msg("An object of class STGrid")
            print_this_msg("Number of counts: ", nb_items(object))
            print_this_msg("Number of features: ", nb_feat(object))
            print_this_msg("Bin size: ", bin_size(object))
            print_this_msg("Number of bins (x axis): ", length(bin_x(object)))
            print_this_msg("Number of bins (y axis): ", length(bin_y(object)))
            print_this_msg("x_min: ", object@x_min)
            print_this_msg("x_max: ", object@x_max)
            print_this_msg("y_min: ", object@y_min)
            print_this_msg("y_max: ", object@y_max)
            print_this_msg(">>> Please, use show_st_methods() to show availables methods <<<")
          })

# -------------------------------------------------------------------------
##    REDEFINE SUMMARY() METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------

#' @title
#' The summary() method of a STGrid object
#' @description
#' The summary method of a STGrid object
#' @param object A STGrid object.
#' @examples
#' example_dataset()
#' summary(Xenium_Mouse_Brain_Coronal_7g)
#' @export summary
setMethod("summary", signature("STGrid"),
          function(object) {
            print_this_msg("An object of class STGrid")
            print_this_msg("Method:", object@method)
            print_this_msg("Bin size:", object@bin_size)
            print_this_msg("Features:", c(utils::head(feat_names(object), 4), "..."))
            print_this_msg("x_min:", object@x_min)
            print_this_msg("x_max:", object@x_max)
            print_this_msg("y_min:", object@y_min)
            print_this_msg("y_max:", object@y_max)
            print_this_msg("path:", object@path)
            print_this_msg("bin_x:", utils::head(bin_x(object), 4), "...")
            print_this_msg("bin_y:", utils::head(bin_y(object), 4), "...")
            print_this_msg("Meta data:", paste0(colnames(object@meta), collapse=" , "))
          })

# -------------------------------------------------------------------------
##    Method for function"[". Subsetting an STGrid obj.
# -------------------------------------------------------------------------
#' @title Subsetting operator of a STGrid object
#' Extract
#' @description
#' The subsetting operator of a STGrid object.
#' The i axis correspond to features (i.e. genes) or bins (x axis). The j axis corresponds to bins
#' @param features (i.e. genes) or bins (x axis) to extract or substract.
#' @param ... See ?'['. Not functionnal here.
#' @param drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. Not functionnal here.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' xen["Ano1",]
#' xen[bin_x(xen)[1:20], bin_y(xen)[1:20]]
#' @keywords internal
setMethod("[", signature(x = "STGrid"),
          function (x, i, j, ..., drop = FALSE) {
            x_is_feat <- FALSE
            x_is_bin <- FALSE

            if (!missing(i)) {
              i <- unique(i)

              if (any(i %in% bin_x(x))) {
                x_is_bin <- TRUE
              } else{
                x_is_bin <- FALSE
                if (any(i %in% feat_names(x))) {
                  x_is_feat <- TRUE
                } else{
                  if (is.numeric(i)) {
                    i <- round(i, 0)
                    if (any(i > nb_feat(x))) {
                      print_this_msg("Numering value i is out of range (should be < nb_feat(x)).",
                                msg_type = "STOP")
                    }
                    i <- feat_names(x)[i]
                    x_is_feat <- TRUE
                  } else{
                    print_this_msg("Check i... Should be a set of features or bins.",
                              msg_type = "STOP")
                  }

                }

              }

              if (!x_is_bin && !x_is_feat)
                print_this_msg("Some features or bin_x or indexes where not found in the object.",
                          msg_type = "STOP")
            }

            if (!missing(j)) {
              if (!all(j %in% bin_y(x)))
                print_this_msg("Some bin_y where not found in the object.", msg_type = "STOP")
            }

            n_coord <- x@coord
            n_bin_mat <- x@bin_mat
            n_meta <- x@meta
            n_bin_x <- x@bin_x
            n_bin_y <- x@bin_y
            n_ripley_k_function <- x@ripley_k_function
            n_moran_index <- x@moran_index
            n_control <- x@control

            if (missing(j)) {
              if (missing(i)) {
                return(x)
              } else {
                if (x_is_feat) {
                  n_coord <- n_coord[n_coord$feature %in% i, ]
                  n_bin_mat <-
                    n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", i)]
                  n_ripley_k_function <-
                    n_ripley_k_function[n_ripley_k_function$feature %in% i, ]
                  n_moran_index <-
                    n_moran_index[n_moran_index$feature %in% i, ]

                } else {
                  n_coord <- n_coord[n_coord$bin_x %in% i, ]
                  test <- n_bin_mat$bin_x %in% i
                  n_bin_mat <- n_bin_mat[test, ]
                  n_meta <- n_meta[test, ,drop=FALSE]
                  feat_left <- unique(n_coord$feature)

                  n_bin_mat <-
                      n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", feat_left)]
                }
              }
            } else{
              if (missing(i)) {
                n_coord <- n_coord[n_coord$bin_y %in% j, ]
                test <- n_bin_mat$bin_y %in% j
                n_bin_mat <- n_bin_mat[test, ]
                n_meta <- n_meta[test, ,drop=FALSE]

              } else{
                if (x_is_feat) {
                  n_coord <- n_coord[n_coord$feature %in% i, ]

                  n_bin_mat <-
                      n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", i)]

                } else {

                  n_coord <- n_coord[n_coord$bin_x %in% i & n_coord$bin_y %in% j, ]
                  test <- n_bin_mat$bin_x %in% i & n_bin_mat$bin_y %in% j
                  n_bin_mat <- n_bin_mat[test, ]
                  n_meta <- n_meta[test, ,drop=FALSE]
                  feat_left <- unique(n_coord$feature)

                  n_bin_mat <-
                      n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", feat_left)]

                }

              }
            }

            STGrid_obj <- methods::new("STGrid",
                              path = x@path)

            STGrid_obj@coord <- n_coord
            STGrid_obj@bin_mat <- n_bin_mat
            STGrid_obj@y_max <- x@y_max
            STGrid_obj@y_min <- x@y_min
            STGrid_obj@x_max <- x@x_max
            STGrid_obj@x_min <- x@x_min
            STGrid_obj@method <- x@method
            STGrid_obj@meta <- n_meta
            STGrid_obj@bin_size <- x@bin_size
            STGrid_obj@bin_x <-
              x@bin_x[x@bin_x %in% n_bin_mat$bin_x]
            STGrid_obj@bin_y <-
              x@bin_y[x@bin_y %in% n_bin_mat$bin_y]
            STGrid_obj@control <- n_control
            STGrid_obj@ripley_k_function <- n_ripley_k_function
            STGrid_obj@moran_index <- n_moran_index

            return(STGrid_obj)

})

# -------------------------------------------------------------------------
##    Method for function "[[". Set/Extract/Replace metadata.
# -------------------------------------------------------------------------
#' @title Subsetting operator "[[" of a STGrid object. Extract metadata.
#' @description
#' Subsetting operator "[[" of a STGrid object. Extract metadata.
#' @param i indices specifying metadata to extract.
#' @param exact See ?"[[".
#' @examples
#' example_dataset()
#' utils::head(Xenium_Mouse_Brain_Coronal_7g[[1]])
#' @keywords internal
#' @export
setMethod("[[", signature(x = "STGrid"),
          function (x, i, exact = TRUE) {
            x@meta[[i]]
          })

#' @title Subsetting operator "$" of a STGrid object. Extract metadata.
#' @description
#' Subsetting operator "$" of a STGrid object. Extract metadata.
#' @param i indices specifying metadata to extract.
#' @param name The metadata to extract.
#' @examples
#' example_dataset()
#' head(Xenium_Mouse_Brain_Coronal_7g$count_sum)
#' head(Xenium_Mouse_Brain_Coronal_7g[[1]])
#' @keywords internal
#' @export
setMethod ("$", "STGrid",
           function (x, name) {
             return(x[[name]])
           })


#' @title This functions makes STGrid objects look like list and allow completion after $.
#' @description
#' This functions makes STGrid objects look like list and allow completion after $.
#' @param x The STGtrid object.
#' @param pattern The pattern to be searched.
#' @keywords internal
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' xen$bla <- "foo"
#' head(xen@meta)
#' @importFrom utils .DollarNames
#' @export
.DollarNames.STGrid <- function (x, pattern = '') {
  grep(pattern,
       colnames(x@meta),
       value = TRUE,
       fixed = TRUE)
}

#' @title Replacing operator "$" of a STGrid object. Extract metadata.
#' @description
#' Subsetting operator "$" of a STGrid object. Extract metadata.
#' @param i indices specifying metadata to extract.
#' @param exact See ?"$".
#' @param value The replacement value.
#' @keywords internal
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g$count_sum <- Xenium_Mouse_Brain_Coronal_7g@meta$count_sum + 1
#' @export
setMethod ("$<-", "STGrid",
           function (x, name, value) {
             if (length(value) != nbin_x(x) * nbin_y(x) && length(value) != 1)
               print_this_msg("The size of the vector needs to be nbin_x(x) * nbin_y(x) or 1.",
                         msg_type = "STOP")
             x@meta[[name]] <- value
             x
           })

#' @title Replacing operator "[[" of a STGrid object. Extract metadata.
#' @description
#' Subsetting operator "[[" of a STGrid object. Extract metadata.
#' @param name The name of the metadata.
#' @param value The replacement value.
#' @keywords internal
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g[[1]] <- Xenium_Mouse_Brain_Coronal_7g[[1]]
#' @export
setMethod (f = '[[<-',
           signature = c('x' = 'STGrid', i = 'numeric', value = 'ANY'),
           definition = function(x, i, ..., value) {
             if (length(value) != nbin_x(x) * nbin_y(x) && length(value) != 1)
               print_this_msg("The size of the vector needs to be nbin_x(x) * nbin_y(x) or 1.",
                              msg_type = "STOP")
             x@meta[[i]] <- value
             return(x)
           })

#' @title Replacing operator "[[" of a STGrid object. Extract metadata.
#' @description
#' Subsetting operator "[[" of a STGrid object. Extract metadata.
#' @param name The name of the metadata.
#' @param value The replacement value.
#' @keywords internal
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g[[1]] <- Xenium_Mouse_Brain_Coronal_7g[[1]]
#' @export
setMethod (f = '[[<-',
           signature = c('x' = 'STGrid', i = 'character', value = 'ANY'),
           definition = function(x, i, ..., value) {
             if (length(value) != nbin_x(x) * nbin_y(x) && length(value) != 1)
               print_this_msg("The size of the vector needs to be nbin_x(x) * nbin_y(x) or 1.",
                              msg_type = "STOP")
             x@meta[[i]] <- value
             return(x)
           })
# -------------------------------------------------------------------------
##    Compute Ripley's K function
# -------------------------------------------------------------------------
#' Estimate Ripley's reduced second moment function for each feature in a spatial grid.
#'
#' This method calculates Ripley's reduced second moment function, K(r), for each feature in a spatial grid.
#'
#' @param object An object of class "STGrid".
#' @param rmax Maximum desired value of the argument r.
#' @param nlarge Efficiency threshold. If the number of points exceeds nlarge, then only the border correction will be computed (by default), using a fast algorithm.
#' @param var.approx Logical. If TRUE, the approximate variance of K(r) under CSR will also be computed.
#' @param ratio Logical. If TRUE, the numerator and denominator of each edge-corrected estimate will also be saved, for use in analyzing replicated point patterns.
#' @param sampling_rate By default a subset of molecules is sampled. Set sampling_rate to 1 to keep all molecules.
#' @param seed The seed for sampling_rate.
#' @param method Whether to to compute Ripley's K-function or an estimate of the L-function.
#' @param verbose Whether to print porgress bar.
#' @return An updated object of class "STGrid" with the Ripley's K function estimates stored in the slot 'ripley_k_function'.
#'
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g <- compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g, verbose=FALSE)
#' @keywords internal
#' @export compute_k_ripley
setGeneric("compute_k_ripley",
           function(object,
                    rmax = 80,
                    nlarge = 2e6,
                    var.approx = FALSE,
                    ratio = FALSE,
                    sampling_rate=0.25,
                    seed=123,
                    method=c("K", "L"),
                    verbose = TRUE)
             standardGeneric("compute_k_ripley"))

#' Estimate Ripley's reduced second moment function for each feature in a spatial grid.
#'
#' This method calculates Ripley's reduced second moment function, K(r), for each feature in a spatial grid.
#'
#' @param object An object of class "STGrid".
#' @param rmax Maximum desired value of the argument r.
#' @param nlarge Efficiency threshold. If the number of points exceeds nlarge, then only the border correction will be computed (by default), using a fast algorithm.
#' @param var.approx Logical. If TRUE, the approximate variance of K(r) under CSR will also be computed.
#' @param ratio Logical. If TRUE, the numerator and denominator of each edge-corrected estimate will also be saved, for use in analyzing replicated point patterns.
#' @param sampling_rate By default a subset of molecules is sampled. Set sampling_rate to 1 to keep all molecules.
#' @param seed The seed for sampling_rate.
#' @param method Whether to to compute Ripley's K-function or an estimate of the L-function.
#' @param verbose Whether to print porgress bar.
#' @return An updated object of class "STGrid" with the Ripley's K function estimates stored in the slot 'ripley_k_function'.
#'
#' @examples
#' example_dataset()
#' Xenium_Mouse_Brain_Coronal_7g <- compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g, verbose=FALSE)
#'
#' @importFrom spatstat.geom ppp owin
#' @importFrom spatstat.explore Kest
#' @export compute_k_ripley
setMethod("compute_k_ripley", signature("STGrid"),
          function(object,
                   rmax = 80,
                   nlarge = 2e6,
                   var.approx = FALSE,
                   ratio = FALSE,
                   sampling_rate=0.25,
                   seed=123,
                   method=c("K", "L"),
                   verbose = TRUE) {

            method <- match.arg(method)

            molecules <- coord(object)[, c("x", "y", "feature")]


            data_out <- data.frame(
              r = NA,
              theo = NA,
              border = NA,
              feature = NA
            )

            x_min <- object@x_min
            x_max <- object@x_max
            y_min <- object@y_min
            y_max <- object@y_max

            if(sampling_rate != 1){
              print_this_msg("Subsetting...", msg_type = "DEBUG")
              print_this_msg("Nb items before:", nrow(molecules), msg_type = "DEBUG")
              set.seed(123)
              subset <- sample(1:nrow(molecules),
                               round(sampling_rate*nrow(molecules), 0),
                               replace = FALSE)
              molecules <- molecules[subset, ]

              print_this_msg("Nb items after:", nrow(molecules), msg_type = "DEBUG")

            }

            rownames(molecules) <- paste0(1:nrow(molecules),
                                          '|',
                                          molecules$feature)


            print_this_msg(">>> Estimating Ripley's reduced second moment function for all features.")

            n_iter <- length(unique(molecules$feature))

            if (verbose) {
              pb <- utils::txtProgressBar(
                min = 0,
                max = n_iter,
                style = 3,
                width = 50,
                char = "="
              )

              n <- 0
            }

            for (g in unique(molecules$feature)) {
              if (verbose) {
                n <- n + 1
                utils::setTxtProgressBar(pb, n)
              }


              list_neighbors <- data.frame(Var1 = NA,
                                           Var2 = NA,
                                           value = NA)

              is_g <- molecules[molecules$feature == g, ]

              X <- spatstat.geom::ppp(is_g$x,
                                      is_g$y,
                                      window = spatstat.geom::owin(c(min(is_g$x),
                                                                     max(is_g$x)),
                                                                   c(min(is_g$y),
                                                                     max(is_g$y))))

              if(method == "K"){
                u <- spatstat.explore::Kest(
                  X,
                  correction=c("border"),
                  rmax = rmax,
                  nlarge = nlarge,
                  var.approx = var.approx,
                  ratio = ratio
                )
              }else{
                u <- spatstat.explore::Lest(
                  X,
                  rmax = rmax,
                  nlarge = nlarge,
                  var.approx = var.approx,
                  ratio = ratio
                )
              }


              u <- as.data.frame(as.matrix(u))
              u$feature <- g

              data_out <- rbind(data_out,
                                u[, c("r", "theo",
                                      "border", "feature")])
            }

            if (verbose)
              close(pb)

            data_out <- stats::na.omit(data_out)

            object@ripley_k_function <- data_out

            return(object)

})


# -------------------------------------------------------------------------
##    Constructor for SpatialTranscriptomicAnalysis class (STGrid)
# -------------------------------------------------------------------------
#' @title Create a Spatial Transcriptomic Grid class (STGrid)
#' @description
#' The load_spatial() function is the entry point of the STarlight package.
#' It will load the molecule coordinates and create a 2D binned grid with default size 25µm.
#'
#' @param path Either a file (if method is set to "coordinates") or a directory (if method
#' is set to "merscope"). If method is set to "coordinates" the file should contain at least 3 columns
#' ("x", "y", "gene"), ("x", "y", "feature") or ("x", "y", "cell").
#' @param method The type of technology/file/directory structure. See details.
#' @param bin_size Numeric value representing the bin size (default to 25). Unit is microns for Xenium and Merscope and pixel for CosMx.
#' @param control A regular expression to identify controls. As the function computes the sum of
#' counts, this will allow to delete these blanks/controls for computation.
#' @param sep The separator when method is set to "coordinates" (default "\\t").
#' @param threads The number of threads (see data.table::fread).
#' @param verbose Whether to display the progress bar.
#' @param constrain Whether to put constrains on a column when method is set to 'coordinates'. E.g: 'list("global_z==0", "fov==0")'. Works currently with 'coordinates'.
#' @param mapping if method='coordinates' is used and non conventional column names are used in the input file,
#'
#' @return An object of class STGrid.
#' @importFrom Seurat ReadVizgen
#' @importFrom Seurat ReadXenium
#' @importFrom Seurat ReadNanostring
#' @importFrom R.utils isZero
#' @importFrom data.table fread
#' @details
#' If method is set to 'coordinates' a flat file with ("x", "y", "gene"), ("x", "y", "feature") or ("x", "y", "cell")
#' expected for path. If method is set to 'merscope_csv' the transcript csv file exported by Merscope should be provided
#' as path (contains 'global_x', 'global_y' and 'gene' column). If method is set to 'xenium' or 'merscope'
#' path is passed to the 'data.dir' argument of Seurat::ReadXenium and Seurat::ReadVizgen respectively (see corresponding docs).
#'
#' @examples
#'   # Coordinate example
#'   fp <- file.path(system.file("extdata", package = "STarlight"), "tyni.txt")
#'   st <- load_spatial(fp, method = "coordinates")
#'   # merscope csv file
#'   fp <- file.path(system.file("extdata", package = "STarlight"), "merscope_122_r1_sub.csv.gz")
#'   st <- load_spatial(fp, method = "merscope_csv", sep=",")
#'   # Xenium
#'   fp <- file.path(system.file("extdata", package = "STarlight"), "xenium_mouse_brain_tx_tiny.csv")
#'   st <- load_spatial(fp,
#'                      method = "coordinates", sep=",",
#'                      mapping = c("x"="x_location", "y"="y_location", "feature"="feature_name"))
#'   # Cosmix
#'   fp <- file.path(system.file("extdata", package = "STarlight"), "Lung5_Rep1_tx_file_tiny.csv.gz")
#'   st <- load_spatial(fp, method = "coordinates", sep=",", mapping=c("x"="x_global_px", "y"="y_global_px", feature="target"))
#'
#' @export load_spatial
load_spatial <- function(path = "",
                         method = c("coordinates",
                                    "merscope",
                                    "merscope_csv",
                                    "xenium",
                                    "cosmx"),
                         bin_size = 25,
                         control = NULL,
                         sep="\t",
                         threads=1,
                         constrain=NULL,
                         verbose = TRUE,
                         mapping=NULL) {
  method <- match.arg(method)

  print_this_msg("Technology is '", method, "'.")
  print_this_msg("Loading data from file:", path)

  # Force the use of R.utils
  out <- R.utils::isZero(0)

  if (method == "merscope") {
    spat_input <-
      Seurat::ReadVizgen(data.dir = path, type = "centroids")
    spat_input <- spat_input$microns[, c("x", "y", "gene")]

    if (is.null(control))
      control <- "^Blank\\-[0-9]+"

  }  else if (method == "xenium") {
    spat_input <-
      Seurat::ReadXenium(data.dir = path, type = "centroids")
    spat_input <- spat_input$microns[, c("x", "y", "gene")]

    if (is.null(control))
      control <-  "(NegControl)|(^BLANK)"

  }else if (method == "cosmx") {
    spat_input <-
      Seurat::ReadNanostring(data.dir = path, type = "centroids")
    spat_input <- spat_input$pixels[, c("x", "y", "gene")]

    if (is.null(control))
      control <-  "^NegPrb"

  }else if(method == "merscope_csv"){
    check_this_file(path, mode = "read")
    spat_input <- as.data.frame(data.table::fread(path,
                                    sep = sep,
                                    head = TRUE,
                                    nThread = threads))

    col_needed <- c("global_x", "global_y", "gene")
    spat_input <- spat_input[, col_needed]

    if (is.null(control))
      control <- "^Blank\\-[0-9]+"

  }else if (method == "coordinates") {
    check_this_file(path, mode = "read")
    spat_input <- as.data.frame(data.table::fread(path,
                                    sep = sep,
                                    head = TRUE, nThread = threads))

    if(!is.null(constrain)){
      if(!is.list(constrain))
        print_this_msg("Argument 'contrain' should be a list", msg_type = "STOP")

      print_this_msg('Evaluating provided contrain.')

      for(this_contrain in constrain){
        print_this_msg("Checking column ('", this_contrain, "').", msg_type = "DEBUG")
        if(!gsub("^(\\w+)\\s?[!=><].*","\\1", this_contrain) %in% colnames(spat_input))
          print_this_msg("Column not found (see 'constrain').", msg_type = "STOP")

        print_this_msg("Evaluating test ('", this_contrain, "').", msg_type = "DEBUG")
        contrain_col <- paste0("spat_input$", this_contrain)
        my_eval <- eval(parse(text=contrain_col))

        print_this_msg("Table size before applying constrain:", nrow(spat_input))
        spat_input <- spat_input[my_eval, ]
        print_this_msg("Table after applying constrain:", nrow(spat_input))
        if(nrow(spat_input) ==0)
          print_this_msg("No line left after 'contrain'...",  msg_type = "STOP")
      }

    }

     if(!is.null(mapping)){
       if(!all(c("x", "y", "feature") %in% names(mapping) )){
         print_this_msg("Please provide a mapping for both x, y and feature.", msg_type = "STOP")
       }

       if(!all(mapping %in% colnames(spat_input) )){
         print_this_msg("Some columns are not part of the data.frame", msg_type = "STOP")
       }

       col_needed <- mapping[c("x", "y", "feature")]

     }else{
       col_needed <- c("x", "y")

       if("cell" %in% colnames(spat_input)){
         col_needed <- c(col_needed, "cell")
       }else if("feature" %in% colnames(spat_input)){
         col_needed <- c(col_needed, "feature")
       }else if("gene" %in% colnames(spat_input)){
         col_needed <- c(col_needed, "gene")
       }else{
         print_this_msg("Could not find a cell/feature/gene columns.", msg_type = "STOP")
       }

       if (any(!col_needed %in% colnames(spat_input))) {
         print_this_msg("Please check the column name of the input file.", msg_type = "STOP")
       }

     }

    spat_input <- spat_input[, col_needed]

    if (is.null(control))
      control <- "^Blank\\-[0-9]+"

  }

  colnames(spat_input) <- c("x","y", "feature")

  bin_matrix <- bin_this_matrix(coord = spat_input,
                                bin_size = bin_size,
                                verbose = verbose)

  count_sum <- sum_of_counts(bin_matrix$spatial_matrix,
                              control = control)

  # create a STGrid object                             ---------

  STGrid_obj <- methods::new("STGrid",
                    path = path)

  STGrid_obj@coord <- bin_matrix$coord
  STGrid_obj@bin_mat <- bin_matrix$spatial_matrix
  STGrid_obj@y_max <- bin_matrix$y_max
  STGrid_obj@y_min <- bin_matrix$y_min
  STGrid_obj@x_max <- bin_matrix$x_max
  STGrid_obj@x_min <- bin_matrix$x_min
  STGrid_obj@path <- path
  STGrid_obj@method <- method
  STGrid_obj@meta <-
    data.frame(row.names = rownames(STGrid_obj@bin_mat),
               count_sum = count_sum)
  STGrid_obj@bin_size <- bin_size
  STGrid_obj@bin_x <- bin_matrix$possible_levels_x
  STGrid_obj@bin_y <- bin_matrix$possible_levels_y
  STGrid_obj@ripley_k_function <- data.frame()
  STGrid_obj@control <- control
  return(STGrid_obj)
}

# -------------------------------------------------------------------------
#      Bin a matrix
# -------------------------------------------------------------------------
#' Compute the Sum of Counts
#'
#' This function calculates the sum of counts for each row in a spatial matrix.
#'
#' @param spatial_matrix The spatial matrix containing molecule counts, typically obtained from a spatial transcriptomics experiment.
#' @param control A regular expression pattern indicating the controls in the dataset. Controls are internal features used for for monitoring experiments.
#'
#' @return A numeric vector containing the sum of counts for each row in the spatial matrix.
#' @export
sum_of_counts <- function(spatial_matrix = NULL,
                          control = NULL) {

  if (is.null(spatial_matrix)) {
    print_this_msg("Please provide a spatial matrix.", msg_type = "STOP")
  }

  print_this_msg("Regexp for controls:", control, msg_type = "DEBUG")

  print_this_msg("Computing sum of counts.")

  pos_bin_xy <- grep("^bin_[xy]$",
                     colnames(spatial_matrix))

  tmp <- spatial_matrix[, -pos_bin_xy, drop=FALSE]

  if(!is.null(control)) {
    pos_ctrl <- grep(control, colnames(tmp), perl = TRUE)
  } else{
    pos_ctrl <- integer(0)
  }

  if (length(pos_ctrl) > 0) {
    print_this_msg("Found the following controls:",
              paste0(utils::head(colnames(tmp)[pos_ctrl], 3),
                     collapse = ", "),
              "...")
    count_sum <- rowSums(tmp[, -pos_ctrl, drop=FALSE])

  } else{

    count_sum <- rowSums(tmp)
  }

  return(count_sum)

}

# -------------------------------------------------------------------------
#      Bin a matrix
# -------------------------------------------------------------------------
#' Create a 2D binned grid (rasterization) from molecule coordinates
#'
#' This function takes molecule coordinates as input (with columns named 'x', 'y', 'feature') and creates a 2D
#' binned grid by dividing the x and y axes into bins. It returns the counts of the number of molecules in
#' each bin along with additional information about the binning process.
#'
#' @param coord A data frame containing molecule coordinates with columns 'x', 'y', and 'feature'. Defaults to NULL.
#' @param bin_size An integer specifying the size of each bin. Defaults to 25.
#' @param verbose Whether to display the progress bar.
#' @return A list containing the binned spatial matrix, updated molecule coordinates, and information about the binning process.
#' @examples
#' example_dataset()
#' res <- bin_this_matrix(Xenium_Mouse_Brain_Coronal_7g@coord)
#'
#' @keywords internal
#' @importFrom data.table rbindlist
#' @export bin_this_matrix
bin_this_matrix <- function(coord = NULL,
                            bin_size = 25,
                            verbose = TRUE) {
  print_this_msg("Binning a matrix...")

  x_min <- min(coord$x)
  x_max <- max(coord$x)
  y_min <- min(coord$y)
  y_max <- max(coord$y)

  if (x_max - x_min < bin_size) {
    print_this_msg(
      "Can't create bin with x_max =",
      x_max,
      ", x_min =",
      x_min,
      "and bin_size =",
      bin_size,
      msg_type = "STOP"
    )
  }

  if (y_max - y_min < bin_size) {
    print_this_msg(
      "Can't create bin with x_max =",
      y_max,
      ", x_min =",
      y_min,
      "and bin_size =",
      bin_size,
      msg_type = "STOP"
    )
  }

  x_lim <- seq(
    from = x_min,
    to = ifelse(x_max %% bin_size != 0,
                x_max + bin_size, x_max),
    by = bin_size
  )


  y_lim <- seq(
    from = y_min,
    to = ifelse(y_max %% bin_size != 0,
                y_max + bin_size, y_max),
    by = bin_size
  )

  possible_levels_x <- levels(cut(
    0,
    breaks = x_lim,
    include.lowest = TRUE,
    right = FALSE
  ))

  possible_levels_y <- levels(cut(
    0,
    breaks = y_lim,
    include.lowest = TRUE,
    right = FALSE
  ))


  spatial_matrix <- data.frame(bin_x = sort(rep(
    possible_levels_x,
    length(possible_levels_y)
  )),
  bin_y = rep(possible_levels_y,
              length(possible_levels_x)))

  rownames(spatial_matrix) <-
    paste(as.character(spatial_matrix$bin_x),
          as.character(spatial_matrix$bin_y),
          sep = "~")


  coord$bin_x <- NA
  coord$bin_y <- NA

  print_this_msg("Looping over genes...")

  if (verbose) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = length(unique(coord$feature)),
      style = 3,
      width = 50,
      char = "="
    )

    n_loop <- 0
  }

  coord_as_list <- split(coord, coord$feature)

  for (goi in unique(coord$feature)) {

    x_molec <- cut(
      coord_as_list[[goi]]$x,
      breaks = x_lim,
      include.lowest = TRUE,
      right = FALSE
    )

    coord_as_list[[goi]]$bin_x <- as.character(x_molec)

    y_molec <- cut(
      coord_as_list[[goi]]$y,
      breaks = y_lim,
      include.lowest = TRUE,
      right = FALSE
    )

    coord_as_list[[goi]]$bin_y <- as.character(y_molec)

    nb_items <- table(paste(x_molec, y_molec, sep = "~"))
    spatial_matrix[, goi] <- 0
    spatial_matrix[names(nb_items), goi] <- nb_items

    if (verbose) {
      n_loop <- n_loop + 1
      utils::setTxtProgressBar(pb, n_loop)
    }
  }

  cat("\n")

  if (verbose)
    close(pb)

  print_this_msg("Ordering")
  spatial_matrix <-
    spatial_matrix[, order(colnames(spatial_matrix))]

  print_this_msg("Merging...")

  coord <- data.table::rbindlist(coord_as_list)

  print_this_msg("Returning")
  return(
    list(
      spatial_matrix = spatial_matrix,
      coord = coord,
      x_max = x_max,
      x_min = x_min,
      y_min = y_min,
      y_max = y_max,
      possible_levels_x = possible_levels_x,
      possible_levels_y = possible_levels_y
    )
  )
}

# -------------------------------------------------------------------------
#      Re-Bin a matrix
# -------------------------------------------------------------------------
#' @title Re-bin a STGrid object.
#' @description Re-bin a STGrid object.
#' @param object The STGrid object
#' @param bin_size The size of the bin.
#' @param verbose Whether to be display progress bar.
#' @examples
#' example_dataset()
#' res <- re_bin(Xenium_Mouse_Brain_Coronal_7g, bin_size=10)
#'
#' @keywords internal
#' @export re_bin
setGeneric("re_bin",
           function(object,
                    bin_size,
                    verbose = TRUE)
             standardGeneric("re_bin"))


#' @title Re-bin a STGrid object.
#' @description Re-bin a STGrid object.
#' @param object The STGrid object.
#' @param bin_size The size of the bin.
#' @param verbose Whether to be display progress bar.
#' @examples
#' example_dataset()
#' res <- re_bin(Xenium_Mouse_Brain_Coronal_7g, bin_size=10)
#' @export re_bin
setMethod("re_bin", signature(object = "STGrid"),
          function(object,
                   bin_size,
                   verbose = TRUE) {
            if (object@bin_size == bin_size) {
              print_this_msg("The bin_size is unchanged.")
              return(object)
            }

            bin_matrix <- bin_this_matrix(coord = object@coord,
                                          bin_size = bin_size,
                                          verbose = verbose)

            print_this_msg("Re-computing sum of counts .")
            count_sum <- sum_of_counts(bin_matrix$spatial_matrix,
                                        control=object@control)

            # create a STGrid object                             ---------
            print_this_msg("Creating an STGrid object")
            print_this_msg("Note that meta data will be lost (except 'count_sum').")

            STGrid_obj <- methods::new("STGrid",
                              path = object@path)

            STGrid_obj@coord <- bin_matrix$coord
            STGrid_obj@bin_mat <- bin_matrix$spatial_matrix
            STGrid_obj@y_max <- bin_matrix$y_max
            STGrid_obj@y_min <- bin_matrix$y_min
            STGrid_obj@x_max <- bin_matrix$x_max
            STGrid_obj@x_min <- bin_matrix$x_min
            STGrid_obj@path <- object@path
            STGrid_obj@method <- object@method
            STGrid_obj@meta <- object@meta

            STGrid_obj@meta <- data.frame(row.names = rownames(STGrid_obj@bin_mat),
                                          count_sum = count_sum)
            STGrid_obj@bin_size <- bin_size
            STGrid_obj@bin_x <- bin_matrix$possible_levels_x
            STGrid_obj@bin_y <- bin_matrix$possible_levels_y
            STGrid_obj@ripley_k_function <- object@ripley_k_function
            STGrid_obj@control <- object@control

            return(STGrid_obj)


          })


# -------------------------------------------------------------------------
#      Return the x/y coordinates of features from a STGrid object
# -------------------------------------------------------------------------
#' @title The x/y coordinates of features from a STGrid object.
#' @description Return the x/y coordinates of features from a STGrid object
#' @param object The STGrid object
#' @param feat_list The list of features.
#' @param as.factor Whether the 'feature' column should be returned as an ordered factor.
#' @examples
#' example_dataset()
#' head(get_coord(Xenium_Mouse_Brain_Coronal_7g))
#' @export
#' @keywords internal
setGeneric("get_coord",
           function(object,
                    feat_list = character(),
                    as.factor = TRUE)
             standardGeneric("get_coord"))

#' @title The x/y coordinates of features from a STGrid object.
#' @description Return the x/y coordinates of features from a STGrid object
#' @param object The STGrid object
#' @param feat_list The list of features.
#' @param as.factor Whether the 'feature' column should be returned as an ordered factor
#' @examples
#' example_dataset()
#' head(get_coord(Xenium_Mouse_Brain_Coronal_7g))
#' @export
setMethod("get_coord", "STGrid",
          function(object,
                   feat_list = character(),
                   as.factor = TRUE) {
            if (!check_features_exist(object, feat_list)) {
              print_this_msg("Some features were not found...",
                        msg_type = "STOP")
            }

            if (length(feat_list) == 0) {
              feat_list <- feat_names(object)
            }

            coord <- object@coord
            coord <- coord[coord$feature %in% feat_list, ]

            if (as.factor)
              coord$feature <- factor(as.character(coord$feature),
                                   levels = feat_list,
                                   ordered = TRUE)

            return(coord)

          })

# -------------------------------------------------------------------------
#      Create a tree from an STGrid object
# -------------------------------------------------------------------------
#' @title Create a tree from an STGrid object
#' @description Create a tree from an STGrid object
#' @param object The STGrid object
#' @param method The agglomeration method to be used. See stats::hclust().
#' @param layout The layout to be used See ggtree::ggtree().
#' @param dist_method The method for distance computation. See stats::cor().
#' @param branch_length A variable for scaling branch, if 'none' draw cladogram. See ggtree::ggtree().
#' @param class_nb An integer indicating the desired number of groups.
#' @param class_name A set of names for the classes.
#' @param class_label Whether to add class labels.
#' @param offset Offset of labels from the clade
#' @param offset.text Offset of text labels from the bar.
#' @param lab_fontsize Label font size.
#' @param lab_barsize Bar label size.
#' @param geom_label Label fomat.
#' @param no_legend Whether to discard legend.
#' @param size The size of the labels.
#' @param colors A set of colors for the classes.
#' @param return_list Logical. Whether to return a list with the tree,
#' associated clusters, distance and correlation matrices.
#' @examples
#' example_dataset()
#' p <- hc_tree(Xenium_Mouse_Brain_Coronal_7g)
#' print(p)
#' @keywords internal
setGeneric("hc_tree",
           function(object,
                    method = c("complete", "ward.D", "ward.D2", "single",
                               "average", "mcquitty", "median", "centroid"),
                    layout = c("circular", "rectangular",
                               "slanted", "fan", "unrooted",
                               "time-scale", "two-dimensional"),
                    dist_method = "pearson",
                    branch_length = "none",
                    class_nb = 1,
                    class_name = NULL,
                    size = 2.25,
                    colors=NULL,
                    class_label=TRUE,
                    offset=5,
                    offset.text=1.5,
                    lab_fontsize=3.88,
                    lab_barsize=2,
                    geom_label=c('text', 'label', 'shadowtext'),
                    no_legend=FALSE,
                    return_list=FALSE)
             standardGeneric("hc_tree"))

#' @title Create a tree from an STGrid object
#' @description Create a tree from an STGrid object
#' @param object The STGrid object
#' @param method The agglomeration method to be used. See stats::hclust().
#' @param layout The layout to be used See ggtree::ggtree().
#' @param dist_method The method for distance computation. See stats::cor().
#' @param branch_length A variable for scaling branch, if 'none' draw cladogram. See ggtree::ggtree().
#' @param class_nb An integer indicating the desired number of groups.
#' @param class_name A set of names for the classes.
#' @param size The size of the labels.
#' @param colors A set of colors for the classes.
#' @param class_label Whether to add class labels.
#' @param offset Offset of labels from the clade
#' @param offset.text Offset of text labels from the bar.
#' @param lab_fontsize Label font size.
#' @param lab_barsize Bar label size.
#' @param geom_label Label fomat.
#' @param no_legend Whether to discard legend.
#' @param return_list Logical. Whether to return a list with the tree,
#' associated clusters, distance and correlation matrices.
#' @importFrom ggtree ggtree geom_hilight geom_tippoint geom_tiplab MRCA geom_cladelab
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggplot2 aes scale_color_viridis_d scale_fill_manual
#' @importFrom stringr str_pad
#' @importFrom tidytree MRCA
#' @importFrom stringr str_pad
#' @examples
#' example_dataset()
#' p <- hc_tree(Xenium_Mouse_Brain_Coronal_7g, class_nb = 4, class_name = letters[1:4])
#' print(p)
#' # Get classes:
#' p$tree_classes
#' @export
setMethod("hc_tree", "STGrid",

          function(object,
                   method = c("complete", "ward.D", "ward.D2", "single",
                              "average", "mcquitty", "median", "centroid"),
                   layout = c("circular", "rectangular",
                              "slanted", "fan", "unrooted",
                              "time-scale", "two-dimensional"),
                   dist_method = "pearson",
                   branch_length = "none",
                   class_nb = 1,
                   class_name = NULL,
                   size = 2.25,
                   colors=NULL,
                   class_label=TRUE,
                   offset=5,
                   offset.text=1.5,
                   lab_fontsize=3.88,
                   lab_barsize=2,
                   geom_label=c('text', 'label', 'shadowtext'),
                   no_legend=FALSE,
                   return_list=FALSE
                   ) {

            method <- match.arg(method)
            layout <- match.arg(layout)
            geom_label <- match.arg(geom_label)

            print_this_msg("Using", method ,"method", msg_type="INFO")
            print_this_msg("Using", layout, "layout", msg_type="INFO")

            if (class_nb > 0) {
              if(is.null(class_name)){

                class_name <- paste0("M", stringr::str_pad(1:class_nb,
                                                           2,
                                                           pad = "0"))


              }else{
                if (length(class_name) != class_nb) {
                  print_this_msg("Please set the right number of class names",
                                 "see 'class_name' and 'class_nb'. ",
                                 msg_type = "STOP")
                }
              }
            } else{
              print_this_msg("The class_nb argument should be an integer > than 0...",
                        msg_type = "STOP")
            }

            print_this_msg("Retrieving binned matrix.", msg_type="DEBUG")

            bin_mat <-
              bin_mat(object, del_bin = TRUE)

            print_this_msg("Computing hierarchical clustering.", msg_type="DEBUG")
            cor_mat <- stats::cor(bin_mat,
                                  method = dist_method)
            dist_mat <- (1 - cor_mat) / 2
            hc_clust <- stats::hclust(stats::as.dist(dist_mat),
                               method = method)

            print_this_msg("Preparing colors.", msg_type="DEBUG")

            if(is.null(colors) & class_nb > 0){
              gg_color_hue <- function(n) {
                hues = seq(15, 375, length = n + 1)
                grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
              }
              colors <- gg_color_hue(class_nb)

            }
            print_this_msg("Building ggplot diagram.", msg_type="DEBUG")

            suppressWarnings(
            p <- ggtree::ggtree(hc_clust,
                                layout = layout,
                                branch.length = branch_length)
            )

            print_this_msg("Cutting tree.", msg_type="DEBUG")

            tree_classes <- stats::cutree(hc_clust, k = class_nb)
            groups <- split(names(tree_classes), tree_classes)

            print_this_msg("Retrieving classes.", msg_type="DEBUG")

            suppressMessages({
            clades <-
              sapply(groups, function(n)
                tidytree::MRCA(p, n))
            })

            print_this_msg("Retrieving annotations", msg_type="DEBUG")
            annotation <- data.frame(id = clades[1:class_nb],
                                     Class = class_name)


            #p$data$cell_type <- "ND"
            #pmath <- match(p$data$label, names(ginfo))
            #p$data$cell_type[!is.na(pmath)] <- names(ginfo)[pmath[!is.na(pmath)]]

            id <- Class <- label <- NULL
            p <- p + ggtree::geom_hilight(
              data = annotation,
              extend = 0.2,
              mapping = ggplot2::aes(node = id,
                            fill = Class),
              alpha = 0.3
            ) +
              ggtree::geom_tiplab(
                ggplot2::aes(label = label),
                offset = 1,
                size = size,
                color = 'black'
              )

            if(class_nb != 0){
              p <- p + ggplot2::scale_fill_manual(values=colors)


              # ggtree::geom_tippoint(
              # ggplot2::aes(color = cell_type),
              # size = 1,
              # inherit.aes = TRUE,
              #show.legend = FALSE
              #) +
              nm <- names(tree_classes)
              tree_classes <- paste0("MOD", stringr::str_pad(tree_classes,
                                                             nchar(as.character(class_nb)),
                                                             pad = "0"))
              names(tree_classes) <- nm

              p[["tree_classes"]] <- split(names(tree_classes), tree_classes)

            }else{
              p[["tree_classes"]] <- list()
            }

            if(class_label){
              p <- p + ggtree::geom_cladelab(data=annotation,
                                     mapping = ggplot2::aes(node=id,
                                                   label=Class,
                                                   color=Class),
                                     geom=geom_label,
                                     offset=offset,
                                     offset.text=offset.text,
                                     angle='auto',
                                     barsize=lab_barsize,
                                     fontsize=lab_fontsize)
              p <- p + ggplot2::scale_color_manual(values=colors)
            }

            if(no_legend)
              p <- p + ggplot2::theme(legend.position="none")

            if(return_list){
              return(list(hc_clust=hc_clust,
                          tree_classes=split(names(tree_classes), tree_classes),
                          annotation=annotation,
                          p=p,
                          cor_mat=cor_mat,
                          dist_mat=dist_mat))
            }else{
              return(p + st_gg_theming())
            }



          })


# -------------------------------------------------------------------------
#      Get the min / max value from ripley's k function
# -------------------------------------------------------------------------
#' Returns features ordered by Ripley's K-function
#'
#' This method returns features ordered by Ripley's K-function results stored in
#' an \code{STGrid}.
#'
#' @param object An \code{STGrid} object.
#'
#' @importFrom dplyr group_by filter arrange
#' @keywords internal
#' @export
setGeneric("order_feat_by_ripley",
           function(object)
             standardGeneric("order_feat_by_ripley"))

#' Returns features ordered by Ripley's K-function
#'
#' This method returns features ordered by Ripley's K-function results stored in
#' an \code{STGrid}.
#'
#' @param object An \code{STGrid} object.
#'
#' @importFrom dplyr group_by filter arrange
#' @import magrittr
#' @export
setMethod("order_feat_by_ripley", "STGrid",
          function(object){
            feature <- border <- r <- NULL
            ripk <- ripley_k_function(object)
            voi <- ripk %>%
              dplyr::group_by(feature) %>%
              dplyr::filter(border == max(border)) %>%
              dplyr::filter(r == max(r)) %>%
              dplyr::arrange(dplyr::desc(border))
            return(voi$feature)
})

# -------------------------------------------------------------------------
#      Check consistency of a list of st_grid
# -------------------------------------------------------------------------
#' Check Validity of STGrid List and Feature List
#'
#' Validates a list of STGrid objects and an optional list of features.
#' Ensures that the list contains at least one STGrid object, all objects in the
#' list are of type STGrid, and all specified features exist in the STGrid objects.
#'
#' @param st_list A list of objects to be validated, where each object should be
#'        of class 'STGrid'.
#' @param feat_list An optional vector of feature names to be checked across
#'        all STGrid objects in 'st_list'. If 'NULL', only the STGrid object
#'        validation is performed.
#'
#' @return Invisible NULL. The function is called for its side effect of
#'         stopping execution with an error message if any validation fails.
#'
#' @examples
#' example_dataset()
#' check_st_list(list(Xenium_Mouse_Brain_Coronal_7g, Xenium_Mouse_Brain_Coronal_7g))
#' @keywords internal
#' @export
#'
check_st_list <- function(st_list,
                          feat_list=NULL){
  if (length(st_list) < 1) {
    print_this_msg("Need at least one experiment !!", msg_type = "STOP")
  }

  if (any(unlist(lapply(lapply(st_list, class), "[", 1)) != "STGrid")) {
    print_this_msg("Object should be of type STGrid", msg_type = "STOP")
  }

  for (i in 1:length(st_list)) {
    if (!check_features_exist(st_list[[i]], feat_list)) {
      print_this_msg("The feature was not found in the object.", msg_type = "STOP")
    }
  }
}

# -------------------------------------------------------------------------
#      Create an STGrid from a data.frame
# -------------------------------------------------------------------------
#' @title Create an STGrid from a data.frame
#' @description
#' Create an STGrid from a data.frame
#'
#' @param this_df A data.frame
#' @param mapping Please provide a named vector indicating which column map to "x", "y", and "feature" in the STGrid object.
#' @param bin_size The size of the bin.
#' @return An STGrid object.
#' @examples
#' d <- data.frame(foo=runif(1000, 0, 1000), bar=runif(1000, 0, 1000), bla=sample(letters, 1000, replace=TRUE))
#' st <- stgrid_from_data_frame(d, mapping=c("x"="foo", "y"="bar", "feature"="bla"))
#' st
#' @importFrom utils write.table
#' @export
stgrid_from_data_frame <- function(this_df=NULL,
                                   mapping=NULL,
                                   bin_size=25){

  if(is.null(this_df) | !inherits(this_df, "data.frame")){
    print_this_msg("Please provide a data.frame object.", msg_type = "STOP")
  }

  if(is.null(mapping)){
    print_this_msg("Please provide a column mapping.")
  }

  if(!all(c("x", "y", "feature") %in% names(mapping) )){
    print_this_msg("Please provide a mapping for x, y and feature.", msg_type = "STOP")
  }

  if(!all(mapping %in% colnames(this_df) )){
    print_this_msg("Some columns are not part of the data.frame", msg_type = "STOP")
  }

  this_df <- this_df[, mapping]
  colnames(this_df) <- names(mapping)

  tmp_file <- tempfile()
  utils::write.table(file=tmp_file, x = this_df, sep="\t", quote = FALSE, col.names = NA)

  STGrid_obj <- load_spatial(path=tmp_file, method = "coordinates", bin_size = bin_size)

  unlink(tmp_file)

  return(STGrid_obj)

}


# -------------------------------------------------------------------------
#      Compute module scores
# -------------------------------------------------------------------------
#' @title From a list of feature (e.g gene) modules, compute a score (e.g. mean value) accross bins.
#'
#' @description
#' From a list of feature (e.g gene) modules, compute a score (e.g. mean value) accross bins.
#'
#' @param object An \code{STGrid} object.
#' @param modules A list of modules.
#' @param rename Force the renaming of the modules to module_1, module_2...
#' @param fun The function to be applied (e.g. mean, median, sd...)
#' @keywords internal
#' @export
setGeneric("compute_module_score",
           function(object,
                    modules=NULL,
                    rename=FALSE,
                    fun="mean"
           )
             standardGeneric("compute_module_score"))

#' @title From a list of feature (e.g gene) modules, compute a score (e.g. mean value) accross bins.
#'
#' @description
#' From a list of feature (e.g gene) modules, compute a score (e.g. mean value) accross bins.
#'
#' @param object An \code{STGrid} object.
#' @param modules A list of modules.
#' @param rename Force the renaming of the modules to module_1, module_2...
#' @param fun The function to be applied (e.g. mean, median, sd...)
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' p <- hc_tree(xen, class_nb = 4, class_name = letters[1:4])
#' print(p)
#' # Get classes:
#' xen <- compute_module_score(xen, modules=p$tree_classes)
#' @export
setMethod("compute_module_score", "STGrid",
          function(object,
                   modules=NULL,
                   rename=FALSE,
                   fun="mean"
                   ){

            print_this_msg("Computing module score...")

            if(is.null(modules))
              print_this_msg("Need at list one module", msg_type = "STOP")

            if(!is.list(modules))
              print_this_msg("'modules' argument should be a list", msg_type = "STOP")

            if(is.null(names(modules)) | rename == TRUE)
              names(modules) <- paste0("module_", 1:length(modules))

            print_this_msg("Iterating over modules...")
            for(mod_nm in names(modules)){
              print_this_msg("Processing module: ", mod_nm, msg_type = "DEBUG")
              check_st_list(list(object), feat_list = modules[[mod_nm]])
              if(fun=="mean"){
                row_m <- rowMeans(bin_mat(object)[, modules[[mod_nm]], drop=FALSE])
                object@meta[[mod_nm]] <-row_m
              }else{
                object@meta[[mod_nm]] <- apply(bin_mat(object)[, modules[[mod_nm]], drop=FALSE],
                                               1,
                                               eval(parse(text =fun)))
              }

            }

      print_this_msg("Returning object...")

      return(object)
})

# -------------------------------------------------------------------------
#      Convert bin_mat slot to x/y matrix (i.e. dcast bin_mat)
#      (e.g. to be displayed using image())
# -------------------------------------------------------------------------
#' @title Convert the bin_mat slot of an STGrid object into an x/y matrix (e.g. to be displayed using image())
#' @description
#' This function converts the bin_mat slot of an STGrid object into an x/y matrix (e.g. to be displayed using image())
#' @param object The STGrid object
#' @param feature The (meta)feature of interest.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' image(as_matrix(xen))
#' @export
setGeneric("as_matrix",
           function(object,
                    feature=NULL)
             standardGeneric("as_matrix"))

#' @title Convert the bin_mat slot of an STGrid object into an x/y matrix (e.g. to be displayed using image())
#' @description
#' This function converts the bin_mat slot of an STGrid object into an x/y matrix (e.g. to be displayed using image())
#' @param object The STGrid object
#' @param feature The (meta)feature of interest.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' image(as_matrix(xen))
#' @export
setMethod("as_matrix", "STGrid",
          function(object,
                   feature=NULL) {

            if(length(feature) > 1)
              print_this_msg("Only a single feature is accepted...")

            if(!check_features_exist(object, feature))
              print_this_msg("Please provide a feature")

            check_st_list(list(object), feature)

            bin_mat <- bin_mat(object,
                               feat_list=feature,
                               as_factor = TRUE)

            spatial_mat <- reshape2::dcast(bin_mat,
                                           bin_x~bin_y,
                                           value.var=feature)

            rownames(spatial_mat) <- spatial_mat[,1]
            spatial_mat <- spatial_mat[,-1]


            return(as.matrix(spatial_mat[bin_x(object), bin_y(object)]))

          }
)


# -------------------------------------------------------------------------
#      Find contiguous cells/bins
# -------------------------------------------------------------------------
#' @title Find contiguous cells using Queen criterion (internal).
#' @description
#' Find contiguous cells in a matrix using Queen criterion. The matrix is first binarized
#' (any value greater than zero is set to 1, any value lower than 1 is set to 0).
#' @param mat A numeric matrix.
#' @param threshold A threshold for binarization.
#' @return A list of matrices.
#' @examples
#' set.seed(123)
#' m <- matrix(sample(0:2, replace=TRUE, size = 26*26, prob = c(0.98, 0.01, 0.01)), nc=26)
#' rownames(m) <- letters
#' colnames(m) <- LETTERS
#' image(m)
#' image(find_contiguous(m))
#' image(find_contiguous(m, ns=2)[[1]])
#' image(find_contiguous(m, threshold=2))
#' @export
#' @keywords internal
find_contiguous <- function(mat,
                            threshold=1) {

  print_this_msg("Preparing matrix", msg_type = "DEBUG")


  mat[mat >= threshold] <- threshold
  mat[mat < threshold] <- 0
  mat[mat == threshold] <- 1
  mat <- as.matrix(mat)

  nr <- nrow(mat)
  nc <- ncol(mat)

  flagged_mat <- mat

  for (i in 1:nr) {
    for (j in 1:nc) {
      if (mat[i, j] == 1) {
        for (di in -1:1) {
          for (dj in -1:1) {
            ni <- i + di
            nj <- j + dj
            if (ni > 0 && ni <= nr && nj > 0 && nj <= nc) {
              flagged_mat[ni, nj] <- 1
            }
          }
        }
      }
    }
  }
  return(flagged_mat)
}

# -------------------------------------------------------------------------
#      Find the neighbors of a Connected Component
# -------------------------------------------------------------------------

#' @title Find the neighbors of a Connected Component.
#' @description
#' Given an STGrid object and a Connected Component (see connected_components()) find the neighboring bins.
#' @param object An STGrid object containing the spatial expression data.
#' @param feature A character string specifying the feature of interest.
#' @return The input STGrid object with additional fields (in meta slot) indicating the satellite bins for the specified feature at each neighborhood size.
#' @examples
#' library(patchwork)
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' p1 <- spatial_image(xen, feat="Necab2")
#' xen <- connected_components(xen, feat_list = "Necab2", threshold = 5, min_size = 30)
#' head(meta(xen))
#' p2 <- spatial_image(xen, feat="Necab2_cpt", as_factor=TRUE, colors = c("black", rainbow(9)))
#' hull_1 <- create_hull(xen, feat="Necab2_cpt", color="red", linew=0.3)
#' p3 <- p2 + hull_1
#' p4 <- spatial_image(xen, feat="Nwd2") + hull_1
#' (p1 | p2 ) / (p3 | p4)
#' @export
#' @keywords internal
setGeneric("cc_neighborhood",
           function(object,
                    feature=NULL)
             standardGeneric("cc_neighborhood"))


#' @title Find the neighbors of a Connected Component.
#' @description
#' Given an STGrid object and a Connected Component (see connected_components()) find the neighboring bins.
#' @param object An STGrid object containing the spatial expression data.
#' @param feature A character string specifying the feature of interest.
#' @return The input STGrid object with additional fields (in meta slot) indicating the satellite bins for the specified feature at each neighborhood size.
#' @examples
#' library(patchwork)
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' p1 <- spatial_image(xen, feat="Necab2")
#' xen <- connected_components(xen, feat_list = "Necab2", threshold = 5, min_size = 30)
#' head(meta(xen))
#' p2 <- spatial_image(xen, feat="Necab2_cpt", as_factor=TRUE, colors = c("black", rainbow(9)))
#' hull_1 <- create_hull(xen, feat="Necab2_cpt", color="red", linew=0.3)
#' p3 <- p2 + hull_1
#' p4 <- spatial_image(xen, feat="Nwd2") + hull_1
#' (p1 | p2 ) / (p3 | p4)
#' @export
setMethod("cc_neighborhood",
          "STGrid",
          function(object,
                   feature=NULL){

            if(is.null(feature) | length(feature) > 1)
              print_this_msg("Please provide one feature...", msg_type = "STOP")

            bin_mat <- bin_mat(object,
                               feat_list=feature,
                               as_factor = TRUE)

            if(!is.numeric(bin_mat[[feature]]))
              print_this_msg("The feature values should be numerical !", msg_type = "STOP")

            if(!all(bin_mat[[feature]]-floor(bin_mat[[feature]])==0))
              print_this_msg("The feature values should be integers !", msg_type = "STOP")

            print_this_msg("Binarizing...")

            bin_mat[[feature]][bin_mat[[feature]] >= 1] <- 1
            bin_mat[[feature]][bin_mat[[feature]] < 1] <- 0

            spatial_mat <- reshape2::dcast(bin_mat, bin_x~bin_y, value.var=feature)
            rownames(spatial_mat) <- spatial_mat[,1]
            spatial_mat <- spatial_mat[,-1]

            print_this_msg("Looking for contiguous bins...")

            isc <- find_contiguous(spatial_mat,
                                   threshold=1)


            tmp <- reshape2::melt(t(isc))
            colnames(tmp) <- c("bin_y", "bin_x", "value")
            rownames(tmp) <- paste0(tmp$bin_x, "~", tmp$bin_y)
            isc <- tmp[rownames(object@meta), ]


            object[[paste0(feature, "_ngb")]] <- isc$value


            return(object)
          })

# -------------------------------------------------------------------------
#      Connected Components Labelling
# -------------------------------------------------------------------------

#' @title Given a Feature and an STGrid object, Label and Store Connected Components.
#' @description
#' This function identifies and labels connected components in a binary matrix representation of features of an STGrid object.
#' The connected components are determined using the specified threshold and a minimum size for the components.
#' @param object An STGrid object containing the spatial expression data.
#' @param feat_list A character vector specifying the features of interest.
#' @param threshold A numeric value specifying the threshold for binarizing the matrix. If NULL, the matrix should already be binary.
#' @param min_size An integer specifying the minimum size of the connected components to retain. Default is 9.
#' @return The input STGrid object with additional fields for each feature's connected components.
#' @importFrom mgc ConnCompLabel
#' @importFrom reshape2 melt
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' xen <- connected_components(xen, feat_list = feat_names(xen)[1:5], threshold = 1, min_size = 9)
#' head(meta(xen))
#' spatial_image(xen, feat=meta_names(xen)[-1], ncol=2, as_factor=TRUE, colors = c("black", rainbow(60)))
#' @export
#' @keywords internal
setGeneric("connected_components",
           function(object=NULL,
                    feat_list=NULL,
                    threshold=NULL,
                    min_size = 9)
             standardGeneric("connected_components"))

#' @title Given a Feature and an STGrid object, Label and Store Connected Components.
#' @description
#' This function identifies and labels connected components in a binary matrix representation of features of an STGrid object.
#' The connected components are determined using the specified threshold and a minimum size for the components.
#' @param object An STGrid object containing the spatial expression data.
#' @param feat_list A character vector specifying the features of interest.
#' @param threshold A numeric value specifying the threshold for binarizing the matrix. If NULL, the matrix should already be binary.
#' @param min_size An integer specifying the minimum size of the connected components to retain. Default is 9.
#' @return The input STGrid object with additional fields for each feature's connected components.
#' @importFrom mgc ConnCompLabel
#' @importFrom reshape2 melt
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' xen <- connected_components(xen, feat_list = feat_names(xen)[1:5], threshold = 1, min_size = 9)
#' head(meta(xen))
#' spatial_image(xen, feat=meta_names(xen)[-1], ncol=2, as_factor=TRUE, colors = c("black", rainbow(60)))
#' @export
setMethod("connected_components",
          "STGrid",
          function(object=NULL,
                   feat_list=NULL,
                   threshold=NULL,
                   min_size = 9){

  if(is.null(feat_list))
    print_this_msg("please provide at least a single feature...", msg_type = "STOP")

  if(!check_features_exist(object, feat_list))
    print_this_msg("Please provide an available feature", msg_type = "STOP")

  check_st_list(list(object), feat_list)

  for(feature in feat_list){

    mat <- as_matrix(object, feature)

    if(!all(unique(mat) %in% c(0,1))){
      if(is.null(threshold)){
        print_this_msg("This should be a binary matrix. Maybe use 'threshold'.", msg_type = "STOP")
      }else{
        mat[mat >= threshold] <- threshold
        mat[mat < threshold] <- 0
        mat[mat == threshold] <- 1
      }

    }else{
      print_this_msg("This is already a binarized matrix. Skipping thresholding...")
    }

    components <- mgc::ConnCompLabel(mat)

    components <- reshape2::melt(components)
    rownames(components) <- paste0(components$Var1, "~", components$Var2)
    components <- components[rownames(object@meta),]

    value_as_factor <- as.factor(components$value)
    cpt_size <- sort(table(value_as_factor), decreasing = TRUE)

    max_size_cpt <- as.numeric(names(cpt_size)[1])

    for(i in 1:length(names(cpt_size))){

      if(cpt_size[i] < min_size){
        print_this_msg("Deleting compartment ", names(cpt_size[i]), "(see min_size).",
                       msg_type = "DEBUG")
        components$value[components$value == as.numeric(names(cpt_size)[i])] <- max_size_cpt
      }

    }

    cpt_name <- paste0(feature, "_cpt")
    object[[cpt_name]] <- components$value
    object[[cpt_name]] <- as.factor(object[[cpt_name]])
    levels(object[[cpt_name]]) <- 0:length(levels(object[[cpt_name]]))
    object[[cpt_name]] <- as.numeric(object[[cpt_name]]) - 1

  }

  return(object)
})


# -------------------------------------------------------------------------
#      Compute Moran's Index
# -------------------------------------------------------------------------

setGeneric("compute_moran_index",
           function(object=NULL,
                    queen_or_rook = c("queen", "rook"),
                    scaling=TRUE)
             standardGeneric("compute_moran_index"))

#' @title Compute Moran's Index.
#' @description
#' Compute Moran's Index.
#' @param queen_or_rook Method for neighborhood computation. Default to 'rook'. May also be queen.
#' @param  scaling	Whether to scale so that the sum of counts equals to 100 for each feature. Default TRUE.
#' @importFrom spdep cell2nb moran nb2listw Szero
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' xen <- moran_index(re_bin(xen, bin_size = 400))
#' print(xen@moran_index)
#' @export
setMethod("compute_moran_index",
          "STGrid",
          function(object=NULL,
                   queen_or_rook = c("queen", "rook"),
                   scaling=TRUE){

            print_this_msg("Starting computation of Moran's index...Please be patient")

            queen_or_rook <- match.arg(queen_or_rook)

            nbin_x <- nbin_x(object)
            nbin_y <- nbin_y(object)

            if(nbin_y * nbin_x > 10000){
              print_this_msg("The number of bins (x * y) is greater than 1e4.")
              print_this_msg("Computing Moran's Index may be time consuming.")
              print_this_msg("Please consider increasing bin size or selecting a sub-region.")
            }

            xygrid <- expand.grid(x = 1:nbin_x,
                                  y = 1:nbin_y)

            print_this_msg("Computing neighborhood.")

            neighborhood <- spdep::cell2nb(nrow = nbin_y,
                                  ncol = nbin_x,
                                  type = queen_or_rook,
                                  legacy = FALSE)

            print_this_msg("Calling adespatial::moran.randtest().")
            nb_as_list <- spdep::nb2listw(neighborhood)

            # moran_index <- moran.randtest(bin_mat(object = object,
            #                                   del_bin = TRUE,
            #                                   feat_list = feat_names(object)),
            #                           nb_as_list)

            bin_mat <- bin_mat(object = object,
                               del_bin = TRUE,
                               feat_list = feat_names(object), )

            moran_index <- data.frame(row.names=feat_names(object))
            tmp_i <- vector()
            tmp_k <- vector()

            for(i in 1:length(feat_names(object))){
               bm <- bin_mat[,i]

               if(scaling){
                bm <- bin_mat[,i] / sum(bin_mat[,i]) * 100
              }

              tmp <- spdep::moran(bm,
                                     nb_as_list,
                                     length(neighborhood),
                                     spdep::Szero(nb_as_list))
              tmp_i[i] <- tmp$I
              tmp_k[i] <- tmp$K
            }

            moran_index$I <- tmp_i
            moran_index$K <- tmp_k


            object@moran_index <- moran_index

            return(object)
})

# -------------------------------------------------------------------------
#      Get the min / max value from Moran's Index
# -------------------------------------------------------------------------
#' Returns features ordered by Moran's Index
#'
#' This method returns features ordered by Moran's Index values stored in
#' an \code{STGrid}.
#' @param object An \code{STGrid} object.
#' @keywords internal
#' @export
setGeneric("order_feat_by_moran_index",
           function(object)
             standardGeneric("order_feat_by_moran_index"))

#' Returns features ordered by Moran's Index
#' This method returns features ordered by Moran's Index values stored in
#' an \code{STGrid}.
#' @param object An \code{STGrid} object.
#' @export
setMethod("order_feat_by_moran_index", "STGrid",
          function(object){
            features <- feat_names(object)
            if(nrow(object@moran_index) == 0){
              print_this_msg("No moran_index slot available in the object",
                        msg_type = "WARNING")
              return(features)
            }else{
              return(features[order(object@moran_index$I)])
            }
})
