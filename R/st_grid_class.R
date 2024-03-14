# -------------------------------------------------------------------------
## Spatial transcriptomic anlysis class (STGrid) definition
# -------------------------------------------------------------------------

#' Spatial transcriptomic anlysis class (STGrid).
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
    control = character()
  )
)


# -------------------------------------------------------------------------
#      NCOL/NROW/DIM METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------
#' @title
#' ncol.STGrid
#' @description
#' The number of columns of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
ncol.STGrid <- function (x) {
  length(x@bin_x)
}

#' @title
#' nrow.STGrid
#' @description
#' The number of rows of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
nrow.STGrid <- function (x) {
  length((x@bin_y))
}

#' @title Dimension of a STGrid object.
#' dim
#' @description
#' The dimension of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
setMethod("dim", signature(x = "STGrid"),
          function(x)
            dim(x@bin_mat))

#' @title Column names of a STGrid object.
#' @description
#' The column names of a STGrid object.
#' @param x The STGrid object
#' @export col_names
#' @keywords internal
setGeneric(
  name = "col_names",
  def = function(x)
    standardGeneric("col_names")
)

#' @title Column names of a STGrid object.
#' @description
#' The column names of a STGrid object.
#' @param x The STGrid object
#' @export col_names
setMethod(
  f = "col_names",
  signature = "STGrid",
  definition = function(x)
    x@bin_x
)

#' @title Row names of a STGrid object.
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' @export row_names
#' @keywords internal
setGeneric("row_names",
           function(x)
             standardGeneric("row_names"))

#' @title Row names of a STGrid object.
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' @export row_names
setMethod("row_names", "STGrid",
          function(x)
            x@bin_y)

#' @title X bins of a STGrid object.
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export bin_x
setGeneric("bin_x",
           function(x)
             standardGeneric("bin_x"))

#' @title X bins of a STGrid object.
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_x
setMethod("bin_x", "STGrid",
          function(x)
            x@bin_x)

#' @title Y bins of a STGrid object.
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export bin_y
setGeneric("bin_y",
           function(x)
             standardGeneric("bin_y"))

#' @title Y bins of a STGrid object.
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_y
setMethod("bin_y", "STGrid",
          function(x)
            x@bin_y)

#' @title Number of bins (x axis) of a STGrid object.
#' @description
#' Number of bins (x axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export
setGeneric("nbin_x",
           function(x)
             standardGeneric("nbin_x"))

#' @title Number of bins (x axis) of a STGrid object.
#' @description
#' Number of bins (x axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export
setMethod("nbin_x", "STGrid",
          function(x)
            length(bin_x(x)))

#' @title Number of bins (y axis) of a STGrid object.
#' @description
#' Number of bins (y axis) of a STGrid object.
#' @param x The STGrid object
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
#' @export
setMethod("nbin_y", "STGrid",
          function(x)
            length(bin_y(x)))
#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param object The STGrid object
#' @export coord
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @keywords internal
setGeneric("coord",
           function(object,
                    as_factor = FALSE)
             standardGeneric("coord"))

#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @export coord
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
#' @keywords internal
#' @export
setGeneric("bin_mat",
           function(object,
                    as_factor = FALSE,
                    feat_list = character(),
                    melt_tab = FALSE,
                    del_bin = FALSE)
             standardGeneric("bin_mat"))

#' @title The binned matrix stored in a STGrid object.
#' @description
#' The binned matrix stored in a STGrid object.
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @param feat_list Whether to subset to some features.
#' @param melt_tab Whether to melt.
#' @param del_bin Whether to delete the bin_x, bin_y columns.
#' @export
setMethod("bin_mat", "STGrid",
          function(object,
                   as_factor = FALSE,
                   feat_list = character(),
                   melt_tab = FALSE,
                   del_bin = FALSE) {
            if (length(feat_list) == 0) {
              feat_list <- feat_names(object)
            }

            this_bin_mat <- object@bin_mat
            if (del_bin) {
              this_bin_mat <- this_bin_mat[, feat_list]
            } else{
              this_bin_mat <- this_bin_mat[, c("bin_x", "bin_y", feat_list)]
            }


            if (as_factor) {
              this_bin_mat$bin_x <- factor(this_bin_mat$bin_x,
                                           levels = bin_x(object),
                                           ordered = TRUE)

              this_bin_mat$bin_y <- factor(this_bin_mat$bin_y,
                                           levels = bin_y(object),
                                           ordered = TRUE)
            }


            if (melt_tab) {
              print_msg("Melting data.frame...", msg_type = "DEBUG")
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


            return(this_bin_mat)
          })

#' @title The number of molecules stored in a STGrid object
#' @description
#' The number of molecules stored in a STGrid object
#' @param x The STGrid object
#' @export nb_molec
#' @keywords internal
setGeneric("nb_molec",
           function(x)
             standardGeneric("nb_molec"))

#' @title The number of molecules stored in a STGrid object
#' @description
#' The number of molecules stored in a STGrid object
#' @param x The STGrid object
#' @export nb_molec
setMethod("nb_molec", "STGrid",
          function(x)
            nrow(x@coord))


#' @title The number of features stored in a STGrid object
#' @description
#' The number of features stored in a STGrid object
#' @param object The STGrid object
#' @keywords internal
#' @export nb_feat
setGeneric("nb_feat",
           function(object) {
             standardGeneric("nb_feat")
           })


#' @title The number of features stored in a STGrid object
#' @description
#' The number of features stored in a STGrid object
#' @param object The STGrid object
#' @export nb_feat
setMethod("nb_feat", signature(object = "STGrid"),
          function(object) {
            length(feat_names(object))
          })


#' @title The features stored in a STGrid object
#' @description
#' The features stored in a STGrid object
#' @param object The STGrid object
#' @param del_control Whether to delete controls.
#' @keywords internal
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
#' @export
setMethod("feat_names", signature(object = "STGrid"),
          function(object,
                   del_control = FALSE) {

            fn <- grep(
                "^bin_[xy]$",
                colnames(object@bin_mat),
                invert = TRUE,
                perl = TRUE,
                val = TRUE
            )


            if (del_control) {
              print(object@control)
              fn <-
                fn[!fn %in% grep(object@control, fn, val = TRUE)]
            }

            fn
          })


#' @title The size of the bins stored in an STGrid object
#' @description
#' The size of the bins stored in an STGrid object
#' @param x The STGrid object
#' @export bin_size
#' @keywords internal
setGeneric("bin_size",
           function(x)
             standardGeneric("bin_size"))

#' @title The size of the bins stored in an STGrid object
#' @description
#' The size of the bins stored in an STGrid object
#' @param x The STGrid object
#' @export bin_size
setMethod("bin_size", "STGrid",
          function(x)
            x@bin_size)

#' @title Remove control (Blank-*) features from a STGrid object.
#' @description
#'  Remove control (Blank-*) features.
#' @param object The STGrid object
#' @param regexp A regexp to search for internal controls.
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
#' @export
setMethod("rm_controls", "STGrid",
          function(object,
                   regexp = NULL) {
            if (is.null(regexp)) {
              regexp <- object@control
            }
            pos <- -grep(regexp, feat_names(object))
            if(length(pos) > 0){
              return(object[pos, ])
            }else{
              return(object)
            }

          })

#' @title  Get Ripley's K function slot from a STGrid object.
#' @description
#'   Get Ripley's K function slot from a STGrid object.
#' @param x The STGrid object
#' @export rm_controls
#' @keywords internal
setGeneric("ripley_k_function",
           function(object)
             standardGeneric("ripley_k_function"))

#' @title Get Ripley's K function slot from a STGrid object.
#' @description
#'  Get Ripley's K function slot from a STGrid object.
#' @param x The STGrid object
#' @export bin_size
#' @keywords internal
#' @export ripley_k_function
setMethod("ripley_k_function", "STGrid",
          function(object)
            object@ripley_k_function)


# -------------------------------------------------------------------------
##    REDEFINE SHOW() METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------

#' @title
#' The show method of a STGrid object
#' @description
#' The show method of a STGrid object
#' @param object A STGrid object.
#' @keywords internal
#' @export
setMethod("show", signature("STGrid"),
          function(object) {
            print_msg("An object of class STGrid")
            print_msg("Memory used: ", object.size(object))
            print_msg("Number of counts: ", nb_molec(object))
            print_msg("Number of features: ", nb_feat(object))
            print_msg("Bin size: ", bin_size(object))
            print_msg("Number of bins (x axis): ", length(bin_x(object)))
            print_msg("Number of bins (y axis): ", length(bin_y(object)))
            print_msg("x_min: ", object@x_min)
            print_msg("x_max: ", object@x_max)
            print_msg("y_min: ", object@y_min)
            print_msg("y_max: ", object@y_max)
            print_msg(">>> Please, use show_methods() to show availables methods <<<")
          })

# -------------------------------------------------------------------------
##    REDEFINE SUMMARY() METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------

#' @title
#' The summary() method of a STGrid object
#' @description
#' The summary method of a STGrid object
#' @param object A STGrid object.
#' @export summary
setMethod("summary", signature("STGrid"),
          function(object) {
            mx_col <- min(ncol(object), 4)
            print_msg("An object of class STGrid")
            print_msg("Method:", object@method)
            print_msg("Bin size:", object@bin_size)
            print_msg("Features:", c(head(feat_names(object), 4), "..."))
            print_msg("x_min:", object@x_min)
            print_msg("x_max:", object@x_max)
            print_msg("y_min:", object@y_min)
            print_msg("y_max:", object@y_max)
            print_msg("path:", object@path)
            print_msg("bin_x:", head(bin_x(object), 4), "...")
            print_msg("bin_y:", head(bin_y(object), 4), "...")
            print_msg("Meta data:", paste0(colnames(object@meta), collapse=" , "))
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
                      print_msg("Numering value i is out of range (should be < nb_feat(x)).",
                                msg_type = "STOP")
                    }
                    i <- feat_names(x)[i]
                    x_is_feat <- TRUE
                  } else{
                    print_msg("Check i... Should be a set of features or bins.",
                              msg_type = "STOP")
                  }

                }

              }

              if (!x_is_bin && !x_is_feat)
                print_msg("Some features or bin_x or indexes where not found in the object.",
                          msg_type = "STOP")
            }

            if (!missing(j)) {
              if (!all(j %in% bin_y(x)))
                print_msg("Some bin_y where not found in the object.", msg_type = "STOP")
            }

            n_coord <- x@coord
            n_bin_mat <- x@bin_mat
            n_bin_x <- x@bin_x
            n_bin_y <- x@bin_y
            n_ripley_k_function <- x@ripley_k_function
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

                } else {
                  n_coord <- n_coord[n_coord$bin_x %in% i, ]
                  n_bin_mat <- n_bin_mat[n_bin_mat$bin_x %in% i, ]
                  feat_left <- unique(n_coord$feature)

                  n_bin_mat <-
                      n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", feat_left)]

                }
              }
            } else{
              if (missing(i)) {
                n_coord <- n_coord[n_coord$bin_y %in% j, ]
                n_bin_mat <- n_bin_mat[n_bin_mat$bin_y %in% j, ]

              } else{
                if (x_is_feat) {
                  n_coord <- n_coord[n_coord$feature %in% i, ]

                  n_bin_mat <-
                      n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", i)]



                } else {
                  n_coord <- n_coord[n_coord$bin_x %in% i, ]
                  n_bin_mat <- n_bin_mat[n_bin_mat$bin_x %in% i, ]
                  feat_left <- unique(n_coord$feature)

                  n_bin_mat <-
                      n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", feat_left)]

                }

                n_coord <- n_coord[n_coord$bin_y %in% j, ]
                n_bin_mat <- n_bin_mat[n_bin_mat$bin_y %in% j, ]
              }
            }

            STGrid_obj <- new("STGrid",
                              path = x@path)

            STGrid_obj@coord <- n_coord
            STGrid_obj@bin_mat <- n_bin_mat
            STGrid_obj@y_max <- x@y_max
            STGrid_obj@y_min <- x@y_min
            STGrid_obj@x_max <- x@x_max
            STGrid_obj@x_min <- x@x_min
            STGrid_obj@method <- x@method
            STGrid_obj@meta <- x@meta
            STGrid_obj@bin_size <- x@bin_size
            STGrid_obj@bin_x <-
              x@bin_x[x@bin_x %in% n_bin_mat$bin_x]
            STGrid_obj@bin_y <-
              x@bin_y[x@bin_y %in% n_bin_mat$bin_y]
            STGrid_obj@meta <- x@meta
            STGrid_obj@control <- n_control

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
#' @keywords internal
#' @export
setMethod ("$", "STGrid",
           function (x, name) {
             return(x[[name]])
           })



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
#' @export
setMethod ("$<-", "STGrid",
           function (x, name, value) {
             if (length(value) != nbin_x(x) * nbin_y(x) && length(value) != 1)
               print_msg("The size of the vector needs to be nbin_x(x) * nbin_y(x) or 1.",
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
#' @export
setMethod ("[[<-", "STGrid",
           function (x, name, value) {
             x@meta[[name]] <- value
             x
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
#' @param verbose Whether to print porgress bar.
#' @return An updated object of class "STGrid" with the Ripley's K function estimates stored in the slot 'ripley_k_function'.
#'
#' @examples
#' # Example usage:
#' # grid_object <- load_spatial(...)
#' # grid_object <- compute_k_ripley(grid_object)
#' @keywords internal
setGeneric("compute_k_ripley",
           function(object,
                    rmax = 200,
                    nlarge = 1e6,
                    var.approx = FALSE,
                    ratio = FALSE,
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
#' @param verbose Whether to print porgress bar.
#' @return An updated object of class "STGrid" with the Ripley's K function estimates stored in the slot 'ripley_k_function'.
#'
#' @examples
#' # Example usage:
#' # grid_object <- load_spatial(...)
#' # grid_object <- compute_k_ripley(grid_object)
#'
#' @importFrom spatstat.geom ppp owin
#' @importFrom spatstat.explore Kest
#' @export compute_k_ripley
setMethod("compute_k_ripley", signature("STGrid"),
          function(object,
                   rmax = 200,
                   nlarge = 1e6,
                   var.approx = FALSE,
                   ratio = FALSE,
                   verbose = TRUE) {

            molecules <- coord(object)[, c("x", "y", "feature")]

            data_out <- data.frame(
              r = NA,
              theo = NA,
              border = NA,
              trans = NA,
              iso = NA,
              feature = NA
            )

            x_min <- object@x_min
            x_max <- object@x_max
            y_min <- object@y_min
            y_max <- object@y_max

            rownames(molecules) <- paste0(1:nrow(molecules),
                                          '|',
                                          molecules$feature)


            print_msg(">>> Estimating Ripley's reduced second moment function for all features.")

            n_iter <- length(unique(molecules$feature))

            if (verbose) {
              pb <- txtProgressBar(
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
                setTxtProgressBar(pb, n)
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

              u <- spatstat.explore::Kest(
                X,
                rmax = rmax,
                nlarge = nlarge,
                var.approx = var.approx,
                ratio = ratio
              )

              u <- as.data.frame(as.matrix(u))
              u$feature <- g

              data_out <- rbind(data_out, u)
            }

            if (verbose)
              close(pb)

            data_out <- na.omit(data_out)

            object@ripley_k_function <- data_out

            return(object)

          })


# -------------------------------------------------------------------------
##    Constructor for SpatialTranscriptomicAnalysis class (STGrid)
# -------------------------------------------------------------------------
#' @title Create a Spatial Transcriptomic Grid class (STGrid)
#' @description
#' The load_spatial() function is the entry point of the stcompr package.
#' It will load the molecule coordinates and create a 2D binned grid with default size 25µm.
#'
#' @param path Either a file (if method is set to "coordinates") or a directory (if method
#' is set to "merscope"). If method is set to "coordinates" the file should contain 3 columns
#' ("x", "y", "gene"), ("x", "y", "feature") or ("x", "y", "cell").
#' @param method The type of technology.
#' @param bin_size Numeric value representing the bin size (default to 25).
#' @param control A regular expression to identify controls. As the function computes the sum of
#' counts, this will allow to delete these blanks/controls for computation.
#' @param verbose Whether to display the progress bar.
#' @return An object of class STGrid.
#'
#' @export load_spatial
load_spatial <- function(path = "",
                         method = c("coordinates",
                                    "merscope",
                                    "xenium"),
                         bin_size = 25,
                         control = NULL,
                         verbose = TRUE) {
  method <- match.arg(method)

  print_msg("Technology is '", method, "'.")
  print_msg("Loading data from file:", path)

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

  }  else if (method == "coordinates") {
    check_file(path, mode = "read")
    spat_input <- as.data.frame(data.table::fread(path,
                                                  sep = "\t",
                                                  head = TRUE))
    col_needed <- c("x", "y")

    if("cell" %in% colnames(spat_input)){
      col_needed <- c(col_needed, "cell")
    }else if("feature" %in% colnames(spat_input)){
      col_needed <- c(col_needed, "feature")
    }else if("gene" %in% colnames(spat_input)){
      col_needed <- c(col_needed, "gene")
    }else{
      print_msg("Could not find a cell/feature/gene columns.", msg_type = "STOP")
    }

    if (any(!col_needed %in% colnames(spat_input))) {
      print_msg("Please check the column name of the input file.", msg_type = "STOP")
    }

    spat_input <- spat_input[, col_needed]

    if (is.null(control))
      control <- "^Blank\\-[0-9]+"

  }

  colnames(spat_input)[3] <- "feature"

  bin_matrix <- bin_this_matrix(coord = spat_input,
                                bin_size = bin_size,
                                verbose = verbose)

  sum_of_cts <- sum_of_counts(bin_matrix$spatial_matrix,
                              control = control)

  # create a STGrid object                             ---------

  STGrid_obj <- new("STGrid",
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
               count_sums = sum_of_cts)
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
    print_msg("Please provide a spatial matrix.", msg_type = "STOP")
  }

  print_msg("Regexp for controls:", control, msg_type = "DEBUG")

  print_msg("Computing sum of counts.")

  pos_bin_xy <- grep("^bin_[xy]$",
                     colnames(spatial_matrix))

  tmp <- spatial_matrix[, -pos_bin_xy, drop=FALSE]

  if(!is.null(control)) {
    pos_ctrl <- grep(control, colnames(tmp), perl = TRUE)
  } else{
    pos_ctrl <- integer(0)
  }

  if (length(pos_ctrl) > 0) {
    print_msg("Found the following controls:",
              paste0(head(colnames(tmp)[pos_ctrl], 3),
                     collapse = ", "),
              "...")
    sum_of_cts <- rowSums(tmp[, -pos_ctrl, drop=FALSE])

  } else{

    sum_of_cts <- rowSums(tmp)
  }

  return(sum_of_cts)

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
#' # Create a binned grid from molecule coordinates
#' bin_mat(coord = my_coord, bin_size = 25)
#' @keywords internal
#' @export bin_this_matrix
bin_this_matrix <- function(coord = NULL,
                            bin_size = 25,
                            verbose = TRUE) {
  print_msg("Binning a matrix...")

  x_min <- min(coord$x)
  x_max <- max(coord$x)
  y_min <- min(coord$y)
  y_max <- max(coord$y)

  if (x_max - x_min < bin_size) {
    print_msg(
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
    print_msg(
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

  if (verbose) {
    pb <- txtProgressBar(
      min = 0,
      max = length(unique(coord$feature)),
      style = 3,
      width = 50,
      char = "="
    )

    n_loop <- 0
  }

  for (goi in unique(coord$feature)) {
    x_molec <- cut(
      coord$x[coord$feature == goi],
      breaks = x_lim,
      include.lowest = TRUE,
      right = FALSE
    )

    coord$bin_x[coord$feature == goi] <- as.character(x_molec)

    y_molec <- cut(
      coord$y[coord$feature == goi],
      breaks = y_lim,
      include.lowest = TRUE,
      right = FALSE
    )
    coord$bin_y[coord$feature == goi] <- as.character(y_molec)

    nb_molec <- table(paste(x_molec, y_molec, sep = "~"))
    spatial_matrix[, goi] <- 0
    spatial_matrix[names(nb_molec), goi] <- nb_molec

    if (verbose) {
      n_loop <- n_loop + 1
      setTxtProgressBar(pb, n_loop)
    }
  }

  cat("\n")

  spatial_matrix <-
    spatial_matrix[, order(colnames(spatial_matrix))]

  if (verbose)
    close(pb)

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
#' @export re_bin
#' @keywords internal
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
#' @export re_bin
setMethod("re_bin", signature(object = "STGrid"),
          function(object,
                   bin_size,
                   verbose = TRUE) {
            if (object@bin_size == bin_size) {
              print_msg("The bin_size is unchanged.")
              return(object)
            }

            bin_matrix <- bin_this_matrix(coord = object@coord,
                                          bin_size = bin_size,
                                          verbose = verbose)

            print_msg("Re-computing sum of counts .")
            sum_of_cts <- sum_of_counts(bin_matrix$spatial_matrix,
                                        control=object@control)

            # create a STGrid object                             ---------
            print_msg("Creating an STGrid object")
            print_msg("Note that meta data will be lost (except 'sum_of_cts').")

            STGrid_obj <- new("STGrid",
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
                                          count_sums = sum_of_cts)
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
#' @export
setMethod("get_coord", "STGrid",
          function(object,
                   feat_list = character(),
                   as.factor = TRUE) {
            if (!all(feat_list %in% feat_names(object))) {
              print_msg("Some features were not found...",
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
#' @param size The size of the labels.
#' @export
#' @keywords internal
setGeneric("hc_tree",
           function(object,
                    method = "complete",
                    layout = "circular",
                    dist_method = "pearson",
                    branch_length = "none",
                    class_nb = 1,
                    class_name = NULL,
                    size = 2.25)
             standardGeneric("hc_tree"))

#' @title Create a tree from an STGrid object
#' @description Create a tree from an STGrid object
#' @param object The STGrid object
#' @param method The agglomeration method to be used. See stats::hclust().
#' @param layout The layout to be used See ggtree::ggtree().
#' @param dist_method The method for distance computation. See stats::cor().
#' @param branch_length A variable for scaling branch, if 'none' draw cladogram. See ggtree::ggtree().
#' @param class_nb An integer indicating the desired number of groups.
#' @param size The size of the labels.
#' @importFrom ggtree ggtree geom_hilight geom_tippoint geom_tiplab MRCA
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggplot2 aes scale_color_viridis_d
#' @importFrom ggsci scale_fill_jco
#' @export
#' @keywords internal
setMethod("hc_tree", "STGrid",

          function(object,
                   method = "complete",
                   layout = "circular",
                   dist_method = "pearson",
                   branch_length = "none",
                   class_nb = 1,
                   class_name = "All",
                   size = 2.25) {
            if (class_nb > 0) {
              if (length(class_name) != class_nb) {
                print_msg("Please set the right number of class names.",
                          msg_type = "STOP")
              }
            } else{
              print_msg("The class_nb argument should be an integer > than 0...",
                        msg_type = "STOP")
            }

            bin_mat <-
              bin_mat(object, del_bin = TRUE)

            hc_clust <- hclust(as.dist((1 - cor(bin_mat,
                                                method = dist_method)) / 2),
                               method = method)

            p <- ggtree::ggtree(hc_clust,
                                layout = layout,
                                branch.length = branch_length)

            tree_classes <- cutree(hc_clust, k = class_nb)
            groups <- split(names(tree_classes), tree_classes)
            clades <-
              sapply(groups, function(n)
                tidytree::MRCA(p, n))
            annotation <- data.frame(id = clades[1:class_nb],
                                     Class = class_name)


            #p$data$cell_type <- "ND"
            #pmath <- match(p$data$label, names(ginfo))
            #p$data$cell_type[!is.na(pmath)] <- names(ginfo)[pmath[!is.na(pmath)]]

            p <- p + ggtree::geom_hilight(
              data = annotation,
              extend = 0.2,
              mapping = aes(node = id, fill = Class),
              alpha = 0.3
            ) +
              ggsci::scale_fill_jco() +
              ggtree::geom_tiplab(
                ggplot2::aes(label = label),
                offset = 1,
                size = 2.25,
                color = 'black'
              ) +
              ggnewscale::new_scale_fill() +
              ggplot2::scale_color_viridis_d()

            # ggtree::geom_tippoint(
            # ggplot2::aes(color = cell_type),
            # size = 1,
            # inherit.aes = TRUE,
            #show.legend = FALSE
            #) +
            return(p)
          })


# -------------------------------------------------------------------------
#      Get the min / max value from ripley's k function
# -------------------------------------------------------------------------

setGeneric("order_feat_by_ripley",
           function(object)
             standardGeneric("order_feat_by_ripley"))

setMethod("order_feat_by_ripley", "STGrid",
          function(object){
            ripk <- ripley_k_function(object)
            voi <- ripk %>%
              dplyr::group_by(feature) %>%
              dplyr::filter(border == max(border)) %>%
              dplyr::filter(r == max(r)) %>%
              dplyr::arrange(desc(border))
            return(voi$feature)
})

