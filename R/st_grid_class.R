# -------------------------------------------------------------------------
## Spatial transcriptomic anlysis class (STGrid) definition
# -------------------------------------------------------------------------

#' Spatial transcriptomic anlysis class (STGrid).
#'
#' A class for representing spatial transcriptomic analysis results using (currently) the "merscope" method.
#' This class store transcript molecule coordinates as a grid.
#'
#' @slot coordinate A data frames representing spatial coordinates (empty by default).
#' May accept several layers in the future.
#' @slot bin_mat A data frames representing the matrix (empty by default).
#' May accept several layers in the future.
#' @slot y_max Numeric value representing the maximum y-coordinate (default is 0).
#' @slot y_min Numeric value representing the minimum y-coordinate (default is 0).
#' @slot x_min Numeric value representing the minimum x-coordinate (default is 0).
#' @slot x_max Numeric value representing the maximum x-coordinate (default is 0).
#' @slot path Character value representing the file path to the original data.
#' @slot meta A list containing meta-information.
#' @slot bin_size Numeric value representing the bin size.
#' @slot bin_x Character vector containing the names of bins/windows along the x-axis.
#' @slot bin_y Character vector containing the names of bins/windows along the y-axis.
#'
#' @export
setClass("STGrid",
         slots = c(
           coord = "data.frame",
           bin_mat = "data.frame",
           y_max = "numeric",
           y_min = "numeric",
           x_min = "numeric",
           x_max = "numeric",
           path = "character",
           method="character",
           meta = "data.frame",
           bin_size = "numeric",
           bin_x = "character",
           bin_y = "character",
           ripley_k_function="data.frame",
           control="character"
         ),
         prototype = list(
           coord = data.frame(),
           bin_mat = data.frame(),
           y_max = 0,
           y_min = 0,
           x_min = 0,
           x_max = 0,
           path = character(),
           method=character(),
           meta = data.frame(),
           bin_size = 0,
           bin_x = character(),
           bin_y = character(),
           ripley_k_function=data.frame(),
           control=character()
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
  length((x@bin_x))
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
setMethod("dim",signature(x="STGrid"),
          function(x)
            dim(x@bin_mat)
)

#' @title Column names of a STGrid object.
#' @description
#' The column names of a STGrid object.
#' @param x The STGrid object
#' @export col_names
#' @keywords internal
setGeneric(name="col_names",
           def=function(x)
             standardGeneric("col_names")
)

#' @title Column names of a STGrid object.
#' @description
#' The column names of a STGrid object.
#' @param x The STGrid object
#' @export col_names
setMethod(f="col_names",
          signature="STGrid",
          definition=function(x) x@bin_x
)

#' @title Row names of a STGrid object.
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' @export row_names
#' @keywords internal
setGeneric("row_names",
           function(x)
             standardGeneric("row_names")
)

#' @title Row names of a STGrid object.
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' @export row_names
setMethod("row_names", "STGrid",
          function(x)
            x@bin_y
)

#' @title X bins of a STGrid object.
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export bin_x
setGeneric("bin_x",
           function(x)
             standardGeneric("bin_x")
)

#' @title X bins of a STGrid object.
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_x
setMethod("bin_x", "STGrid",
          function(x)
            x@bin_x
)

#' @title Y bins of a STGrid object.
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export bin_y
setGeneric("bin_y",
           function(x)
             standardGeneric("bin_y")
)

#' @title Y bins of a STGrid object.
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_y
setMethod("bin_y", "STGrid",
          function(x)
            x@bin_y
)

#' @title Number of bins (x axis) of a STGrid object.
#' @description
#' Number of bins (x axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export
setGeneric("nbin_x",
           function(x)
             standardGeneric("nbin_x")
)

#' @title Number of bins (x axis) of a STGrid object.
#' @description
#' Number of bins (x axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export
setMethod("nbin_x", "STGrid",
           function(x)
             length(bin_x(x))

)

#' @title Number of bins (y axis) of a STGrid object.
#' @description
#' Number of bins (y axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export
setGeneric("nbin_y",
           function(x)
             standardGeneric("nbin_y")
)

#' @title Number of bins (y axis) of a STGrid object.
#' @description
#'  Number of bins (y axis) of a STGrid object.
#' @param x The STGrid object
#' @keywords internal
#' @export
setMethod("nbin_y", "STGrid",
           function(x)
             length(bin_y(x))

)
#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param object The STGrid object
#' @export coord
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @keywords internal
setGeneric("coord",
           function(object,
                    as_factor=FALSE)
             standardGeneric("coord")
)

#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @export coord
setMethod("coord", "STGrid",
          function(object,
                   as_factor=FALSE){
            if(!as_factor){
              return(object@coord)
            }else{
              tmp <- object@coord
              tmp$bin_x <- factor(tmp$bin_x,
                                  levels = bin_x(object),
                                  ordered=TRUE)
              tmp$bin_y <- factor(tmp$bin_y,
                                  levels = bin_y(object),
                                  ordered=TRUE)
              return(tmp)
            }
          }
)

#' @title The binned matrix stored in a STGrid object
#' @description
#' The binned matrix stored in a STGrid object
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @param gene_list Whether to subset to some genes.
#' @param melt Whether to melt.
#' @keywords internal
#' @export
setGeneric("bin_mat",
           function(object,
                    as_factor=FALSE,
                    gene_list=character(),
                    melt=FALSE)
             standardGeneric("bin_mat")
)

#' @title The binned matrix stored in a STGrid object.
#' @description
#' The binned matrix stored in a STGrid object.
#' @param object The STGrid object
#' @param as_factor Should bin_x and bin_y columns be returned as ordered factor ?
#' @param gene_list Whether to subset to some genes.
#' @param melt Whether to melt.
#' @export
setMethod("bin_mat", "STGrid",
          function(object,
                   as_factor=FALSE,
                   gene_list=character(),
                   melt=FALSE){

            if(length(gene_list) == 0){
              gene_list <- gene_names(object, all_genes = TRUE)
            }

            this_bin_mat <- object@bin_mat
            this_bin_mat <- this_bin_mat[, c("bin_x", "bin_y", gene_list)]

            if(!as_factor){
              return(this_bin_mat)
            }else{
                this_bin_mat$bin_x <- factor(this_bin_mat$bin_x,
                                            levels = bin_x(object),
                                            ordered=TRUE)

                this_bin_mat$bin_y <- factor(this_bin_mat$bin_y,
                                            levels = bin_y(object),
                                            ordered=TRUE)
            }

            if(melt){
              this_bin_mat <- reshape2::melt(this_bin_mat, id.vars=c("bin_x", "bin_y"))
              colnames(this_bin_mat) <- c("bin_x", "bin_y", "gene", "value")
            }

              return(this_bin_mat)
}
)

#' @title The number of molecules stored in a STGrid object
#' @description
#' The number of molecules stored in a STGrid object
#' @param x The STGrid object
#' @export nb_molec
#' @keywords internal
setGeneric("nb_molec",
           function(x)
             standardGeneric("nb_molec")
)

#' @title The number of molecules stored in a STGrid object
#' @description
#' The number of molecules stored in a STGrid object
#' @param x The STGrid object
#' @export nb_molec
setMethod("nb_molec", "STGrid",
          function(x)
            nrow(x@coord)
)


#' @title The number of genes stored in a STGrid object
#' @description
#' The number of genes stored in a STGrid object
#' @param x The STGrid object
#' @keywords internal
#' @export
setGeneric("nb_genes",
             function(object){
               standardGeneric("nb_genes")
             }
)


#' @title The number of genes stored in a STGrid object
#' @description
#' The number of genes stored in a STGrid object
#' @param object The STGrid object
#' @export
setMethod("nb_genes", signature(object = "STGrid"),
          function(object){
            length(gene_names(object))
          }
)


#' @title The genes stored in a STGrid object
#' @description
#' The genes stored in a STGrid object
#' @param object The STGrid object
#' @keywords internal
#' @export
setGeneric("gene_names",
           function(object, all_genes=FALSE)
                  standardGeneric("gene_names")

)

#' @title The genes stored in a STGrid object
#' @description
#' The genes stored in a STGrid object
#' @param object The STGrid object
#' @param all_genes Whether "all_genes" should be included. Default to FALSE.
#' @export
setMethod("gene_names", signature(object = "STGrid"),
           function(object, all_genes=FALSE){

             if(all_genes){
               grep("^bin_[xy]$",
                    colnames(object@bin_mat),
                    invert = TRUE,
                    val=TRUE)
             }else{
               grep("(^bin_[xy]$)|(^all_genes$)",
                    colnames(object@bin_mat),
                    invert = TRUE,
                    perl=TRUE,
                    val=TRUE)
             }
           }
)


#' @title The size of the bins stored in an STGrid object
#' @description
#' The size of the bins stored in an STGrid object
#' @param x The STGrid object
#' @export bin_size
#' @keywords internal
setGeneric("bin_size",
           function(x)
             standardGeneric("bin_size")
)

#' @title The size of the bins stored in an STGrid object
#' @description
#' The size of the bins stored in an STGrid object
#' @param x The STGrid object
#' @export bin_size
setMethod("bin_size", "STGrid",
          function(x)
            x@bin_size
)

#' @title Remove control (Blank-*) "genes" from a STGrid object.
#' @description
#'  Remove control (Blank-*) "genes".
#' @param x The STGrid object
#' @keywords internal
setGeneric("rm_controls",
           function(x,
                    regexp="^Blank\\-.*")
             standardGeneric("rm_controls")
)

#' @title Remove control (Blank-*) "genes" from a STGrid object.
#' @description
#'  Remove control (Blank-*) "genes".
#' @param x The STGrid object
#' @export rm_controls
setMethod("rm_controls", "STGrid",
          function(x,
                   regexp="^Blank\\-.*")
            return(x[-grep(regexp, gene_names(x)), ])
)

#' @title  Get Ripley's K function slot from a STGrid object.
#' @description
#'   Get Ripley's K function slot from a STGrid object.
#' @param x The STGrid object
#' @export rm_controls
#' @keywords internal
setGeneric("ripley_k_function",
           function(object)
             standardGeneric("ripley_k_function")
)

#' @title Get Ripley's K function slot from a STGrid object.
#' @description
#'  Get Ripley's K function slot from a STGrid object.
#' @param x The STGrid object
#' @export bin_size
#' @keywords internal
#' @export ripley_k_function
setMethod("ripley_k_function", "STGrid",
          function(object)
            object@ripley_k_function
)


# -------------------------------------------------------------------------
##    REDEFINE SHOW() METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------

#' @title
#' The show method of a STGrid object
#' @description
#' The show method of a STGrid object
#' @param object A STGrid object.
#' @export show
#' @keywords internal
setMethod(
  "show", signature("STGrid"),
  function(object) {
    print_msg("An object of class STGrid")
    print_msg("Memory used: ", object.size(object))
    print_msg("Number of molecules: ", nb_molec(object))
    print_msg("Number of genes: ", nb_genes(object))
    print_msg("Bin size: ", bin_size(object))
    print_msg("Number of bins (x axis): ", length(bin_x(object)))
    print_msg("Number of bins (y axis): ", length(bin_y(object)))
    print_msg("x_min: ", object@x_min)
    print_msg("x_max: ", object@x_max)
    print_msg("y_min: ", object@y_min)
    print_msg("y_max: ", object@y_max)
    print_msg(">>> Please, use show_methods() to show availables methods <<<")
  }
)

# -------------------------------------------------------------------------
##    REDEFINE SUMMARY() METHOD FOR CLASS OBJECT : STGrid
# -------------------------------------------------------------------------

#' @title
#' The summary() method of a STGrid object
#' @description
#' The summary method of a STGrid object
#' @param object A STGrid object.
#' @export summary
setMethod(
  "summary", signature("STGrid"),
  function(object) {
    print_msg("An object of class STGrid")
    print_msg("Method:", object@method)
    print_msg("Bin size:", object@bin_size)
    print_msg("Coord:")
    cat("\n")
    print(head(object@coord, 3))
    cat("\n")
    print_msg("bin_mat:")
    cat("\n")
    print(head(object@bin_mat, 3))
    cat("\n")
    print_msg("x_min:", object@x_min)
    print_msg("x_max:", object@x_max)
    print_msg("y_min:", object@y_min)
    print_msg("y_max:", object@y_max)
    print_msg("path:", object@path)
    print_msg("bin_x:", head(bin_x(object)), "...")
    print_msg("bin_y:", head(bin_y(object)), "...")
    print_msg("Meta data:", object@meta)
  }
)

# -------------------------------------------------------------------------
##    Method for function"[". Subsetting an STGrid obj.
# -------------------------------------------------------------------------
#' @title Subsetting operator of a STGrid object
#' Extract
#' @description
#' The subsetting operator of a STGrid object.
#' The i axis correspond to genes. j axis is not implemented yet.
#' @param i indices specifying genes to extract or substract.
#' @param ... See ?'['. Not functionnal here.
#' @param drop For matrices and arrays. If TRUE the result is coerced to the lowest possible dimension. Not functionnal here.
#' @keywords internal
setMethod("[", signature(x = "STGrid"),
          function (x, i, j, ..., drop = FALSE) {


            x_is_gene <- FALSE
            x_is_bin <- FALSE
            x_is_num <- FALSE

            if(!missing(i)){

              i <- unique(i)

              if(any(i %in% bin_x(x))){
                x_is_bin <- TRUE
              }else{
                x_is_bin <- FALSE
                if(any(i %in% gene_names(x, all_genes=TRUE))){
                  x_is_gene <- TRUE
                }else{
                  x_is_gene <- FALSE
                  if(is.numeric(i)){
                    if(any(i > nb_genes(x))){
                      print_msg("Numering value i is out of range (should be < nb_genes(x)).", msg_type = "STOP")
                    }
                    x_is_num <- TRUE
                  }

                }

              }

              if(!x_is_bin && !x_is_gene && !x_is_num)
                print_msg("Some gene or bin_x or indexes where not found in the object.", msg_type = "STOP")
            }

            if(!missing(j)){
              if(!all(j %in% bin_y(x)))
                print_msg("Some bin_y where not found in the object.", msg_type = "STOP")
            }

            n_coord <- x@coord
            n_bin_mat <- x@bin_mat
            n_bin_x <- x@bin_x
            n_bin_y <- x@bin_y
            n_ripley_k_function <- x@ripley_k_function
            n_control <- x@control

            if(missing(j)) {

              if (missing(i)) {
                return(x)
              } else {
                if(x_is_gene){

                  n_coord <- n_coord[n_coord$gene %in% i, ]
                  n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", i)]
                  n_ripley_k_function <- n_ripley_k_function[n_ripley_k_function$gene %in% i, ]

                }else if(x_is_bin){

                  n_coord <- n_coord[n_coord$bin_x %in% i, ]
                  n_bin_mat <- n_bin_mat[n_bin_mat$bin_x %in% i, ]
                  gene_left <- unique(n_coord$gene)

                  if("all_genes" %in% colnames(n_bin_mat)){
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", "all_genes", gene_left)]
                  }else{
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", gene_left)]
                  }


                }else{

                  g <- gene_names(x, all_genes=TRUE)[i]
                  n_coord <- n_coord[n_coord$gene %in% g, ]

                  if("all_genes" %in% colnames(n_bin_mat)){
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", "all_genes", g)]
                  }else{
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", g)]
                  }

                  n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", g)]
                  n_ripley_k_function <- n_ripley_k_function[n_ripley_k_function$gene %in% g,]

                }
              }
            } else{
              if (missing(i)) {

                n_coord <- n_coord[n_coord$bin_y %in% j, ]
                n_bin_mat <- n_bin_mat[n_bin_mat$bin_y %in% j, ]

              }else{
                if(x_is_gene){
                  n_coord <- n_coord[n_coord$gene %in% i, ]

                  if("all_genes" %in% colnames(n_bin_mat)){
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", "all_genes", i)]
                  }else{
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", i)]
                  }


                }else if(x_is_bin){

                  n_coord <- n_coord[n_coord$bin_x %in% i, ]
                  n_bin_mat <- n_bin_mat[n_bin_mat$bin_x %in% i, ]
                  gene_left <- unique(n_coord$gene)

                  if("all_genes" %in% colnames(n_bin_mat)){
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", "all_genes", gene_left)]
                  }else{
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", gene_left)]
                  }



                }else{

                  g <- gene_names(x, all_genes=TRUE)[i]
                  n_coord <- n_coord[n_coord$gene %in% g, ]

                  if("all_genes" %in% colnames(n_bin_mat)){
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", "all_genes", g)]
                  }else{
                    n_bin_mat <- n_bin_mat[, colnames(n_bin_mat) %in% c("bin_x", "bin_y", g)]
                  }

                  n_ripley_k_function <- n_ripley_k_function[n_ripley_k_function$gene %in% g,]
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
            STGrid_obj@bin_x <- x@bin_x[x@bin_x %in% n_bin_mat$bin_x]
            STGrid_obj@bin_y <- x@bin_y[x@bin_y %in% n_bin_mat$bin_y]
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
             grep(
               pattern,
               colnames(x@meta),
               value = TRUE,
               fixed = TRUE
             )
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
             if(length(value) != nbin_x(x) * nbin_y(x) && length(value) != 1)
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
#' Estimate Ripley's reduced second moment function for each gene in a spatial grid.
#'
#' This method calculates Ripley's reduced second moment function, K(r), for each gene in a spatial grid.
#'
#' @param object An object of class "STGrid".
#' @param rmax Maximum desired value of the argument r.
#' @param nlarge Efficiency threshold. If the number of points exceeds nlarge, then only the border correction will be computed (by default), using a fast algorithm.
#' @param var.approx Logical. If TRUE, the approximate variance of K(r) under CSR will also be computed.
#' @param ratio Logical. If TRUE, the numerator and denominator of each edge-corrected estimate will also be saved, for use in analyzing replicated point patterns.
#'
#' @return An updated object of class "STGrid" with the Ripley's K function estimates stored in the slot 'ripley_k_function'.
#'
#' @examples
#' # Example usage:
#' # grid_object <- load_spatial(...)
#' # grid_object <- compute_k_ripley(grid_object)
#' @keywords internal
setGeneric("compute_k_ripley",
           function(object,
                    rmax=200,
                    nlarge=1e6,
                    var.approx=FALSE,
                    ratio=FALSE)
             standardGeneric("compute_k_ripley")
)

#' Estimate Ripley's reduced second moment function for each gene in a spatial grid.
#'
#' This method calculates Ripley's reduced second moment function, K(r), for each gene in a spatial grid.
#'
#' @param object An object of class "STGrid".
#' @param rmax Maximum desired value of the argument r.
#' @param nlarge Efficiency threshold. If the number of points exceeds nlarge, then only the border correction will be computed (by default), using a fast algorithm.
#' @param var.approx Logical. If TRUE, the approximate variance of K(r) under CSR will also be computed.
#' @param ratio Logical. If TRUE, the numerator and denominator of each edge-corrected estimate will also be saved, for use in analyzing replicated point patterns.
#'
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
setMethod(
  "compute_k_ripley", signature("STGrid"),
  function(object,
           rmax=200,
           nlarge=1e6,
           var.approx=FALSE,
           ratio=FALSE) {

    molecules <- coord(object)[, c("x", "y", "gene")]

    data_out <- data.frame(r=NA,
                           theo=NA,
                           border=NA,
                           trans=NA,
                           iso=NA,
                           gene=NA)

    x_min <- object@x_min
    x_max <- object@x_max
    y_min <- object@y_min
    y_max <- object@y_max

    rownames(molecules) <- paste0(1:nrow(molecules),
                                  '|',
                                  molecules$gene)


    print_msg(">>> Estimating Ripley's reduced second moment function for all genes.")

    n_iter <- length(unique(molecules$gene))

    # Initializes the progress bar
    pb <- txtProgressBar(min = 0,
                         max = n_iter,
                         style = 3,
                         width = 50,
                         char = "=")   # Character used to create the bar

    n <- 0
    for(g in unique(molecules$gene)){
      n <- n + 1
      setTxtProgressBar(pb, n)

      list_neighbors <- data.frame(Var1 = NA,
                                   Var2 = NA,
                                   value = NA)

      is_g <- molecules[molecules$gene == g, ]

      X <- spatstat.geom::ppp(is_g$x,
                              is_g$y,
                              window=spatstat.geom::owin(c(min(is_g$x),
                                                           max(is_g$x)),
                                                         c(min(is_g$y),
                                                           max(is_g$y))))

      u <- spatstat.explore::Kest(X,
                                  rmax=rmax,
                                  nlarge=nlarge,
                                  var.approx=var.approx,
                                  ratio=ratio)

      u <- as.data.frame(as.matrix(u))
      u$gene <- g

      data_out <- rbind(data_out, u)
    }
    close(pb)

    data_out <- na.omit(data_out)

    object@ripley_k_function <- data_out

    return(object)

  })


# -------------------------------------------------------------------------
##    Constructor for SpatialTranscriptomicAnalysis class (STGrid)
# -------------------------------------------------------------------------

#' Constructor for SpatialTranscriptomicAnalysis class (STGrid)
#'
#' @param path Either a file (if method is set to "coordinates") or a directory (if method
#' is set to "merscope"). If method is set to "coordinates" the file should contain 3 columns
#' ("x", "y", "gene").
#' @param method The type of technology ("merscope" in the only supported method at the moment).
#' @param bin_size Numeric value representing the bin size (default to 25).
#' @param control A regular expression to identify controls. As the function computes the sum of
#' expression levels ("all_genes" features), this will allow to delete these blanks/controls for computation.
#' @return An object of class STGrid.
#'
#' @export load_spatial
load_spatial <- function(path = "",
                         method=c("coordinates",
                                  "merscope"),
                         bin_size=25,
                         control="^Blank\\-[0-9]+") {

  method <- match.arg(method)


  print_msg("Technology is '", method, "'.")
  print_msg("Loading data from file:", path)

  if(method == "merscope"){

    spat_input <- Seurat::ReadVizgen(data.dir = path, type = "centroids")
    spat_input <- spat_input$microns

  }else if(method=="coordinates"){

    check_file(path, mode="read")
    spat_input <- as.data.frame(data.table::fread(path,
                                                  sep="\t",
                                                  head=TRUE))

    if(any(!c("x", "y", "gene") %in% colnames(spat_input))){
      print_msg("Please check the column name of the input file.", msg_type = "STOP")
    }

    spat_input <- spat_input[, c("x", "y", "gene")]

  }

  bin_matrix <- bin_this_matrix(coord=spat_input,
                                bin_size=bin_size,
                                control=control)

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
  STGrid_obj@meta <- data.frame(row.names = rownames(STGrid_obj@bin_mat))
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
#' Create a 2D binned grid (rasterization) from molecule coordinates
#'
#' This function takes molecule coordinates as input (with columns named 'x', 'y', 'gene') and creates a 2D
#' binned grid by dividing the x and y axes into bins. It returns the counts of the number of molecules in
#' each bin along with additional information about the binning process.
#'
#' @param coord A data frame containing molecule coordinates with columns 'x', 'y', and 'gene'. Defaults to NULL.
#' @param bin_size An integer specifying the size of each bin. Defaults to 25.
#' @return A list containing the binned spatial matrix, updated molecule coordinates, and information about the binning process.
#' @examples
#' # Create a binned grid from molecule coordinates
#' bin_mat(coord = my_coord, bin_size = 25)
#' @keywords internal
#' @export bin_this_matrix
bin_this_matrix <- function(coord=NULL,
                            bin_size=25,
                            control=NULL){

  print_msg("Binning a matrix...")

  x_min <- min(coord$x)
  x_max <- max(coord$x)
  y_min <- min(coord$y)
  y_max <- max(coord$y)

  if(x_max - x_min < bin_size){
    print_msg("Can't create bin with x_max =",
              x_max,
              ", x_min =",
              x_min,
              "and bin_size =",
              bin_size,
              msg_type = "STOP")
  }

  if(y_max - y_min < bin_size){
    print_msg("Can't create bin with x_max =",
              y_max,
              ", x_min =",
              y_min,
              "and bin_size =",
              bin_size,
              msg_type = "STOP")
  }

  x_lim <- seq(from = x_min,
               to = ifelse(x_max %% bin_size != 0,
                           x_max + bin_size, x_max),
               by = bin_size)


  y_lim <- seq(from = y_min,
               to = ifelse(y_max %% bin_size != 0,
                           y_max + bin_size, y_max),
               by = bin_size)

  possible_levels_x <- levels(cut(0, breaks = x_lim,
                                  include.lowest = TRUE,
                                  right=FALSE))

  possible_levels_y <- levels(cut(0, breaks = y_lim,
                                  include.lowest = TRUE,
                                  right=FALSE))


  spatial_matrix <- data.frame(bin_x = sort(rep(possible_levels_x,
                                                length(possible_levels_y))),
                               bin_y = rep(possible_levels_y,
                                           length(possible_levels_x)))

  rownames(spatial_matrix) <- paste(as.character(spatial_matrix$bin_x),
                                    as.character(spatial_matrix$bin_y), sep="~")


  coord$bin_x <- NA
  coord$bin_y <- NA


  pb <- txtProgressBar(min = 0,
                       max = length(unique(coord$gene)),
                       style = 3,
                       width = 50,
                       char = "=")

  n_loop <- 0

  for(goi in unique(coord$gene)){

    x_molec <- cut(coord$x[coord$gene == goi],
                   breaks = x_lim, include.lowest = TRUE, right=FALSE)

    coord$bin_x[coord$gene == goi] <- as.character(x_molec)

    y_molec <- cut(coord$y[coord$gene == goi],
                   breaks = y_lim, include.lowest = TRUE, right=FALSE)
    coord$bin_y[coord$gene == goi] <- as.character(y_molec)

    nb_molec <- table(paste(x_molec, y_molec, sep="~"))
    spatial_matrix[, goi] <- 0
    spatial_matrix[names(nb_molec), goi] <- nb_molec

    n_loop <- n_loop + 1
    setTxtProgressBar(pb, n_loop)

  }

  cat("\n")

  print_msg("Compute sum of counts ('all_genes').")


  pos_bin_xy <- grep("^bin_[xy]$",
                     colnames(spatial_matrix))

  spatial_matrix$all_genes <- 0

  tmp <- spatial_matrix[, -pos_bin_xy]

  pos_ctrl <- grep(control, colnames(tmp), perl=TRUE)

  if(length(pos_ctrl) > 0){
    print_msg("Found the following controls:",
              paste0(head(colnames(tmp)[pos_ctrl], 3),
                     collapse = ", "),
              "...")
    all_genes <- rowSums(tmp[, -pos_ctrl])
  }else{
    print_msg("No control found")
    all_genes <- rowSums(tmp)
  }

  spatial_matrix$all_genes <- all_genes

  spatial_matrix <- spatial_matrix[, order(colnames(spatial_matrix))]

  return(list(spatial_matrix=spatial_matrix,
              coord=coord,
              x_max=x_max,
              x_min=x_min,
              y_min=y_min,
              y_max=y_max,
              possible_levels_x=possible_levels_x,
              possible_levels_y=possible_levels_y
  ))
}

# -------------------------------------------------------------------------
#      Re-Bin a matrix
# -------------------------------------------------------------------------
#' @title Re-bin a STGrid object.
#' @description Re-bin a STGrid object.
#' @param object The STGrid object
#' @param bin_size The size of the bin.
#' @export re_bin
#' @keywords internal
setGeneric("re_bin",
           function(object,
                    bin_size)
             standardGeneric("re_bin")
)


#' @title Re-bin a STGrid object.
#' @description Re-bin a STGrid object.
#' @param x The STGrid object.
#' @param bin_size The size of the bin.
#' @export re_bin
setMethod("re_bin",signature(object="STGrid"),
          function(object,
                   bin_size){

            if(object@bin_size == bin_size){
              print_msg("The bin_size is unchanged.")
              return(object)
            }

            bin_matrix <- bin_this_matrix(coord=object@coord,
                                          bin_size=bin_size,
                                          control=object@control)

            # create a STGrid object                             ---------
            print_msg("Creating an STGrid object")

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
            STGrid_obj@bin_size <- bin_size
            STGrid_obj@bin_x <- bin_matrix$possible_levels_x
            STGrid_obj@bin_y <- bin_matrix$possible_levels_y
            STGrid_obj@ripley_k_function <- object@ripley_k_function
            STGrid_obj@control <- object@control

            return(STGrid_obj)

          }

)


# -------------------------------------------------------------------------
#      Return the x/y coordinates of genes from a STGrid object
# -------------------------------------------------------------------------
#' @title The x/y coordinates of genes from a STGrid object.
#' @description Return the x/y coordinates of genes from a STGrid object
#' @param object The STGrid object
#' @param gene_list The list of genes
#' @param as.factor Whether the gene column should be returned as an ordered factor.
#' @export
#' @keywords internal
setGeneric("get_gn_coord",
           function(object,
                    gene_list=character(),
                    as.factor=TRUE)
             standardGeneric("get_gn_coord")
)

#' @title The x/y coordinates of genes from a STGrid object.
#' @description Return the x/y coordinates of genes from a STGrid object
#' @param object The STGrid object
#' @param gene_list The list of genes
#' @param as.factor Whether the gene column should be returned as an ordered factor
#' @export
setMethod("get_gn_coord", "STGrid",
           function(object,
                    gene_list=character(),
                    as.factor=TRUE){

             if(!all(gene_list %in% gene_names(object))){
               print_msg("Some genes were not found...",
                         msg_type = "STOP")
             }

             if(length(gene_list)==0){
               gene_list <- gene_names(object)
             }

             coord <- object@coord
             coord <- coord[coord$gene %in% gene_list, ]

             if(as.factor)
               coord$gene <- factor(as.character(coord$gene),
                                       levels=gene_list,
                                       ordered=TRUE)

             return(coord)

      }
)

