# -------------------------------------------------------------------------
## Spatial transcriptomic anlysis class (STGrid) definition
# -------------------------------------------------------------------------

#' Spatial transcriptomic anlysis class (STGrid).
#'
#' A class for representing spatial transcriptomic analysis results using (currently) the "merscope" method.
#' This class store transcript molecule coordinates as a grid.
#'
#' @slot coordinate A list of data frames representing spatial coordinates (empty by default).
#' May accept several layers in the future.
#' @slot bin_mat A list of data frames representing binned matrices (empty by default).
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
           coord = "list",
           bin_mat = "list",
           y_max = "numeric",
           y_min = "numeric",
           x_min = "numeric",
           x_max = "numeric",
           path = "character",
           method="character",
           meta = "list",
           bin_size = "numeric",
           bin_x = "character",
           bin_y = "character"
         ),
         prototype = list(
           coord = list(),
           bin_mat = list(),
           y_max = 0,
           y_min = 0,
           x_min = 0,
           x_max = 0,
           path = character(),
           method=character(),
           meta = list(),
           bin_size = 0,
           bin_x = character(),
           bin_y = character()
         )
)

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
#'
#' @return An object of class STGrid.
#'
#' @export load_spatial
load_spatial <- function(path = "",
                         method=c("coordinates", "merscope"),
                         bin_size=25) {

  method <- match.arg(method)


  print_msg("Technology is '", method, "'.")
  print_msg("Loading data from file", path)

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


  print_msg("Binning...")

  x_min <- min(spat_input$x)
  x_max <- max(spat_input$x)
  y_min <- min(spat_input$y)
  y_max <- max(spat_input$y)

  if(x_max - x_min < bin_size){
    print_msg("Can't create bin with x_max =",
              x_max,
              ", x_min =",
              x_min,
              "and bin_size =",
              bin_size)
  }

  if(y_max - y_min < bin_size){
    print_msg("Can't create bin with x_max =",
              y_max,
              ", x_min =",
              y_min,
              "and bin_size =",
              bin_size)
  }

  x_lim <- seq(from = x_min,
               to = x_max,
               by = bin_size)

  y_lim <- seq(from = y_min,
               to = y_max,
               by = bin_size)

  possible_levels_x <- levels(cut(0, breaks = x_lim, include.lowest = TRUE, right=FALSE))
  possible_levels_y <- levels(cut(0, breaks = y_lim, include.lowest = TRUE, right=FALSE))

  spatial_matrix <- data.frame(bin_x = sort(rep(possible_levels_x,
                                                length(possible_levels_y))),
                               bin_y = rep(possible_levels_y,
                                           length(possible_levels_x)))

  rownames(spatial_matrix) <- paste(as.character(spatial_matrix$bin_x),
                                    as.character(spatial_matrix$bin_y), sep="~")


  spat_input$bin_x <- NA
  spat_input$bin_y <- NA

  for(goi in unique(spat_input$gene)){

    x_molec <- cut(spat_input$x[spat_input$gene == goi],
                   breaks = x_lim, include.lowest = TRUE, right=FALSE)
    spat_input$bin_x[spat_input$gene == goi] <- as.character(x_molec)

    y_molec <- cut(spat_input$y[spat_input$gene == goi],
                   breaks = y_lim, include.lowest = TRUE, right=FALSE)
    spat_input$bin_y[spat_input$gene == goi] <- as.character(y_molec)

    nb_molec <- table(paste(x_molec, y_molec, sep="~"))
    spatial_matrix[names(nb_molec), goi] <- nb_molec

  }


  spatial_matrix <- spatial_matrix[, order(colnames(spatial_matrix))]


  # create a STGrid object                             ---------

  STGrid_obj <- new("STGrid",
                    path = path)

  STGrid_obj@coord <- list("layer_1"=spat_input)
  STGrid_obj@bin_mat <- list("layer_1"=spatial_matrix)
  STGrid_obj@y_max <- y_max
  STGrid_obj@y_min <- y_min
  STGrid_obj@x_max <- x_max
  STGrid_obj@x_min <- x_min
  STGrid_obj@path <- path
  STGrid_obj@method <- method
  STGrid_obj@meta <- list()
  STGrid_obj@bin_size <- bin_size
  STGrid_obj@bin_x <- possible_levels_x
  STGrid_obj@bin_y <- possible_levels_y

  return(STGrid_obj)
}

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
            dim(x@bin_mat[[1]])
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
#' row_names
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
#' row_names
#' @description
#' The row names of a STGrid object.
#' @param x The STGrid object
#' @export row_names
setMethod("row_names", "STGrid",
          function(x)
            x@bin_y
)

#' @title X bins of a STGrid object.
#' row_names
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_x
#' @keywords internal
setGeneric("bin_x",
           function(x)
             standardGeneric("bin_x")
)

#' @title X bins of a STGrid object.
#' row_names
#' @description
#' X bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_x
setMethod("bin_x", "STGrid",
          function(x)
            x@bin_x
)

#' @title Y bins of a STGrid object.
#' row_names
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_y
#' @keywords internal
setGeneric("bin_y",
           function(x)
             standardGeneric("bin_y")
)

#' @title Y bins of a STGrid object.
#' row_names
#' @description
#' Y bins of a STGrid object.
#' @param x The STGrid object
#' @export bin_y
setMethod("bin_y", "STGrid",
          function(x)
            x@bin_y
)

#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param x The STGrid object
#' @export coord
#' @keywords internal
setGeneric("coord",
           function(x)
             standardGeneric("coord")
)

#' @title The coordinates stored in a STGrid object
#' @description
#' The coordinates stored in a STGrid object
#' @param x The STGrid object
#' @export coord
setMethod("coord", "STGrid",
          function(x)
            x@coord$layer_1
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
            nrow(x@coord$layer_1)
)

#' @title The number of genes stored in a STGrid object
#' @description
#' The number of genes stored in a STGrid object
#' @param x The STGrid object
#' @export nb_genes
#' @keywords internal
setGeneric("nb_genes",
           function(x)
             standardGeneric("nb_genes")
)

#' @title The number of genes stored in a STGrid object
#' @description
#' The number of genes stored in a STGrid object
#' @param x The STGrid object
#' @export nb_genes
setGeneric("nb_genes",
           function(x)
             length(unique(x@coord$layer_1$gene))
)

#' @title The genes stored in a STGrid object
#' @description
#' The genes stored in a STGrid object
#' @param x The STGrid object
#' @export gene_names
#' @keywords internal
setGeneric("gene_names",
           function(x)
             standardGeneric("gene_names")
)

#' @title The genes stored in a STGrid object
#' @description
#' The genes stored in a STGrid object
#' @param x The STGrid object
#' @export gene_names
setGeneric("gene_names",
           function(x)
             unique(x@coord$layer_1$gene)
)

setMethod("nb_molec", "STGrid",
          function(x)
            nrow(x@coord$layer_1)
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
#' @keywords internal
setMethod(
  "summary", signature("STGrid"),
  function(object) {
    print_msg("An object of class STGrid")
    print_msg("Method:", object@method)
    print_msg("Bin size:", object@bin_size)
    print_msg("Coord:")
    lapply(object@coord, function(x) print(head(x, 2)))
    print_msg("bin_mat:")
    lapply(object@bin_mat, function(x) print(head(x, 2)))
    print_msg("x_min:", object@x_min)
    print_msg("x_max:", object@x_max)
    print_msg("y_min:", object@y_min)
    print_msg("y_max:", object@y_max)
    print_msg("path:", object@path)
    print_msg("bin_x:", head(bin_x(object)), "...")
    print_msg("bin_y:", head(bin_y(object)), "...")
    print_msg("Meta data:", object@metadata)
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

            i <- unique(i)
            x_is_gene <- FALSE
            x_is_bin <- FALSE

            if(!missing(i)){

              if(any(i %in% bin_x(x))){
                x_is_bin <- TRUE
              }else{
                x_is_bin <- FALSE
                if(any(i %in% gene_names(x))){
                  x_is_gene <- TRUE
                }else{
                  x_is_gene <- FALSE
                }

              }

              if(!x_is_bin && !x_is_gene)
                print_msg("Some gene or bin_x where not found in the object.", msg_type = "STOP")
            }

            if(!missing(j)){
              if(!all(j %in% bin_y(x)))
                print_msg("Some bin_y where not found in the object.", msg_type = "STOP")
            }

            n_coord <- x@coord
            n_bin_mat <- x@bin_mat
            n_bin_x <- x@bin_x
            n_bin_y <- x@bin_y

            if(missing(j)) {

              if (missing(i)) {
                return(x)
              } else {
                if(x_is_gene){
                  n_coord <- lapply(n_coord,
                                    function(x, y) x[x$gene %in% y, ], i)
                  n_bin_mat <- lapply(n_bin_mat,
                                      function(x, y) x[, colnames(x$gene) %in% c("bin_x", "bin_y", y)], i)
                }else{

                  n_coord <- lapply(n_coord,
                                    function(x, y) x[x$bin_x %in% y, ], i)

                  n_bin_mat <- lapply(n_bin_mat,
                                      function(x, y) x[x$bin_x %in% y, ], i)
                }
              }
            } else{
              if (missing(i)) {
                n_coord <- lapply(n_coord,
                                  function(x, y) x[x$bin_y %in% y, ], j)
                n_bin_mat <- lapply(n_bin_mat,
                                    function(x, y) x[x$bin_y %in% y, ], j)
              }else{
                if(x_is_gene){
                  n_coord <- lapply(n_coord,
                                    function(x, y) x[x$gene %in% y, ], i)
                  n_bin_mat <- lapply(n_bin_mat,
                                      function(x, y) x[, colnames(x$gene) %in% c("bin_x", "bin_y", y)], i)
                }else{
                  n_coord <- lapply(n_coord,
                                    function(x, y) x[x$bin_x %in% y, ], i)
                  n_bin_mat <- lapply(n_bin_mat,
                                      function(x, y) x[x$bin_x %in% y, ], i)
                }

                n_coord <- lapply(n_coord,
                                  function(x, y) x[x$bin_y %in% y, ], j)
                n_bin_mat <- lapply(n_bin_mat,
                                    function(x, y) x[x$bin_y %in% y, ], j)
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
            STGrid_obj@bin_x <- x@bin_x[x@bin_x %in% n_bin_mat$layer_1$bin_x]
            STGrid_obj@bin_y <- x@bin_y[x@bin_y %in% n_bin_mat$layer_1$bin_y]

            return(STGrid_obj)

})




