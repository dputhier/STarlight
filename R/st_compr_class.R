# -------------------------------------------------------------------------
## Spatial transcriptomic comparison class (STCompR) definition
# -------------------------------------------------------------------------

#' Spatial transcriptomic comparison class (STCompR).
#'
#' This class represents the result of a spatial transcriptomic comparison between
#' two STGrid objects. It stores information such as conditions, methods, raw counts,
#' scaling factors, and statistical test results.
#'
#' @slot conditions A character vector specifying the conditions being compared.
#' @slot method A character vector specifying the technology used for each condition.
#' @slot raw_counts A list storing raw count data for each condition.
#' @slot scaling_factor A numeric vector representing the scaling factors used during normalization (to be divided).
#' @slot stat_test A data frame containing statistical test results, including log2 ratios,
#'                odds ratios, and p-values for each gene.
#' @slot pseudo_count The value for the pseudo_count.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' st_grid_1 <- create_STGrid(...)
#' st_grid_2 <- create_STGrid(...)
#' st_compr_result <- stcompr(st_grid_1, st_grid_2, name_1 = "Condition1", name_2 = "Condition2")
#' st_compr_object <- new("STCompR", conditions = c("Condition1", "Condition2"),
#'                        method = c("Method1", "Method2"), raw_counts = list(...),
#'                        scaling_factor = c(1.2, 0.8), stat_test = data.frame(...))
#'
setClass("STCompR",
         slots = c(
           conditions = "character",
           method="character",
           raw_counts = "list",
           scaling_factor = "numeric",
           stat_test = "data.frame",
           pseudo_count= "numeric",
           neighborhood="list"
         ),
         prototype = list(
           conditions = "",
           method="",
           raw_counts = list(),
           scaling_factor = numeric(),
           stat_test = data.frame(),
           pseudo_count= 1,
           neighborhood=list()
         )
)

# -------------------------------------------------------------------------
##  Estimate Size Factors using DESeq2 RLE Method
# -------------------------------------------------------------------------
#' @title Estimate Size Factors using DESeq2 RLE Method
#'
#' @description
#'  Estimate size factors for normalization using the Relative Log Expression (RLE) method
#' as in the DESeq2 package.
#'
#' This function computes size factors for normalization based on the DESeq2 RLE method.
#' Size factors are calculated as the median of the ratios of counts to their
#' geometric mean across all samples.
#'
#' @param cds A DESeqDataSet object containing count data.
#'
#' @return A numeric vector of size factors for normalization.
#'
#' @details The DESeq2 RLE method estimates size factors for each sample based on
#' the relative log expression (RLE). The geometric mean of counts for each feature
#' across all samples is computed. Zero counts are set to NA to avoid division by zero.
#' Each count is then divided by its corresponding geometric mean, and the median
#' of these ratios across all samples is computed to obtain the size factor for each sample.
#'
#' @keywords internal
estimSf <- function(cts) {

  geomMean <- function(x) prod(x)^(1/length(x))

  # Compute the geometric mean over the line
  gm.mean  <-  apply(cts, 1, geomMean)

  # Zero values are set to NA (avoid subsequent division by 0)
  gm.mean[gm.mean == 0] <- NA

  # Divide each line by its corresponding geometric mean
  cts <- sweep(cts, 1, gm.mean, FUN="/")

  # Compute the median over the columns
  med <- apply(cts, 2, median, na.rm=TRUE)

  # Force sf to have geometric mean of 1
  med <- med/exp(mean(log(med)))

  # Return the scaling factor
  return(med)
}

# -------------------------------------------------------------------------
##    Constructor for SpatialTranscriptomic comparison class (STCompR)
# -------------------------------------------------------------------------

#' @title Create a STCompR compare class to compare two STGrid objects.
#' @description
#' Create a STCompR compare class to compare two STGrid objects.
#'
#' @param object_1 A STGrid object.
#' @param object_2 A STGrid object.
#' @param name_1 A character string specifying the name for Condition 1.
#' @param name_2 A character string specifying the name for Condition 2.
#' @param stat_test A character string specifying the statistical test to be used.
#'                 Currently, only "Fisher" is supported.
#' @param norm_method A character string specifying the normalization method to be used.
#'                   Currently, only "RLE" is supported.
#' @param pseudo_count A numeric value specifying the pseudo-count to be added during normalization.
#'
#' @return An object of class 'STCompR'.
#'
#' @examples
#' # Example usage:
#' st_grid_1 <- load_spatial(...)
#' st_grid_2 <- load_spatial(...)
#' st_compr_result <- stcompr(st_grid_1, st_grid_2, name_1 = "Condition1", name_2 = "Condition2")
#' @export stcompr
stcompr <- function(object_1,
                    object_2,
                    name_1="Condition_1",
                    name_2="Condition_2",
                    stat_test=c("Fisher"),
                    norm_method=c("RLE"),
                    pseudo_count=1) {

  stat_test <- match.arg(stat_test)
  norm_method <- match.arg(norm_method)

  check_var(name_1)
  check_var(name_2)

  if(class(object_1)[1] != "STGrid" | class(object_2)[1] != "STGrid"){
    print_msg("Please provide STGrid objects.")
  }

  gn_1 <- gene_names(object_1)
  gn_2 <- gene_names(object_1)

  if(!all(gn_1 %in% gn_2) | !all(gn_2 %in% gn_1)){
    print_msg("Objects do not contain the same genes.")
  }

  tb_1 <- data.frame(row.names = names(table(coord(object_1)$gene)),
                     counts=as.numeric(table(coord(object_1)$gene)))

  tb_2 <- data.frame(row.names = names(table(coord(object_2)$gene)),
                     counts=as.numeric(table(coord(object_2)$gene)))

  raw_counts <- cbind(tb_1, tb_2)
  colnames(raw_counts) <- c(name_1, name_2)

  print_msg("Normalizing...")
  scaling_factor <- estimSf(raw_counts)

  norm_counts <- sweep(raw_counts, MARGIN = 2, STATS =  scaling_factor, FUN="/")

  norm_counts <- norm_counts + pseudo_count

  print_msg("Computing log2 ratio...")
  ratio <- norm_counts[, name_2]/norm_counts[, name_1]
  log2_ratio <- log2(ratio)

  p_values <- c()
  odd_ratio <- c()

  print_msg("Computing Fisher tests...")

  for(g in 1:nrow(norm_counts)){

    a <- norm_counts[g, name_2]
    b <- norm_counts[g, name_1]
    c <- sum(norm_counts[-g,name_2])
    d <- sum(norm_counts[-g, name_1])

    m <- matrix(c(round(a, 0), round(b, 0), round(c, 0) , round(d, 0)), byrow = TRUE, nc=2)
    ft <- fisher.test(m)
    p_values[g] <- ft$p.value
    odd_ratio[g] <- ft$estimate["odds ratio"]
  }

  stats <- data.frame(row.names = rownames(norm_counts),
                      log2_ratio=log2_ratio,
                      odd_ratio=odd_ratio,
                      p_values=p_values)

  stats$p_values[stats$p_values < 1e-320] <- 1e-320

  print_msg("Adjusting p-values...")

  for(corr in c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")){
    stats[, paste0("padj_", corr)] <- p.adjust(stats$p_values, method=corr)
  }

  print_msg("Preparing neighborhood analysis...")

  bin_mat_1 <- bin_mat(object_1, all_genes = FALSE, del_bin = TRUE)
  bin_mat_2 <- bin_mat(object_2, all_genes = FALSE, del_bin = TRUE)

  bin_mat_2 <- bin_mat_2[, colnames(bin_mat_1)]

  spatial_matrix_ratio_1 <- matrix(NA,
                                     nrow= ncol(bin_mat_1),
                                     ncol= ncol(bin_mat_1))

  colnames(spatial_matrix_ratio_1) <- colnames(bin_mat_1)
  rownames(spatial_matrix_ratio_1) <- colnames(bin_mat_1)

  spatial_matrix_ratio_2 <- spatial_matrix_ratio_1

  bin_mat_1 <- bin_mat_1[rowSums(bin_mat_1) > 0, ]
  bin_mat_2 <- bin_mat_2[rowSums(bin_mat_2) > 0, ]

  bin_mat_1 <- sweep(bin_mat_1,
                     MARGIN = 2,
                     STATS = colSums(bin_mat_1, na.rm = TRUE), FUN = "/")

  bin_mat_2 <- sweep(bin_mat_2,
                     MARGIN = 2,
                     STATS = colSums(bin_mat_2, na.rm = TRUE), FUN = "/")


  spatial_matrix_ratio_1 <- as.matrix(dist(t(bin_mat_1), method = "manhattan"))
  spatial_matrix_ratio_2 <- as.matrix(dist(t(bin_mat_2), method = "manhattan"))

  print_msg("Preparing an STCompR object... ")

  STCompR <- new("STCompR")
  STCompR@neighborhood <- list(spatial_matrix_ratio_1,
                               spatial_matrix_ratio_2)
  names(STCompR@neighborhood) <- c(name_1, name_2)
  STCompR@conditions <- c(name_1, name_2)
  STCompR@method <- c(object_1@method, object_2@method)
  names(STCompR@method) <- STCompR@conditions
  STCompR@scaling_factor <- scaling_factor
  STCompR@stat_test <- stats
  STCompR@raw_counts <- raw_counts
  STCompR@pseudo_count <- pseudo_count

  return(STCompR)

}

# -------------------------------------------------------------------------
##    Generate a heatmap comparison for an 'STCompR' object
# -------------------------------------------------------------------------
#' @title Generate a heatmap comparison for an 'STCompR' object
#'
#' @description
#' Generate a heatmap to visualize neighborhood changes between conditions.
#' The changes are computed as the difference in manhattan distances between the two conditions
#' stored in the object.
#'
#' @param object An object of class 'STCompR'.
#' @param what A character string specifying what to plot: "changes" (the difference in manhattan distances between conditions 2 and 1),
#' "changes_2" (the difference in manhattan distances between conditions 1 and 2), "c1" (the manhattan distances of conditions 1), or "c2"
#' (the manhattan distances of conditions 2).
#' @param hclust_method A character string specifying the agglomerative criteria for hierarchical clustering ("ward.D", "ward.D2",
#' "single", "complete", "average", "mcquitty", "median", "centroid").
#' @param dist_method A character string specifying the distance method to be used ("euclidean",
#' "maximum", "manhattan", "canberra", "binary","minkowski").
#' @param del_gene A character vector specifying genes to be excluded from the heatmap.
#' @param only_genes A character vector specifying genes to be included in the heatmap.
#' @param filter A numeric value specifying the threshold for filtering out low values in the heatmap.
#' @param size A numeric value specifying the size of text in the heatmap.
#'
#' @return A heatmap plot comparing conditions in the 'STCompR' object.
#'
#' @examples
#' # Example usage:
#' heatmap_cmp(object = my_STCompR_object, hclust_method = "ward.D", dist_method = "euclidean")
#' @keywords internal
#' @export
setGeneric("heatmap_cmp",
           function(object,
                    what=c("changes", "changes_2",  "c1", "c2"),
                    hclust_method=c("ward.D", "ward.D2",
                                    "single", "complete", "average",
                                    "mcquitty", "median", "centroid"),
                    dist_method=c("euclidean",
                                  "maximum",
                                  "manhattan",
                                  "canberra",
                                  "binary",
                                  "minkowski"),
                    del_gene=NULL,
                    only_genes=NULL,
                    filter=0.2,
                    size=6)
             standardGeneric("heatmap_cmp")
)


#' @title Generate a heatmap comparison for an 'STCompR' object
#'
#' @description
#' Generate a heatmap to visualize neighborhood changes between conditions.
#' The changes are computed as the difference in manhattan distances between the two conditions
#' stored in the object.
#'
#' @param object An object of class 'STCompR'.
#' @param what A character string specifying what to plot: "changes" (the difference in manhattan distances between conditions 2 and 1),
#' "changes_2" (the difference in manhattan distances between conditions 1 and 2), "c1" (the manhattan distances of conditions 1), or "c2"
#' (the manhattan distances of conditions 2).
#' @param hclust_method A character string specifying the agglomerative criteria for hierarchical clustering ("ward.D", "ward.D2",
#' "single", "complete", "average", "mcquitty", "median", "centroid").
#' @param dist_method A character string specifying the distance method to be used ("euclidean",
#' "maximum", "manhattan", "canberra", "binary","minkowski").
#' @param del_gene A character vector specifying genes to be excluded from the heatmap.
#' @param only_genes A character vector specifying genes to be included in the heatmap.
#' @param filter A numeric value specifying the threshold for filtering out low values in the heatmap.
#' @param size A numeric value specifying the size of text in the heatmap.
#'
#' @return A heatmap plot comparing conditions in the 'STCompR' object.
#'
#' @examples
#' # Example usage:
#' heatmap_cmp(object = my_STCompR_object, hclust_method = "ward.D", dist_method = "euclidean")
#' @export
setMethod(
  "heatmap_cmp", signature("STCompR"),
  function(object,
           what=c("changes", "changes_2",  "c1", "c2"),
           hclust_method=c("ward.D", "ward.D2",
                           "single", "complete", "average",
                           "mcquitty", "median", "centroid"),
           dist_method=c("euclidean",
                         "maximum",
                         "manhattan",
                         "canberra",
                         "binary",
                         "minkowski"),
           del_gene=NULL,
           only_genes=NULL,
           filter=0.2,
           size=6) {

    what <- match.arg(what)
    dist_method <- match.arg(dist_method)
    hclust_method <- match.arg(hclust_method)
  if(what=="changes"){
    print_msg("Parameter 'what' is set to 'changes'...")
    obj <- object@neighborhood[[2]] - object@neighborhood[[1]]
    legend_name <- "Changes"
  } else if(what == "changes_2"){
    print_msg("Parameter 'what' is set to 'changes_2'...")
    obj <- object@neighborhood[[1]] - object@neighborhood[[2]]
    legend_name <- "Changes"
  }else if(what=="c1"){
    print_msg("Parameter 'what' is set to 'c1'...")
    obj <- object@neighborhood[[1]]
    legend_name <- "Dist"
  }else {
    print_msg("Parameter 'what' is set to 'c2'...")
    obj <- object@neighborhood[[2]]
    legend_name <- "Dist"
  }

  if(!is.null(del_gene)){
    del_gene <- unique(del_gene)
    if(all(del_gene %in% gene_names(object))){
      obj <- obj[!rownames(obj) %in% del_gene,
                 !colnames(obj) %in% del_gene]
    }else{
      print_msg("Some gene to delete were not found in the object", msg_type = "STOP")
    }
  }

    if(!is.null(only_genes)){
      only_genes <- unique(only_genes)
      if(all(del_gene %in% gene_names(object))){
        obj <- obj[rownames(obj) %in% only_genes,
                   colnames(obj) %in% only_genes]
      }else{
        print_msg("Some gene to delete were not found in the object", msg_type = "STOP")
      }
    }

  if(ncol(obj) == 0){
    print_msg("No more gene left", msg_type = "STOP")
  }

  if(!is.null(filter)){
    TF <- abs(obj) > filter
    obj <- obj[rowSums(TF) > 0 , colSums(TF) > 0]
  }


  print_msg("Calling ggheatmap...")
  p <- ggheatmap(obj,
                 dist_method=dist_method,
                 hclust_method = hclust_method,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 legendName=legend_name,
                 scale = "none",
                 color=RColorBrewer::brewer.pal(7, "Spectral"),
                 border = "darkgrey",
                 tree_color_rows = "darkgrey",
                 tree_color_cols = "darkgrey",
                 annotation_rows = NULL,
                 annotation_cols = NULL,
                 annotation_color = NULL
  ) %>% ggheatmap_theme(1,
                        theme=list(theme(axis.text.x = element_text(angle = 90,
                                                                   size = size, hjust=1),
                                        axis.text.y = element_text(size=size))))
  print(p)

}
)

# -------------------------------------------------------------------------
##    REDEFINE SHOW() METHOD FOR CLASS OBJECT : STCompR
# -------------------------------------------------------------------------

#' @title The show method of a STCompR object
#' @description
#' The show method of a STCompR object displays information about the object,
#' including conditions, memory usage, and the number of genes.
#'
#' @param object A STCompR object.
#' @export show
#' @keywords internal
#'
#' @examples
#' # Example usage:
#' show(st_compr_result)
#'
setMethod(
  "show", signature("STCompR"),
  function(object) {
    print_msg("An object of class STCompR")
    for(i in 1:length(object@conditions)){
      print_msg("Condition", i, object@conditions[i])
    }
    print_msg("Memory used: ", object.size(object))
    print_msg("Number of genes: ", nrow(object@raw_counts))
    print_msg(">>> Please, use show_methods(class = 'STCompR') to show availables methods <<<")
  }
)


# -------------------------------------------------------------------------
##    Some basic functions
# -------------------------------------------------------------------------

if(!isGeneric("nb_genes")){
  setGeneric("nb_genes",
            function(object)
              standardGeneric("nb_genes")
  )
}

#' @title The number of genes stored in a STCompR object
#' @description
#' The number of genes stored in a STCompR object
#' @param object The STCompR object
#' @export nb_genes
setMethod("nb_genes", signature(object = "STCompR"),
           function(object)
             nrow(object@stat_test)
)


if(!isGeneric("gene_names")){
  setGeneric("gene_names",
             function(object)
               standardGeneric("gene_names")
  )
}
#' @title The genes stored in a STCompR object
#' @description
#' The genes stored in a STCompR object
#' @param x The STCompR object
#' @export
setMethod("gene_names", signature(object="STCompR"),
           function(object)
             rownames(object@stat_test)
)

# -------------------------------------------------------------------------
##    The normalized counts and test results of a STCompR object.
# -------------------------------------------------------------------------

#' @title The normalized counts and test results of a STCompR object.
#' @description
#' The normalized counts and test results of a STCompR object.
#' @param x The STCompR object
#' @param transform Whether the count are transformed (the pseudo count defined for the object is added).
#' @param count_only Returns only counts.
#' @param melted_count Returns a melted table with counts.
#' @param normalized Counts are normalized.
#' @param genes The genes that should be returned. Default to all (NULL).
#' @keywords internal
setGeneric("stat_test",
           function(object,
                    transform=c("None", "log2", "log10", "log"),
                    count_only=FALSE,
                    melted_count=FALSE,
                    normalized=TRUE,
                    genes=NULL)
             standardGeneric("stat_test")
)

#' @title The normalized counts and test results of a STCompR object.
#' @description
#' The normalized counts and test results of a STCompR object.
#' @param x The STCompR object
#' @param transform Whether the count are transformed (the pseudo count defined for the object is added).
#' @param count_only Returns only counts.
#' @param melted_count Returns a melted table with counts.
#' @param normalized Counts are normalized.
#' @param genes The genes that should be returned. Default to all (NULL).
#' @export stat_test
setMethod("stat_test", "STCompR",
           function(object,
                    transform=c("None", "log2", "log10", "log"),
                    count_only=FALSE,
                    melted_count=FALSE,
                    normalized=TRUE,
                    genes=NULL){

             transform <- match.arg(transform)

             counts <- object@raw_counts

             if(!is.null(genes)){
               if(!all(genes %in% gene_names(object))){
                 print_msg("Some genes where not found.", msg_type = "STOP")
               }
             }else{
               genes <- gene_names(object)
             }

             if(normalized)
               counts <- sweep(counts, MARGIN = 2, STATS =  object@scaling_factor, FUN="/")

              if(transform == "log2"){
               counts <- log2(counts + object@pseudo_count)
             }else if(transform == "log10"){
               counts <- log10(counts + object@pseudo_count)
             }else if(transform == "log"){
               counts <- log(counts + object@pseudo_count)
             }

             if(count_only){
               d <- counts
               if(melted_count){
                 d <- reshape2::melt(as.matrix(d), id.vars=integer())
                 colnames(d) <- c("Genes", "Conditions", "Counts")
               }
             }else{
               d <- cbind(counts,
                          object@stat_test)
             }

             if(melted_count){
               return(d[d$Genes %in% genes, ])
             }else{
               return(d[genes,])
             }

}
)

