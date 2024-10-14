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
#' @slot raw_counts A data.frame storing raw count data for each condition.
#' @slot norm_counts A data.frame storing normalized counts for each condition.
#' @slot scaling_factor A numeric vector representing the scaling factors used during normalization (to be divided).
#' @slot stat_test A data frame containing statistical test results, including log2 ratios,
#'                odds ratios, and p-values for each feature.
#' @slot pseudo_count The value for the pseudo_count.
#' @slot neighborhood a list of feature-feature neighborhood for the two conditions to be compared (not functional at the moment).
#' @slot neighborhood_changes A feature-feature matrix indicative of neighborhood changes (not functional at the moment).
#'
#' @export
#'
#' @examples
#' # Example usage:
#' example_dataset("11284296/files/cmp_xen")
#' cmp_xen
#'
setClass("STCompR",
         slots = c(
           conditions = "character",
           method="character",
           raw_counts = "data.frame",
           norm_counts = "data.frame",
           scaling_factor = "numeric",
           stat_test = "data.frame",
           pseudo_count= "numeric",
           neighborhood="list",
           neighborhood_changes="data.frame"
         ),
         prototype = list(
           conditions = "",
           method="",
           raw_counts = data.frame(),
           norm_counts = data.frame(),
           scaling_factor = numeric(),
           stat_test = data.frame(),
           pseudo_count= 1,
           neighborhood=list(),
           neighborhood_changes=data.frame()
         )
)


# -------------------------------------------------------------------------
##    Constructor for SpatialTranscriptomic comparison class (STCompR)
# -------------------------------------------------------------------------

#' @title Create a STCompR compare class to compare two sets of STGrid objects.
#' @description
#' Create a STCompR compare class to compare two sets of STGrid objects. By default differential
#' expression is performed with DESeq2. If only one sample of each classe is provided, a Fisher's exact
#' test is performed.
#'
#' @param object_1 A STGrid object.
#' @param object_2 A STGrid object.
#' @param name_1 A string to name classe 1.
#' @param name_2 A string to name classe 2.
#' @param fit_type See  estimateDispersions() in DESeq2 library.
#' @param p_adj_method The pvalue correction method. See stats::p.adjust.
#' @param pseudo_count A pseudo-count value to be stored in the object.
#' @param cooksCutoff See DESeq2::results().
#' @return An object of class 'STCompR'.
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts estimateDispersions nbinomWaldTest results sizeFactors
#' @importFrom stats p.adjust
#' @examples
#' # Example usage:
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r2 <- xen[x_bins, y_bins]
#' cmp <- stcompr(list(xen_r1=xen_r1), list(xen_r2=xen_r2))
#'
#' @export stcompr
stcompr <- function(object_1,
                    object_2,
                    name_1="Condition_1",
                    name_2="Condition_2",
                    fit_type=c("local", "mean",
                               "parametric",  "glmGamPoi"),
                    p_adj_method=c("BH", "holm",
                                   "hochberg", "hommel",
                                   "bonferroni", "BY",
                                   "fdr", "none"),
                    pseudo_count=1,
                    cooksCutoff=FALSE) {

  fit_type <- match.arg(fit_type)
  p_adj_method <- match.arg(p_adj_method)

  print_this_msg("Checking objects.")

  if(is.list(object_1) & is.list(object_2)){

    print_this_msg("Object were provided as lists", msg_type = "DEBUG")

    check_st_list(object_1)
    check_st_list(object_2)

  }else if(inherits(object_1, "STGrid") & inherits(object_2, "STGrid")){

    print_this_msg("Objects were provided as STGrid")
    print_this_msg("Converting to lists.")

    object_1 <- list(object_1)
    object_2 <- list(object_2)

    check_st_list(object_1)
    check_st_list(object_2)

  }else{
    print_this_msg("Please provide two STGrid objects or two lists of STGrid.", msg_type = "STOP")
  }

  two_sample_comparison <- FALSE
  if(length(object_1) == 1 & length(object_2) == 1)
    two_sample_comparison <- TRUE

  print_this_msg("Retrieving feature names.", msg_type = "DEBUG")

  gn_1 <- lapply(object_1, feat_names)
  gn_2 <- lapply(object_2, feat_names)

  gn_1 <- Reduce(intersect, gn_1)
  gn_2 <- Reduce(intersect, gn_2)

  print_this_msg("Checking feature names.")

  if(length(gn_1) == 0){
    print_this_msg("No common genes found in condition 1.", msg_type = "STOP")
  }

  if(length(gn_2) == 0){
    print_this_msg("No common genes found in condition 2.", msg_type = "STOP")
  }

  inter_gn <- intersect(gn_1, gn_2)

  if(length(inter_gn) == 0){
    print_this_msg("No shared features between conditions")
  }

  print_this_msg("Retrieving counts.")

  tb_1 <- lapply(object_1, table_st)

  tb_2 <- lapply(object_2, table_st)

  print_this_msg("Ordering rows", msg_type = "DEBUG")

  for(i in 1:length(tb_1)){
    tb_1[[i]] <- tb_1[[i]][inter_gn]
  }

  for(i in 1:length(tb_2)){

    tb_2[[i]] <- tb_2[[i]][inter_gn]
  }

  print_this_msg("Merging counts.")

  tb_1 <- do.call("cbind", tb_1)
  tb_2 <- do.call("cbind", tb_2)

  if(is.null(names(object_1)) | any(names(object_1) == "")){
    print_this_msg("Naming samples.", msg_type = "DEBUG")
    colnames(tb_1) <- paste0("Condition_1_smp_", 1:ncol(tb_1))
  }else{
    colnames(tb_1) <- names(object_1)
  }

  if(is.null(names(object_2)) | any(names(object_2) == "")){
    colnames(tb_2) <- paste0("Condition_2_smp_", 1:ncol(tb_2))
  }else{
    colnames(tb_2) <- names(object_2)
  }

  raw_counts <- cbind(tb_1, tb_2)
  rownames(raw_counts) <- inter_gn

  exp_design <- data.frame(condition=as.factor(c(rep("condition_1", ncol(tb_1)),
                                       rep("condition_2", ncol(tb_2)))),
                             row.names = colnames(raw_counts))

  print_this_msg("Normalizing counts...")

  dds0 <- DESeq2::DESeqDataSetFromMatrix(countData = raw_counts,
                                         colData = exp_design,
                                         design = ~condition)

  dds.norm <-  DESeq2::estimateSizeFactors(dds0)
  norm_counts <- DESeq2::counts(dds.norm, normalized=TRUE)

  if(!two_sample_comparison){
    print_this_msg("Calling DESeq2")
    dds.disp <- DESeq2::estimateDispersions(dds.norm, fitType=fit_type)
    wald.test <- DESeq2::nbinomWaldTest(dds.disp)
    res_deseq2 <- DESeq2::results(wald.test,
                                  pAdjustMethod="BH", cooksCutoff=cooksCutoff)
    res_deseq2 <- res_deseq2[rownames(raw_counts), ]
    stats <- data.frame(row.names = rownames(raw_counts),
                        log2_ratio=res_deseq2$log2FoldChange,
                        p_values=res_deseq2$pvalue,
                        padj=res_deseq2$padj)

  }else{

    print_this_msg("Only two samples. Performing Fisher's exact test.")

    stats <- data.frame(row.names = rownames(raw_counts),
                        log2_ratio=rep(NA, length(rownames(raw_counts))),
                        p_values=rep(NA, length(rownames(raw_counts))),
                        padj=rep(NA, length(rownames(raw_counts))))


    for(g in 1:nrow(norm_counts)){

      a <- norm_counts[g, 2]
      b <- norm_counts[g, 1]
      c <- sum(norm_counts[-g,2])
      d <- sum(norm_counts[-g, 1])

      m <- matrix(c(round(a, 0),
                    round(b, 0),
                    round(c, 0) ,
                    round(d, 0)),
                  byrow = TRUE, ncol=2)

      ft <- stats::fisher.test(m)
      stats$p_values[g] <- ft$p.value
    }

    ratio <- (norm_counts[, 2] + pseudo_count) / (norm_counts[, 1] + + pseudo_count)
    stats$log2_ratio <- log2(ratio)
    stats$p_values[stats$p_values < 1e-320] <- 1e-320
    stats$padj <- stats::p.adjust(stats$p_values, method=p_adj_method)
    stats$padj[stats$padj < 1e-320] <- 1e-320
  }

  # print_this_msg("Preparing neighborhood analysis...")
  #
  # bin_mat_1 <- bin_mat(object_1, del_bin = TRUE)
  # bin_mat_2 <- bin_mat(object_2, del_bin = TRUE)
  # bin_mat_1 <- bin_mat_1[, inter_gn]
  # bin_mat_2 <- bin_mat_2[, inter_gn]
  #
  # bin_mat_1 <- bin_mat_1[rowSums(bin_mat_1) > 0, ]
  # bin_mat_2 <- bin_mat_2[rowSums(bin_mat_2) > 0, ]
  #
  # bin_mat_1 <- sweep(bin_mat_1,
  #                    MARGIN = 2,
  #                    STATS = colSums(bin_mat_1, na.rm = TRUE), FUN = "/")
  #
  # bin_mat_2 <- sweep(bin_mat_2,
  #                    MARGIN = 2,
  #                    STATS = colSums(bin_mat_2, na.rm = TRUE), FUN = "/")
  #
  # print_this_msg("Computing distances (matrix 1)...")
  # spatial_matrix_ratio_1 <- as.matrix(stats::dist(t(bin_mat_1),
  #                                          method = "manhattan"))
  #
  # print_this_msg("Computing distances (matrix 2)...")
  # spatial_matrix_ratio_2 <- as.matrix(stats::dist(t(bin_mat_2),
  #                                          method = "manhattan"))

  # print_this_msg("Computing contrasts... ")
  # neighborhood_changes <- as.data.frame(spatial_matrix_ratio_2 - spatial_matrix_ratio_1)

  print_this_msg("Preparing an STCompR object... ")

  STCompR <- methods::new("STCompR")
  # STCompR@neighborhood <- list(spatial_matrix_ratio_1,
  #                              spatial_matrix_ratio_2)
  #
  # STCompR@neighborhood_changes <- neighborhood_changes
  # names(STCompR@neighborhood) <- c(name_1, name_2)
  STCompR@conditions <- c(rep(name_1, length(object_1)),
                          rep(name_2, length(object_2)))
  names(STCompR@conditions) <- colnames(raw_counts)
  STCompR@method <- c(unlist(lapply(object_1, function(x) x@method)),
                      unlist(lapply(object_2, function(x) x@method)))
  names(STCompR@method) <- colnames(raw_counts)
  STCompR@scaling_factor <- DESeq2::sizeFactors(dds.norm)
  STCompR@stat_test <- stats
  STCompR@raw_counts <- as.data.frame(raw_counts)
  STCompR@norm_counts <- as.data.frame(norm_counts)
  STCompR@pseudo_count <- pseudo_count

  return(STCompR)

}


# -------------------------------------------------------------------------
##    Generate a heatmap comparison for an 'STCompR' object
# -------------------------------------------------------------------------
#### @title Generate a heatmap comparison for an 'STCompR' object
####
#### @description
#### Generate a heatmap to visualize neighborhood changes between conditions.
#### The changes are computed as the difference in manhattan distances between the two conditions
#### stored in the object.
####
#### @param object An object of class 'STCompR'.
#### @param what A character string specifying what to plot: "changes" (the difference in pairwise Manhattan distances between conditions 2 and 1),
#### "changes_2" (the difference in pairwise Manhattan distances between conditions 1 and 2), "c1" (the pairwise Manhattan distances of features in condition 1), or "c2"
#### (the pairwise Manhattan distances of features in condition 2).
#### @param hclust_method A character string specifying the agglomerative criteria for hierarchical clustering ("ward.D", "ward.D2",
#### "single", "complete", "average", "mcquitty", "median", "centroid").
#### @param dist_method A character string specifying the distance method to be used ("euclidean",
#### "maximum", "manhattan", "canberra", "binary","minkowski").
#### @param del_feat A character vector specifying features to be excluded from the heatmap.
#### @param only_feat A character vector specifying features to be included in the heatmap.
#### @param filter_method If 'cv', select row/col by checking whether the variation coefficient of absolute changes is greater than 'filter'.
#### If 'diff' select a row/col if at least one absolute value greater than 'filter' is observed.
#### @param filter A numeric value as threshold.
#### @param size A numeric value specifying the size of text in the heatmap.
#### @param title A title for the diagram.
####
#### @return A heatmap plot comparing conditions in the 'STCompR' object.
####
#### @examples
#### # Example usage:
#### example_dataset("11284296/files/cmp_xen")
#### heatmap_cmp(object = cmp_xen, hclust_method = "ward.D", dist_method = "euclidean")
#### @keywords internal
#### @export
# setGeneric("heatmap_cmp",
#            function(object,
#                     what=c("changes", "changes_2",  "c1", "c2"),
#                     hclust_method=c("ward.D", "ward.D2",
#                                     "single", "complete", "average",
#                                     "mcquitty", "median", "centroid"),
#                     dist_method=c("euclidean",
#                                   "maximum",
#                                   "manhattan",
#                                   "canberra",
#                                   "binary",
#                                   "minkowski"),
#                     del_feat=NULL,
#                     only_feat=NULL,
#                     filter_method=c("cv", "diff"),
#                     filter=0.5,
#                     size=6,
#                     title=NULL)
#              standardGeneric("heatmap_cmp")
# )


#### @title Generate a heatmap comparison for an 'STCompR' object
####
#### @description
#### Generate a heatmap to visualize neighborhood changes between conditions.
#### The changes are computed as the difference in manhattan distances between the two conditions
#### stored in the object.
####
#### @param object An object of class 'STCompR'.
#### @param what A character string specifying what to plot: "changes" (the difference in pairwise Manhattan distances between conditions 2 and 1),
#### "changes_2" (the difference in pairwise Manhattan distances between conditions 1 and 2), "c1" (the pairwise Manhattan distances of features in condition 1), or "c2"
#### (the pairwise Manhattan distances of features in condition 2).
#### @param hclust_method A character string specifying the agglomerative criteria for hierarchical clustering ("ward.D", "ward.D2",
#### "single", "complete", "average", "mcquitty", "median", "centroid").
#### @param dist_method A character string specifying the distance method to be used ("euclidean",
#### "maximum", "manhattan", "canberra", "binary","minkowski").
#### @param del_feat A character vector specifying feature to be excluded from the heatmap.
#### @param only_feat A character vector specifying features to be included in the heatmap.
#### @param filter A numeric value specifying the threshold for filtering out low absolute values in the heatmap.
#### @param size A numeric value specifying the size of text in the heatmap.
#### @param title A title for the diagram.
####
#### @return A heatmap plot comparing conditions in the 'STCompR' object.
####
#### @examples
#### # Example usage:
#### example_dataset("11284296/files/cmp_xen")
#### heatmap_cmp(object = cmp_xen, hclust_method = "ward.D", dist_method = "euclidean")
#### @importFrom ggheatmap ggheatmap ggheatmap_theme
#### @importFrom RColorBrewer brewer.pal
#### @import magrittr
#### @export
# setMethod(
#   "heatmap_cmp", signature("STCompR"),
#   function(object,
#            what=c("changes", "changes_2",  "c1", "c2"),
#            hclust_method=c("ward.D", "ward.D2",
#                            "single", "complete", "average",
#                            "mcquitty", "median", "centroid"),
#            dist_method=c("euclidean",
#                          "maximum",
#                          "manhattan",
#                          "canberra",
#                          "binary",
#                          "minkowski"),
#            del_feat=NULL,
#            only_feat=NULL,
#            filter_method=c("diff", "cv"),
#            filter=0.2,
#            size=6,
#            title=NULL) {
#
#     what <- match.arg(what)
#     dist_method <- match.arg(dist_method)
#     hclust_method <- match.arg(hclust_method)
#     filter_method <- match.arg(filter_method)
#
#   if(what=="changes"){
#
#     print_this_msg("Parameter 'what' is set to 'changes'...")
#     obj <- as.matrix(object@neighborhood_changes)
#     legend_name <- "Changes"
#     if(is.null(title))
#       title <- paste0("d(", object@conditions[2], ") - d(",object@conditions[1],")" )
#
#   } else if(what == "changes_2"){
#
#     print_this_msg("Parameter 'what' is set to 'changes_2'...")
#     obj <- object@neighborhood[[1]] - object@neighborhood[[2]]
#     legend_name <- "Changes"
#     if(is.null(title))
#       title <- paste0("d(", object@conditions[1], ") - d(",object@conditions[2],")" )
#
#   }else if(what=="c1"){
#
#     print_this_msg("Parameter 'what' is set to 'c1'...")
#     obj <- object@neighborhood[[1]]
#     legend_name <- "Dist"
#     if(is.null(title))
#       title <- paste0("d(", object@conditions[1], ")")
#
#   }else {
#
#     print_this_msg("Parameter 'what' is set to 'c2'...")
#     obj <- object@neighborhood[[2]]
#     legend_name <- "Dist"
#     title <- paste0("d(", object@conditions[2], ")")
#
#   }
#
#   if(!is.null(del_feat)){
#     del_feat <- unique(del_feat)
#     if(all(del_feat %in% feat_names(object))){
#       obj <- obj[!rownames(obj) %in% del_feat,
#                  !colnames(obj) %in% del_feat]
#     }else{
#       print_this_msg("Some feature to delete were not found in the object", msg_type = "STOP")
#     }
#   }
#
#   if(!is.null(only_feat)){
#       only_feat <- unique(only_feat)
#       if(all(del_feat %in% feat_names(object))){
#         obj <- obj[rownames(obj) %in% only_feat,
#                    colnames(obj) %in% only_feat]
#       }else{
#         print_this_msg("Some feature to delete were not found in the object", msg_type = "STOP")
#       }
#   }
#
#   if(nrow(obj) == 0){
#     print_this_msg("No more feature left", msg_type = "STOP")
#   }
#
#   if(filter_method == "diff"){
#     if(!is.null(filter)){
#       TF <- abs(obj) > filter
#       obj <- obj[rowSums(TF) > 0 , colSums(TF) > 0]
#
#       if(is.vector(obj))
#         print_this_msg("No feature left after filtering. Adapt 'filter' argument please.", msg_type = "STOP")
#     }else{
#       print_this_msg("'filter' is set to NULL, no filtering...")
#     }
#   }else if(filter_method=="cv"){
#     if(!is.null(filter)){
#
#       sd_row <- apply(abs(obj), 1, sd)
#       mean_row <- apply(abs(obj), 1, mean)
#       cv_row <- sd_row/mean_row
#
#       sd_col <- apply(abs(obj), 2, sd)
#       mean_col <- apply(abs(obj), 2, mean)
#       cv_col <- sd_col/mean_col
#
#       obj <- obj[cv_row > filter , cv_col > filter]
#
#       if(is.vector(obj))
#         print_this_msg("No feature left after filtering. Adapt 'filter' argument please.", msg_type = "STOP")
#     }else{
#       print_this_msg("'filter' is set to NULL, no filtering...")
#     }
#   }
#
#
#
#   print_this_msg("Calling ggheatmap...")
#   p <- ggheatmap::ggheatmap(obj,
#                  dist_method=dist_method,
#                  hclust_method = hclust_method,
#                  cluster_rows = TRUE,
#                  cluster_cols = TRUE,
#                  legendName=legend_name,
#                  scale = "none",
#                  color=RColorBrewer::brewer.pal(7, "Spectral"),
#                  border = "darkgrey",
#                  tree_color_rows = "darkgrey",
#                  tree_color_cols = "darkgrey",
#                  annotation_rows = NULL,
#                  annotation_cols = NULL,
#                  annotation_color = NULL
#   )
#
#   p$plotlist[[1]] <- p$plotlist[[1]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
#                                                                                           size = size,
#                                                                                           hjust=1),
#                                                       axis.text.y = ggplot2::element_text(size=size))
#
#   p$plotlist[[3]] <- p$plotlist[[3]] + ggplot2::ggtitle(title)
#
#   p
#
# }
# )

# -------------------------------------------------------------------------
##    REDEFINE SHOW() METHOD FOR CLASS OBJECT : STCompR
# -------------------------------------------------------------------------

#' @title The show method of a STCompR object
#' @description
#' The show method of a STCompR object displays information about the object,
#' including conditions, memory usage, and the number of features
#' @param object A STCompR object.
#' @keywords internal
#' @examples
#' example_dataset("11284296/files/cmp_xen")
#' show(cmp_xen)
#' @import methods
#' @export
setMethod("show", signature("STCompR"),
  function(object) {
    print_this_msg("An object of class STCompR")
    for(i in 1:length(object@conditions)){
      print_this_msg(names(object@conditions)[i], object@conditions[i])
    }
    print_this_msg("Number of features: ", nrow(object@raw_counts))
    print_this_msg(">>> Please, use show_st_methods(class = 'STCompR') to show availables methods <<<")
  }
)


# -------------------------------------------------------------------------
##    Some basic functions
# -------------------------------------------------------------------------

#' @title The number of features stored in a STCompR object
#' @description
#' The number of features stored in a STCompR object
#' @param object The STCompR object
#' @examples
#' example_dataset("11284296/files/cmp_xen")
#' nb_feat(cmp_xen)
#' @export
setMethod("nb_feat", signature(object = "STCompR"),
           function(object)
             nrow(object@stat_test)
)

#' @title The features stored in a STCompR object
#' @description
#' The features stored in a STCompR object
#' @param object The STCompR object
#' @keywords internal
#' @examples
#' example_dataset("11284296/files/cmp_xen")
#' feat_names(cmp_xen)
#' @export
#' @name feat_names
if(!isGeneric("feat_names")){
  setGeneric("feat_names",
             function(object)
               standardGeneric("feat_names")
  )
}

#' @title The features stored in a STCompR object
#' @description
#' The features stored in a STCompR object
#' @param object The STCompR object
#' @examples
#' example_dataset("11284296/files/cmp_xen")
#' feat_names(cmp_xen)
#' @export
setMethod("feat_names", signature(object="STCompR"),
           function(object)
             rownames(object@stat_test)
)

# -------------------------------------------------------------------------
##    The normalized counts and test results of a STCompR object.
# -------------------------------------------------------------------------

#' @title The normalized counts and test results of a STCompR object.
#' @description
#' The normalized counts and test results of a STCompR object.
#' @param object The STCompR object
#' @param transform Whether the count are transformed (the pseudo count defined for the object is added).
#' @param count_only Returns only counts.
#' @param melted_count Returns a melted table with counts.
#' @param normalized Counts are normalized.
#' @param features The features that should be returned. Default to all (NULL).
#' @examples
#' example_dataset("11284296/files/cmp_xen")
#' stat_test(cmp_xen)
#' @keywords internal
setGeneric("stat_test",
           function(object,
                    transform=c("None", "log2", "log10", "log"),
                    count_only=FALSE,
                    melted_count=FALSE,
                    normalized=TRUE,
                    features=NULL)
             standardGeneric("stat_test")
)

#' @title The normalized counts and test results of a STCompR object.
#' @description
#' The normalized counts and test results of a STCompR object.
#' @param object The STCompR object
#' @param transform Whether the count are transformed (the pseudo count defined for the object is added).
#' @param count_only Returns only counts.
#' @param melted_count Returns a melted table when count_only is TRUE.
#' @param normalized Counts are normalized.
#' @param features The features that should be returned. Default to all (NULL).
#' @examples
#' example_dataset("11284296/files/cmp_xen")
#' stat_test(cmp_xen)
#' @export stat_test
setMethod(f="stat_test",
          signature=signature(object="STCompR"),
           function(object,
                    transform=c("None", "log2", "log10", "log"),
                    count_only=FALSE,
                    melted_count=FALSE,
                    normalized=TRUE,
                    features=NULL){

             transform <- match.arg(transform)

             counts <- object@raw_counts

             if(!is.null(features)){
               if(!all(features %in% feat_names(object))){
                 print_this_msg("Some features where not found.", msg_type = "STOP")
               }
             }else{
               features <- feat_names(object)
             }

             counts <- counts[features, ]

             if(normalized)
               counts <- object@norm_counts

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
                 colnames(d) <- c("Features", "Samples", "Counts")

               }
             }else{
               d <- cbind(counts,
                          object@stat_test)
             }

            return(d)

})

