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
#' @slot scaling_factor A numeric vector representing the scaling factors used during normalization (to be divided).
#' @slot stat_test A data frame containing statistical test results, including log2 ratios,
#'                odds ratios, and p-values for each feature.
#' @slot pseudo_count The value for the pseudo_count.
#' @slot neighborhood a list of feature-feature neighborhood for the two conditions to be compared.
#' @slot neighborhood_changes A feature-feature matrix indicative of neighborhood changes.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' example_dataset("10819270/files/cmp_xen")
#' cmp_xen
#'
setClass("STCompR",
         slots = c(
           conditions = "character",
           method="character",
           raw_counts = "data.frame",
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
           scaling_factor = numeric(),
           stat_test = data.frame(),
           pseudo_count= 1,
           neighborhood=list(),
           neighborhood_changes=data.frame()
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
#' @param cts A matrix with count data.
#'
#' @return A numeric vector of size factors for normalization.
#'
#' @details The DESeq2 RLE method estimates size factors for each sample based on
#' the relative log expression (RLE). The geometric mean of counts for each feature
#' across all samples is computed. Zero counts are set to NA to avoid division by zero.
#' Each count is then divided by its corresponding geometric mean, and the median
#' of these ratios across all samples is computed to obtain the size factor for each sample.
#'
#' @export
estimSf <- function(cts) {

  geomMean <- function(x) prod(x)^(1/length(x))

  # Compute the geometric mean over the line
  gm.mean  <-  apply(cts, 1, geomMean)

  # Zero values are set to NA (avoid subsequent division by 0)
  gm.mean[gm.mean == 0] <- NA

  # Divide each line by its corresponding geometric mean
  cts <- sweep(cts, 1, gm.mean, FUN="/")

  # Compute the median over the columns
  med <- apply(cts, 2, stats::median, na.rm=TRUE)

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
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r2 <- xen[x_bins, y_bins]
#' cmp <- stcompr(xen_r1, xen_r2)
#'
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

  check_this_var(name_1)
  check_this_var(name_2)

  print_this_msg("Checking objects.")

  if(class(object_1)[1] != "STGrid" | class(object_2)[1] != "STGrid"){
    print_this_msg("Please provide STGrid objects.")
  }

  gn_1 <- feat_names(object_1)
  gn_2 <- feat_names(object_2)

  print_this_msg("Checking feature names.")

  inter_gn <- intersect(gn_1, gn_2)

  if(length(inter_gn) == 0){
    print_this_msg("No shared features between objects")
  }

  print_this_msg("Retrieving counts.")

  tb_1 <- data.frame(row.names = names(table(coord(object_1)$feature)),
                     counts=as.numeric(table(coord(object_1)$feature)))

  tb_2 <- data.frame(row.names = names(table(coord(object_2)$feature)),
                     counts=as.numeric(table(coord(object_2)$feature)))

  print_this_msg("Merging counts.")

  raw_counts <- cbind(tb_1[inter_gn, ], tb_2[inter_gn, ])
  rownames(raw_counts) <- inter_gn
  colnames(raw_counts) <- c(name_1, name_2)

  print_this_msg("Normalizing...")
  scaling_factor <- estimSf(raw_counts)

  norm_counts <- sweep(raw_counts, MARGIN = 2, STATS =  scaling_factor, FUN="/")

  norm_counts <- norm_counts + pseudo_count

  print_this_msg("Computing log2 ratio...")
  ratio <- norm_counts[, name_2]/norm_counts[, name_1]
  log2_ratio <- log2(ratio)

  p_values <- c()
  odd_ratio <- c()

  print_this_msg("Computing Fisher tests...")

  for(g in 1:nrow(norm_counts)){

    a <- norm_counts[g, name_2]
    b <- norm_counts[g, name_1]
    c <- sum(norm_counts[-g,name_2])
    d <- sum(norm_counts[-g, name_1])

    m <- matrix(c(round(a, 0),
                  round(b, 0),
                  round(c, 0) ,
                  round(d, 0)),
                byrow = TRUE, ncol=2)

    ft <- stats::fisher.test(m)
    p_values[g] <- ft$p.value
    odd_ratio[g] <- ft$estimate["odds ratio"]
  }

  stats <- data.frame(row.names = rownames(norm_counts),
                      log2_ratio=log2_ratio,
                      odd_ratio=odd_ratio,
                      p_values=p_values)

  stats$p_values[stats$p_values < 1e-320] <- 1e-320

  print_this_msg("Adjusting p-values...")

  for(corr in c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")){
    stats[, paste0("padj_", corr)] <- stats::p.adjust(stats$p_values, method=corr)
  }

  print_this_msg("Preparing neighborhood analysis...")

  bin_mat_1 <- bin_mat(object_1, del_bin = TRUE)
  bin_mat_2 <- bin_mat(object_2, del_bin = TRUE)
  bin_mat_1 <- bin_mat_1[, inter_gn]
  bin_mat_2 <- bin_mat_2[, inter_gn]

  bin_mat_1 <- bin_mat_1[rowSums(bin_mat_1) > 0, ]
  bin_mat_2 <- bin_mat_2[rowSums(bin_mat_2) > 0, ]

  bin_mat_1 <- sweep(bin_mat_1,
                     MARGIN = 2,
                     STATS = colSums(bin_mat_1, na.rm = TRUE), FUN = "/")

  bin_mat_2 <- sweep(bin_mat_2,
                     MARGIN = 2,
                     STATS = colSums(bin_mat_2, na.rm = TRUE), FUN = "/")


  spatial_matrix_ratio_1 <- as.matrix(stats::dist(t(bin_mat_1),
                                           method = "manhattan"))

  spatial_matrix_ratio_2 <- as.matrix(stats::dist(t(bin_mat_2),
                                           method = "manhattan"))

  print_this_msg("Preparing an STCompR object... ")

  STCompR <- methods::new("STCompR")
  STCompR@neighborhood <- list(spatial_matrix_ratio_1,
                               spatial_matrix_ratio_2)

  STCompR@neighborhood_changes <- as.data.frame(spatial_matrix_ratio_2 - spatial_matrix_ratio_1)

  names(STCompR@neighborhood) <- c(name_1, name_2)
  STCompR@conditions <- c(name_1, name_2)
  STCompR@method <- c(object_1@method, object_2@method)
  names(STCompR@method) <- STCompR@conditions
  STCompR@scaling_factor <- scaling_factor
  STCompR@stat_test <- stats
  STCompR@raw_counts <- as.data.frame(raw_counts)
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
#' @param what A character string specifying what to plot: "changes" (the difference in pairwise Manhattan distances between conditions 2 and 1),
#' "changes_2" (the difference in pairwise Manhattan distances between conditions 1 and 2), "c1" (the pairwise Manhattan distances of features in condition 1), or "c2"
#' (the pairwise Manhattan distances of features in condition 2).
#' @param hclust_method A character string specifying the agglomerative criteria for hierarchical clustering ("ward.D", "ward.D2",
#' "single", "complete", "average", "mcquitty", "median", "centroid").
#' @param dist_method A character string specifying the distance method to be used ("euclidean",
#' "maximum", "manhattan", "canberra", "binary","minkowski").
#' @param del_feat A character vector specifying features to be excluded from the heatmap.
#' @param only_feat A character vector specifying features to be included in the heatmap.
#' @param filter A numeric value specifying the threshold for filtering out low values in the heatmap.
#' @param size A numeric value specifying the size of text in the heatmap.
#'
#' @return A heatmap plot comparing conditions in the 'STCompR' object.
#'
#' @examples
#' # Example usage:
#' example_dataset("10819270/files/cmp_xen")
#' heatmap_cmp(object = cmp_xen, hclust_method = "ward.D", dist_method = "euclidean")
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
                    del_feat=NULL,
                    only_feat=NULL,
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
#' @param what A character string specifying what to plot: "changes" (the difference in pairwise Manhattan distances between conditions 2 and 1),
#' "changes_2" (the difference in pairwise Manhattan distances between conditions 1 and 2), "c1" (the pairwise Manhattan distances of features in condition 1), or "c2"
#' (the pairwise Manhattan distances of features in condition 2).
#' @param hclust_method A character string specifying the agglomerative criteria for hierarchical clustering ("ward.D", "ward.D2",
#' "single", "complete", "average", "mcquitty", "median", "centroid").
#' @param dist_method A character string specifying the distance method to be used ("euclidean",
#' "maximum", "manhattan", "canberra", "binary","minkowski").
#' @param del_feat A character vector specifying feature to be excluded from the heatmap.
#' @param only_feat A character vector specifying features to be included in the heatmap.
#' @param filter A numeric value specifying the threshold for filtering out low values in the heatmap.
#' @param size A numeric value specifying the size of text in the heatmap.
#'
#' @return A heatmap plot comparing conditions in the 'STCompR' object.
#'
#' @examples
#' # Example usage:
#' example_dataset("10819270/files/cmp_xen")
#' heatmap_cmp(object = cmp_xen, hclust_method = "ward.D", dist_method = "euclidean")
#' @importFrom ggheatmap ggheatmap ggheatmap_theme
#' @importFrom RColorBrewer brewer.pal
#' @import magrittr
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
           del_feat=NULL,
           only_feat=NULL,
           filter=0.2,
           size=6) {

    what <- match.arg(what)
    dist_method <- match.arg(dist_method)
    hclust_method <- match.arg(hclust_method)

  if(what=="changes"){

    print_this_msg("Parameter 'what' is set to 'changes'...")
    obj <- as.matrix(object@neighborhood_changes)
    legend_name <- "Changes"

  } else if(what == "changes_2"){

    print_this_msg("Parameter 'what' is set to 'changes_2'...")
    obj <- object@neighborhood[[1]] - object@neighborhood[[2]]
    legend_name <- "Changes"

  }else if(what=="c1"){

    print_this_msg("Parameter 'what' is set to 'c1'...")
    obj <- object@neighborhood[[1]]
    legend_name <- "Dist"

  }else {

    print_this_msg("Parameter 'what' is set to 'c2'...")
    obj <- object@neighborhood[[2]]
    legend_name <- "Dist"

  }

  if(!is.null(del_feat)){
    del_feat <- unique(del_feat)
    if(all(del_feat %in% feat_names(object))){
      obj <- obj[!rownames(obj) %in% del_feat,
                 !colnames(obj) %in% del_feat]
    }else{
      print_this_msg("Some feature to delete were not found in the object", msg_type = "STOP")
    }
  }

  if(!is.null(only_feat)){
      only_feat <- unique(only_feat)
      if(all(del_feat %in% feat_names(object))){
        obj <- obj[rownames(obj) %in% only_feat,
                   colnames(obj) %in% only_feat]
      }else{
        print_this_msg("Some feature to delete were not found in the object", msg_type = "STOP")
      }
  }

  if(nrow(obj) == 0){
    print_this_msg("No more feature left", msg_type = "STOP")
  }

  if(!is.null(filter)){
    TF <- abs(obj) > filter
    obj <- obj[rowSums(TF) > 0 , colSums(TF) > 0]
  }


  print_this_msg("Calling ggheatmap...")
  p <- ggheatmap::ggheatmap(obj,
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
  )
  #%>% ggheatmap::ggheatmap_theme(1,
  #                     theme=list(ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
  #                                                                 size = size, hjust=1),
  #                                      axis.text.y = ggplot2::element_text(size=size))))
  p

}
)

# -------------------------------------------------------------------------
##    REDEFINE SHOW() METHOD FOR CLASS OBJECT : STCompR
# -------------------------------------------------------------------------

#' @title The show method of a STCompR object
#' @description
#' The show method of a STCompR object displays information about the object,
#' including conditions, memory usage, and the number of features
#'
#' @param object A STCompR object.
#' @keywords internal
#' @examples
#' example_dataset("10819270/files/cmp_xen")
#' show(cmp_xen)
#' @import methods
#' @export
setMethod("show", signature("STCompR"),
  function(object) {
    print_this_msg("An object of class STCompR")
    for(i in 1:length(object@conditions)){
      print_this_msg("Condition", i, object@conditions[i])
    }
    print_this_msg("Memory used: ", utils::object.size(object))
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
#' example_dataset("10819270/files/cmp_xen")
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
#' example_dataset("10819270/files/cmp_xen")
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
#' example_dataset("10819270/files/cmp_xen")
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
#' example_dataset("10819270/files/cmp_xen")
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
#' @param melted_count Returns a melted table with counts.
#' @param normalized Counts are normalized.
#' @param features The features that should be returned. Default to all (NULL).
#' @examples
#' example_dataset("10819270/files/cmp_xen")
#' stat_test(cmp_xen)
#' @export stat_test
setMethod("stat_test",
          "STCompR",
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
                 colnames(d) <- c("Features", "Conditions", "Counts")
               }
             }else{
               d <- cbind(counts,
                          object@stat_test)
             }

             if(melted_count){
               return(d[d$Features %in% features, ])
             }else{
               return(d[features,])
             }

}
)

