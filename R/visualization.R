# -------------------------------------------------------------------------
##    Spatial image
# -------------------------------------------------------------------------
#' Color-coded representation of the object (e.g. molecules) density
#'
#' This function displays a color-coded representation of the object (e.g. molecules) density observed in a spatial transcriptomics experiment.
#'
#' @param object The STGrid object.
#' @param features The name of the features for which the spatial image will be created.
#' @param saturation The ceiling level for the feature expression values. Defaults to 1 (no ceiling).
#' @param scale Logical value indicating whether to scale the feature expression values. Defaults to TRUE.
#' @param colors The colors to use for gradient fill in the spatial image. Defaults to a set of colors.
#' @param coord_fixed Logical value indicating whether to keep the aspect ratio fixed. Defaults to TRUE.
#' @param overlay_feature The feature to overlay on the spatial image. Defaults to NULL.
#' @param colors_overlay The colors to use for gradient fill in the overlay feature. Defaults to a set of colors.
#' @param grid_by Whether to overlay a grid with horizontal and vertical lines at particular interval.
#' No overlay if NULL otherwise the size of the interval (e.g. 20).
#' @param color_grid A color for the grid.
#' @param size The size of the overlayed points.
#' @param logb The basis for the log transformation. Default to 10. If NULL no log transformation.
#' @param pseudo_count a value for the pseudo count used for log transformation (default to 1).
#' @param ncol The number of columns for the facets.
#' @keywords internal
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' spatial_image(xen,
#'              features="Chat",
#'              grid_by=50)
#' bx <- bin_x(xen)[175:nbin_x(xen)]
#' by <- bin_y(xen)[100:nbin_y(xen)]
#' sub <- xen[bx, by]
#' spatial_image(sub,
#'              features=c("Chat", "Ano1"))
#' spatial_image(sub,
#'              features=c("Chat", "Ano1"),
#'              overlay_feature="Nwd2",
#'              size=0.05)
#' @export
setGeneric("spatial_image",
           function(object = NULL,
                    features = NULL,
                    saturation = 1,
                    scale = TRUE,
                    colors = c("black", "#33FF00", "#FFFF00",
                               "#FF0000", "#CC00FF"),
                    coord_fixed = TRUE,
                    overlay_feature = NULL,
                    colors_overlay = c("#DEEBF7", "#9ECAE1", "#3182BD", "blue"),
                    grid_by = NULL,
                    color_grid = "white",
                    size = 0.5,
                    logb = 10,
                    pseudo_count = 1,
                    ncol=4)
           standardGeneric("spatial_image"))

#' Color-coded representation of the object (e.g. molecules) density
#'
#' This function displays a color-coded representation of the object (e.g. molecules) density observed in a spatial transcriptomics experiment.
#'
#' @param object The STGrid object.
#' @param features The name of the features for which the spatial image will be created.
#' @param saturation The ceiling level for the feature expression values. Defaults to 1 (no ceiling).
#' @param scale Logical value indicating whether to scale the feature expression values. Defaults to TRUE.
#' @param colors The colors to use for gradient fill in the spatial image. Defaults to a set of colors.
#' @param coord_fixed Logical value indicating whether to keep the aspect ratio fixed. Defaults to TRUE.
#' @param overlay_feature The feature to overlay on the spatial image. Defaults to NULL.
#' @param colors_overlay The colors to use for gradient fill in the overlay feature. Defaults to a set of colors.
#' @param grid_by Whether to overlay a grid with horizontal and vertical lines at particular interval.
#' No overlay if NULL otherwise the size of the interval (e.g. 20).
#' @param color_grid A color for the grid.
#' @param size The size of the overlayed points.
#' @param logb The basis for the log transformation. Default to 10. If NULL no log transformation.
#' @param pseudo_count a value for the pseudo count used for log transformation (default to 1).
#' @param ncol The number of columns for the facets.
#' @importFrom ggplot2 aes coord_fixed facet_wrap geom_point geom_tile scale_fill_gradientn theme xlab ylab element_blank element_rect element_text scale_x_discrete scale_y_discrete
#' @importFrom reshape2 melt
#' @importFrom viridis inferno
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' spatial_image(xen,
#'              features="Chat",
#'              grid_by=50)
#' bx <- bin_x(xen)[175:nbin_x(xen)]
#' by <- bin_y(xen)[100:nbin_y(xen)]
#' sub <- xen[bx, by]
#' spatial_image(sub,
#'              features=c("Chat", "Ano1"))
#' spatial_image(sub,
#'              features=c("Chat", "Ano1"),
#'              overlay_feature="Nwd2",
#'              size=0.05)
#' @export
setMethod("spatial_image",
          signature(object = "STGrid"),
          function(object = NULL,
                   features = NULL,
                   saturation = 1,
                   scale = TRUE,
                   colors = c("black", "#33FF00", "#FFFF00",
                              "#FF0000", "#CC00FF"),
                   coord_fixed = TRUE,
                   overlay_feature = NULL,
                   colors_overlay = c("#DEEBF7", "#9ECAE1", "#3182BD", "blue"),
                   grid_by = NULL,
                   color_grid = "white",
                   size = 0.5,
                   logb = 10,
                   pseudo_count = 1,
                   ncol=4) {

            check_this_var(grid_by, null_accepted = TRUE, type = "int")

            if (saturation > 1 | saturation < 0)
              print_this_msg("Saturation should be between 0 and 1.",
                             msg_type = "STOP")

            print_this_msg("Checking arguments", msg_type = "DEBUG")

            if (is.null(object))
              print_this_msg("Please provide an STGrid object.",
                             msg_type = "STOP")

            if (!is.null(overlay_feature)) {
              if (!overlay_feature %in% c(feat_names(object), meta_names(object)))
                print_this_msg("The feature to overlay was not found in the object.",
                               msg_type = "STOP")

            }

            print_this_msg("Checking features", msg_type = "DEBUG")

            if (is.null(features))
              print_this_msg("Please provide a feature name (see feature arguments).",
                             msg_type = "STOP")


            if (!all(features %in% c(feat_names(object), colnames(object@meta)))) {
              print_this_msg("The feature was not found in the object.", msg_type = "STOP")
            }

            if (any(colnames(object@meta) %in% features)) {
              print_this_msg("Using a feature from meta slot.", msg_type = "DEBUG")
              spatial_matrix <- object@bin_mat[, c("bin_x", "bin_y",
                                                   setdiff(features, colnames(object@meta)))]
              for(i in features[features %in% colnames(object@meta)]){
                spatial_matrix[,i] <- object@meta[, i, drop=FALSE]
              }


            } else{
              spatial_matrix <- object@bin_mat[, c("bin_x", "bin_y", features)]
            }

            tmp <- spatial_matrix[, features, drop = FALSE]

            for (i in 1:ncol(tmp)) {
              if (saturation < 1) {
                q_sat <- stats::quantile(tmp[, i][tmp[, i] != 0], saturation)
                tmp[tmp[, i] > q_sat, i] <- q_sat
              }

              if (!is.null(logb)) {
                tmp[, i] <- log(tmp[, i] + pseudo_count, base = logb)
              }

              if (scale) {
                tmp[, i] <- (tmp[, i] - min(tmp[, i])) / (max(tmp[, i]) - min(tmp[, i]))
              }
            }


            spatial_matrix[, features] <- tmp

            spatial_matrix$bin_x <- factor(spatial_matrix$bin_x,
                                           levels = bin_x(object),
                                           ordered = TRUE)

            spatial_matrix$bin_y <- factor(spatial_matrix$bin_y,
                                           levels = bin_y(object),
                                           ordered = TRUE)

            spatial_matrix_melted <-
              reshape2::melt(spatial_matrix, id.vars = c("bin_x", "bin_y"))
            spatial_matrix_melted$variable <-
              factor(spatial_matrix_melted$variable,
                     levels = features,
                     ordered = TRUE)

            bin_x <- bin_y <- value <- .data <- NULL
            p <- ggplot2::ggplot(data = spatial_matrix_melted,
                                 mapping = ggplot2::aes(x = bin_x,
                                                        y = bin_y,
                                                        fill = value)) +
              ggplot2::geom_tile() +
              ggplot2::xlab("")  +
              ggplot2::ylab("") +
              ggplot2::theme(
                axis.ticks = ggplot2::element_blank(),
                axis.text = ggplot2::element_blank(),
                strip.background = ggplot2::element_rect(fill =
                                                           "gray30"),
                strip.text = ggplot2::element_text(color =
                                                     "white")
              ) +
              ggplot2::scale_fill_gradientn(colours = colors) +
              ggplot2::facet_wrap( ~ variable, ncol = ncol)

            if (coord_fixed)
              p <- p + ggplot2::coord_fixed()


            if (!is.null(overlay_feature)) {
              if(overlay_feature %in% feat_names(object)){
                over <- bin_mat(object,
                                as_factor = TRUE)[, c("bin_x", "bin_y", overlay_feature)]
                over <- over[over[[overlay_feature]] != 0, ]
              }else{
                over <- bin_mat(object,
                                as_factor = TRUE)[, c("bin_x", "bin_y")]
                over[[overlay_feature]] <- object@meta[[overlay_feature]]
                over <- over[over[[overlay_feature]] != 0, ]
              }


              p <-
                p + ggplot2::geom_point(
                  data = over,
                  mapping = ggplot2::aes(
                    x = bin_x,
                    y = bin_y,
                    color = log10(.data[[overlay_feature]])
                  ),
                  size = size,
                  inherit.aes = FALSE
                ) +
                ggplot2::scale_color_gradientn(colors = colors_overlay)
            }

            if (!is.null(grid_by)) {
              lev_bin_x <- levels(spatial_matrix$bin_x)
              lev_bin_y <- levels(spatial_matrix$bin_y)

              x_seq <- seq(from = 1,
                           to = length(lev_bin_x),
                           by = grid_by)
              y_seq <- seq(from = 1,
                           to = length(lev_bin_y),
                           by = grid_by)

              label_x <-
                stats::setNames(lev_bin_x, 1:length(lev_bin_x))
              names(label_x)[-x_seq] <- ""
              label_y <-
                stats::setNames(lev_bin_y, 1:length(lev_bin_y))
              names(label_y)[-y_seq] <- ""

              p <-
                p + geom_vline(
                  data = data.frame(bin_x = levels(spatial_matrix$bin_x)[x_seq]),
                  mapping = aes(xintercept = bin_x),
                  color = color_grid
                ) +
                geom_hline(
                  data = data.frame(bin_y = levels(spatial_matrix$bin_y)[y_seq]),
                  mapping = aes(yintercept = bin_y),
                  color = color_grid
                ) +
                ggplot2::theme(axis.text = element_text(size = 6, angle = 45)) +
                ggplot2::scale_x_discrete("x", labels = names(label_x)) +
                ggplot2::scale_y_discrete("y", labels = names(label_y))
            }

            return(p)

          })

# -------------------------------------------------------------------------
##    The spatial_plot function
# -------------------------------------------------------------------------
#' Plot x/y coordinates of molecules
#'
#' Plot x/y coordinates of molecules of a spatial transcriptomics experiment.
#'
#' @param object The STGrid object.
#' @param feat_list The name of the features for which the spatial image will be created.
#' @param colors The colors to use for features in the spatial image.
#' @param size The size of the points
#' @param coord_fixed Logical value indicating whether to keep the aspect ratio fixed. Defaults to TRUE.
#' @export
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' spatial_plot(xen,
#'              feat_list=c("Chat", "Nwd2", "Ano1"),
#'              size=0.05)
#' @keywords internal
setGeneric("spatial_plot",
           function(object = NULL,
                    feat_list = NULL,
                    colors = NULL,
                    size = 0.1,
                    coord_fixed = TRUE)
             standardGeneric("spatial_plot"))

#' Plot x/y coordinates of molecules
#'
#' Plot x/y coordinates of molecules of a spatial transcriptomics experiment.
#'
#' @param object The STGrid object.
#' @param feat_list The name of the features for which the spatial image will be created.
#' @param colors The colors to use for features in the spatial image.
#' @param size The size of the points
#' @param coord_fixed Logical value indicating whether to keep the aspect ratio fixed. Defaults to TRUE.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' spatial_plot(xen,
#'              feat_list=c("Chat", "Nwd2", "Ano1"),
#'              size=0.05)
#' @importFrom ggplot2 aes geom_point scale_color_manual xlab ylab theme_minimal theme coord_fixed geom_point scale_color_manual aes
#' @export
setMethod("spatial_plot", "STGrid",
          function(object = NULL,
                   feat_list = NULL,
                   colors = NULL,
                   size = 0.1,
                   coord_fixed = TRUE) {
            if (is.null(object))
              print_this_msg("Please provide an STGrid object.",
                             msg_type = "STOP")
            if (!is.null(colors)) {
              if (length(colors) < length(feat_list))
                print_this_msg("More colors are needed...", msg_type = "STOP")
            }

            if (is.null(feat_list))
              print_this_msg("Please provide feature names (see feat_list arguments).",
                             msg_type = "STOP")

            if (!all(feat_list %in% feat_names(object)))
              print_this_msg("One or several features were not found in the object.", msg_type = "STOP")

            coord <-
              get_coord(object, feat_list = feat_list, as.factor = TRUE)

            x <- y <- feature <- NULL
            p <- ggplot2::ggplot(data = coord,
                                 mapping = ggplot2::aes(x = x,
                                                        y = y,
                                                        color = feature)) +
              ggplot2::geom_point(size = size) +
              ggplot2::xlab("x")  +
              ggplot2::ylab("y") +
              ggplot2::theme_minimal() +
              ggplot2::theme(axis.text = element_text(size = 6))

            if (!is.null(colors))
              p <- p + ggplot2::scale_color_manual(values = colors)

            if (coord_fixed)
              p <- p + ggplot2::coord_fixed()

            return(p)

          })

# -------------------------------------------------------------------------
##     cmp_bar_plot()
# -------------------------------------------------------------------------
#' @title Create a barplot to show counts for selected features.
#' @description
#' Create a barplot to show counts for selected features.
#' @param object A STGrid object.
#' @param features The list of features (NULL for all of them).
#' @param normalized  Whether counts should be normalized.
#' @param transform Whether the count should be transformed (the pseudo count defined for the object is added).
#' @param colors A set of colors.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r2 <- xen[x_bins, y_bins]
#' cmp <- stcompr(xen_r1, xen_r2)
#' cmp_bar_plot(cmp,
#'              features=c("Chat", "Nwd2", "Ano1"))
#' @keywords internal
setGeneric("cmp_bar_plot",
           function(object,
                    features = utils::head(feat_names(object)),
                    normalized = FALSE,
                    transform = c("None", "log2", "log10", "log"),
                    colors = c("#3074BB", "#BE5B52"))
             standardGeneric("cmp_bar_plot"))

#' @title Create a barplot to show counts for selected features.
#' @description
#' Create a barplot to show counts for selected features.
#' @param object A STCompR object.
#' @param features The list of features (NULL for all of them).
#' @param normalized  Whether counts should be normalized.
#' @param transform Whether the count should be transformed (the pseudo count defined for the object is added).
#' @param colors A set of colors.
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r2 <- xen[x_bins, y_bins]
#' cmp <- stcompr(xen_r1, xen_r2)
#' cmp_bar_plot(cmp,
#'              features=c("Chat", "Nwd2", "Ano1"))
#' @export
setMethod("cmp_bar_plot", signature("STCompR"),
          function(object,
                   features = utils::head(feat_names(object)),
                   normalized = FALSE,
                   transform = c("None", "log2", "log10", "log"),
                   colors = c("#3074BB", "#BE5B52")) {
            if (is.null(features)) {
              print_this_msg("Please provide some features...", msg_type = "STOP")
            }

            transform <- match.arg(transform)

            counts <- stat_test(
              object,
              normalized = normalized,
              count_only = TRUE,
              melted_count = TRUE,
              transform = transform,
              features = features
            )

            counts$Conditions <- object@conditions[as.character(counts$Samples)]
            counts$Conditions <- as.factor(counts$Condition)

            ylabel <- ifelse(
              transform %in% c("log2", "log10", "log"),
              paste0(transform, "(Feature counts)"),
              "Feature counts"
            )

            counts$Features <- factor(counts$Features,
                                      levels = features,
                                      ordered = TRUE)

            Samples <- Counts <- Conditions <- NULL

            ggplot2::ggplot(data = counts,
                            mapping = ggplot2::aes(x = Samples,
                                                   y = Counts,
                                                   fill = Conditions)) +
              ggplot2::geom_col(color = "black",
                                linewidth = 0,
                                position = "dodge") +
              ggplot2::theme_bw() +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                  size = 10,
                  angle = 45,
                  vjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(size = 8),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank(),
                panel.grid.minor.y = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank()
              ) +
              ggplot2::ylab(ylabel) +
              ggplot2::scale_fill_manual(values = colors) +
              ggplot2::facet_wrap(~Features)

})

# -------------------------------------------------------------------------
##    Compare counts across multiple st_grid_objects
# -------------------------------------------------------------------------
#' @title Compare feature counts across multiple STGrid objects.
#' @description
#' This function compares feature counts across multiple STGrid objects.
#'
#' @param ... A set of STGrid objects.
#' @param features A character vector specifying the features (e.g., genes) to compare. Default is NULL.
#' @param normalized Logical indicating whether the counts should be normalized. Default is FALSE.
#' @param type Type of plot to generate. Currently only barplot is supported.
#' @param names Optional character vector specifying names for each STGrid object. If NULL, default names will be assigned.
#' @param transform Transformation method for the counts. Options are "None" (default), "log2", "log10", or "log".
#' @param fill_color Optional character vector specifying fill colors for different conditions in the plot.
#' @param border_color Color for the border of bars in bar plots. Default is "black".
#'
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual theme_bw element_text
#'
#' @return A ggplot object representing the comparison of counts across different conditions or samples.
#'
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_r2 <- xen[x_bins, y_bins]
#'
#' cmp_counts_st(xen, xen_r1, xen_r2, features = c("Chat", "Ano1"), normalized = FALSE)
#' cmp_counts_st(xen, xen_r1, xen_r2,
#'              features = c("Chat", "Ano1"),
#'               normalized = TRUE, fill_color=rainbow(3))
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors counts
#' @export
cmp_counts_st <- function(...,
         features = NULL,
         normalized = FALSE,
         type=c("barplot"),
         names=NULL,
         transform = c("None",
                       "log2",
                       "log10",
                       "log"),
         fill_color = NULL,
         border_color="black") {

  type <- match.arg(type)
  transform <- match.arg(transform)

  st_list <- list(...)
  print_this_msg("Found ", length(st_list), "STGrid objects.",
                 msg_type = "INFO")

  if(is.list(st_list[[1]]))
    st_list <- st_list[[1]]

  if(is.null(features))
    print_this_msg("Provide at least a feature...",
                   msg_type = "STOP")

  check_st_list(st_list, feat_list = features)

  all_features <- Reduce(intersect, lapply(st_list, feat_names))

  if(length(all_features) == 0){
      print_this_msg("No shared features between objects...",
                     msg_type = "STOP")
  }

  if (is.null(names)) {
    names <- paste("Condition_", 1:length(st_list), sep = "")
  } else{
    if (length(names) != length(st_list))
      print_this_msg("The number of names should be same as the number of objects.",
                     msg_type = "STOP")
  }


  print_this_msg("Naming objects.", msg_type = "DEBUG")

  names(st_list) <- names

  print_this_msg("Getting 'coords' slot.", msg_type = "DEBUG")

  st_list <- lapply(
    st_list,
    coord
  )

  print_this_msg("Subsetting STGrid objects.", msg_type = "DEBUG")

  st_list <- lapply(
    st_list,
    "[[",
    "feature"
  )

  print_this_msg("Counting...", msg_type = "DEBUG")

  st_list <- lapply(
    st_list,
    table
  )

  print_this_msg("Merging...", msg_type = "DEBUG")

  st_list <- do.call("cbind", st_list)

  if(normalized){
    exp_design <- data.frame(condition=as.factor(rep("condition_1", ncol(st_list))),
                             row.names = colnames(st_list))

    dds0 <- DESeq2::DESeqDataSetFromMatrix(countData = st_list,
                                           colData = exp_design,
                                           design = ~1)

    dds.norm <-  DESeq2::estimateSizeFactors(dds0)
    st_list <- DESeq2::counts(dds.norm, normalized=TRUE)

  }

  st_list <- st_list[features, , drop=FALSE]

  count_per_gene <- reshape2::melt(as.matrix(st_list))

  colnames(count_per_gene) <- c("Gene", "Conditions", "value")
  count_per_gene$Conditions <- factor(count_per_gene$Conditions,
                                      levels=names, ordered = TRUE)

  count_per_gene$Gene <- factor(count_per_gene$Gene,
                                      levels=features, ordered = TRUE)

  if(!is.null(fill_color)){
    if(length(fill_color) < length(unique(count_per_gene$Conditions))){
      print_this_msg("Not enough color provided. Using default palette...", msg_type = "WARNING")
      fill_color <- NULL
    }
  }

  if(transform == "log2"){
    y_label <- "Log2(counts)"
    count_per_gene$value <- log2(count_per_gene$value)
  }else if(transform ==  "log10"){
    y_label <- "Log10(counts)"
    count_per_gene$value <- log10(count_per_gene$value)
  }else if(transform ==  "log"){
    y_label <- "Log(counts)"
    count_per_gene$value <- log(count_per_gene$value)
  }else{
    y_label <- "Counts"
  }

  Gene <- value <- Conditions <- NULL

  if(type=="barplot"){

    p <- ggplot2::ggplot(data = count_per_gene,
                         mapping = ggplot2::aes(x = Gene,
                                                y = value,
                                                fill = Conditions)) +
      ggplot2::geom_col(position="dodge2")

  }else if(type=="radar"){
      print_this_msg("Not implemented yet. Sry", msg="STOP")
  }

  if(!is.null(fill_color))
    p <- p + ggplot2::scale_fill_manual(values=fill_color)

  p <- p + ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 10,
        angle = 45,
        vjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )

    p <- p + ggplot2::ylab(y_label) +
      ggplot2::xlab("Conditions")

  return(p)

}

# -------------------------------------------------------------------------
##    Compare count distributions across multiple st_grid_objects
# -------------------------------------------------------------------------
#' @title Compare distribution of features across various STGrid objects.
#'
#' @description
#'  This function compares distributions of features across various STGrid objects using different plot types, including histograms, density plots, boxplots, or boxjitter plots.
#'
#' @param ... One or more objects of class "STGrid" to be compared.
#' @param normalized Logical, indicating whether counts should be normalized. Defaults to FALSE.
#' @param names A character vector specifying the names for each object. Defaults to NULL, where default names are assigned.
#' @param type The type of plot to be generated. Can be one of "hist" (histogram), "density" (density plot), "boxplot", or "boxjitter". Defaults to "hist".
#' @param transform The transformation to be applied to the counts. Can be one of "None", "log2", "log10", or "log". Defaults to "None".
#' @param border_color The color of the border for plot elements. Defaults to "#333333".
#' @param fill_color A character vector specifying the fill colors for the plot. If NULL, default palette is used. Defaults to NULL.
#' @param ncol Number of columns for facet wrap. Defaults to 2.
#' @importFrom ggplot2 scale_fill_manual aes facet_wrap geom_boxplot geom_density geom_histogram scale_fill_manual theme_bw theme element_text element_blank
#' @importFrom ggpol geom_boxjitter
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[180:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[100:nbin_y(xen)]
#' xen_r1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[60:100]
#' y_bins <-  bin_y(xen)[100:nbin_y(xen)]
#' xen_r2 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[20:60]
#' y_bins <-  bin_y(xen)[20:60]
#' xen_r3 <- xen[x_bins, y_bins]
#' dist_st(xen_r1, xen_r2, xen_r3,
#'         fill_color=c("red", "black", "green"),
#'         type="density",
#'         transform="log10")
#' dist_st(xen_r1, xen_r2, xen_r3,
#'         fill_color=c("red", "black", "green"),
#'         type="boxjitter",
#'         transform="log2")
#' @export
dist_st <- function(...,
                    normalized = FALSE,
                    names=NULL,
                    type=c("boxjitter", "hist", "density", "boxplot"),
                    transform = c("None", "log2", "log10", "log"),
                    border_color = "#333333",
                    fill_color = NULL,
                    ncol=2) {

  print_this_msg("Checking STGrid objects", msg_type = "DEBUG")

  type <- match.arg(type)
  transform <- match.arg(transform)
  st_list <- list(...)

  if(is.list(st_list[[1]]))
    st_list <- st_list[[1]]

  check_st_list(st_list, feat_list = NULL)

  if(length(border_color) > 1)
    border_color <- border_color[1]


    features <- Reduce(intersect, lapply(st_list, feat_names))
    if(length(features) == 0)
      print_this_msg("No shared features between objects...",
                     msg_type = "STOP")


  check_st_list(st_list)

  if (is.null(names)) {
    names <- paste("Condition_", 1:length(st_list), sep = "")
  } else{
    if (length(names) != length(st_list))
      print_this_msg("The number of names should be same as the number of objects.",
                     msg_type = "STOP")
  }

  print_this_msg("Subsetting STGrid objects.", msg_type = "DEBUG")

  names(st_list) <- names

  get_feat <- function(x){coord(x)$feature}

  st_list <- lapply(
    st_list,
    get_feat
  )

  st_list <- lapply(
    st_list,
    table
  )

  st_list <- do.call("cbind", st_list)

  if(normalized){
    if(ncol(st_list) > 1){
      print_this_msg("Normalizing...")
      exp_design <- data.frame(condition=as.factor(rep("condition_1", ncol(st_list))),
                               row.names = colnames(st_list))

      dds0 <- DESeq2::DESeqDataSetFromMatrix(countData = st_list,
                                             colData = exp_design,
                                             design = ~1)

      dds.norm <-  DESeq2::estimateSizeFactors(dds0)
      st_list <- DESeq2::counts(dds.norm, normalized=TRUE)

    }

  }
  count_per_gene <- reshape2::melt(st_list)

  colnames(count_per_gene) <- c("Gene", "Conditions", "value")

  if(!is.null(fill_color)){
    if(length(fill_color) < length(unique(count_per_gene$Conditions))){
      print_this_msg("Not enough color provided. Using default palette...", msg_type = "WARNING")
      fill_color <- NULL
    }
  }

  if(transform == "log2"){
    x_label <- "Log2(counts)"
    count_per_gene$value <- log2(count_per_gene$value)
  }else if(transform ==  "log10"){
    x_label <- "Log10(counts)"
    count_per_gene$value <- log10(count_per_gene$value)
  }else if(transform ==  "log"){
    x_label <- "Log(counts)"
    count_per_gene$value <- log(count_per_gene$value)
  }else{
    x_label <- "Counts"
  }

  if(type=="hist"){

      p <- ggplot2::ggplot(data = count_per_gene,
                           mapping = ggplot2::aes(x = value,
                                                  fill=Conditions)) +
        ggplot2::geom_histogram(color = border_color)

  }else if(type=="boxplot"){
    p <- ggplot2::ggplot(data = count_per_gene,
                         mapping = ggplot2::aes(x = Conditions,
                                                y=value,
                                                fill=Conditions)) +
      ggplot2::geom_boxplot(color = border_color)

  }else if(type=="density"){
    p <- ggplot2::ggplot(data = count_per_gene,
                         mapping = ggplot2::aes(x = value,
                                                fill=Conditions)) +
      ggplot2::geom_density(color = border_color)


  Conditions <- value <- NULL

  }else if(type=="boxjitter"){
    p <- ggplot2::ggplot(data = count_per_gene,
                         mapping = ggplot2::aes(x = Conditions,
                                                y = value,
                                                fill = Conditions)) +
      ggpol::geom_boxjitter(color = border_color)
  }

  if(!is.null(fill_color))
    p <- p + ggplot2::scale_fill_manual(values=fill_color)

  p <- p + ggplot2::theme_bw() +
       ggplot2::theme(
        axis.text.x = ggplot2::element_text(
        size = 10,
        angle = 45,
        vjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )


  if(type %in% c("hist", "density")){
    p <- p + ggplot2::facet_wrap(~Conditions,
                                 ncol = ncol) +
      ggplot2::ylab("Number") +
      ggplot2::xlab(x_label)
  }else{
    p <- p + ggplot2::ylab(x_label) +
      ggplot2::xlab("Conditions")
  }

  return(p)

}

# -------------------------------------------------------------------------
##    Boxplot / jitter
# -------------------------------------------------------------------------
#' @title Create a boxplot/jitter plot to show molecule counts distribution.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#'
#' @param object A STCompR object.
#' @param normalized A logical value indicating whether to use normalized counts.
#' @param transform A character string specifying the transformation to be applied to molecule counts.
#'                  Options include "None", "log2", "log10", and "log".
#' @param colors A vector of colors for the boxplot/jitter plot.
#' @param ... Additional arguments to be passed to the underlying ggpol::geom_boxjitter function.
#'
#' @return A ggplot object displaying the molecule counts distribution.
#' @examples
#' example_dataset("11210787/files/cmp_xen")
#' cmp_xen
#' cmp_boxplot(cmp_xen, normalized = TRUE, transform = "log2", colors = c("blue", "red"))
#'
#' @keywords internal
setGeneric("cmp_boxplot",
           function(object,
                    normalized = TRUE,
                    transform = c("None", "log2", "log10", "log"),
                    colors = c("#3074BB", "#BE5B52"),
                    ...)
             standardGeneric("cmp_boxplot"))

#' @title Create a boxplot/jitter plot to show molecule counts distribution.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#'
#' @param object A STCompR object.
#' @param normalized A logical value indicating whether to use normalized counts.
#' @param transform A character string specifying the transformation to be applied to molecule counts.
#'                  Options include "None", "log2", "log10", and "log".
#' @param colors A vector of colors for the boxplot/jitter plot.
#' @param ... Additional arguments to be passed to the underlying ggpol::geom_boxjitter function.
#'
#' @return A ggplot object displaying the molecule counts distribution.
#' @export cmp_boxplot
#' @examples
#' example_dataset("11210787/files/cmp_xen")
#' cmp_xen
#' cmp_boxplot(cmp_xen, normalized = TRUE, transform = "log2", colors = c("blue", "red"))
#' @importFrom ggplot2 ggplot aes theme_bw ylab scale_fill_manual
#' @importFrom ggsci pal_npg
#' @importFrom ggpol geom_boxjitter
setMethod("cmp_boxplot", signature("STCompR"),
          function(object,
                   normalized = TRUE,
                   transform = c("None", "log2", "log10", "log"),
                   colors = c("#3074BB", "#BE5B52"),
                   ...) {
            transform <- match.arg(transform)

            counts <- stat_test(
              object,
              normalized = normalized,
              count_only = TRUE,
              melted_count = TRUE,
              transform = transform
            )

            ylabel <- ifelse(
              transform %in% c("log2", "log10", "log"),
              paste0(transform, "(Molecule counts)"),
              "Molecule counts"
            )

            Conditions <- Counts <- Conditions <- NULL

            ggplot2::ggplot(
              data = counts,
              mapping = ggplot2::aes(x = Conditions, y = Counts, fill = Conditions)
            ) +
              ggpol::geom_boxjitter(...) +
              ggplot2::theme_bw() +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(
                  size = 10,
                  angle = 45,
                  vjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(size = 12),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank(),
                panel.grid.minor.y = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank()
              ) +
              ggplot2::ylab(ylabel) +
              ggplot2::scale_fill_manual(values = colors)

          })


# -------------------------------------------------------------------------
##    Volcano
# -------------------------------------------------------------------------
#' @title Create a volcano plot to compare molecule counts between 2 conditions.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#' @param object A STCompR object.
#' @param y_axis Character vector specifying the y-axis variable.
#' Options are "p_values" or "padj".
#' @param x_lim Numeric vector specifying the x-axis limits. Default is c(-2.5, 2.5).
#' @param colors Color palette for the plot. Default is RColorBrewer::brewer.pal(7, "Spectral").
#' @param text_y_lim Numeric value specifying the threshold for feature labels on the y-axis. Genes with -log10(p-value) less than this threshold will not be labeled. Default is 100.
#' @param text_x_lim Numeric value specifying the threshold for feature labels on the x-axis. Genes with absolute x-values less than this threshold will not be labeled. Default is 2.
#' @param text_size Numeric value specifying the size of text labels in the plot. Default is 5.
#' @param title A title for the diagram.
#' @param max.overlaps The maximum number of label overlaps.
#' @keywords internal
#' @examples
#' example_dataset("11210787/files/cmp_xen")
#' cmp_xen
#' cmp_volcano(cmp_xen)
#' @keywords internal
setGeneric("cmp_volcano",
           function(object,
                    y_axis = c("p_values",
                               "padj"),
                    x_lim = c(-2.5, 2.5),
                    colors = rev(RColorBrewer::brewer.pal(7, "Spectral")),
                    text_y_lim = 100,
                    text_x_lim = 1,
                    text_size = 4,
                    title=NULL,
                    max.overlaps=10)
           standardGeneric("cmp_volcano"))

#' @title Create a volcano plot to compare molecule counts between 2 conditions.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#' @param object A STCompR object.
#' @param y_axis Character vector specifying the y-axis variable.
#' Options are "p_values" or "padj".
#' @param x_lim Numeric vector specifying the x-axis limits. Default is c(-2.5, 2.5).
#' @param colors Color palette for the plot. Default is RColorBrewer::brewer.pal(7, "Spectral").
#' @param text_y_lim Numeric value specifying the threshold for feature labels on the y-axis. Genes with -log10(p-value) less than this threshold will not be labeled. Default is 100.
#' @param text_x_lim Numeric value specifying the threshold for feature labels on the x-axis. Genes with absolute x-values less than this threshold will not be labeled. Default is 2.
#' @param text_size Numeric value specifying the size of text labels in the plot. Default is 5.
#' @param title A title for the diagram.
#' @param max.overlaps The maximum number of label overlaps.
#' @keywords internal
#' @examples
#' example_dataset("11210787/files/cmp_xen")
#' cmp_xen
#' cmp_volcano(cmp_xen)
#' @importFrom ggplot2 aes geom_vline geom_hline geom_point theme_bw xlab ylab expand_limits scale_fill_gradientn ggtitle
#' @importFrom  ggrepel geom_text_repel
#' @importFrom stringr str_to_title
#' @export cmp_volcano
setMethod("cmp_volcano", signature("STCompR"),
          function(object,
                   y_axis = c("p_values",
                              "padj"),
                   x_lim = c(-2.5, 2.5),
                   colors = rev(RColorBrewer::brewer.pal(7, "Spectral")),
                   text_y_lim = 100,
                   text_x_lim = 1,
                   text_size = 4,
                   title=NULL,
                   max.overlaps=10) {

            y_axis <- match.arg(y_axis)

            if (y_axis != "p_values")
              y_axis <- paste0("padj")

            uniq_conditions <- unique(object@conditions)
            x_axis_lab <- paste0("log2(",
                                 uniq_conditions[2],
                                   "/",
                                 uniq_conditions[1],
                                   ")")


            if(is.null(title)){
              title <- paste0("Comparison of ",
                              object@conditions[1],
                              " vs ",
                              object@conditions[2])
            }

            volc_data <- stat_test(
              object,
              count_only = FALSE,
              melted_count = FALSE,
              features = NULL
            )

            volc_data <- volc_data[, c("log2_ratio", y_axis)]

            colnames(volc_data) <- c("x", "y")

            counts <- stat_test(
              object,
              count_only = TRUE,
              melted_count = FALSE,
              features = NULL
            )

            volc_data$mean_counts <- rowMeans(counts)
            volc_data$feature <- rownames(volc_data)
            volc_data$feature[-log10(volc_data$y) < text_y_lim] <- NA
            volc_data$feature[abs(volc_data$x) < text_x_lim] <- NA

            x <- y <- mean_counts <- feature <- NULL


            ggplot2::ggplot(
              data = volc_data,
              mapping = ggplot2::aes(
                x = x,
                y = -log10(y),
                fill = x,
                size = mean_counts
              )
            ) +
              ggplot2::geom_vline(xintercept = 0) +
              ggplot2::geom_hline(yintercept = 0) +
              ggplot2::geom_point(shape = 21,
                                  color = "black",
                                  stroke = 0.2) +
              ggplot2::theme_bw() +
              ggrepel::geom_text_repel(
                data = stats::na.omit(volc_data),
                mapping = ggplot2::aes(label = feature),
                size = text_size,
                color = "black",
                max.overlaps=max.overlaps
              ) +
              ggplot2::ggtitle(title) +
              ggplot2::xlab(x_axis_lab) +
              ggplot2::theme(
                axis.text.x = ggplot2::element_text(size = 10,  vjust = 0.5),
                axis.text.y = ggplot2::element_text(size = 8),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor.x = ggplot2::element_blank(),
                panel.grid.minor.y = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank()
              ) +
              ggplot2::ylab(paste0("-log10(", y_axis, ")")) +
              ggplot2::expand_limits(x = x_lim) +
              ggplot2::scale_fill_gradientn(colors = colors,
                                            name = "Log2 ratio")

})

# -------------------------------------------------------------------------
##    Ripley's K function
# -------------------------------------------------------------------------
#' @title Plot the results from the Ripley's k function stored in the STGrid object.
#' @description
#' Plot the results from the Ripley's k function stored in the STGrid object.
#' @param object A STCompR object.
#' @param correction kind of correction should be displayed.
#' @param max_feat_label The maximum number of features to display.
#' The features with the highest value (whatever the radius)
#' are displayed
#' @param color The colors for the features to be displayed.
#' @param size The size of the labels.
#' @keywords internal
setGeneric("plot_rip_k",
           function(object,
                    correction = c("border",
                                   "isotropic",
                                   "Ripley",
                                   "translate"),
                    max_feat_label = 8,
                    color = NULL,
                    size = 4)
             standardGeneric("plot_rip_k"))

#' @title Call the Ripley's k function.
#' @description
#' Call the Ripley's k function.
#' @param object A STCompR object.
#' @param correction kind of correction should be displayed.
#' @param max_feat_label The maximum number of features to display.
#' The features with the highest correction value (whatever the radius)
#' are displayed.
#' @param color The colors for the features to be displayed.
#' @param size The size of the labels.
#' @import magrittr
#' @importFrom dplyr group_by filter arrange
#' @examples
#' example_dataset()
#' plot_rip_k(compute_k_ripley(Xenium_Mouse_Brain_Coronal_7g[c("Ano1", "Chat", "Ebf3")]))
#'
#' @export plot_rip_k
setMethod("plot_rip_k", signature("STGrid"),
          function(object,
                   correction = c("border",
                                  "isotropic",
                                  "Ripley",
                                  "translate"),
                   max_feat_label = 8,
                   color = NULL,
                   size = 4) {
            gg_color_hue <- function(n) {
              hues <-  seq(15, 375, length = n + 1)
              grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
            }


            if (is.null(color)) {
              color <- gg_color_hue(max_feat_label)
            }


            if (nrow(ripley_k_function(object)) == 0) {
              print_this_msg("Please run ripley_k_function() first.", msg_type = "STOP")
            }

            correction <- match.arg(correction)
            ripk <- ripley_k_function(object)

            if (max_feat_label > nrow(ripk))
              print_this_msg("Too much selected gehes...", msg_type = "STOP")

            voi <- ripk %>%
              dplyr::group_by(feature) %>%
              dplyr::filter(border == max(border)) %>%
              dplyr::filter(r == max(r)) %>%
              dplyr::arrange(dplyr::desc(border)) %>%
              utils::head(n = max_feat_label)

            voi <- voi[!duplicated(voi$feature), ]

            goi <- unique(voi$feature)

            ripk_sub <- ripk[ripk$feature %in% goi,]

            r <- border <- feature <- NULL
            p <- ggplot2::ggplot(data = ripk) +
              ggplot2::geom_line(mapping = ggplot2::aes(x = r,
                                                        y = border,
                                                        group = feature),
                                 color = "black") +
              ggplot2::theme_bw() +
              ggplot2::theme(
                legend.text = ggplot2::element_text(size = 4),
                legend.position = "none",
                panel.grid.minor =  ggplot2::element_blank()
              ) +
              ggplot2::geom_line(
                data = ripk_sub,
                mapping = ggplot2::aes(
                  x = r,
                  y = border,
                  group = feature,
                  color = feature
                ),
                inherit.aes = FALSE
              ) +
              ggrepel::geom_label_repel(
                data = voi,
                mapping = ggplot2::aes(
                  x = r,
                  y = border,
                  label = feature,
                  color = feature
                ),
                inherit.aes = FALSE,
                size = size,
                force = 20
              ) +
              ggplot2::ylab(paste0("Ripley's K function (correction=", correction, ")")) +
              ggplot2::scale_color_manual(values = color)

            return(p)
          })


# -------------------------------------------------------------------------
##    Compare Images
# -------------------------------------------------------------------------
#' Compare Images
#'
#' This function compares and visualizes multiple spatial images obtained from STGrid objects. It displays a color-coded representation of the feature (e.g., molecules) density observed in spatial transcriptomics experiments for different conditions.
#'
#' @param ... STGrid objects to be compared.
#' @param feat_list A list of features for which the spatial images will be created.
#' @param names Names for the conditions or experiments. Defaults to NULL.
#' @param colors The colors to use for gradient fill in the spatial images. Defaults to viridis::inferno(10).
#' @param saturation The ceiling level for the feature expression values. Defaults to 1 (no ceiling).
#' @param coord_fixed Logical value indicating whether to keep the aspect ratio fixed. Defaults to TRUE.
#' @param scale Logical value indicating whether to scale the feature expression values. Defaults to TRUE.
#' @param logb The basis for the log transformation. Default to 10. If NULL, no log transformation.
#' @param pseudo_count A value for the pseudo count used for log transformation (default to 1).
#' @param condition_vs_feat Logical indicating whether to facet by condition vs. feature (TRUE) or by feature vs. condition (FALSE). Defaults to TRUE.
#'
#' @importFrom ggplot2 aes geom_tile scale_fill_gradientn theme xlab ylab element_blank element_rect element_text
#' @importFrom ggh4x facet_grid2
#' @importFrom grid unit
#' @examples
#' example_dataset()
#' xen <- Xenium_Mouse_Brain_Coronal_7g
#' x_bins <-  bin_x(xen)[181:nbin_x(xen)]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_1 <- xen[x_bins, y_bins]
#' x_bins <-  bin_x(xen)[61:101]
#' y_bins <-  bin_y(xen)[101:nbin_y(xen)]
#' xen_2 <- xen[x_bins, y_bins]
#' cmp_images(xen_1, xen_2, feat_list = c("Nwd2", "Kctd8", "Necab2", "Nrp2"))
#' @export
cmp_images <- function(...,
                       feat_list = NULL,
                       names = NULL,
                       colors = c("black", "#33FF00", "#FFFF00",
                                  "#FF0000", "#CC00FF"),
                       saturation = 1,
                       coord_fixed = TRUE,
                       scale = TRUE,
                       logb = 10,
                       pseudo_count = 1,
                       condition_vs_feat = TRUE) {

  if (saturation > 1 | saturation < 0)
    print_this_msg("Saturation should be between 0 and 1.",
                   msg_type = "STOP")

  if (is.null(feat_list))
    print_this_msg("Please provide a feature list.", msg_type = "STOP")

  print_this_msg("Checking STGrid objects", msg_type = "DEBUG")

  st_list <- list(...)
  if(is.list(st_list[[1]]))
    st_list <- st_list[[1]]

  check_st_list(st_list, feat_list)

  if (is.null(names)) {
    names <- paste("Condition_", 1:length(st_list), sep = "")
  } else{
    if (length(names) != length(st_list))
      print_this_msg("The number of names should be same as the number of objects.",
                     msg_type = "STOP")
  }


  print_this_msg("Subsetting STGrid objects.", msg_type = "DEBUG")

  # Save STGrid version
  st_list_grid <- st_list

  st_list <- lapply(
    st_list,
    bin_mat,
    melt_tab = TRUE,
    as_factor = TRUE,
    feat_list = feat_list
  )

  # Lot of manipulation to try to fix a weird bug...
  for (i in 1:length(st_list)) {
    st_list[[i]]$bin_x <- factor(st_list[[i]]$bin_x,
                                 levels=bin_x(st_list_grid[[i]]),
                                 ordered = TRUE)
    st_list[[i]]$bin_x <- as.numeric(st_list[[i]]$bin_x)
    st_list[[i]]$bin_x <- as.factor(st_list[[i]]$bin_x)

    st_list[[i]]$bin_y <- factor(st_list[[i]]$bin_y,
                                 levels=bin_y(st_list_grid[[i]]),
                                 ordered = TRUE)
    st_list[[i]]$bin_y <- as.numeric(st_list[[i]]$bin_y)
    st_list[[i]]$bin_y <- as.factor(st_list[[i]]$bin_y)
  }

  print_this_msg("Preparing data.", msg_type = "DEBUG")

  for (i in 1:length(st_list)) {
    for (j in feat_list) {
      tmp <- st_list[[i]]$value[st_list[[i]]$feature == j]

      if (saturation < 1) {
        print_this_msg("Ceiling.", msg_type = "DEBUG")
        q_sat <- stats::quantile(tmp[tmp != 0], saturation)
        tmp[tmp > q_sat] <- q_sat
      }

      if (!is.null(logb)) {
        print_this_msg("Transforming in log base ", logb, ".", msg_type = "DEBUG")
        tmp <- log(tmp + pseudo_count, base = logb)
      }

      if (scale) {
        print_this_msg("Rescaling", msg_type = "DEBUG")
        tmp <- (tmp - min(tmp)) / (max(tmp) - min(tmp))
      }

      st_list[[i]]$value[st_list[[i]]$feature == j] <- tmp

    }

  }

  for (i in 1:length(st_list)) {
    st_list[[i]]$condition <- names[i]
  }


  print_this_msg("Merging data.", msg_type = "DEBUG")

  st_list <- do.call(rbind, st_list)

  print_this_msg("Converting columns 'condition' to ordered factor.", msg_type = "DEBUG")

  st_list$condition <- factor(st_list$condition,
                              levels = names,
                              ordered = TRUE)

  print_this_msg("Converting columns 'gene' to ordered factor.", msg_type = "DEBUG")

  st_list$gene <- factor(st_list$feature,
                         levels = feat_list,
                         ordered = TRUE)

  print_this_msg("Building diagram", msg_type = "DEBUG")

  value <- bin_x <- bin_y <- NULL

  st_list$bin_x <- factor(st_list$bin_x,
                          levels = levels(st_list$bin_x),
                          ordered = TRUE)

  st_list$bin_y <- factor(st_list$bin_y,
                          levels = levels(st_list$bin_y),
                          ordered = TRUE)

  st_list$feature <- factor(st_list$feature,
                            levels = feat_list,
                            ordered = TRUE)

  st_list$condition <- factor(st_list$condition,
                            levels = names,
                            ordered = TRUE)

  p <- ggplot2::ggplot(data = st_list,
                       mapping = ggplot2::aes(x = bin_x,
                                              y = bin_y,
                                              fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::xlab("")  +
    ggplot2::ylab("") +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "gray30"),
      strip.text = ggplot2::element_text(color = "white")
    ) +
    ggplot2::scale_fill_gradientn(colours = colors)

  if (condition_vs_feat) {
    p <-
      p + ggh4x::facet_grid2(condition ~ feature,
                             scale = "free",
                             independent = "x") +
      ggplot2::theme(panel.spacing = grid::unit(0.25, "lines"))
  } else{
    p <-
      p + ggh4x::facet_grid2(feature ~ condition,
                             scale = "free",
                             independent = "y") +
      ggplot2::theme(panel.spacing = grid::unit(0.25, "lines"))
  }

  p
}

