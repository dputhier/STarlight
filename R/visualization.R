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
#' @param size The size of the overlayed points.
#' @param logb The basis for the log transformation. Default to 10. If NULL no log transformation.
#' @param pseudo_count a value for the pseudo count used for log transformation (default to 1).
#' @keywords internal
#' @export
setGeneric("spatial_image",
           function(object=NULL,
                    features=NULL,
                    saturation=1,
                    scale=TRUE,
                    colors=viridis::inferno(10),
                    coord_fixed=TRUE,
                    overlay_feature=NULL,
                    colors_overlay=c("#DEEBF7", "#9ECAE1", "#3182BD"),
                    grid_by=NULL,
                    color_grid="white",
                    size=0.5,
                    logb=10,
                    pseudo_count=1)
             standardGeneric("spatial_image")
)

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
#' @param size The size of the overlayed points.
#' @param logb The basis for the log transformation. Default to 10. If NULL no log transformation.
#' @param pseudo_count a value for the pseudo count used for log transformation (default to 1).
#' @importFrom ggplot2 aes coord_fixed facet_wrap geom_tile scale_fill_gradientn theme xlab ylab element_blank element_rect element_text
#' @importFrom reshape2 melt
#' @importFrom viridis inferno
#' @export
setMethod("spatial_image",
          signature(object = "STGrid"),
          function(object=NULL,
                   features=NULL,
                   saturation=1,
                   scale=TRUE,
                   colors=viridis::inferno(10),
                   coord_fixed=TRUE,
                   overlay_feature=NULL,
                   colors_overlay=c("#DEEBF7", "#9ECAE1", "#3182BD"),
                   grid_by=NULL,
                   color_grid="white",
                   size=0.5,
                   logb=10,
                   pseudo_count=1) {

            print_msg("Checking arguments", msg_type = "DEBUG")

            if(is.null(object))
                print_msg("Please provide an STGrid object.",
                        msg_type = "STOP")

            if(!is.null(overlay_feature)){
              if(!overlay_feature %in% feat_names(object))
                print_msg("The feature to overlay was not found in the object.",
                          msg_type = "STOP")

            }

            print_msg("Checking features", msg_type = "DEBUG")

            if(is.null(features))
              print_msg("Please provide a feature name (see feature arguments).", msg_type = "STOP")


            if(!all(features %in% c(feat_names(object), "sum_of_cts"))){
              print_msg("The feature was not found in the object.", msg_type = "STOP")
            }

            if("sum_of_cts" %in% features){

              print_msg("Using sum_of_cts", msg_type = "DEBUG")

              spatial_matrix <- object@bin_mat

              if(length(features) == 1){
                sub_feat <- feat_names(object)
                spatial_matrix <- spatial_matrix[, c("bin_x", "bin_y", sub_feat)]
                tmp <- spatial_matrix[, sub_feat]
                tmp$sum_of_cts <- rowSums(tmp)
                spatial_matrix$sum_of_cts <- tmp$sum_of_cts
                spatial_matrix <- spatial_matrix[, c("bin_x", "bin_y","sum_of_cts")]
                tmp <- tmp[, "sum_of_cts", drop=FALSE]
              }else{
                sub_feat <- setdiff(features, "sum_of_cts")
                spatial_matrix <- spatial_matrix[, c("bin_x", "bin_y", sub_feat)]
                tmp <- spatial_matrix[, sub_feat, drop=FALSE]
                tmp$sum_of_cts <- rowSums(tmp)
                spatial_matrix$sum_of_cts <- tmp$sum_of_cts
              }

            }else{
              spatial_matrix <- object@bin_mat[, c("bin_x", "bin_y", features)]
              tmp <- spatial_matrix[, features, drop=FALSE]
            }

            for(i in 1:ncol(tmp)){
              if(saturation < 1){
                q_sat <- quantile(tmp[,i][tmp[,i] != 0], saturation)
                tmp[tmp[,i] > q_sat, i] <- q_sat
              }

              if(!is.null(logb)){
                tmp[,i] <- log(tmp[,i] + pseudo_count, base = logb)
              }

              if(scale){
                tmp[,i] <- (tmp[,i] - min(tmp[,i])) / (max(tmp[,i]) - min(tmp[,i]))
              }
            }


            spatial_matrix[, features] <- tmp

            spatial_matrix$bin_x <- factor(spatial_matrix$bin_x,
                                           levels = bin_x(object),
                                           ordered=TRUE)

            spatial_matrix$bin_y <- factor(spatial_matrix$bin_y,
                                           levels = bin_y(object),
                                           ordered=TRUE)

            spatial_matrix_melted <- reshape2::melt(spatial_matrix, id.vars=c("bin_x", "bin_y"))
            spatial_matrix_melted$variable <- factor(spatial_matrix_melted$variable,
                                                     levels = features,
                                                     ordered=TRUE)

            p <- ggplot2::ggplot(data=spatial_matrix_melted,
                                 mapping = ggplot2::aes(x=bin_x, y=bin_y, fill= value)) +
              ggplot2::geom_tile() +
              ggplot2::xlab("")  +
              ggplot2::ylab("") +
              ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                             axis.text = ggplot2::element_blank(),
                             strip.background = ggplot2::element_rect(fill="gray30"),
                             strip.text = ggplot2::element_text(color="white")) +
              ggplot2::scale_fill_gradientn(colours=colors) +
              ggplot2::facet_wrap(~variable, nrow = 3, ncol=3)

            if(coord_fixed)
              p <- p + ggplot2::coord_fixed()


            if(!is.null(overlay_feature)){
              over <- bin_mat(object,
                              as_factor = TRUE)[ ,c("bin_x", "bin_y", overlay_feature)]
              over <- over[over[[overlay_feature]] != 0,]

              p <- p + geom_point(data=over, mapping=aes(x=bin_x,
                                                         y=bin_y,
                                                         color=log10(.data[[overlay_feature]])),
                                  size=size,
                                  inherit.aes = FALSE) +
                scale_color_gradientn(colors=colors_overlay)
            }

            if(!is.null(grid_by)){

              lev_bin_x <- levels(spatial_matrix$bin_x)
              lev_bin_y <- levels(spatial_matrix$bin_y)

              x_seq <- seq(from=1, to=length(lev_bin_x), by=grid_by)
              y_seq <- seq(from=1, to=length(lev_bin_y), by=grid_by)

              label_x <- setNames(lev_bin_x, 1:length(lev_bin_x))
              names(label_x)[-x_seq] <- ""
              label_y <- setNames(lev_bin_y, 1:length(lev_bin_y))
              names(label_y)[-y_seq] <- ""

              p <- p + geom_vline(data=data.frame(bin_x=levels(spatial_matrix$bin_x)[x_seq]),
                                  mapping=aes(xintercept=bin_x), color = color_grid) +
                       geom_hline(data=data.frame(bin_y=levels(spatial_matrix$bin_y)[y_seq]),
                                  mapping=aes(yintercept=bin_y), color = color_grid) +
                  theme(axis.text = element_text(size=6, angle = 45)) +
                  scale_x_discrete("bla", labels=names(label_x)) +
                  scale_y_discrete("bla", labels=names(label_y))
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
#' @keywords internal
setGeneric("spatial_plot",
           function(object=NULL,
                    feat_list=NULL,
                    colors=NULL,
                    size=0.1,
                    coord_fixed=TRUE)
             standardGeneric("spatial_plot")
)

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
setMethod("spatial_plot", "STGrid",
           function(object=NULL,
                    feat_list=NULL,
                    colors=NULL,
                    size=0.1,
                    coord_fixed=TRUE){

             if(is.null(object))
               print_msg("Please provide an STGrid object.",
                         msg_type = "STOP")

             if(is.null(feat_list))
               print_msg("Please provide feature names (see feat_list arguments).",
                         msg_type = "STOP")

             if(!all(feat_list %in% feat_names(object)))
               print_msg("One or several features were not found in the object.", msg_type = "STOP")

             coord <- get_coord(object, feat_list = feat_list, as.factor=TRUE)

             p <- ggplot2::ggplot(data=coord,
                                  mapping = ggplot2::aes(x=x,
                                                         y=y,
                                                         color= feature)) +
               ggplot2::geom_point(size=size) +
               ggplot2::xlab("x")  +
               ggplot2::ylab("y") +
               ggplot2::theme_minimal() +
               theme(axis.text = element_text(size=6))

              if(!is.null(colors))
                p <- p +ggplot2::scale_color_manual(values=colors)

             if(coord_fixed)
               p <- p + ggplot2::coord_fixed()

             return(p)

           }
)

# -------------------------------------------------------------------------
##    Molecule counts
# -------------------------------------------------------------------------
#' @title Create a barplot to show counts for selected features.
#' @description
#' Create a barplot to show counts for selected features.
#' @param object A STGrid object.
#' @param features The list of features (NULL for all of them).
#' @param normalized  Whether counts should be normalized.
#' @param transform Whether the count should be transformed (the pseudo count defined for the object is added).
#' @param colors A set of colors.
#' @keywords internal
setGeneric("cmp_bar_plot",
           function(object,
                    features=head(feat_names(object)),
                    normalized=FALSE,
                    transform=c("None", "log2", "log10", "log"),
                    colors=c("#3074BB", "#BE5B52"))
             standardGeneric("cmp_bar_plot")
)

#' @title Create a barplot to show counts for selected features.
#' @description
#' Create a barplot to show counts for selected features.
#' @param object A STGrid object.
#' @param features The list of features (NULL for all of them).
#' @param normalized  Whether counts should be normalized.
#' @param transform Whether the count should be transformed (the pseudo count defined for the object is added).
#' @param colors A set of colors.
#' @export
setMethod(
  "cmp_bar_plot", signature("STCompR"),
    function(object,
             features=head(feat_names(object)),
             normalized=FALSE,
             transform=c("None", "log2", "log10", "log"),
             colors=c("#3074BB", "#BE5B52")) {


    if(is.null(features)){
      print_msg("Please provide some features...", msg_type = "STOP")
    }

    transform <- match.arg(transform)

    counts <- stat_test(object,
                        normalized=normalized,
                        count_only = TRUE,
                        melted_count = TRUE,
                        transform=transform,
                        features=features)

    ylabel <- ifelse(transform %in% c("log2", "log10", "log"),
                     paste0(transform, "(Molecule counts)"),
                     "Molecule counts")
    counts$Features <- factor(counts$Features,
                              levels=features,
                              ordered = TRUE)

    ggplot2::ggplot(data=counts,
                    mapping = ggplot2::aes(x=Features, y=Counts,
                                           fill=Conditions)) +
    ggplot2::geom_col(color="black", linewidth=0, position="dodge") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, angle=45, vjust = 0.5),
          axis.text.y = ggplot2::element_text(size=8),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          panel.grid.minor.y =ggplot2:: element_blank(),
          panel.border = ggplot2::element_blank()) +
    ggplot2::ylab(ylabel) +
    ggplot2::scale_fill_manual(values = colors)

  }
)

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
#' # Example usage:
#' st_grid_1 <- create_STGrid(...)
#' st_grid_2 <- create_STGrid(...)
#' st_compr_result <- stcompr(st_grid_1, st_grid_2, name_1 = "Condition1", name_2 = "Condition2")
#' cmp_boxplot(st_compr_result, normalized = TRUE, transform = "log2", colors = c("blue", "red"))
#'
#' @keywords internal
setGeneric("cmp_boxplot",
           function(object,
                    normalized=TRUE,
                    transform=c("None", "log2", "log10", "log"),
                    colors=c("#3074BB", "#BE5B52"),
                    ...)
             standardGeneric("cmp_boxplot")
)

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
#'
#' @examples
#' # Example usage:
#' st_grid_1 <- create_STGrid(...)
#' st_grid_2 <- create_STGrid(...)
#' st_compr_result <- stcompr(st_grid_1, st_grid_2, name_1 = "Condition1", name_2 = "Condition2")
#' cmp_boxplot(st_compr_result, normalized = TRUE, transform = "log2", colors = c("blue", "red"))
#'
#' @importFrom ggplot2 ggplot aes theme_bw ylab scale_fill_manual
#' @importFrom ggsci pal_npg
#' @importFrom ggpol geom_boxjitter
setMethod(
  "cmp_boxplot", signature("STCompR"),
  function(object,
           normalized=TRUE,
           transform=c("None", "log2", "log10", "log"),
           colors=c("#3074BB", "#BE5B52"),
           ...) {

    transform <- match.arg(transform)

    counts <- stat_test(object,
                        normalized=normalized,
                        count_only = TRUE,
                        melted_count = TRUE,
                        transform=transform)

    ylabel <- ifelse(transform %in% c("log2", "log10", "log"),
                     paste0(transform, "(Molecule counts)"),
                     "Molecule counts")

    ggplot2::ggplot(data=counts,
           mapping = ggplot2::aes(x=Conditions, y=Counts, fill=Conditions)) +
      ggpol::geom_boxjitter(...) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, angle=45, vjust = 0.5),
            axis.text.y = ggplot2::element_text(size=12),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank()) +
      ggplot2::ylab(ylabel) +
      ggplot2::scale_fill_manual(values = colors)

  }
)


# -------------------------------------------------------------------------
##    Volcano
# -------------------------------------------------------------------------
#' @title Create a volcano plot to compare molecule counts between 2 conditions.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#' @param object A STCompR object.
#' @param x_axis Character vector specifying the x-axis variable.
#' Options are "log2_ratio", "log2_odds", or "odd_ratio".
#' @param y_axis Character vector specifying the y-axis variable.
#' Options are "p_values", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", or "fdr".
#' @param x_lim Numeric vector specifying the x-axis limits. Default is c(-2.5, 2.5).
#' @param colors Color palette for the plot. Default is RColorBrewer::brewer.pal(7, "Spectral").
#' @param text_y_lim Numeric value specifying the threshold for feature labels on the y-axis. Features with -log10(p-value) less than this threshold won't be labeled. Default is 100.
#' @param text_x_lim Numeric value specifying the threshold for feature labels on the x-axis. Features with absolute x-values less than this threshold won't not be labeled. Default is 2.
#' @param text_size Numeric value specifying the size of text labels in the plot. Default is 5.
#' @keywords internal
#' @importFrom ggplot2 aes geom_vline geom_hline geom_point theme_bw xlab ylab expand_limits scale_fill_gradientn
#' @importFrom  ggrepel geom_text_repel
#' @keywords internal
setGeneric("cmp_volcano",
           function(object,
                    x_axis=c("log2_ratio", "log2_odds", "odd_ratio"),
                    y_axis=c("p_values", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
                    x_lim=c(-2.5, 2.5),
                    colors=RColorBrewer::brewer.pal(7, "Spectral"),
                    text_y_lim=100,
                    text_x_lim=2,
                    text_size=5)
             standardGeneric("cmp_volcano")
)

#' @title Create a volcano plot to compare molecule counts between 2 conditions.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#' @param object A STCompR object.
#' @param x_axis Character vector specifying the x-axis variable.
#' Options are "log2_ratio", "log2_odds", or "odd_ratio".
#' @param y_axis Character vector specifying the y-axis variable.
#' Options are "p_values", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", or "fdr".
#' @param x_lim Numeric vector specifying the x-axis limits. Default is c(-2.5, 2.5).
#' @param colors Color palette for the plot. Default is RColorBrewer::brewer.pal(7, "Spectral").
#' @param text_y_lim Numeric value specifying the threshold for feature labels on the y-axis. Genes with -log10(p-value) less than this threshold will not be labeled. Default is 100.
#' @param text_x_lim Numeric value specifying the threshold for feature labels on the x-axis. Genes with absolute x-values less than this threshold will not be labeled. Default is 2.
#' @param text_size Numeric value specifying the size of text labels in the plot. Default is 5.
#' @keywords internal
#' @importFrom ggplot2 aes geom_vline geom_hline geom_point theme_bw xlab ylab expand_limits scale_fill_gradientn
#' @importFrom  ggrepel geom_text_repel
#' @export cmp_volcano
setMethod(
  "cmp_volcano", signature("STCompR"),
  function(object,
             x_axis=c("log2_ratio", "log2_odds", "odd_ratio"),
           y_axis=c("p_values", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"),
           x_lim=c(-2.5, 2.5),
           colors=RColorBrewer::brewer.pal(7, "Spectral"),
           text_y_lim=100,
           text_x_lim=1,
           text_size=4) {

    x_axis <- match.arg(x_axis)
    y_axis <- match.arg(y_axis)

    if(y_axis != "p_values")
      y_axis <- paste0("padj_", y_axis)

    volc_data <- stat_test(object,
                          count_only = FALSE,
                          melted_count = FALSE,
                          features=NULL)

    volc_data <- volc_data[, c(x_axis, y_axis)]
    colnames(volc_data) <- c("x", "y")

    counts <- stat_test(object,
                        count_only = TRUE,
                        melted_count = FALSE,
                        features=NULL)

    volc_data$mean_counts <- rowMeans(counts)
    volc_data$feature <- rownames(volc_data)
    volc_data$feature[-log10(volc_data$y) < text_y_lim] <- NA
    volc_data$feature[abs(volc_data$x) < text_x_lim] <- NA

    ggplot2::ggplot(data=volc_data,
           mapping = ggplot2::aes(x=x, y=-log10(y),
                         fill=x,
                         size=mean_counts)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_point(shape=21,
                          color="black",
                          stroke=0.2) +
      ggplot2::theme_bw() +
      ggrepel::geom_text_repel(data=na.omit(volc_data), mapping=ggplot2::aes(label=feature),
                               size=text_size,
                               color="black") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=10,  vjust = 0.5),
            axis.text.y = ggplot2::element_text(size=8),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank()) +
      ggplot2::xlab(paste0(stringr::str_to_title(gsub("_", " ", x_axis)))) +
      ggplot2::ylab(paste0("-log10(", y_axis, ")")) +
      ggplot2::expand_limits(x = x_lim) +
      ggplot2::scale_fill_gradientn(colors=colors,
                                    name =paste0(stringr::str_to_title(gsub("_", " ", x_axis))))

  }
)

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
                    correction=c("border",
                                 "isotropic",
                                 "Ripley",
                                 "translate"),
                    max_feat_label=8,
                    color=NULL,
                    size=4)
             standardGeneric("plot_rip_k")
)

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
#' @export plot_rip_k
setMethod(
  "plot_rip_k", signature("STGrid"),
  function(object,
           correction=c("border",
                        "isotropic",
                        "Ripley",
                        "translate"),
           max_feat_label=8,
           color=NULL,
           size=4
           ) {

    gg_color_hue <- function(n) {
      hues <-  seq(15, 375, length = n + 1)
      grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
    }


    if(is.null(color)){
      color <- gg_color_hue(max_feat_label)
    }


    if(nrow(ripley_k_function(object)) == 0){
      print_msg("Please run ripley_k_function() first.", msg_type = "STOP")
    }

    correction <- match.arg(correction)
    ripk <- ripley_k_function(object)

    if(max_feat_label > nrow(ripk))
      print_msg("Too much selected gehes...", msg_type = "STOP")

    voi <- ripk %>%
      dplyr::group_by(feature) %>%
      dplyr::filter(border == max(border)) %>%
      dplyr::filter(r == max(r)) %>%
      dplyr::arrange(desc(border)) %>%
          head(n=max_feat_label)

    voi <- voi[!duplicated(voi$feature),]

    goi <- unique(voi$feature)

    ripk_sub <- ripk[ripk$feature %in% goi, ]

    p <- ggplot2::ggplot(data= ripk) +
      ggplot2::geom_line(mapping = ggplot2::aes(x=r, y=border, group=feature), color="black") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.text = ggplot2::element_text(size=4),
                      legend.position = "none",
                      panel.grid.minor =  ggplot2::element_blank()) +
      ggplot2::geom_line(data=ripk_sub,
                 mapping=ggplot2::aes(x=r,
                                      y=border,
                                      group=feature,
                                      color=feature),
                 inherit.aes = FALSE) +
      ggrepel::geom_label_repel(data=voi,
                                mapping=ggplot2::aes(x=r, y=border, label=feature, color=feature),
                                inherit.aes = FALSE,
                                size = size,
                                force=20) +
      ggplot2::ylab(paste0("Ripley's K function (correction=", correction, ")")) +
      scale_color_manual(values=color)

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
#' @export
cmp_images <- function(...,
                       feat_list=NULL,
                       names=NULL,
                       colors=viridis::inferno(10),
                       saturation=1,
                       coord_fixed=TRUE,
                       scale=TRUE,
                       logb=10,
                       pseudo_count=1,
                       condition_vs_feat=TRUE){

  if(is.null(feat_list))
    print_msg("Please provide a feature list.", msg_type = "STOP")

  print_msg("Checking STGrid objects", msg_type = "DEBUG")

  st_list <- list(...)
  if(any(unlist(lapply(lapply(st_list, class), "[", 1)) != "STGrid")){
    print_msg("Object should be of type STGrid", msg_type = "STOP")
  }

  if(length(st_list) < 1){
    print_msg("Need at least one experiment !!", msg_type = "STOP")
  }

  if(is.null(names)){
    names <- paste("Condition_", 1:length(st_list), sep="")
  }else{
    if(length(names) != length(st_list))
      print_msg("The number of names should be same as the number of objects.",
                msg_type = "STOP")
  }

  print_msg("Subsetting STGrid objects.", msg_type = "DEBUG")

  for(i in 1:length(st_list)){
    st_list[[i]] <- st_list[[i]][feat_list, ]
  }

  st_list <- lapply(st_list,
                    bin_mat,
                    melt_tab = TRUE,
                    as_factor = TRUE,
                    feat_list =feat_list)

  print_msg("Preparing data.", msg_type = "DEBUG")

  for(i in 1:length(st_list)){
      for(j in feat_list){

          tmp <- st_list[[i]]$value[st_list[[i]]$feature == j]

          if(saturation < 1){
            print_msg("Ceiling.", msg_type = "DEBUG")
            q_sat <- quantile(tmp[tmp != 0], saturation)
            tmp[tmp > q_sat] <- q_sat
          }

          if(!is.null(logb)){
            print_msg("Transforming in log base ", logb, ".", msg_type = "DEBUG")
            tmp <- log(tmp + pseudo_count, base = logb)
          }

          if(scale){
            print_msg("Rescaling", msg_type = "DEBUG")
            tmp <- (tmp - min(tmp)) / (max(tmp) - min(tmp))
          }

          st_list[[i]]$value[st_list[[i]]$feature == j] <- (tmp - min(tmp)) / (max(tmp) - min(tmp))

      }

  }


  for(i in 1:length(st_list)){
    st_list[[i]]$condition <- names[i]
  }

  print_msg("Merging data.", msg_type = "DEBUG")

  st_list <- do.call(rbind, st_list)

  print_msg("Converting columns 'condition' to ordered factor.", msg_type = "DEBUG")

  st_list$condition <- factor(st_list$condition,
                              levels=names,
                              ordered = TRUE)

  print_msg("Converting columns 'gene' to ordered factor.", msg_type = "DEBUG")

  st_list$gene <- factor(st_list$feature,
                              levels=feat_list,
                              ordered = TRUE)

  print_msg("Building diagram", msg_type = "DEBUG")

  p <- ggplot2::ggplot(data=st_list,
                       mapping = ggplot2::aes(x=bin_x,
                                              y=bin_y,
                                              fill= value)) +
    ggplot2::geom_tile() +
    ggplot2::xlab("")  +
    ggplot2::ylab("") +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill="gray30"),
                   strip.text = ggplot2::element_text(color="white")) +
    ggplot2::scale_fill_gradientn(colours=colors)

    if(condition_vs_feat){
      p <- p + ggh4x::facet_grid2(condition~feature, scale="free", independent = "x")
    }else{
      p <- p + ggh4x::facet_grid2(feature~condition, scale="free", independent = "y")
    }

  p
}





