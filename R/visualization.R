# -------------------------------------------------------------------------
##    Spatial image
# -------------------------------------------------------------------------
#' @title Create a spatial image of based on a  STGrid object.
#' @description
#' Create a spatial image of based on a  STGrid object.
#' @param object The STGrid object.
#' @keywords internal
#'
setGeneric("spatial_image",
           function(object=NULL,
                    gene=NULL,
                    saturation=0.75,
                    scale=TRUE,
                    colors=c('#000003', '#410967', '#932567', '#DC5039', '#FBA40A', '#FCFEA4'),
                    ceil=TRUE,
                    coord_fixed=TRUE,
                    overlay_gene=NULL,
                    colors_overlay=c("#DEEBF7", "#9ECAE1", "#3182BD"))
             standardGeneric("spatial_image")
)

#' @title Create a spatial image of based on a  STGrid object.
#' @description
#' Create a spatial image of based on a  STGrid object.
#' @param object The STGrid object.
#' @export spatial_image
setMethod("spatial_image",
          signature(object = "STGrid"),
          function(object=NULL,
                   gene=NULL,
                   saturation=0.75,
                   scale=TRUE,
                   colors=c('#000003', '#410967', '#932567', '#DC5039', '#FBA40A', '#FCFEA4'),
                   ceil=TRUE,
                   coord_fixed=TRUE,
                   overlay_gene=NULL,
                   colors_overlay=c("#DEEBF7", "#9ECAE1", "#3182BD")) {

            if(is.null(object))
                print_msg("Please provide an STGrid object.",
                        msg_type = "STOP")

            if(!is.null(overlay_gene)){
              if(!overlay_gene %in% gene_names(object))
                print_msg("The gene was not found in the object.",
                          msg_type = "STOP")

            }

            if(is.null(gene))
              print_msg("Please provide a gene name (see gene arguments).", msg_type = "STOP")

            if(length(gene) > 1)
              print_msg("Please provide a single gene name (see gene arguments).", msg_type = "STOP")

            if(!gene %in% gene_names(object))
              print_msg("The gene was not found in the object.", msg_type = "STOP")

              spatial_matrix <- object@bin_mat[, c("bin_x", "bin_y", gene)]
              tmp <- spatial_matrix[, gene]
              tmp[is.na(tmp)] <- 0

            if(saturation != 0){
              q_sat <- quantile(tmp[tmp != 0], saturation)
              tmp[tmp > q_sat] <- q_sat
            }

            if(scale){
              if(max(tmp) > 1){
                tmp <- (tmp - min(tmp)) / (max(tmp) - min(tmp))
              }
            }

            spatial_matrix[, gene] <- tmp

            spatial_matrix$bin_x <- factor(spatial_matrix$bin_x,
                                           levels = bin_x(object),
                                           ordered=TRUE)

            spatial_matrix$bin_y <- factor(spatial_matrix$bin_y,
                                           levels = bin_y(object),
                                           ordered=TRUE)

            spatial_matrix_melted <- reshape2::melt(spatial_matrix, id.vars=c("bin_x", "bin_y"))

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


            if(!is.null(overlay_gene)){
              over <- na.omit(bin_mat(object,
                              as_factor = TRUE)[ ,c("bin_x", "bin_y", overlay_gene)])
              p <- p + geom_point(data=over, mapping=aes(x=bin_x,
                                                         y=bin_y,
                                                         color=log10(.data[[overlay_gene]])),
                                  size=0.1,
                                  inherit.aes = FALSE) +
                scale_color_gradientn(colors=colors_overlay)
            }

            return(p)
})


# -------------------------------------------------------------------------
##    Molecule counts
# -------------------------------------------------------------------------
#' @title Create a barplot to show molecule counts of selected genes
#' @description
#' Create a barplot to show molecule counts of selected genes
#' @param object A STGrid object.
#' @param genes The list of genes (NULL for all of them).
#' @param normalized  Whether counts should be normalized.
#' @param transform Whether the count should be transformed (the pseudo count defined for the object is added).
#' @param facet_shared Use ggpol::facet_shared.
#' @keywords internal
setGeneric("cmp_bar_plot",
           function(object,
                    genes=head(gene_names(object)),
                    normalized=FALSE,
                    transform=c("None", "log2", "log10", "log"),
                    colors=c("#3074BB", "#BE5B52"))
             standardGeneric("cmp_bar_plot")
)

#' @title Create a barplot to show molecule counts of selected genes
#' @description
#' Create a barplot to show molecule counts of selected genes
#' @param object A STGrid object.
#' @param genes The list of genes (NULL for all of them).
#' @param normalized  Whether counts should be normalized.
#' @param transform Whether the count should be transformed (the pseudo count defined for the object is added).
#' @export cmp_bar_plot
setMethod(
  "cmp_bar_plot", signature("STCompR"),
    function(object,
             genes=head(gene_names(object)),
             normalized=FALSE,
             transform=c("None", "log2", "log10", "log"),
             colors=c("#3074BB", "#BE5B52")) {


    if(is.null(genes)){
      print_msg("Please provide some genes...", msg_type = "STOP")
    }

    transform <- match.arg(transform)

    counts <- stat_test(object,
                        normalized=normalized,
                        count_only = TRUE,
                        melted_count = TRUE,
                        transform=transform,
                        genes=genes)

    ylabel <- ifelse(transform %in% c("log2", "log10", "log"),
                     paste0(transform, "(Molecule counts)"),
                     "Molecule counts")

    ggplot2::ggplot(data=counts,
                    mapping = ggplot2::aes(x=Genes, y=Counts,
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
#' @import ggsci
#' @import ggpol
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
#' @param text_y_lim Numeric value specifying the threshold for gene labels on the y-axis. Genes with -log10(p-value) less than this threshold will not be labeled. Default is 100.
#' @param text_x_lim Numeric value specifying the threshold for gene labels on the x-axis. Genes with absolute x-values less than this threshold will not be labeled. Default is 2.
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
#' @param text_y_lim Numeric value specifying the threshold for gene labels on the y-axis. Genes with -log10(p-value) less than this threshold will not be labeled. Default is 100.
#' @param text_x_lim Numeric value specifying the threshold for gene labels on the x-axis. Genes with absolute x-values less than this threshold will not be labeled. Default is 2.
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
                          genes=NULL)

    volc_data <- volc_data[, c(x_axis, y_axis)]
    colnames(volc_data) <- c("x", "y")

    counts <- stat_test(object,
                        count_only = TRUE,
                        melted_count = FALSE,
                        genes=NULL)

    volc_data$mean_counts <- rowMeans(counts)
    volc_data$gene <- rownames(volc_data)
    volc_data$gene[-log10(volc_data$y) < text_y_lim] <- NA
    volc_data$gene[abs(volc_data$x) < text_x_lim] <- NA

    ggplot2::ggplot(data=volc_data,
           mapping = ggplot2::aes(x=x, y=-log10(y),
                         fill=x,
                         size=mean_counts)) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_point(shape=21, color="black", stroke=0.2) +
      ggplot2::theme_bw() +
      ggrepel::geom_text_repel(data=na.omit(volc_data), mapping=ggplot2::aes(label=gene),
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
#' @param Which kind of correction should be displayed.
#' @param max_gene_label The maximum number of genes to display.
#' The genes with the highest correction value (whatever the radius)
#' are displayed
#' @keywords internal
setGeneric("plot_rip_k",
           function(object,
                    correction=c("border",
                                 "isotropic",
                                 "Ripley",
                                 "translate"),
                    max_gene_label=10)
             standardGeneric("plot_rip_k")
)

#' @title Call the Ripley's k function.
#' @description
#' Call the Ripley's k function.
#' @param object A STCompR object.
#' @param Which kind of correction should be displayed.
#' @param max_gene_label The maximum number of genes to display.
#' The genes with the highest correction value (whatever the radius)
#' are displayed.
#' @export plot_rip_k
setMethod(
  "plot_rip_k", signature("STGrid"),
  function(object,
           correction=c("border",
                        "isotropic",
                        "Ripley",
                        "translate"),
           max_gene_label=10) {

    if(nrow(ripley_k_function(object)) == 0){
      print_msg("Please run ripley_k_function() first.", msg_type = "STOP")
    }

    correction <- match.arg(correction)
    ripk <- ripley_k_function(object)

    if(max_gene_label > nrow(ripk))
      print_msg("Too much selected gehes...", msg_type = "STOP")

    voi <- ripk %>%
      dplyr::group_by(gene) %>%
      dplyr::filter(border == max(border)) %>%
      dplyr::arrange(desc(border)) %>%
          head(n=max_gene_label)

    goi <- rownames(voi)

    ggplot2::ggplot(data= ripk) +
      ggplot2::geom_line(mapping = ggplot2::aes(x=r, y=border, group=gene), color="black") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.text = ggplot2::element_text(size=4),
                      legend.position = "none",
                      panel.grid.minor =  ggplot2::element_blank()) +
      ggrepel::geom_label_repel(data=voi,
                                mapping=ggplot2::aes(x=r, y=border, label=gene, color=gene),
                                inherit.aes = FALSE,
                                force=20) +
      ggplot2::ylab(paste0("Ripley's K function (correction=", correction, ")"))
})










