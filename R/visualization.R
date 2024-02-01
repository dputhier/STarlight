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
           function(object,
                    gene=NULL,
                    saturation=0.75,
                    scale=TRUE,
                    colors=c('#000003', '#410967', '#932567', '#DC5039', '#FBA40A', '#FCFEA4'),
                    ceil=TRUE,
                    coord_fixed=TRUE)
             standardGeneric("spatial_image")
)

#' @title Create a spatial image of based on a  STGrid object.
#' @description
#' Create a spatial image of based on a  STGrid object.
#' @param object The STGrid object.
#' @export spatial_image
setMethod("spatial_image",
          signature(object = "STGrid"),
          function(object,
                   gene=NULL,
                   saturation=0.75,
                   scale=TRUE,
                   colors=c('#000003', '#410967', '#932567', '#DC5039', '#FBA40A', '#FCFEA4'),
                   ceil=TRUE,
                   coord_fixed=TRUE) {

            if(is.null(gene))
              print_msg("Please provide a gene name (see gene arguments).", msg_type = "STOP")

            if(length(gene) > 1)
              print_msg("Please provide a single gene name (see gene arguments).", msg_type = "STOP")

            if(!gene %in% c(gene_names(object), "all_genes"))
              print_msg("The gene not found in the object.", msg_type = "STOP")

            if(gene != "all_genes"){
              spatial_matrix <- object@bin_mat$layer_1[, c("bin_x", "bin_y", gene)]
              tmp <- spatial_matrix[, gene]
              tmp[is.na(tmp)] <- 0
            }else{
              spatial_matrix <- object@bin_mat$layer_1
              spatial_matrix$all_genes <- 0
              tmp <- spatial_matrix[, setdiff(gene_names(object), c("bin_x", "bin_y"))]
              spatial_matrix <- spatial_matrix[, c("bin_x", "bin_y", gene)]
              tmp[is.na(tmp)] <- 0
              tmp <- rowSums(tmp)
              print("YO")
            }

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
            return(p)
})


# -------------------------------------------------------------------------
##    Molecule counts
# -------------------------------------------------------------------------
#' @title Create a barplot to show molecule counts of selected genes
#' @description
#' Create a barplot to show molecule counts of selected genes
#' @param object The STGrid object.
#' @keywords internal
setGeneric("cmp_mol_counts",
           function(object,
                    genes=head(gene_names(object_1)),
                    name_1="Condition_1",
                    name_2="Condition_2",
                    colors=c("#4DBBD5B2", "#F39B7FB2"))
             standardGeneric("cmp_mol_counts")
)

#' @title Create a barplot to show molecule counts of selected genes
#' @description
#' Create a barplot to show molecule counts of selected genes
#' @param object_1 A STGrid object.
#' @param object_2 A STGrid object.
#' @export cmp_mol_counts
setMethod(
  "cmp_mol_counts", signature("STCompR"),
    function(object,
             genes=head(gene_names(object_1)),
             name_1="Condition_1",
             name_2="Condition_2",
             colors=c("#4DBBD5B2", "#F39B7FB2")) {

    check_var(name_1)
    check_var(name_2)

    if(is.null(genes)){
      print_msg("Please provide some genes...", msg_type = "STOP")
    }

    if(!all(genes %in% gene_names(object)))
       print_msg("Some genes where not found.", msg_type = "STOP")

    merge_coord <- rbind(tb_1, tb_2)
    colnames(merge_coord) <- c("Genes", "Counts", "Conditions")
    merge_coord$Genes <- factor(merge_coord$Genes, levels=genes, ordered = TRUE)

    ggplot2::ggplot(data=merge_coord,
                    mapping = ggplot2::aes(x=Genes, y=Counts,
                                           fill=Conditions)) +
    ggplot2::geom_col(color="black", linewidth=0, position="dodge") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=4, angle=45, vjust = 0.5),
          axis.text.y = ggplot2::element_text(size=8),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          panel.grid.minor.y =ggplot2:: element_blank(),
          panel.border = ggplot2::element_blank()) +
    ggplot2::ylab("Molecule count") +
    ggplot2::scale_fill_manual(values = colors)

  }
)

# -------------------------------------------------------------------------
##    Boxplot / jitter
# -------------------------------------------------------------------------
#' @title Create a boxplot/jitter plot to show molecule counts distribution.
#' @description
#' Create a boxplot/jitter plot to show molecule counts distribution.
#' @param object_1 A STGrid object.
#' @param object_2 A STGrid object.
#' @keywords internal
setGeneric("cmp_boxplot",
           function(object_1,
                    object_2,
                    genes=head(gene_names(object_1)),
                    name_1="Condition_1",
                    name_2="Condition_2",
                    colors=c("#4DBBD5B2", "#F39B7FB2"),
                    ...)
             standardGeneric("cmp_boxplot")
)

#' @title Create a boxplot/jitter plot to show molecule counts.
#' @description
#' Create a boxplot/jitter plot to show molecule counts.
#' @param object_1 A STGrid object.
#' @param object_2 A STGrid object.
#' ... Additiona parameters to pass to ggpol::geom_boxjitter().
#' @export cmp_boxplot
setMethod(
  "cmp_boxplot", signature("STGrid"),
  function(object_1,
           object_2,
           name_1="Condition_1",
           name_2="Condition_2",
           colors=ggsci::pal_npg("nrc", alpha = 1)(2),
           ...) {

    check_var(name_1)
    check_var(name_2)

    tb_1 <- as.data.frame(table(coord(object_1)$gene))
    tb_1$condition <- name_1
    print(tb_1)
    tb_2 <- as.data.frame(table(coord(object_2)$gene))
    tb_2$condition <- name_2


    merge_coord <- rbind(tb_1, tb_2)
    colnames(merge_coord) <- c("Genes", "Counts", "Conditions")

    merge_coord$Genes <- factor(merge_coord$Genes, levels=genes, ordered = TRUE)

    ggplot2::ggplot(data=merge_coord,
           mapping = ggplot2::aes(x=Conditions, y=log10(Counts), fill=Conditions)) +
      ggpol::geom_boxjitter(alpha=0.5, ...) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, angle=45, vjust = 0.5),
            axis.text.y = ggplot2::element_text(size=12),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank()) +
      ggplot2::ylab("log10(Molecule count)") +
      ggplot2::scale_fill_manual(values = colors)

  }
)

