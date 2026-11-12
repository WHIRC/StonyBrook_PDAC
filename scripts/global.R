library(tidyverse)
library(glue)
library(limma)
library(ConsensusClusterPlus)

set.seed(1)

# Read and filter LVL3 matrix -------------------------------------------------
mtx <- read.csv("./input/stonybrook_pooled_lvl3_matrix.csv") %>%
  as.data.frame() %>%
  `rownames<-`(., glue("{.$Protein.Accession};{.$Gene.Symbol}")) %>%
  select(
    -contains(".2"),
    -c("Protein.Accession", "Gene.Symbol")
  ) %>%
  rename_with(~ sub("^X", "", .x))

mtx <- mtx[, order(colnames(mtx))]

# Read metadata ---------------------------------------------------------------
metadata <- read.csv("./input/stonybrook_pooled_workorder.csv") %>%
  `rownames<-`(., .$sample.id) %>%
  select(-file.name) |>
  mutate(
    cell.type = cell.type |> sub(" |\\/", "_", x = _)
  ) |>
  filter(sample.id %in% colnames(mtx)) %>%
  arrange(sample.id)

stopifnot(colnames(mtx) == row.names(metadata))

# Set theme for all plots -----------------------------------------------------
theme_set(
  theme_classic() +
    theme(
      text = element_text(
        family = "Helvetica",
        size = 12
      ),
      panel.background = element_rect(
        fill = "grey98",
        color = "grey90"
      ),
      plot.background = element_rect(
        fill = "white",
        color = "white"
      ),
      legend.background = element_rect(
        fill = "white",
        color = "white",
        linewidth = 0.25
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(
        angle = 0,
        vjust = 0.5,
        hjust = 0.5,
        color = "grey20",
        # face = "bold"
      ),
      axis.title = element_text(
        angle = 0,
        vjust = 1,
        hjust = 0.5,
        color = "grey20",
        # face = "bold"
      ),
      plot.title = element_text(
        hjust = 0,
        color = "black",
        face = "bold"
      ),
      plot.subtitle = element_text(
        hjust = 0,
        color = "black"
      ),
      legend.title = element_text(
        color = "black"
      ),
      legend.text = element_text(
        color = "black"
      ),
      legend.spacing.y = unit(0.05, "cm"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "left",
      strip.background = element_rect(
        linetype = "solid",
        color = "black",
        fill = "grey90",
        linewidth = 1
      ),
      strip.text = element_text(
        face = "bold",
        color = "black"
      ),
    )
)

cell.type.pal <- c(
  DUAL_NEG = "#9ff9fd",
  GATA6 = "#f0c209",
  K17 = "#2f08fe",
  K17_GATA6 = "#6aa32b",
  All_others = "grey"
)

class(cell.type.pal) <- "colors"

#' @description subset input matrix to top n median absolute deviated (MAD) rows
#' @param mtx numeric matrix
#' @param n integer; number of rows in output matrix; default 100
#' @return numeric matrix with n rows
top_n_mad <- function(mtx, n = 100) {
  mtx[
    as.matrix(mtx) %>%
      matrixStats::rowMads() %>%
      order(decreasing = TRUE) %>%
      head(n),
  ] %>%
    return(.)
}

#' @description Recreate ConsensusClusterPlus (CCP) heatmaps but with metadata.
#' CCP heatmaps don't preserve column/row names and thus is born this function.
#' This function uses the exact same clustering/distance algorithm as CCP and
#' plots additional metdata.
#'
#' @params ccp_output -- list; results from ConsensusClusterPlus() function.
#'
#' @params metadata -- dataframe; dataframe of metadata. Rownames MUST be a
#' superset of names(ccp_output[["consensusClass"]])
#'
#' @params metadata_palettes -- named list; must pass a list of vectors, where
#' the class of each vector is "colors" (IMPORTANT!). The name of each vector
#' in this list must be a column name in metadata -- only these will appear on
#' plot. For each represented column in metadata, there must a name of each
#' element in the corresponding vector of metadata_palettes.
#' For example, if mtcars was our metadata and you want to plot "cyl" and
#' "gear" as metadata, then:
#' # create vector of class "colors" for cylinder
#' cylinder <- paletteer::paletter_d(
#'     palette = "ggthemes::colorblind",
#'     n = length(mtcars$cyl)
#' )
#' # make sure all mtcars$cyl values are represented in names()
#' names(cylinder) <- unique(mtcars$cyl)
#' # do the same for gears
#' gear <- paletteer::paletter_d(
#'     palette = "ggsci::default_jama",
#'     n = length(mtcars$gear)
#' )
#' names(gear) <- unique(mtcars$gear)
#' # create metadata_palettes with names matching metadata columns of interest
#' metadata_palettes <- list(
#'     "cyl" = cylinder,
#'     "gear" = gear
#' )
#'
#' @params ccp_palette; named list; same as metadata_palettes but for CCP 
#' clusters. If NULL, function will create palette for you using
#' palette = "ggthemr::light"
#'
#' @params k -- integer; which CCP number of clusters to plot
#' 
#' @params title -- string; title to show on plot, if desired
#' 
#' @param save_plot boolean; whether or not to save plot as png. default FALSE
#' 
#' @param save_options named list; parameters for saving plot. default
#' list(save_where = "results/heatmap.png", height = 6, width = 4, res = 500)
#'
#' @param save_ccp_flat string; file path to which CCP results excel file is
#' saved. Pass NULL if not desired. Default NULL
#' 
#' @returns pheatmap
ccp_heatmap_with_metadata <- function(
    ccp_output = NULL,
    metadata = NULL,
    metadata_palettes = NULL,
    ccp_palette = NULL,
    k = 3,
    title = NULL,
    save_plot = FALSE,
    save_options = list(
        save_where = "results/heatmap.png",
        height = 8,
        width = 8,
        res = 500
    ),
    save_ccp_flat = NULL
) {
    ccp_clusters <- ccp_output[[k]][["consensusClass"]]
    sample_names <- names(ccp_clusters)
    # ccp matrix with row and column names
    ccp_mtx <- ccp_output[[k]][["consensusMatrix"]] %>%
        `rownames<-`(., sample_names) %>%
        `colnames<-`(., sample_names)
    # sample_names[which(!sample_names%in% row.names(metadata))]
    # make sure that all ccp samples are represented in metadata
    if (
        !all( sample_names %in% row.names(metadata) )
    ) {
        stop(
            "All sample names from CCP must be in rownames of metadata"
        )
    } else {
        metadata$CaseID <- row.names(metadata)
    }
    # merge metadata with ccp cluster designations
    metadata <- ccp_clusters %>%
        as.data.frame() %>%
        `colnames<-`(., c("Consensus Cluster")) %>%
        mutate(
            `Consensus Cluster` = sub("^", "Cluster ", `Consensus Cluster`),
            CaseID = row.names(.)
        ) %>%
        left_join(., metadata, by = "CaseID") %>%
        `rownames<-`(., .$CaseID) %>%
        select(
            any_of(
                c("Consensus Cluster", names(metadata_palettes))
            )
        )
    # make sure all metadata are represented in metadata_palettes. If so,
    # subset each metadata_palette to elements in metadata.
    metadata_palettes <- mapply(
        function(plot_column, metadata_palette) {
            if (!all(
                metadata[[plot_column]] %in% names(metadata_palette)
            )) {
                stop("Not all metadata in plot_columns represented in metadata_color_palette")
            }
            metadata_palette[
                names(metadata_palette) %in% metadata[[plot_column]]
            ] %>%
                return(.)
        },
        names(metadata_palettes),
        metadata_palettes,
        SIMPLIFY = FALSE
    )
    if ( length( ccp_clusters %>% unique() ) > 9 ) {
        stop("Too many clusters (>9), must change color palette or use fewer clusters")
    }
    ccp_palette <- paletteer::paletteer_d(
      palette = "ggthemr::light"
    )[1:k] # %>% as.list()
    names(ccp_palette) <- seq_along(ccp_palette) %>%
      sub("^", "Cluster ", .)
    meta_names <- c(names(metadata_palettes), "Consensus Cluster")
    # crazy that append doesn't preserve the class() of each element
    metadata_palettes[[ length(metadata_palettes) + 1 ]] <- ccp_palette
    names(metadata_palettes) <- meta_names
    hc = hclust(as.dist(1 - ccp_mtx), method = "average")
    heatmap_palette <- paletteer::paletteer_c(
        palette = "ggthemes::Classic Blue",
        n = 30
    ) %>% as.vector()
    heatmap <- pheatmap::pheatmap(
        mat = ccp_mtx,
        cluster_rows = hc,
        cluster_cols = hc,
        cellwidth = 4,
        cellheight = 142 / nrow(ccp_mtx),
        border_color = NA,
        color = heatmap_palette,
        treeheight_row = 0,
        treeheight_col = 30,
        cutree_cols = k,
        cutree_rows = k,
        show_rownames = FALSE,
        show_colnames = FALSE,
        labels_col = sub("\\..*", "", colnames(ccp_mtx)),
        angle_col = 90,
        main = title,
        legend = TRUE,
        annotation_col = metadata, 
        annotation_names_col = FALSE,
        annotation_colors = metadata_palettes
    )
    if (!is.null(save_ccp_flat)) {
        library(openxlsx2)
        wb_workbook() %>%
            wb_add_worksheet(sheet = "CCP matrix") %>%
            wb_add_data(
                x = ccp_mtx,
                row_names = TRUE
            ) %>%
            wb_add_worksheet(sheet = "CCP clusters") %>%
            wb_add_data(
                x = ccp_clusters %>%
                    as.data.frame() %>%
                    `colnames<-`(., c("CCP Cluster")) %>%
                    mutate(Sample = row.names(.)) %>%
                    select(Sample, "CCP Cluster"),
                row_names = FALSE
            ) %>%
            wb_save(file = save_ccp_flat)
    }
    if (save_plot) {
      library(grid)

        dir.create(
            path = dirname(save_options[["save_where"]]),
            recursive = TRUE,
            showWarnings = FALSE
        )
        png(
            file = save_options[["save_where"]],
            height = save_options[["height"]],
            width = save_options[["width"]],
            res = save_options[["res"]],
            unit = "in"
        )
        grid.draw(rectGrob(gp=gpar(fill="grey96", lwd=0)))
        grid.draw(heatmap)
        # print(heatmap)
        dev.off()
    }
    return(heatmap)
}

odds_ratio_analysis <- function(
    meta= NULL,
    X = NULL,
    Y = NULL,
    method = "midp",
    heatmap_title,
    save_plot = FALSE,
    save_options = list(
        save_where = "results/or.png",
        height = 5,
        width = 5,
        res = 500
    )
) {
    xs <- meta[[X]] %>% unique() %>% sort()
    ys <- meta[[Y]] %>% unique() %>% sort()
    OR.statistics <- lapply(
        xs, function(x) {
            lapply(
                ys, function(y) {
                    # build contigency table:
                    # ______| x.yes | x.no |
                    # y.yes | ...   | ...  |
                    # y.no  | ...   | ...  |
                    OR <- matrix(
                        c(
                            meta %>% # x.yes y.yes
                                filter(get(X) == x) %>%
                                filter(get(Y) == y) %>%
                                nrow(),
                            meta %>% # x.no y.yes
                                filter(get(X) != x) %>%
                                filter(get(Y) == y) %>%
                                nrow(),
                            meta %>% # x.yes y.no
                                filter(get(X) == x) %>%
                                filter(get(Y) != y) %>%
                                nrow(),
                            meta %>% # x.no y.no
                                filter(get(X) != x) %>%
                                filter(get(Y) != y) %>%
                                nrow()
                        ),
                        nrow = 2,
                        ncol = 2,
                        byrow = TRUE
                    ) %>%
                        # run oddsratio analysis
                        epitools::oddsratio(
                            method = method
                        )
                    (OR$measure["Exposed2", "estimate"] + 0.001) %>%
                        log2() %>%
                        round(2) %>%
                        return(.)
                }
            ) %>%
                # Stack X statistics into one column
                do.call(rbind, .) %>%
                as.data.frame() %>%
                `rownames<-`(., ys) %>%
                `colnames<-`(., c(x))
        }
    ) %>%
        # Bring together all columns
        do.call(cbind, .)
    # Same thing but this time to get p-values
    OR.pvalues <- lapply(
        xs, function(x) {
            lapply(
                ys, function(y) {
                    OR <- matrix(
                        c(
                            meta %>% # x.yes y.yes
                                filter(get(X) == x) %>%
                                filter(get(Y) == y) %>%
                                nrow(),
                            meta %>% # x.no y.yes
                                filter(get(X) != x) %>%
                                filter(get(Y) == y) %>%
                                nrow(),
                            meta %>% # x.yes y.no
                                filter(get(X) == x) %>%
                                filter(get(Y) != y) %>%
                                nrow(),
                            meta %>% # x.no y.no
                                filter(get(X) != x) %>%
                                filter(get(Y) != y) %>%
                                nrow()
                        ),
                        nrow = 2,
                        ncol = 2,
                        byrow = TRUE
                    ) %>%
                        epitools::oddsratio(
                            method = method
                        )
                    OR$p.value["Exposed2", 1] %>%
                        # round(3) %>%
                        gtools::stars.pval() %>%
                        sub(" ", "ns", .) %>%
                        return(.)
                }
            ) %>%
                do.call(rbind, .) %>%
                as.data.frame() %>%
                `rownames<-`(., ys) %>%
                `colnames<-`(., c(x))
        }
    ) %>%
        do.call(cbind, .)
    stopifnot(row.names(OR.statistics) == row.names(OR.pvalues))
    stopifnot(colnames(OR.statistics) == colnames(OR.pvalues))
    # display_matrix <- matrix(
    #     NA,
    #     ncol = ncol(OR.statistics),
    #     nrow = nrow(OR.statistics)
    # )
    # for (r in seq_along(OR.statistics)) {
    #     for (c in seq_along(OR.statistics))
    #         display_matrix[r, c] <- paste(
    #             OR.statistics[r, c], OR.pvalues[r, c]
    #         )
    # }
    heatmap_palette <- paletteer::paletteer_c(
        # palette = "grDevices::Purp", # unidirectional
        palette = "grDevices::Blue-Red 2", # bidirectional
        # palette = "ggthemes::Orange-Blue Diverging", # bidirectional
        n = 30
    ) %>% as.vector()
    hm <- pheatmap::pheatmap(
        mat = OR.statistics,
        color = heatmap_palette,
        border_color = "white",
        cellwidth = 30,
        cellheight = 30,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        angle_col = 90,
        treeheight_row = 0,
        treeheight_col = 0,
        fontsize = 7,
        fontsize_row = 9,
        fontsize_col = 9,
        fontsize_number = 9,
        display_numbers = OR.pvalues,
        number_color = "black",
        main = heatmap_title,
        breaks = seq( # center scale at 0
            -(OR.statistics %>% abs() %>% max()),
            OR.statistics %>% abs() %>% max(),
            length.out = length(heatmap_palette)
        )
    )
    if (save_plot) {
        dir.create(
            path = dirname(save_options[["save_where"]]),
            recursive = TRUE,
            showWarnings = FALSE
        )
        png(
            file = save_options[["save_where"]],
            height = save_options[["height"]],
            width = save_options[["width"]],
            res = save_options[["res"]],
            unit = "in"
        )
        print(hm)
        dev.off()
    }
    list(
        statistic_matrix = OR.statistics,
        pvalue_matrix = OR.pvalues
    ) %>%
        return(.)
}

pca_projection_analysis <- function(
  start,
  end,
  mtx,
  metadata
) {

  pca <- prcomp(t(mtx), scale. = TRUE)

  pca_df <- as.data.frame(pca$x) %>%
    rownames_to_column("sample.id") %>%
    left_join(metadata, by = "sample.id")

  centroids <- pca_df %>%
    group_by(cell.type) %>%
    summarise(across(starts_with("PC"), mean))

  vector_names <- metadata$cell.type %>% unique()
  names(vector_names) <- vector_names

  for (transition in vector_names) {
    for (n_pcs in c(2:10)) {

      vectors <- lapply(
        vector_names,
        function(n) {
          centroids %>%
            filter(cell.type == n) %>%
            select(PC1:paste0("PC", n_pcs)) %>%
            unlist()
        }
      )

      axis <- vectors[[end]] - vectors[[start]]
      h <- vectors[[transition]] - vectors[[start]]

      t_proj <- sum(h * axis) / sum(axis^2)

      proj_point <- vectors[[start]] + t_proj * axis
      orth_dist <- sqrt(sum((vectors[[transition]] - proj_point)^2))
      axis_length <- sqrt(sum(axis^2))

      n_boot <- 5000

      boot_results <- replicate(n_boot, {

        start_boot <- pca_df %>%
          filter(cell.type == start) %>%
          slice_sample(prop = 1, replace = TRUE)

        end_boot <- pca_df %>%
          filter(cell.type == end) %>%
          slice_sample(prop = 1, replace = TRUE)

        transition_boot <- pca_df %>%
          filter(cell.type == transition) %>%
          slice_sample(prop = 1, replace = TRUE)

        pcs <- paste0("PC", 1:n_pcs)
        start_c_b <- colMeans(start_boot[pcs])
        end_c_b <- colMeans(end_boot[pcs])
        transition_c_b <- colMeans(transition_boot[pcs])

        axis_b <- end_c_b - start_c_b
        h_b <- transition_c_b - start_c_b
        t_b <- sum(h_b * axis_b) / sum(axis_b^2)

        proj_b <- start_c_b + t_b * axis_b
        orth_b <- sqrt(sum((transition_c_b - proj_b)^2))
        axis_length_b <- sqrt(sum(axis_b^2))

        c(t_proj = t_b, orth_ratio = orth_b / axis_length_b)

      }, simplify = TRUE)

      boot_df <- as.data.frame(t(boot_results))

      lower_quant <- quantile(boot_df$t_proj, 0.025) %>% round(2)
      upper_quant <- quantile(boot_df$t_proj, 0.975) %>% round(2)

      horse <- tibble(
        `Projected cell.type` = transition,
        `Number of PCs` = n_pcs,
        `Projection point, P` = round(t_proj, 3),
        `95% CI` = glue("[{lower_quant} - {upper_quant}]"),
        `Orthogonal distance, O` = round(orth_dist / axis_length, 3),
      )
      summary_path <- file.path(output, "PC_projection.csv")

      if (!file.exists(summary_path)) {
        write.csv(horse, summary_path, row.names = FALSE)
      } else {
        write.table(
          horse,
          summary_path,
          row.names = FALSE,
          col.names = FALSE,
          sep = ",",
          append = TRUE
        )
      }

      glue(
        "
        ---------------------------------------------------------------------------
        PCs:{n_pcs}
        ---------------------------------------------------------------------------

        Projection point P on {start} -> {end} axis:
        {round(t_proj, 3)}

        Orthogonal distance O as fraction of axis:
        {round(orth_dist / axis_length, 3)}

        Bootstrap 95% CI:
        [{lower_quant} - {upper_quant}]

        "
      ) %>%
      print()

      horse <- ggplot(
        boot_df,
        aes(x = t_proj)
      ) +

        geom_histogram(
          bins = 50,
          fill = "steelblue",
          alpha = 0.5
        ) +

        # Observed value
        geom_vline(
          xintercept = t_proj,
          color = "firebrick",
          linewidth = 1,
          linetype = "solid"
        ) +

        # # Standard deviation
        # geom_vline(
        #   xintercept = mean(boot_df$t_proj) - sd(boot_df$t_proj),
        #   color = "steelblue",
        #   linetype = "dashed",
        #   alpha = 0.7
        # ) +
        # geom_vline(
        #   xintercept = mean(boot_df$t_proj) + sd(boot_df$t_proj),
        #   color = "steelblue", 
        #   linetype = "dashed",
        #   alpha = 0.7
        # ) +
        #
        # # Mean
        # geom_vline(
        #   xintercept = mean(boot_df$t_proj),
        #   color = "steelblue",
        #   linewidth = 0.7,
        #   alpha = 0.7
        # ) +

        # 95% CI
        geom_vline(
          xintercept = quantile(boot_df$t_proj, 0.025),
          color = "gray30",
          linetype = "dashed"
        ) +
        geom_vline(
          xintercept = quantile(boot_df$t_proj, 0.975),
          color = "gray30",
          linetype = "dashed"
        ) +

        # Reference lines
        geom_vline(
          xintercept = 0,
          color = "darkgreen",
          linetype = "dotted"
        ) +
        geom_vline(
          xintercept = 1,
          color = "darkorange",
          linetype = "dotted"
        ) +

        annotate(
          "text",
          x= t_proj,
          y = Inf,
          label = "Observed",
          vjust = -1,
          hjust = -0.1,
          color = "firebrick",
          angle = 270,
          fontface = "bold",
          size = 3
        ) +

        annotate(
          "text",
          x = 0,
          y = Inf,
          label = start,
          vjust = -1,
          hjust = -0.1,
          color = "darkgreen",
          angle = 270,
          fontface = "bold",
          size = 3
        ) +

        annotate(
          "text",
          x = 1,
          y = Inf,
          label = end,
          vjust = -1.2,
          hjust = -0.1,
          color = "darkorange",
          angle = 270,
          fontface = "bold",
          size = 3
        ) +
        # scale_x_continuous(breaks = seq(0, 1.2, by = 0.1)) +

        scale_fill_manual(
          values = cell.type.pal
        ) +

        labs(
          title = "Bootstrap distribution of transitional centroid projection",
          subtitle = glue("Along Benign → Carcinoma axis in {n_pcs}-D PC space"),
          x = "Projection (0 = Benign, 1 = Carcinoma)",
          y = ""
        ) +

        theme(
            axis.ticks = element_line()
        )

      png(
        file = file.path(
          output,
          glue("bootstrap_distribution_{transition}_{n_pcs}.png")
        ),
        width = 4.5, height = 3.5, res = 300, unit = "in"
      )
      print(horse)
      dev.off()
    }
  }

  return(NULL)
}
