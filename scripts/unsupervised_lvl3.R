source("./scripts/global.R")
output <- "./output/unsupervised_analysis/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

cell.type.pal <- cell.type.pal[names(cell.type.pal) != "All_others"]

pca_mtx <- mtx %>% top_n_mad(500)

sdev <- mtx %>%
  prcomp(., rank. = 2) %>%
  .$sdev
variance <- round(100 * (sdev^2 / sum(sdev^2))[1:2], 1)

pca_fit <- prcomp(t(mtx), scale. = TRUE)

pca_mtx <- as.data.frame(pca_fit$x[, 1:2]) %>%
  rownames_to_column("id") %>%
  left_join(
    metadata %>% rownames_to_column("id"),
    by = "id"
  )

# Plot centroid distances as annotated light grey lines. Mark centroids with
# an X or something
png(
  file = file.path(
    output, "PCA_LVL3.png"
  ),
  width = 4.5, height = 3.5, res = 300, unit = "in"
)
pca_mtx %>%
  ggplot(
    aes(
      x = PC1,
      y = PC2,
      color = cell.type,
      fill = cell.type
    )
  ) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.1,
    linewidth = 0.1,
  ) +
  geom_point(
    size = 2.4,
    shape = 1,
    color = "black",
    show.legend = FALSE
  ) +
  geom_point(
    size = 2,
    alpha = 1
  ) +
  scale_color_manual(
    values = cell.type.pal
  ) +
  scale_fill_manual(
    values = cell.type.pal
  ) +
  labs(
    x = paste0("PC1: ", variance[1], "% variance"),
    y = paste0("PC2: ", variance[2], "% variance"),
    title = "PCA, top 500 MAD",
    # title = "PCA, K17 vs GATA6 DE proteins",
    subtitle = glue("{nrow(mtx)} Proteins")
  ) +
  theme(
    legend.title = element_blank(),
    axis.text = element_blank()
  )
dev.off()


annotation.palette <- list(
  cell.type = cell.type.pal
)

# Create heatmap colors
heatmap_palette <- paletteer::paletteer_c(
  palette = "grDevices::Blue-Red 2",
  n = 30
) %>%
  as.vector()


hm_mtx <- mtx %>% top_n_mad(n = 100)
# Row-wise median center matrix, just for visualization purposes
hm_mtx <- t(apply(hm_mtx, 1, function(x) x - median(x, na.rm = TRUE)))
hm_title <- glue("LVL3 matrix, top 500 MAD")

# Create heatmap and save
png(
  file = file.path(
    output, "heatmap_LVL3.png"
  ),
  width = 7, height = 4, res = 500, unit = "in"
)
pheatmap::pheatmap(
  mat = hm_mtx,
  main = hm_title,
  clustering_method = "ward.D2",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  cellwidth = 5,
  cellheight = 150 / nrow(hm_mtx),
  border_color = "white",
  color = heatmap_palette,
  treeheight_row = 0,
  treeheight_col = 20,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cutree_cols = 2,
  angle_col = 90,
  legend = TRUE,
  annotation_col = metadata %>% select(cell.type),
  annotation_names_col = FALSE,
  annotation_colors = annotation.palette,
  breaks = seq(
    # center scale at 0
    -(hm_mtx %>% abs() %>% max(na.rm = TRUE)),
    hm_mtx %>% abs() %>% max(na.rm = TRUE),
    length.out = 30
  )
)
dev.off()


ccp_out <- file.path(
  output, "CCP"
)
dir.create(ccp_out, showWarnings = FALSE, recursive = TRUE)

ccp <- mtx %>%
  top_n_mad(
    n = round(nrow(.) * 0.25, 0)
  ) %>%
  as.matrix() %>%
  ConsensusClusterPlus(
    d = .,
    maxK = 8,
    reps = 10000,
    pItem = 0.9,
    pFeature = 0.9,
    title = ccp_out,
    clusterAlg = "pam",
    distance = "spearman",
    seed = 1,
    plot = "png"
  )

# decide optimal number of k and retrieve respective results
k <- 3
kluster_assignment <- ccp[[k]][["consensusClass"]]

metadata <- metadata %>%
  left_join(
    data.frame(
      sample.id = names(kluster_assignment),
      ccp.cluster = kluster_assignment,
      check.names = FALSE
    ) %>%
      mutate(
        ccp.cluster = sub("^", "Cluster ", ccp.cluster) 
      )
  ) %>%
  `rownames<-`(., .$sample.id)


x <- ccp_heatmap_with_metadata(
  ccp_output = ccp,
  metadata = metadata,
  metadata_palettes = list(
    cell.type = cell.type.pal
  ),
  k = k,
  title = glue("CCP, k={k}"),
  save_plot = TRUE,
  save_options = list(
    save_where = file.path(
      output, "heatmap_CCP_k3.png"
    ),
    height = 3.5,
    width = 7,
    res = 200
  ),
  save_ccp_flat = file.path(
    output, "CCP_flatfile.xlsx"
  )
)

horse <- odds_ratio_analysis(
  meta = metadata,
  Y = "ccp.cluster",
  X = "cell.type",
  method = "small",
  heatmap_title = glue(
    "log2(odds ratio), k = {k}"
  ),
  save_plot = TRUE,
  save_options = list(
    save_where = file.path(
      output, "heatmap_odds-ratio.png"
    ),
    height = 3,
    width = 3,
    res = 500
  )
)
