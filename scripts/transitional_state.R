source("./scripts/global.R")
output <- "./output/transitional_state/"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Testing for dual negative projection beyond K17 (P > 1)

# 1. Unsupervised PCA (all 4,766 proteins and top 500 MAD proteins)
# DUAL_NEG clustered near GATA6 in PC space, not beyond K17. No visual evidence
# of a linear four-point trajectory.

# 2. Centroid projection onto GATA6→K17 axis (GATA6 vs K17 DE proteins)
# DUAL_NEG projected at P ~ 0.5 — halfway between GATA6 and K17, not beyond K17
# as the four-point model would predict.

# 3. Limma differential analysis
# DUAL_NEG vs GATA6 produced only 119 significant proteins at p < 0.01 — the
# weakest pairwise contrast in the entire dataset. If DUAL_NEG were a
# downstream state from K17, you'd expect it to be highly differentiated from
# GATA6.

# 4. Jonckheere-Terpstra monotonicity test
# Three-point ordering (GATA6 → hybrid → K17): 286 monotonic proteins, DUAL_NEG
# projects at P ~ 0.5 on those proteins
# Four-point ordering (GATA6 → hybrid → K17 → DUAL_NEG): only 111 monotonic
# proteins, and the projection is circular by construction


# -----------------------------------------------------------------------------
cell.type.pal <- cell.type.pal[
  !names(cell.type.pal) %in% c("All_others") ]


# # -----------------------------------------------------------------------------
# # Determine which proteins are monotonic as disease progresses along apriori
# # defined groups.
# # We will subset the matrix to these proteins before running PCA projection
# # analysis.
# #
# # Jonckheere-Terpstra test:
# # statology.org/what-is-the-jonckheere-terpstra-test-definition-example/
#
# # Define ordered group factor
# groups <- c("GATA6", "K17_GATA6", "K17")
#
# # Subset matrix to groups only
# metadata <- metadata %>%
#   filter(cell.type %in% groups)
#
# mtx <- mtx[, rownames(metadata)]
#
# group_order <- factor(
#   metadata$cell.type,
#   levels = groups,
#   ordered = TRUE
# )
#
# # Test each protein
# jt_results <- apply(mtx, 1, function(protein_values) {
#   clinfun::jonckheere.test(protein_values, group_order)$p.value
# })
#
# # Identify significant proteins
# monotonic_proteins <- names(jt_results[jt_results < 0.01])
#
# # Write out significant proteins
# monotonic_proteins %>%
#   sub(".*;", "", .) %>%
#   write.csv(
#     x = .,
#     file = file.path(output, "monotonic_proteins.csv"),
#     row.names = FALSE
#   )
#
# source("./scripts/global.R")
#
# mtx <- mtx %>%
#   select(any_of(metadata$sample.id)) %>%
#   filter(row.names(.) %in% monotonic_proteins)


# # # -----------------------------------------------------------------------------
# Filter to only proteins idenfitied as significant in limma GATA6 vs K17
genes <- read.csv("output/supervised_analysis/limma_GATA6_vs_K17.csv") %>%
  filter(P.Value < 0.01) %>%
  pull(Gene)
mtx <- mtx %>%
  select(any_of(metadata$sample.id)) %>%
  mutate(gene = row.names(.) %>% sub(".*;", "", .)) %>%
  filter(gene %in% genes) %>%
  select(-gene)

# -----------------------------------------------------------------------------
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

    start <- "GATA6"
    end <- "K17"

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
        title = "Bootstrap distribution",
        subtitle = glue("Along Classical → Basal axis; {n_pcs}-D PCs"),
        x = "Projection (0 = Classical, 1 = Basal)",
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
      width = 4.5, height = 4.5, res = 300, unit = "in"
    )
    print(horse)
    dev.off()
  }
}


# -----------------------------------------------------------------------------
vectors <- lapply(
  list(
    gata = "GATA6",
    k17 = "K17",
    hybrid = "K17_GATA6"
  ),
  function(n) {
    centroids %>%
      filter(cell.type == n) %>%
      select(PC1:paste0("PC", 2)) %>%
      unlist()
  }
)

# Project hybrid onto GATA6 -> K17 vector
v <- vectors[['k17']] - vectors[['gata']]
h <- vectors[['hybrid']] - vectors[['gata']]

# Scalar projection: how far along the GATA6->K17 axis is the hybrid?
t_proj <- sum(h * v) / sum(v * v)

# Orthogonal distance from the axis
proj_point <- vectors[['gata']] + t_proj * v
orth_dist <- sqrt(sum((vectors[['hybrid']] - proj_point)^2))

axis_length <- sqrt(sum(v * v))


png(
  file = file.path(
    output, "PCA_transitional_state_dual_positive.png"
  ),
  width = 5, height = 3.5, res = 300, unit = "in"
)
pca_df %>%
  ggplot(
    aes(
      x     = PC1,
      y     = PC2,
      color = cell.type,
      fill  = cell.type
    )
  ) +

  geom_point(
    size  = 2,
    alpha = 0.5
  ) +

  stat_ellipse(
    geom      = "polygon",
    alpha     = 0.5,
    linewidth = 0.05,
  ) +

  # GATA6 --> K17 vector
  geom_segment(
    aes(
      x    = vectors[['gata']]['PC1'],
      xend = vectors[['k17']]['PC1'],
      y    = vectors[['gata']]['PC2'],
      yend = vectors[['k17']]['PC2']
    ),
    color       = cell.type.pal[["K17"]],
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +

  # GATA6 --> Hybrid vector
  geom_segment(
    aes(
      x    = vectors[['gata']]['PC1'],
      xend = vectors[['hybrid']]['PC1'],
      y    = vectors[['gata']]['PC2'],
      yend = vectors[['hybrid']]['PC2']
    ),
    color       = cell.type.pal[["K17_GATA6"]],
    linewidth   = 0.8,
    inherit.aes = FALSE
  ) +

  # Perpendicular drop from hybrid centroid to axis
  geom_segment(
    aes(
      x    = vectors[['hybrid']]['PC1'],
      xend = proj_point['PC1'],
      y    = vectors[['hybrid']]['PC2'],
      yend = proj_point['PC2']
    ),
    color       = "black",
    linetype    = "dotted",
    linewidth   = 0.6,
    inherit.aes = FALSE
  ) +

  # Projected point on axis
  geom_point(
    aes(
      x = proj_point['PC1'],
      y = proj_point['PC2']
    ),
    size        = 3,
    color       = "black",
    inherit.aes = FALSE
  ) +

  # Projected point label
  annotate(
    "text",
    x        = proj_point['PC1'] + 2,
    y        = proj_point['PC2'] + 2,
    label    = "P",
    size     = 5,
    color    = "black",
    fontface = "bold"
  ) +

  geom_point(
    aes(
      x = vectors[['gata']]['PC1'],
      y = vectors[['gata']]['PC2']
    ),
    shape       = 23,
    fill        = cell.type.pal[["GATA6"]],
    color       = 'black',
    size        = 4,
    show.legend = FALSE
  ) +

  # K16 centroid
  geom_point(
    aes(
      x = vectors[['k17']]['PC1'],
      y = vectors[['k17']]['PC2']
    ),
    shape       = 23,
    fill        = cell.type.pal[["K17"]],
    color       = 'black',
    size        = 4,
    show.legend = FALSE
  ) +

  # Hybrid centroid
  geom_point(
    aes(
      x = vectors[['hybrid']]['PC1'],
      y = vectors[['hybrid']]['PC2']
    ),
    shape       = 23,
    fill        = cell.type.pal[["K17_GATA6"]],
    color       = 'black',
    size        = 4,
    show.legend = FALSE
  ) +

  scale_color_manual(
    values = cell.type.pal
  ) +

  scale_fill_manual(
    values = cell.type.pal
  ) +

  coord_equal(
    xlim = c(-30, 15),
    ylim = c(-15, 15)
  ) +

  labs(
    title    = glue("Hybrid as a transitional state"),
    subtitle = glue("Projection of hybrid centroid onto GATA6 --> K17 axis"),
    x        = "PC1",
    y        = "PC2"
  ) +

  theme(
    legend.title = element_blank(),
    axis.text = element_blank()
  )
dev.off()
