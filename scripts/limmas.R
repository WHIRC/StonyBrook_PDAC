source("./scripts/global.R")
source("./scripts/limma_analysis_function.R")

pairs <- list(
  c("DUAL_NEG", "GATA6"),
  c("DUAL_NEG", "K17_GATA6"),
  c("DUAL_NEG", "K17"),
  c("GATA6", "K17_GATA6"),
  c("GATA6", "K17"),
  c("K17_GATA6", "K17")
)

limmas <- list()

for (pair in pairs) {

  this.meta <- metadata %>% filter(cell.type %in% pair)

  this.mtx <- mtx[, row.names(this.meta)]
  stopifnot(colnames(this.mtx) == row.names(this.meta))

  numerator <- pair[1]
  denominator <- pair[2]

  L <- limma_analysis(
    mtx = this.mtx,
    meta = this.meta,
    variable = "cell.type",
    numerator = numerator,
    denominator = denominator,
    gsea_gmt_files = "./input/msigdb_v2023.2.Hs_GMTs",
    output_dir = file.path("output", "supervised_analysis")
  )
  L <- L[["limma.results"]] %>%
    mutate(X = glue("{Accession}|{Gene}")) %>%
    column_to_rownames("X")
  limmas[[length(limmas) + 1]] <- L

}

# Hybrid vs all others --------------------------------------------------------
this.meta <- metadata %>%
  mutate(cell.type = ifelse(
    cell.type == "K17_GATA6",
    "K17_GATA6",
    "All_others"
  ))

this.mtx <- mtx[, row.names(this.meta)]
stopifnot(colnames(this.mtx) == row.names(this.meta))

numerator <- "K17_GATA6"
denominator <- "All_others"

L <- limma_analysis(
  mtx = this.mtx,
  meta = this.meta,
  variable = "cell.type",
  numerator = numerator,
  denominator = denominator,
  gsea_gmt_files = "./input/msigdb_v2023.2.Hs_GMTs",
  output_dir = file.path("output", "supervised_analysis")
)
L <- L[["limma.results"]] %>%
  mutate(X = glue("{Accession}|{Gene}")) %>%
  column_to_rownames("X")
limmas[[length(limmas) + 1]] <- L

l.names <- lapply(
  pairs, function(P) {paste(P, collapse = ".vs.")}
)

l.names[[length(l.names)+1]] <- "K17_GATA6.vs.All_others"

names(limmas) <- l.names

# which vs all others ---------------------------------------------------------
for (which in c("GATA6", "K17", "DUAL_NEG")) {

  this.meta <- metadata %>%
    mutate(cell.type = ifelse(
      cell.type == which,
      which,
      "All_others"
    ))

  this.mtx <- mtx[, row.names(this.meta)]
  stopifnot(colnames(this.mtx) == row.names(this.meta))

  numerator <- which
  denominator <- "All_others"

  L <- limma_analysis(
    mtx = this.mtx,
    meta = this.meta,
    variable = "cell.type",
    numerator = numerator,
    denominator = denominator,
    gsea_gmt_files = "./input/msigdb_v2023.2.Hs_GMTs",
    output_dir = file.path("output", "supervised_analysis")
  )
  L <- L[["limma.results"]] %>%
    mutate(X = glue("{Accession}|{Gene}")) %>%
    column_to_rownames("X")
  limmas[[length(limmas) + 1]] <- L

  l.names[[length(l.names)+1]] <- glue("{which}.vs.All_others")

  names(limmas) <- l.names

}

# -----------------------------------------------------------------------------
library(ComplexUpset)

sig_lists <- lapply(limmas, function(x) {
  rownames(x[x$P.Value < 0.01, ])
})

# Get union of all sig genes
all_proteins <- unique(unlist(sig_lists))

# Build binary presence/absence df
upset_df <- as.data.frame(
  lapply(sig_lists, function(x) as.integer(all_proteins %in% x)),
  row.names = all_proteins
)

upset_df %>%
  tibble::rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "comparison", values_to = "sig") %>%
  filter(sig == 1) %>%
  group_by(Protein) %>%
  summarise(Intersection = paste(sort(comparison), collapse = " & ")) %>%
  mutate(
    Protein.Accession = sub("\\|.*", "", Protein),
    Gene.Symbol = sub(".*\\|", "", Protein)
  ) %>%
  select(-Protein) %>%
  write.csv(
    x = .,
    file = file.path(
      "output",
      "supervised_analysis",
      "upset_flat_file.csv"
    ),
    row.names = FALSE
  )

upset_plot <- upset(
  upset_df,
  intersect = names(upset_df),
  keep_empty_groups = TRUE,
  # max_degree = 2,
  min_degree = 2,
  min_size = 3,
  # set_size = FALSE,
  sort_sets = "descending", # Keep ordering of Comparison_names
  height_ratio = 0.62, # Inverse golden ratio of matrix:bar plot
  base_annotations = list( # Change y axis label
    "Intersection size" = (
      intersection_size(
        text_colors = c( # Change all numbers to black
          on_background = "black",
          on_bar = "white"
        ),
        text = list(
          angle = 90,
          vjust = 0.5,
          hjust = 0
        )
      ) +
        ylab("Number of coaltered proteins")
    )
  ),
  themes = upset_default_themes(
    # panel.background = element_rect(fill = "grey96"),
    axis.text = element_text(
      size = 11,
      face = "bold",
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
)

png(
  filename = file.path(
    "output", "supervised_analysis", "upset_plot_limmas.png"
  ),
  width = 10, height = 5, units = "in", res = 500
)
upset_plot
dev.off()



# -----------------------------------------------------------------------------
# For each limma, extract proteins upregulated in numerator and denominator
group_sig_proteins <- list()

for (comparison_name in names(limmas)) {
  limma_data <- limmas[[comparison_name]]
  
  # Parse groups from comparison name
  groups <- str_split(comparison_name, "\\.vs\\.")[[1]]
  numerator_group <- groups[1]
  denominator_group <- groups[2]
  
  # Proteins upregulated in numerator (logFC > 0)
  numerator_proteins <- limma_data %>%
    filter(P.Value < 0.01) %>%
    filter(across(starts_with("logFC"), ~ . > 0)) %>%
    pull(Accession)
  
  # Proteins upregulated in denominator (logFC < 0)
  denominator_proteins <- limma_data %>%
    filter(P.Value < 0.01) %>%
    filter(across(starts_with("logFC"), ~ . < 0)) %>%
    pull(Accession)
  
  # Accumulate into group lists
  group_sig_proteins[[numerator_group]] <- union(
    group_sig_proteins[[numerator_group]], numerator_proteins
  )
  group_sig_proteins[[denominator_group]] <- union(
    group_sig_proteins[[denominator_group]], denominator_proteins
  )
}

group_sig_proteins[["All_others"]] <- NULL

# Build presence/absence matrix
all_proteins <- unique(unlist(group_sig_proteins))

upset_df <- as.data.frame(
  lapply(group_sig_proteins, function(protein_set) as.integer(all_proteins %in% protein_set)),
  row.names = all_proteins
)

upset_df %>%
  tibble::rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "cell.type", values_to = "sig") %>%
  filter(sig == 1) %>%
  group_by(Protein) %>%
  summarise(Intersection = paste(sort(cell.type), collapse = " & ")) %>%
  mutate(
    Protein.Accession = sub("\\|.*", "", Protein),
    Gene.Symbol = sub(".*\\|", "", Protein)
  ) %>%
  select(-Protein) %>%
  write.csv(
    x = .,
    file = file.path(
      "output",
      "supervised_analysis",
      "upset_presence_absence.csv"
    ),
    row.names = FALSE
  )

upset_plot <- upset(
  upset_df,
  intersect = names(upset_df),
  keep_empty_groups = TRUE,
  # max_degree = 2,
  min_degree = 2,
  min_size = 3,
  # set_size = FALSE,
  sort_sets = "descending", # Keep ordering of Comparison_names
  height_ratio = 0.62, # Inverse golden ratio of matrix:bar plot
  base_annotations = list( # Change y axis label
    "Intersection size" = (
      intersection_size(
        text_colors = c( # Change all numbers to black
          on_background = "black",
          on_bar = "white"
        ),
        text = list(
          angle = 90,
          vjust = 0.5,
          hjust = 0
        )
      ) +
        ylab("Number of coaltered proteins")
    )
  ),
  themes = upset_default_themes(
    # panel.background = element_rect(fill = "grey96"),
    axis.text = element_text(
      size = 11,
      face = "bold",
      angle = 0,
      vjust = 0.5,
      hjust = 0.5,
      color = "black"
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
)

png(
  filename = file.path(
    "output", "supervised_analysis", "upset_plot_presence_absence.png"
  ),
  width = 10, height = 5, units = "in", res = 500
)
upset_plot
dev.off()
