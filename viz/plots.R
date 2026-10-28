# Pictures from steps 10 and 11.
# Sourced by viz/app.R. From the project root:
#   Rscript viz/plots.R
# writes heatmap, per_sample, same_its, tree_circular and pca (PDF and PNG) to results/viz/.

if (!exists("project_root", mode = "function")) {
  source(file.path(if (file.exists("tree.R")) "." else "viz", "tree.R"))
}

read_tsv <- function(project, path) {
  file <- file.path(project, path)
  if (!file.exists(file)) stop("Not found: ", path, ". This run has not reached that step yet.")
  read.delim(file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

draw_heatmap <- function(project = project_root(), top_n = 20) {
  org <- read_tsv(project, "results/10_organisms/organisms.tsv")
  org <- org[order(org$position), ]
  org <- head(org, top_n)
  samples <- grep("^barcode", names(org), value = TRUE)

  long <- do.call(rbind, lapply(samples, function(sample) {
    data.frame(organism = org$organism, sample = sample, reads = org[[sample]],
               stringsAsFactors = FALSE)
  }))
  long$organism <- factor(long$organism, levels = rev(org$organism))
  long$sample <- factor(long$sample, levels = samples)
  long$reads_present <- ifelse(long$reads == 0, NA, long$reads)

  ggplot(long, aes(sample, organism, fill = reads_present)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_viridis_c(trans = "log10", name = "Reads", na.value = "grey92") +
    scale_x_discrete(position = "top") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 0),
          panel.grid = element_blank()) +
    labs(x = NULL, y = NULL,
         caption = "Most abundant organisms. Grey means no reads in that sample. Color is on a log scale.")
}

draw_per_sample <- function(project = project_root()) {
  counts <- read_tsv(project, "results/10_organisms/per_sample.tsv")
  shared <- counts$organisms - counts$only_in_this_sample
  long <- rbind(
    data.frame(sample = counts$sample, n = shared, part = "Also in other samples"),
    data.frame(sample = counts$sample, n = counts$only_in_this_sample, part = "Only in this sample")
  )
  long$sample <- factor(long$sample, levels = counts$sample)
  long$part <- factor(long$part, levels = c("Also in other samples", "Only in this sample"))

  ggplot(long, aes(sample, n, fill = part)) +
    geom_col(width = 0.72) +
    geom_text(data = counts, aes(sample, organisms, label = organisms),
              inherit.aes = FALSE, vjust = -0.4, size = 3.4) +
    scale_fill_manual(values = c("Also in other samples" = "#0072B2",
                                 "Only in this sample" = "#D55E00"),
                      name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1),
          legend.position = "bottom",
          panel.grid.major.x = element_blank()) +
    labs(x = NULL, y = "Organisms",
         caption = paste("Fungal organisms at high or low confidence.",
                         "Samples have different numbers of reads, so the heights are not a direct comparison."))
}

draw_same_its <- function(project = project_root()) {
  same <- read_tsv(project, "results/11_alignments/summary.tsv")
  same$pct <- as.numeric(same$pct_reads_at_least_98.5)
  same$organism <- factor(same$organism, levels = rev(same$organism))
  same$group <- ifelse(same$sequences < 2, "Only one sequence",
                       ifelse(same$pct >= 98.5, "At least 98.5% of reads",
                              "Split across ITS variants"))
  same$group <- factor(same$group, levels = c("At least 98.5% of reads",
                                              "Split across ITS variants",
                                              "Only one sequence"))

  ggplot(same, aes(organism, pct, fill = group)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sprintf("%.1f", pct)), hjust = -0.15, size = 3.3) +
    scale_fill_manual(values = c("At least 98.5% of reads" = "#0072B2",
                                 "Split across ITS variants" = "#D55E00",
                                 "Only one sequence" = "grey70"),
                      name = NULL) +
    scale_y_continuous(limits = c(0, 118), expand = c(0, 0), breaks = seq(0, 100, 25)) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
    labs(x = NULL, y = "Reads within 98.5% of the representative (%)",
         caption = paste("Each consensus compared with its organism's representative.",
                         "A grey bar has nothing to compare: that organism has one sequence."))
}

draw_circular <- function(project = project_root(), label_size = 2.6) {
  org <- read_tsv(project, "results/10_organisms/organisms.tsv")
  tree_file <- file.path(project, "results/12_tree/tree.treefile")
  if (!file.exists(tree_file)) stop("Tree not found. This run has not reached step 12 yet.")
  tree <- read.tree(tree_file)
  org$label <- tip_key(org$organism)
  samples <- grep("^barcode", names(org), value = TRUE)
  rows <- org[match(tree$tip.label, org$label), ]

  # Zero reads stay empty cells, not the darkest color.
  reads <- as.data.frame(lapply(rows[samples], function(x) ifelse(x > 0, x, NA)))
  rownames(reads) <- tree$tip.label
  names(reads) <- sub("^barcode", "", samples)

  tips <- data.frame(label = tree$tip.label, organism = rows$organism,
                     phylum = rows$phylum, stringsAsFactors = FALSE)
  depth <- max(node.depth.edgelength(tree))
  ring_offset <- 0.04 * depth
  ring_width <- 0.55

  plot <- ggtree(tree, layout = "fan", open.angle = 18, linewidth = 0.35) %<+% tips +
    geom_tippoint(aes(color = phylum), size = 1.8) +
    scale_color_manual(values = PHYLUM_COLORS, name = "Phylum", na.value = "grey50")

  plot <- gheatmap(plot, reads, offset = ring_offset, width = ring_width,
                   color = "white", colnames = FALSE) +
    scale_fill_viridis_c(trans = "log10", name = "Reads", na.value = "grey92")

  plot +
    geom_tiplab(aes(label = organism), size = label_size, align = TRUE, linetype = NA,
                offset = ring_offset + ring_width * depth + 0.04 * depth) +
    xlim(0, depth * 2.6) +
    theme(legend.position = "right",
          plot.caption = element_text(size = 9, color = "grey30", hjust = 0)) +
    labs(caption = paste0("Step 12 tree with reads per barcode as rings, inner to outer: ",
                          paste(samples, collapse = ", "),
                          ".\nGrey is no reads. Color is on a log scale."))
}

draw_pca <- function(project = project_root(), loadings = 6) {
  org <- read_tsv(project, "results/10_organisms/organisms.tsv")
  samples <- grep("^barcode", names(org), value = TRUE)
  counts <- t(as.matrix(org[samples]))
  colnames(counts) <- org$organism

  # Hellinger: square root of each organism's share of the sample.
  # Read depth differs between barcodes, so raw counts would mostly measure depth.
  hellinger <- sqrt(counts / rowSums(counts))
  pca <- prcomp(hellinger, center = TRUE, scale. = FALSE)
  explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  scores <- data.frame(sample = rownames(counts), PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                       organisms = rowSums(counts > 0), stringsAsFactors = FALSE)

  arrows <- as.data.frame(pca$rotation[, 1:2])
  arrows$organism <- rownames(arrows)
  arrows <- arrows[order(rowSums(arrows[, 1:2]^2), decreasing = TRUE), ][seq_len(loadings), ]
  stretch <- 0.8 * max(abs(scores[, c("PC1", "PC2")])) / max(abs(arrows[, c("PC1", "PC2")]))
  arrows$x <- arrows$PC1 * stretch
  arrows$y <- arrows$PC2 * stretch

  ggplot(scores, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, color = "grey85") +
    geom_vline(xintercept = 0, color = "grey85") +
    geom_segment(data = arrows, aes(x = 0, y = 0, xend = x, yend = y),
                 arrow = arrow(length = unit(0.18, "cm")), color = "grey55") +
    geom_text(data = arrows, aes(x * 1.05, y * 1.08, label = organism,
                                 hjust = ifelse(x < 0, 1, 0)),
              color = "grey35", size = 3, fontface = "italic") +
    scale_x_continuous(expand = expansion(mult = 0.3)) +
    geom_point(aes(size = organisms), color = "#0072B2", alpha = 0.85) +
    geom_text(aes(label = sample), vjust = -1.1, size = 3.6) +
    scale_size_continuous(range = c(3, 8), name = "Organisms") +
    coord_equal() +
    theme_minimal(base_size = 12) +
    labs(x = sprintf("PC1 (%.1f%%)", explained[1]),
         y = sprintf("PC2 (%.1f%%)", explained[2]),
         caption = paste0("PCA of Hellinger-transformed organism shares (", ncol(counts),
                          " organisms). Arrows: the ", loadings,
                          " organisms that separate the samples most."))
}

if (sys.nframe() == 0) {
  root <- project_root()
  pictures <- list(
    heatmap = draw_heatmap(root),
    per_sample = draw_per_sample(root),
    same_its = draw_same_its(root),
    tree_circular = draw_circular(root),
    pca = draw_pca(root)
  )
  out_dir <- file.path(root, "results/viz")
  for (name in names(pictures)) {
    pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
    save_tree(pictures[[name]], pdf_path, sub("\\.pdf$", ".png", pdf_path))
    message("Wrote ", pdf_path)
  }
}
