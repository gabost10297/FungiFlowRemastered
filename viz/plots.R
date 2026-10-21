# Pictures from steps 10 and 11.
# Sourced by viz/app.R. From the project root:
#   Rscript viz/plots.R
# writes results/viz/heatmap.pdf, same_its.pdf and per_sample.pdf (and PNGs).

if (!exists("project_root", mode = "function")) {
  source(file.path(if (file.exists("tree.R")) "." else "viz", "tree.R"))
}

read_tsv <- function(project, path) {
  read.delim(file.path(project, path), sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
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

if (sys.nframe() == 0) {
  root <- project_root()
  pictures <- list(
    heatmap = draw_heatmap(root),
    per_sample = draw_per_sample(root),
    same_its = draw_same_its(root)
  )
  out_dir <- file.path(root, "results/viz")
  for (name in names(pictures)) {
    pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
    save_tree(pictures[[name]], pdf_path, sub("\\.pdf$", ".png", pdf_path))
    message("Wrote ", pdf_path)
  }
}
