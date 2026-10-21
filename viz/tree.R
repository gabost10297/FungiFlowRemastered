# Draw the step 12 tree.
#
# Tip names in the Newick file are the organism names with spaces turned into
# "_" and dots removed (the same rule as scripts/12_tree.sh). organisms.tsv
# supplies the readable name, phylum, class and read count.
#
# Run from the project root, after viz/install.R (or inside the fungiflow-viz image):
#   Rscript viz/tree.R
# writes results/viz/tree.pdf and results/viz/tree.png.

library(ape)
library(ggtree)
library(ggplot2)

PHYLUM_COLORS <- c(Ascomycota = "#0072B2", Basidiomycota = "#D55E00")

project_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]))
    return(normalizePath(file.path(dirname(script), "..")))
  }
  if (file.exists("viz/tree.R")) return(normalizePath(getwd()))
  if (file.exists("tree.R") && file.exists("app.R")) {
    return(normalizePath(file.path(getwd(), "..")))
  }
  normalizePath(getwd())
}

tip_key <- function(organism) {
  gsub("\\.", "", gsub(" ", "_", organism))
}

model_name <- function(project) {
  path <- file.path(project, "results/12_tree/tree.iqtree")
  if (!file.exists(path)) return("IQ-TREE")
  line <- grep("Best-fit model according to BIC:", readLines(path, warn = FALSE), value = TRUE)
  if (!length(line)) return("IQ-TREE")
  sub(".*BIC:\\s*", "", line[[1]])
}

draw_tree <- function(project = project_root(), bootstrap = 70,
                      color_by = "phylum", label_size = 3.2) {
  tree_file <- file.path(project, "results/12_tree/tree.treefile")
  org_file <- file.path(project, "results/10_organisms/organisms.tsv")
  if (!file.exists(tree_file)) stop("Tree not found: ", tree_file, ". Run script 12 first.")
  if (!file.exists(org_file)) stop("Organism table not found: ", org_file, ". Run script 10 first.")
  if (!color_by %in% c("phylum", "class", "reads")) {
    stop("color_by must be phylum, class or reads.")
  }

  tree <- read.tree(tree_file)
  org <- read.delim(org_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  org$label <- tip_key(org$organism)

  tips <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
  keep <- c("label", "organism", "phylum", "class", "reads")
  tips <- merge(tips, org[, keep], by = "label", all.x = TRUE)
  missing <- tips$label[is.na(tips$organism)]
  if (length(missing)) {
    message("No organism row for: ", paste(missing, collapse = ", "))
    tips$organism[is.na(tips$organism)] <- gsub("_", " ", missing)
  }
  tips$reads <- as.numeric(tips$reads)

  # A column the color scale can use for every option.
  tips$color_value <- if (color_by == "reads") pmax(tips$reads, 1) else tips[[color_by]]

  support <- "Bootstrap support from 1000 ultrafast replicates"
  caption <- paste0(
    support, ". Values below ", bootstrap, " are grey dots. ",
    "Scale: substitutions per site. Model: ", model_name(project), "."
  )

  plot <- ggtree(tree, ladderize = TRUE, linewidth = 0.4) %<+% tips +
    geom_tippoint(aes(color = color_value), size = 2.4) +
    geom_tiplab(aes(label = organism, color = color_value),
                size = label_size, offset = 0.012, show.legend = FALSE) +
    geom_nodelab(
      aes(subset = !is.na(as.numeric(label)) & as.numeric(label) >= !!bootstrap,
          label = label),
      size = 2.5, hjust = 1.15, vjust = -0.35, color = "grey20"
    ) +
    geom_point2(
      aes(subset = !isTip & !is.na(as.numeric(label)) & as.numeric(label) < !!bootstrap),
      color = "grey70", size = 1.8
    ) +
    hexpand(0.35) +
    theme_tree2() +
    theme(
      legend.position = "right",
      plot.caption = element_text(size = 9, color = "grey30", hjust = 0)
    ) +
    labs(caption = caption)

  if (color_by == "reads") {
    plot <- plot + scale_color_viridis_c(trans = "log10", name = "Reads")
  } else if (color_by == "phylum") {
    plot <- plot + scale_color_manual(values = PHYLUM_COLORS, name = "Phylum", na.value = "grey50")
  } else {
    plot <- plot + scale_color_discrete(name = "Class")
  }
  plot
}

save_tree <- function(plot, pdf_path, png_path = NULL) {
  dir.create(dirname(pdf_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(pdf_path, plot, width = 12, height = 9, device = "pdf")
  if (!is.null(png_path)) {
    ggsave(png_path, plot, width = 12, height = 9, dpi = 130, device = "png")
  }
  invisible(pdf_path)
}

if (sys.nframe() == 0) {
  root <- project_root()
  picture <- draw_tree(root)
  out <- file.path(root, "results/viz/tree.pdf")
  png <- sub("\\.pdf$", ".png", out)
  save_tree(picture, out, png)
  message("Wrote ", out, " and ", png)
}
