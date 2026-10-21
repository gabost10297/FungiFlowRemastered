# Install the R packages for viz/tree.R and viz/app.R.
# Usage, from the project root: Rscript viz/install.R

# ggplot2 4 removed check_linewidth, which ggtree 3.14 still imports.
cran <- "https://packagemanager.posit.co/cran/2025-08-01"
install.packages(c("shiny", "ggplot2", "ape", "BiocManager"), repos = cran)
BiocManager::install(c("treeio", "ggtree"), ask = FALSE, update = FALSE)
if (!requireNamespace("ggtree", quietly = TRUE)) {
  stop("ggtree did not install. See the messages above.")
}
