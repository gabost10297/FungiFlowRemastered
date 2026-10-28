FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# ggplot2 4 removed check_linewidth, which ggtree 3.14 (Bioconductor 3.20) still imports.
# The 2025-08-01 snapshot is the last CRAN freeze before that break.
RUN R -q -e 'options(repos = c(CRAN = "https://p3m.dev/cran/__linux__/noble/2025-08-01")); install.packages(c("shiny", "ggplot2", "ape", "BiocManager"))' \
    && R -q -e 'options(repos = c(CRAN = "https://p3m.dev/cran/__linux__/noble/2025-08-01")); BiocManager::install(c("treeio", "ggtree"), ask = FALSE, update = FALSE); if (!requireNamespace("ggtree", quietly = TRUE)) quit(status = 1)'

# The app starts pipeline containers through the host's Docker engine.
COPY --from=docker:27-cli /usr/local/bin/docker /usr/local/bin/docker

WORKDIR /data
EXPOSE 3838
CMD ["R", "-q", "-e", "shiny::runApp('viz', host = '0.0.0.0', port = 3838)"]
