# syntax=docker/dockerfile:1
# biotapy's R golden-file image (contracts/r-golden-parity). Bumping a pin here is its own commit,
# followed by regenerating every file (playbooks/regenerate-golden-files).
FROM rocker/r-ver:4.5.3
# rocker pins CRAN to a dated Posit Package Manager snapshot, so install.packages() is reproducible
# and installs Ubuntu binaries.
ENV DEBIAN_FRONTEND=noninteractive
# Bioconductor packages compile from source on Linux; these are the headers phyloseq's dependency tree needs.
RUN apt-get update && apt-get install -y --no-install-recommends \
    zlib1g-dev libbz2-dev liblzma-dev libxml2-dev libglpk-dev libgmp-dev libcurl4-openssl-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'install.packages("BiocManager")' \
    && Rscript -e 'BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)'
RUN Rscript -e 'BiocManager::install("phyloseq", version = "3.22", ask = FALSE, update = FALSE)'
RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE))'
WORKDIR /work
CMD ["Rscript", "tests/r/export_golden.R"]
