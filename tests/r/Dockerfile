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
# picante (CRAN, same P3M snapshot): the Faith PD golden file. vegan and ape come with phyloseq.
RUN Rscript -e 'install.packages("picante")'
# philr (Bioconductor, GPL-3): the pp.philr golden file; biotapy never calls it. Its ggtree -> treeio -> fs
# chain installs fs as a P3M binary that links libuv at load time.
RUN apt-get update && apt-get install -y --no-install-recommends libuv1 && rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'BiocManager::install("philr", version = "3.22", ask = FALSE, update = FALSE)'
# MicrobiomeStat (CRAN, GPL-3): the da.linda golden file; biotapy reimplements linda and never calls it.
RUN Rscript -e 'install.packages("MicrobiomeStat")'
# ANCOMBC (Bioconductor, Artistic-2.0): the da.ancombc2 golden file; biotapy calls scikit-bio's ancombc2 instead.
# ANCOMBC 2.12.0 imports CVXR::solve, which CVXR 1.8 (the CRAN snapshot's version) removed: pin CVXR 1.0-15, the last
# release before it, built from the CRAN archive after its dependencies. The gsl binary among ANCOMBC's dependencies
# links libgsl at load time.
RUN apt-get update && apt-get install -y --no-install-recommends libgsl27 && rm -rf /var/lib/apt/lists/*
RUN Rscript -e 'install.packages(c("Rmpfr", "gmp", "ECOSolveR", "scs", "osqp", "bit64", "cli", "Rcpp", "RcppEigen"))' \
    && Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/CVXR/CVXR_1.0-15.tar.gz", repos = NULL, type = "source")'
RUN Rscript -e 'BiocManager::install("ANCOMBC", version = "3.22", ask = FALSE, update = FALSE)'
# ALDEx2 (Bioconductor, GPL-3): the da.aldex2 golden file and the bridge's r tests; biotapy calls it through rpy2.
RUN Rscript -e 'BiocManager::install("ALDEx2", version = "3.22", ask = FALSE, update = FALSE)'
RUN Rscript -e 'stopifnot(requireNamespace("phyloseq", quietly = TRUE), requireNamespace("Biostrings", quietly = TRUE), requireNamespace("picante", quietly = TRUE), requireNamespace("philr", quietly = TRUE), requireNamespace("MicrobiomeStat", quietly = TRUE), requireNamespace("ANCOMBC", quietly = TRUE), requireNamespace("ALDEx2", quietly = TRUE))'
WORKDIR /work
CMD ["Rscript", "tests/r/export_golden.R"]
