FROM rocker/r-base:4.3.1

# Set OCI labels
LABEL org.opencontainers.image.authors="Pasha Mazin <pavel.mazin@sanger.ac.uk>"
LABEL org.opencontainers.image.title="reticulate-free single cell format conversion"
LABEL org.opencontainers.image.description="Package to load scanpy h5ad into R as list, SingleCellExperiment or Seurat object."
LABEL org.opencontainers.image.source="https://github.com/cellgeni/schard"
LABEL org.opencontainers.image.licenses="MIT"

RUN apt-get update && \
    apt-get install -y make zlib1g-dev libcurl4-openssl-dev libssl-dev \
     libfontconfig1-dev libfreetype6-dev libfribidi-dev libharfbuzz-dev \
     libjpeg-dev libpng-dev libtiff-dev pandoc libicu-dev libxml2-dev git \
     libgit2-dev zlib1g-dev libicu-dev libcurl4-openssl-dev \
     libglpk-dev libgmp3-dev libpng-dev python3
        
RUN Rscript -e 'install.packages(c("devtools","BiocManager"))'
RUN Rscript -e 'devtools::install_github("satijalab/seurat", "seurat5", quiet = TRUE)'
RUN Rscript -e 'BiocManager::install(c("rhdf5", "SingleCellExperiment"))'
RUN Rscript -e 'devtools::install_github("cellgeni/schard")'
