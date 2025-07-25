# Use R version 4.4.1 
FROM rocker/r-ver:4.4.1

# Install system-level dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libglpk-dev \
    libhdf5-dev \
    && apt-get clean

# Install remotes and BiocManager
RUN R -e "install.packages(c('remotes', 'BiocManager'))"

# Install exact CRAN versions
RUN R -e "remotes::install_version('Seurat', version = '5.3.0')"
RUN R -e "remotes::install_version('ggplot2', version = '3.5.2')"
RUN R -e "remotes::install_version('dplyr', version = '1.1.4')"
RUN R -e "remotes::install_version('patchwork', version = '1.3.0')"
RUN R -e "remotes::install_version('cowplot', version = '1.1.3')"
RUN R -e "remotes::install_version('tidyverse', version = '2.0.0')"

# Install Harmony (Bioconductor)
RUN R -e "BiocManager::install('harmony', version = '3.18')"

# Copy your R scripts
WORKDIR /app
COPY scripts/ ./scripts

# Default command (add all the scripts)
#CMD ["Rscript", "scripts/run_all.R"]
