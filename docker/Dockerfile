FROM r-base:4.3.2

RUN apt-get update && apt-get -yq dist-upgrade && \
  apt-get install -yq build-essential  libcurl4-openssl-dev libxml2-dev  libgsl-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev && \
  apt-get clean && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages("BiocManager", repos="http://cloud.r-project.org/" )' 
RUN Rscript -e 'library(BiocManager)'
RUN Rscript -e 'BiocManager::install("DEWSeq")'
RUN Rscript -e 'library(DEWSeq)'
RUN Rscript -e 'BiocManager::install("IHW")'
RUN Rscript -e 'library(IHW)'
RUN Rscript -e 'install.packages("tidyverse", repos="http://cloud.r-project.org/" )' 
RUN Rscript -e 'library(tidyverse)'
