FROM bioconductor/bioconductor_docker:devel

MAINTAINER alan.ocallaghan@outlook.com
LABEL authors="alan.ocallaghan@@outlook.com" \
    description="Docker image containing the Bioconductor package in a Bioconductor-devel container."

WORKDIR /home/build/package

COPY . /home/build/package 

ENV R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

RUN Rscript -e "devtools::install('.', dependencies=TRUE, repos = BiocManager::repositories(), build_vignettes = TRUE)"
