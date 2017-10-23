FROM rocker/geospatial:latest
MAINTAINER "Adriaan Dokter" amd427@cornell.edu

COPY docker_install.R /opt/
RUN Rscript /opt/docker_install.R && rm /opt/docker_install.R
COPY mosaic_aws.R basemap.RData /opt/
