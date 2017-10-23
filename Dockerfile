FROM rocker/geospatial:latest
MAINTAINER "Adriaan Dokter" amd427@cornell.edu

COPY docker_R_packages.R install.R \
  && 
RUN R CMD BATCH install.R \
  && cat install.Rout
