FROM rocker/r-base:3.6.3
RUN apt-get update && apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev/unstable

RUN mkdir /ASimulatoR /input /output
COPY ./ /ASimulatoR
WORKDIR /ASimulatoR

ENV RENV_VERSION 0.9.3-71
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org')); \
  remotes::install_github('rstudio/renv@${RENV_VERSION}'); \
  renv::consent(provided = TRUE); \
  renv::restore(); \
  devtools::install(quick = T)"

# usage: docker run --user $(id -u):$(id -g) -v input_host:/input -v ouput_host:/output_container image
ENTRYPOINT Rscript /input/runASS.R /input /output
