# Example shiny app docker file


# get shiny server and R from the rocker project
FROM rocker/shiny:4.0.5

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev
  

# install R packages required 
# Change the packages list to suit your needs
RUN R -e 'install.packages(c(\
              "shiny", \
              "shinydashboard", \
              "ggplot2", \
              "plotly", \
              "ggcorrplot", \
              "shinythemes", \

              "tibble", \
              "dplyr", \
              "colorRamps", \
              "DT", \
              "seqRFLP", \
              "seqinr", \
              "xlsx", \
              "openxlsx", \
              "BiocManager" \
            ), \
            repos="https://packagemanager.rstudio.com/cran/__linux__/focal/2021-04-23"\
          )'

RUN R -e 'BiocManager::install("msa")'
RUN R -e 'BiocManager::install("Biostrings")'

# copy the app directory into the image
COPY ./shiny-app/* /srv/shiny-server/

# run app
CMD ["/usr/bin/shiny-server"]
