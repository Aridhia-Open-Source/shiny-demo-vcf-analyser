#######################
##### DEPENDECIES #####
#######################


# Packages needed to run the app
packages <- c("shiny", "dplyr", "DBI", "RPostgreSQL", "sqldf", "GenomicRanges", "shinyBS", "ellipsis", "shinydashboard")


# Install packages if not already installed
package_install <- function(x){
  
  for (i in x){
    # Check if the package is installed
    if (!require(i, character.only = TRUE)){
      install.packages(i, dependencies = TRUE)
    }
  }
  
}


package_install(packages)


if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install("ggbio")
  BiocManager::install("IRanges")
}

