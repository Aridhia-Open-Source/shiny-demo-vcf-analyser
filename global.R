##################
##### GLOBAL #####
##################

# Load libraries
library(shiny)
library(readr)
library(shinyBS)
library(shinydashboard)
library(dplyr)
library(DBI)
library(DT)
library(ggvis)
library(qqman)


# Read Data
variants <- read.csv("./data/vcf_clean.csv")
vcf_header <- read.csv("./data/vcf_header_clean.csv")
clinvar_variant_summary <- read.csv("./data/clinvar_variant_summary_clean.csv")

# Source everything on the code folder
for (file in list.files('code', full.names = TRUE)){
  source(file, local = TRUE)
}

