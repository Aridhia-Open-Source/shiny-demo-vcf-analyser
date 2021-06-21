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
library(RPostgreSQL)
library(sqldf)
library(ggbio)
library(IRanges)
library(DT)
library(GenomicRanges)
library(ggvis)
library(qqman)


# Set color=
getOption("biovizBase")$cytobandColor

# Load data from biovizBase package
data(hg19IdeogramCyto, package = "biovizBase")
data(ideoCyto, package = "biovizBase")
data(darned_hg19_subset500, package = "biovizBase")
biovizBase::isIdeogram(ideoCyto$hg19)



# Read Data
sample_data <- read.csv("./data/vcf_sample_data.csv")
variants <- read.csv("./data/vcf_data.csv")
vcf_header <- read.csv("./data/vcf_header.csv")
clinvar_variant_summary <- read.csv("./data/clinvar_variant_summary.csv")

# Source everything on the code folder
for (file in list.files('code', full.names = TRUE)){
  source(file, local = TRUE)
}

