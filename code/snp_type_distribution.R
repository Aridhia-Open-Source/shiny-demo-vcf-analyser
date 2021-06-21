#################################
##### SNP TYPE DISTRIBUTION #####
#################################

# This app is designed to show the proportions of SNP types across chromosomes.
# It is part of a suite of apps being used to perform QC and exploratory data analysis.
# The data is presented as a heatmap.


snp_type_distributions <- 
  tabItem(tabName = "snp_type_distribution", 
                         h1("SNP Type Distribution"),
                         br(),
                         fluidRow(
                           column(3,
                                  
                                  box(title = h2("Inputs"), solidHeader = FALSE, width = 12,
                                      selectInput("heatBins", 
                                                  label = h3("Chromsome"), 
                                                  choices = list("Chromosome 1" = "chr1", "Chromosome 2" = "chr2",
                                                                 "Chromosome 3" = "chr3", "Chromosome 4" = "chr4",
                                                                 "Chromosome 5" = "chr5", "Chromosome 6" = "chr6",
                                                                 "Chromosome 7" = "chr7", "Chromosome 8" = "chr8",
                                                                 "Chromosome 9" = "chr9", "Chromosome 10" = "chr10",
                                                                 "Chromosome 11" = "chr11", "Chromosome 12" = "chr12",
                                                                 "Chromosome 13" = "chr13", "Chromosome 14" = "chr14",
                                                                 "Chromosome 15" = "chr15", "Chromosome 16" = "chr16",
                                                                 "Chromosome 17" = "chr17", "Chromosome 18" = "chr18",
                                                                 "Chromosome 19" = "chr19", "Chromosome 20" = "chr20",
                                                                 "Chromosome 21" = "chr21", "Chromosome X" = "chrX", "Chromosome Y" = "chrY"
                                                  ),
                                                  selected = "chr1", multiple = FALSE), 
                                      sliderInput("heatPosition", label = h3("Genomic range"), min = 0, 
                                                  max = 260000000, value = c(1, 260000000))
                                      
                                  )
                                  
                           ),
                           column(9,
                                  ggvisOutput("ggvis_output_heatmap")
                           ),
                         ), #End of fluidRow
                         
                         
                         
) # End of tabItem