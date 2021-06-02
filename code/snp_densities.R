#################################
###### SNP DENSITIES ##########
###############################

# This app is designed to sumarise SNP distributions. 
# The user selects the sample they wish tor browse, and can also filter on the chromosome.
# It provides some first line QC of the data


snp_densities <- tabItem(tabName = "snp_densities", 
                         h1("SNP Densities"),
                         br(),
                           fluidRow(
                             column(3,
                                    
                                    box(title = h2("Inputs"), solidHeader = FALSE, width = 12,
                                        selectInput("bins", 
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
                                        sliderInput("position", label = h3("Genomic range"), min = 0, 
                                                    max = 260000000, value = c(1, 260000000))
                                        
                                        )

                             ),
                             column(9,
                                    p('Density plot of SNP locations'),
                                    ggvisOutput("ggvis_output_density")
                                    ),
                           ), #End of fluidRow
                           
                         
                         
                         ) # End of tabItem