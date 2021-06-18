#################################
###### SNP DENSITIES ##########
###############################


snp_densities <- tabItem(tabName = "snp_densities", 
                         h1("SNP Densities"),
                         br(),
                           fluidRow(
                             column(3,
                                    
                                    box(title = h2("Inputs"), solidHeader = FALSE, width = 12,
                                        selectInput("bins_snp", 
                                                    label = h3("Filter by Chromsome:"), 
                                                    choices = c(1:23),
                                                    selected = 1, multiple = FALSE), 
                                        sliderInput("position_snp", label = h3("Genomic range"), min = 0, 
                                                    max = 260000000, value = c(1, 260000000))
                                        
                                        )

                             ),
                             column(9,
                                    p('Density plot of SNP locations'),
                                    ggvisOutput("ggvis_output_density_snp")
                                    ),
                           ), #End of fluidRow
                           
                         
                         
                         ) # End of tabItem