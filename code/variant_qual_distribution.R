####################################
##### VARIANT QUAL DISTRIBUION #####
####################################

# This app shows the distribution across chromosomes of variants failing QC, based on the QUAL field in the VCF.



variant_qual_distribution <- 
  tabItem(tabName = "variant_qual_distribution", 
          h1("Variant Quality Distribution"),
          br(),
          fluidRow(
            column(3,
                   
                   box(title = h2("Inputs"), solidHeader = FALSE, width = 12,
                       selectInput("QBins", 
                                   label = h3("Chromsome"), 
                                   choices = c(1:25),
                                   selected = 1:25, multiple = TRUE)
                       
                   )
                   
            ),
            column(9,
                   plotOutput("manhattan")
            ),
          ), #End of fluidRow
          
          
          
  ) # End of tabItem