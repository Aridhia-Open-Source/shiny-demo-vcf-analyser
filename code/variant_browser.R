browser_tab <- 
  
  tabItem(tabName = "variant_browser", 
          h1("Explore the variants and view annotation information"),
          br(),
          fluidRow(
            column(2,
                shinydashboard::box(title = h2("Inputs"), solidHeader = FALSE, width = 12,
                       # Drop down menu to select sample ID
                       selectInput(
                         inputId='selected_sample_browser',
                         label='Sample Id',
                         choices=as.character(unique(variants[,c('sample_id')])),
                         selected = 'OTHER'
                       ),
                       selectInput(
                         inputId = "chromosome_browser",
                         label = "Filter by chromosome",
                         choices= c("All", as.character(unique(variants[,c('CHROM')]))),
                         selected = 'All'
                         ),
                       uiOutput("coordinates_browser")
                  ), # End of box
            ), # End of column1
            column(6,
                   # Table that prints filtered data by sample ID
                   DTOutput("variants_table"),
                   verbatimTextOutput('selected_idx'),
                   div(id = "demo", class = "collapse", tableOutput('v_header'))
            ), # End of column2
            column(4,
                   shinydashboard::box(
                     title = h2("Selection Information"), solidHeader = FALSE, width = 12,
                     h4(htmlOutput('v_coords')),
                     h4(htmlOutput('v_ensemble_link')), 
                     h4(htmlOutput('v_diff_small')),
                     h4(textOutput('v_quality')),
                     h4(textOutput('v_clinvar_count')),
                     h4(htmlOutput("v_collapse_button"))
                   ), #End of box2
                   br(),
                   ), # End of column3
            
          ) # End of fluid Row
    
    
    ) # End of tab
             

      