######################
####### SERVER #######
######################

#donors <- unique(donors_sample$icgc_donor_id)


server <- function(input, output, session) {
  
  ###########################
  ##### VARIANT BROWSER #####
  ###########################
  
  start.loc <- 0
  end.loc <- 150000
  chr1.gr <- GenomicRanges::GRanges("X", IRanges(start.loc, end.loc))
  
  # Filter data by the sample selector
  sample <- reactive({
    
    # Filter by chromosome
    if (input$chromosome_browser != "All") {
      variants[variants$CHROM == input$chromosome_browser, c('sample_id','CHROM','POS','ID','REF','ALT','QUAL','FILTER')]
    } else {
      variants[variants$sample_id == input$selected_sample_browser, c('sample_id','CHROM','POS','ID','REF','ALT','QUAL','FILTER')]
      
    }
    
  })
  
  output$coordinates_browser <- renderUI({
    sliderInput(
      inputId = "coordinates_browser",
      label = "Filter by position",
      min = min(sample()$POS),
      max = max(sample()$POS),
      value = c(min(sample()$POS), max(sample()$POS)),
      round = TRUE
    )
  })
  
  available_sample <- reactive({

      sample()[sample()$POS >= input$coordinates_browser[1] & sample()$POS <= input$coordinates_browser[2], c('sample_id','CHROM','POS','ID','REF','ALT','QUAL','FILTER')]
    
  })
  
  available_header <- reactive({
    
    vcf_header[vcf_header$sample_id==input$selected_sample_browser,c('sample_id','category','id','number','type','description','assembly','length')]
    
  })
  
  output$variants_table <- DT::renderDataTable({available_sample()},
                                               selection = 'single', server = TRUE)
  
 selected_id <- reactive ({
   req(input$variants_table_rows_selected)
   variants[input$variants_table_rows_selected, c('CHROM', 'POS', 'REF', 'ALT')]
  })
  
  output$selected_idx <- renderPrint({
    selected_id()
  })
  
  add_chr <- reactive({  
    
    if ((grepl("chr", variants[input$variants_table_rows_selected, c('CHROM')], ignore.case = T) == F & is.null(input$variants_table_rows_selected) == F) && (!is.null(input$variants_table_rows_selected))){
      paste0("chr", variants[input$variants_table_rows_selected, c('CHROM')])
    }
    
  })
  
  output$v_collapse_button = renderText(
    "<button type='button' class='btn btn-primary' data-toggle='collapse' data-target='#demo'>View VCF Header</button>"
  )
  
  # Selected Coordinates 
  output$v_coords = renderText({
    ifelse(is.null(input$variants_table_rows_selected), 'Select a Row to display coordinates',
           paste0(
             'Coordinates: ', 
             variants[input$variants_table_rows_selected, c('CHROM')],
             ':',
             variants[input$variants_table_rows_selected, c('POS')],
             ' (', (input$variants_table_rows_selected), ')',
             '</span>'
           )
    )
  })
  
  # Link to ensemble
  output$v_ensemble_link = renderText(
    
    ifelse(is.null(input$variants_table_rows_selected), 'No variants selected.',
           paste0(
             '<a target="_new" href="http://grch37.ensembl.org/Homo_sapiens/Location/View?r=',
             
             substring(variants[input$variants_table_rows_selected, c('CHROM')], 4, 4),
             ':',
             variants[input$variants_table_rows_selected, c('POS')],
             '-',
             variants[input$variants_table_rows_selected, c('POS')],
             '">View on Ensembl</a>'
           )        
    )
  )
  
  
  output$v_header <- renderTable({
    available_header()
  })
  
  
  # TODO - color coding
  output$v_diff <- renderText(paste0(
    variants[input$variants_table_rows_selected, c('REF')], 
    '<span style="color:blue;">', 
    variants[input$variants_table_rows_selected, c('ALT')], 
    '</span>'));
  
  output$v_diff_small = renderText(paste0(
    'REF: ',
    ifelse(is.null(input$variants_table_rows_selected), '-', variants[input$variants_table_rows_selected, c('REF')]), 
    ' | <span style="color:blue;">ALT: ', 
    ifelse(is.null(input$variants_table_rows_selected), '-', variants[input$variants_table_rows_selected, c('ALT')]), 
    '</span>'));
  
  #TODO would be nice to render the icon too, conditionally?
  output$v_quality = renderText(paste0('Quality: ', ifelse(is.null(input$variants_table_rows_selected), 'No row Selected', variants[input$variants_table_rows_selected, c('QUAL')])))
  
  #Anything in clinvar?
  cvc <- reactive({
    req(input$variants_table_rows_selected)
    
    chr <- substring(variants[input$variants_table_rows_selected, c('CHROM')], 4)
    pos <- variants[input$variants_table_rows_selected, c('POS')]
    
    t <-  clinvar_variant_summary[clinvar_variant_summary$chromosome == chr & clinvar_variant_summary$start >= pos & clinvar_variant_summary$start <= pos, ]
    rs1 <- nrow(t)
    
    if(nrow(t) == 0) {
      'None'
    }
    else {
      nrow(t)
    }
  }) ;
  
  # Number of clinical variants on the position
  output$v_clinvar_count <- renderText(paste0('Clinvar variants: ', cvc()))
  
  # Add an ideogram
  ideogram <- reactive({ 
    if (input$chromosome_browser == "All") {
      p.ideo <- Ideogram(genome = "hg19", aspect.ratio=1/10, size=2)
      p.ideo
    } else {
      
      chr <- variants[variants$CHROM == input$chromosome_browser, 'CHROM']
      p.ideo <- Ideogram(genome = "hg19", aspect.ratio=1/10, size=2)
      #start <- input$coordinates[1]
      #end <- input$coordinates[2]
      #gr <- GRanges(chr, IRanges(start = start, end = end))
      #print(gr)
      p.ideo
    }
    
    
  });
  
  output$v_ideogram <- renderPlot(    
    
    ideogram()
  );
  
  info_field <- reactive ({   
    
    validate(
      need(is.numeric(input$variants_table_rows_selected), "Select a variant to view annotations.")
    )
    
    info_string <- variants[input$variants_table_rows_selected, c('INFO')]
    
    info_vec <- c(strsplit(info_string, ";")[[1]][1:length(strsplit(info_string, ";")[[1]])]) 
    info_table <- data.frame("Key" = gsub("(.*)(=)(.*)", "\\1", info_vec), "Value" = gsub("(.*)(=)(.*)", "\\3", info_vec))
    info_table      
  })
  
  
  output$v_info <- renderTable({
    
    info_field()
    
  },
  include.rownames=FALSE)
  
  ############################
  ###### SNP DENSITIES #######
  ############################
  
  snps <- reactive(
    subset(variants, CHROM == input$bins_snp & nchar(as.character(REF)) == 1 & nchar(as.character(ALT)) == 1 & as.numeric(POS) > input$position_snp[1] & as.numeric(POS) < input$position_snp[2])
  )
  
  snps %>% ggvis(~POS, fill:="red") %>%
    layer_densities() %>%
    add_axis("x", title_offset = 60, title = "Chromosome Position") %>%
    add_axis("y", title_offset = 110, title = "Density") %>%
    set_options(width = "auto", resizable=FALSE) %>%
    bind_shiny("ggvis_output_density")
  
  
  #################################
  ##### SNP TYPE DISTRIBUTION #####
  #################################
  
  mut_counts <- reactive(
    subset(variants, CHROM == input$heatBins & nchar(as.character(REF)) < 2 & nchar(as.character(ALT)) < 2 & as.numeric(POS) > as.numeric(input$heatPosition[1]) & as.numeric(POS) < as.numeric(input$heatPosition[2]))
    %>%
      group_by(ALT,REF) %>% summarize(count = n(), .groups = 'drop')     
  )
  
  mut_counts%>%
    ggvis(~ALT, ~REF, fill=~count)%>%
    layer_rects(width = band(), height = band()) %>%
    
    layer_text(
      x = prop("x", ~ALT, scale = "xcenter"),
      y = prop("y", ~REF, scale = "ycenter"),
      text:=~count, fontSize := 20, fill:="white", baseline:="middle", align:="center") %>%
    
    scale_nominal("x", padding = 0, points = FALSE) %>%
    scale_nominal("y", padding = 0, points = FALSE) %>%
    add_legend("fill", title="Mutations rates (SNPs)") %>%
    scale_nominal("x", name = "xcenter", padding = 1, points = TRUE) %>%
    scale_nominal("y", name = "ycenter", padding = 1, points = TRUE) %>%
    set_options(width = "auto", resizable=FALSE) %>%
    bind_shiny("ggvis_output_heatmap")
  
  ####################################
  ##### VARIANT QUAL DISTRIBUION #####
  ####################################
  
  variants1 <- variants
  
  variants1$CHROM <- gsub("chr","",variants1$CHROM)
  variants1$CHROM <- gsub("X","23",variants1$CHROM)
  variants1$CHROM <- gsub("Y","24",variants1$CHROM)
  variants1$CHROM <- gsub("MT","25",variants1$CHROM)
  variants1$CHROM <- as.numeric(as.character(variants1$CHROM))
  
  qual <- reactive(
    if (input$QBins == 'All'){
      variants1
    } else{
      subset(variants1, CHROM == input$QBins)
    }
    
  )
  
  output$manhattan <- renderPlot(      
    manhattan(qual(), snp = "ID", bp = "POS", chr = "CHROM", p = "QUAL", 
              suggestiveline = 33, genomewideline = FALSE, chrlabs = c(1:21, "X", "Y"), 
              main = "Manhattan Plot of VCF QUAL scores", col = c("blue4", "orange3"), 
              logp = FALSE, ylab = "QUAL" , cex = 0.3)
  )
  
  
}