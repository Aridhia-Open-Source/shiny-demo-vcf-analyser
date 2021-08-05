documentation_tab <- function() {
  tabPanel("Help",
           fluidPage(width = 12,
                     fluidRow(column(
                       6,
                       h3("VCF Analyser"), 
                       p("This mini-app allows you to perform basic exploratory analysis of VCF files. Variant Calling Format (VCF) files are used for storing gene sequence variations, and are composed by a header, that provides metadata, and a body, which is tab separated into 8 mandatory columns (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO) plus other optional columns."),
                       tags$ol(
                         tags$li(strong("Variant Browser"), " shows a table displaying the body of the file. This table can be filtered by sample, chromosome and position. The user can also select one variant by clicking on the row, then information about the variant and the link to Ensemble will fill the right-side column. It also allows to view the header of the VCF file."), 
                         tags$li(strong("SNP Densities"), 
                                 " displays a graph showing the SNP distributions, the user can filter the data by choosing a chromosome and a specific genomic range."),
                         tags$li(strong("SNP Type Distribution "), 
                                 "shows a heatmap with the number of SNPs by base change. The user can filter the data by chromosome and genomic range."),
                         tags$li(strong("Quality Distribuion "), "shows a Manhattan plot of the quality scores for each variant. The user can filter by chromosome.")
                       ),
                       br(),
                       p("The data used in this demo app was extracted from ClinVar FTP.")
                     ),
                     column(
                       6,
                       h3("Walkthrough video"),
                       tags$video(src="vcf.mp4", type = "video/mp4", width="100%", height = "350", frameborder = "0", controls = NA),
                       p(class = "nb", "NB: This mini-app is for provided for demonstration purposes, is unsupported and is utilised at user's 
                       risk. If you plan to use this mini-app to inform your study, please review the code and ensure you are 
                       comfortable with the calculations made before proceeding. ")
                       
                     )
                     )
                     
                     
                     
                     
           ))
}