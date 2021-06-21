##################
####### UI #######
##################

# Sidebar ----------------------------------------------------------------------
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Variant Browser", tabName = 'variant_browser', icon = icon("search")),
    menuItem("SNP Densities", tabName = "snp_densities", icon = icon("chart-area")),
    menuItem("SNP Type Distribution", tabName = "snp_type_distribution", icon = icon("chart-bar")),
    menuItem("Variants Quality Distribution", tabName = "variant_qual_distribution", icon = icon("chart-bar"))
  )
)


# Body -------------------------------------------------------------------------

# Define UI for application that plots random distributions 

body <- dashboardBody(
  includeCSS("./www/styles.css"),
  
  # Browser tab title
  titlePanel(
    windowTitle = "VCF Analyser",
    title = tags$head(tags$link(rel = "shortcut icon", href="favicon.ico"))
  ),
  
  # Tabs
  tabItems(
    # Variant Browser tab
    browser_tab,
    snp_densities,
    snp_type_distributions,
    variant_qual_distribution
  )
)


# Defining the page
dashboardPage(
  dashboardHeader(disable = TRUE),
  sidebar,
  body
)



