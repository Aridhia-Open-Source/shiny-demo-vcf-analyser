##################
####### UI #######
##################

# Sidebar ----------------------------------------------------------------------
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Variant Browser", tabName = 'variant_browser', icon = icon("search")),
    menuItem("SNP Densities", tabName = "snp_densities", icon = icon("line-chart"))
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
    snp_densities
  )
)


# Defining the page
dashboardPage(
  dashboardHeader(disable = TRUE),
  sidebar,
  body
)



