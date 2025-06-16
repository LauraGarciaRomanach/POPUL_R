library(shiny)
library(shinyBS)
library(bslib)
library(shinythemes)
library(ggplot2)
library(plotrix)
library(networkD3)
library(DT)
library(dplyr)
library(tidyverse)
library(cowplot)
library(bs4Dash)
library(shinyWidgets)
library(htmltools)
library(grid)
library(gridExtra)
# library(extrafont)
library(showtext)
library(sysfonts)
library(InteractiveComplexHeatmap)#Needs citation!
library(ComplexHeatmap) #Needs citation!
# # library(circlize) #Needs citation!
library(RColorBrewer)
library(htmlwidgets)
library(shinyjs)
library(shinydashboard)
library(svglite)

#paquete pacman: pacman load

#Functions
create_label <- function(name, image, width) {
  paste0("<img src = '", image, "' width = '", width, "' /><br>", name)
}

plot_expression <- function(data, labels, data_expr) {
  ggplot(data, aes(x = `Month/Treatment`, y = mean, ymin = mean - se,
                   ymax = mean + se)) +
    geom_line(aes(group = 1), color = "#009E73",
              linewidth = 1.7) +
    geom_ribbon(aes(group = 1, y = mean, ymin = mean - se, ymax = mean + se), fill = "#009E73", alpha = 0.2) +
    labs(x = NULL,
         y = "Expression") +
         # title = values$gene) +
    scale_x_discrete(name = NULL,
                     labels = labels)+
    ylim(-0.3, round(max(data_expr) + 3)) +
    theme_classic() +
    theme(plot.title = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 0, b = 20, l = 0)),
          axis.title = element_text(size = 15,color = "black"),
          axis.title.y = element_text(margin = margin(r = 20)),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = ggtext::element_markdown(size = 13, color = "black"))
}

greenhouse_images <- c(rep("www/AB_0.png", 7), rep("www/Leaf.png", 5))
greenhouse_widths <- c(rep("12.5", 7), rep("33.5", 5))
greenhouse_names <- c("SD15", "CT2", "CT8", "CT10", "LD1", "LD2", "LD3", 
                      "LD4", "SD1", "SD2", "SD3", "SD10")

labels_gh <- setNames(
  mapply(create_label, greenhouse_names, greenhouse_images, greenhouse_widths),
  greenhouse_names
)


outdoors_images <- c(rep("www/AB_0.png", 8), rep("www/Leaf.png", 3))
outdoors_widths <- c(rep("12.5", 8), rep("33.5", 3))
outdoors_names <- c("SEP", "OCT", "DEC", "JAN", "FEB", "MAR", "APR", 
                    "MAY", "JUN", "JUL", "AUG")

labels_out <- setNames(
  mapply(create_label, outdoors_names, outdoors_images, outdoors_widths),
  outdoors_names
)

data <- readRDS("expression_data/expression_data_mini.rds") 
nodes <- readRDS("network/nodes_mini.rds")
edges <- readRDS("network/edges_mini.rds")
annot_col <- read_rds("annot_col.rds")
legend_colors <- readRDS("legend_colors.rds")
GO_modules_mini <- readRDS("GO_modules_mini.rds")
GOIs <- read_delim("GOI_mini_app.txt", col_names = FALSE)

TF <- nodes %>% 
  filter(`Gene name` %in% data$Gene) #Make sure that only the families for our GOIs appear

TF_list <- sort(discard(as.vector(unique(TF$Family)), is.na))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

custom_theme <- bslib::bs_theme(
  version = 3,           # Using Bootstrap 3
  bootswatch = "flatly", # Using the Flatly theme from Bootswatch
  fg = "black",          # Customize foreground color
  primary = "#4A8B94",
  info = "#006F8E",
  link = "#006F8E", # Customize primary color
  bg = "#f0f8f7"         # Light background color
)

#https://stackoverflow.com/questions/14452465/how-to-create-textarea-as-input-in-a-shiny-webapp-in-r
textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

vals <- unique(data$Location)
location_map <- gg_color_hue(length(vals))
names(location_map) <- vals

month_palette <- rev(colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(11))
treatment_palette <- rev(colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(12))

annot_colors = list(
  Location = location_map,
  Tissue   = c(Bud = "chocolate", Leaf = "darkgreen"),
  `Month/Treatment` = c(
    SEP = month_palette[1],
    OCT = month_palette[2],
    DEC = month_palette[3],
    JAN = month_palette[4],
    FEB = month_palette[5],
    MAR = month_palette[6],
    APR = month_palette[7],
    MAY = month_palette[8],
    JUN = month_palette[9],
    JUL = month_palette[10],
    AUG = month_palette[11],
    SD15 = treatment_palette[1],
    CT2 = treatment_palette[2],
    CT8 = treatment_palette[3],
    CT10 = treatment_palette[4],
    LD1 = treatment_palette[5],
    LD2 = treatment_palette[6],
    LD3 = treatment_palette[7],
    LD4 = treatment_palette[8],
    SD1 = treatment_palette[9],
    SD2 = treatment_palette[10],
    SD3 = treatment_palette[11],
    SD10 = treatment_palette[12]
  )
)


color_scale_js <- paste0("d3.scaleOrdinal().domain([", 
                         paste0("'", legend_colors$Number, "'", collapse = ", "), 
                         "]).range([", 
                         paste0("'", legend_colors$Color, "'", collapse = ", "), 
                         "])")

# ".selectize-dropdown {position: static}",
ui <- fluidPage(
  useShinyjs(),
  # theme = shinytheme("flatly"),
  theme = custom_theme, 
  tags$style(
   
    HTML("
/* Increase the font size of the tab titles */
    .nav-tabs > li > a {
      font-size: 18px; /* Adjust this value to change the size */
    }
    
  /* Custom styling for Bootstrap 3 panels */
  .panel-heading {
    background-color: #C1D6E0 !important; /* Custom background color for header */
  }
  .panel {
        border: 1px solid #C1D6E0; /* Custom border color for the panel */
        background-color: #FFFFFF; /* Custom background color for the panel body */
        font-size: 16px; /* Custom text size for panel body */
    color: #333333; /* Custom text color for panel body */
  }

  .panel-title {
    font-size: 20px; /* Custom font size for panel title */
    color: #000000;
    font-weight: bold;
  }
  .title-text {
    color: #2C4B54; /* Custom color for the title text */
    font-size: 30px; /* Custom font size for the title */
    font-weight: bold; /* Bold font for emphasis */
    margin-bottom: 30px;
    margin-top: 25px;
  }

 #  /* Custom styles for selectizeInput placeholder */
 #    .selectize-input {
 #      font-size: 16px !important; /* Increase the font size of the input field */
 #    }
 #    .selectize-dropdown { font-size: 16px; line-height: 16px; }
 #
 
 # /* Change the text color of the selected option in the input box */
 #  .selectize-input.single.has-items {
 #    color: black !important; /* Set the text color to black when an item is selected */
 #  }

  /* Change the text color in the dropdown options to black */
  .selectize-dropdown-content .option {
    color:black !important; /* Ensure default text color for dropdown options */
  }

  # /* Change the text color on hover */
  # .selectize-dropdown-content .option:hover {
  #   background-color: #D3D3D3
  #   color: black !important; /* Text color on hover */
  # }
  # 
#   .selectize-input.single.has-items .selectize-dropdown-content .option:hover {
#   background-color: #D3D3D3 !important; /* Keep the hover color grey */
#   color: black !important; /* Ensure text color stays black */
# }
 
  # /* Change the text color when the input is focused */
  # .selectize-input:focus {
  #   color: black !important; /* Focused text color */
  # }
     .logo {
      position: absolute;
      top: 10px;
      right: 10px;
      height: 110px; 
    }
  
  /* Ensure no column overlap */
  .col-sm-9 {
    float: right;
  }
  .col-sm-3 {
    float: left;
  }
  
  .dataTables_paginate .pagination .paginate_button a {
  color: black !important; /* Change text color to black */
}

# /* Style active pagination link */
# .dataTables_paginate .pagination .paginate_button.active a {
#   color: black !important; /* Keep the active page text color black */
# }
# 
# /* Style disabled pagination link */
# .dataTables_paginate .pagination .paginate_button.disabled a {
#   color: grey !important; /* Set color for disabled page links (optional) */
# }
  
")
    
  ),
  div(
    class = "title-container",
    h1(class= "title-text",
       p("POPUL-R: A transcriptional roadmap of the yearly growth cycle in", em("Populus"), "trees")),
    tags$img(src = "POPUL_R_symbol.png", class = "logo")
  ), 
  tabsetPanel(
    tabPanel("Introduction",
             fluidRow(  
               column(
                 width = 6,
                 bsCollapse(
                   id = "welcome", open = "Welcome",  
                   bsCollapsePanel(
                     "Welcome",
                     htmlOutput("Welcome"),
                     style = "primary"
                   )
                 )
               ),
               column(
                 width = 6,
                 bsCollapse(
                   id = "info", open = "How it works",  
                   bsCollapsePanel(
                     "How it works",
                     htmlOutput("Summary"),
                     style = "primary"
                   )
                 )
               )
             ),
             fluidRow(  
               column(
                 width = 6,
                 bsCollapse(
                   id = "tutorial", open = "Tutorial video",  
                   bsCollapsePanel(
                     "Tutorial video",
                     div(class = "video-container",
                         tags$video(src = "www/Tutorial.mp4", type = "video/mp4", controls = TRUE, width = "650", style = "max-width: 100%;"), 
                         style = "margin: 0 auto;"  
                     ),
                     style = "primary"
                   )
                 )
               )
             )
    ),
    tabPanel("Materials and methods",
             fluidRow(
               column(
                 width = 12,
                 bsCollapse(
                   id = "methods", open = "",
                   bsCollapsePanel(
                     "",
                     style = "primary",
                     # style = "height: 500px, width = 1100px",
                     htmlOutput("genes_methods")
                   )
                 )
               ))),
    
    tabPanel("Gene profile",
             fluidRow(
               column(
                 width = 12,
                 bsCollapse(
                   id = "instructions", open = "NULL",
                   bsCollapsePanel(
                     "Instructions",
                     style = "primary",
                     # style = "height: 500px, width = 1100px",
                     htmlOutput("genes_instructions")
                   )
                 )
               )
               #,?
              
             ),
             fluidRow(
               column(
                 width = 3,
                 style = "position: sticky; top: 20px; height: 100vh; overflow-y: auto;",  # Sidebar styling for stickiness
                 bsCollapse(
                   id = "gene_selection", open = "Gene Selection",
                   bsCollapsePanel(
                     "Gene Selection",
                     radioButtons(inputId="List_gen",label="1. Choose or paste a list of genes", selected=character(0),
                                  choices=c("Choose a transcription factor family" = "fam",
                                            "Paste a list of genes" = "paste")),
                     conditionalPanel('input.List_gen === "fam"', selectizeInput(inputId = "TFfamil", label = "2. Pick a TF family",
                                                                                 choices = TF_list, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1))), 
                     conditionalPanel('input.List_gen === "paste"', 
                                      tags$p("Attention: since this is the mini version of POPUL-R, only the genes in the following file can be used", style = "color: red; font-weight: bold; margin-top: 10px;"),
                                      downloadButton("download_gois", "GOIs.txt", style = "background-color: #f0f0f0; color: black; border: 2px solid black; padding: 10px 20px; font-size: 14px; margin-bottom: 10px;"),
                                      textareaInput(id="Genes_list_Tab2","2. Paste a list of genes","",rows=20),
                                      div(HTML("<i>One gene per line, no separator</i>"),style = "margin-bottom:15px")),
                                     
                     actionButton("show_heatmap", "Generate expression profiles", status = "info"),
                     tags$script(HTML("
    // JavaScript function to open a specific bsCollapsePanel
    $(document).on('click', '#show_heatmap', function() {
      $('#expression_single .panel-collapse').collapse('show');  
      $('#expression_mult .panel-collapse').collapse('show');  
    });
  ")),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("error_message"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("error_message_2"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("threshold_message"),
                     div(style = "margin-bottom: 40px;"), 
                     sliderInput("weight_thr", "Edge Weight Threshold", min = 0.3, max = 1, value = 0.3, step = 0.05),
                     actionButton("submit_table", "Get first neighbours", status = "info"),
                     tags$script(HTML("
    // JavaScript function to open a specific bsCollapsePanel
    $(document).on('click', '#submit_table', function() {
      $('#neigh .panel-collapse').collapse('show');  // Opens the panel
    });
  ")),
                     div(style = "margin-bottom: 40px;"), 
                     actionButton("plot_network", "Plot Network", status = "info"),
                     tags$script(HTML("
    // JavaScript function to open a specific bsCollapsePanel
    $(document).on('click', '#plot_network', function() {
      $('#network .panel-collapse').collapse('show');  // Opens the panel
    });
  "))
                     ))),
               column(
                 width = 9,
                 bsCollapse(
                   id = "expression_single", open = NULL,
                   bsCollapsePanel(
                     "Combined expression profile",
                     style = "primary",
                     plotOutput("expression_plots"),
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("plot_single_footnote"),
                     div(style = "display: flex; justify-content: flex-start; margin-top: 20px;",
                         # downloadButton("download_plot_png", "Expression profile.png", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"),
                         downloadButton("download_plot_svg", "Expression profile.svg", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"), 
                         downloadButton("download_plot_pdf", "Expression profile.pdf", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px; "),
                         downloadButton("download_expression_data", "Expression data.txt", class = "btn-info", style = "width: 200px; padding: 5px 10px;")))
                 )
               ),
               
               column(
                 width = 9,
                 bsCollapse(
                   id = "expression_mult", open = NULL,
                   bsCollapsePanel(
                     "Heatmap",
                     style = "primary; height: 600px; overflow-y: auto; padding: 15px;",
                     # InteractiveComplexHeatmapOutput("heatmap_output", height1 = 880, height2 = 880, width1 = 500, width2 = 500),
                     plotOutput("expression_heatmap", height = "600px"),
                     htmlOutput("plot_footnote"),
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("plot_footnote"),
                     div(style = "display: flex; justify-content: flex-start; margin-top: 20px;",
                         # downloadButton("download_plot_png", "Expression profile.png", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"),
                         downloadButton("download_heatmap_svg", "Expression profile.svg", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"), 
                         downloadButton("download_heatmap_pdf", "Expression profile.pdf", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px; "),
                         downloadButton("download_data_heatmap", "Expression data.txt", class = "btn-info", style = "width: 200px; padding: 5px 10px;")))
                 )
               ),
               column(
                 width = 9,
                 bsCollapse(
                   id = "neigh", open = NULL,
                   bsCollapsePanel(
                     "First neighbours",
                     style = "primary",
                     div(
                       style = "display: flex; flex-direction: column; align-items: flex-start; overflow-x: auto; width: 100%;",
                       htmlOutput("node_table_footnote"),
                       DT::dataTableOutput("node_table"),
                       div(style = "display: flex; justify-content: flex-start; margin-top: 20px;", # Add some margin to create space between the table and button
                           downloadButton("download_table", "Nodes", class = "btn-info", style = "width: 100px; margin-right: 10px; padding: 5px 10px;"),
                           downloadButton("download_edges", "Edges", class = "btn-info", style = "width: 100px; padding: 5px 10px;"))
                     ),
                     tags$style(HTML("
                              #download_btn_wrapper {
                              margin-top: 20px;  /* Minimum space */
                              }
                              #my_table table {
                              margin-bottom: 20px;  /* Minimum space below table */
                              }"))
                   ))),
               column(
                 width = 9,
                 bsCollapse(
                   id = "network", open = NULL,
                   bsCollapsePanel(
                     "Co-expression network",
                     style = "primary",
                     fluidRow(
                       column(
                         width = 8,  # Adjust the width of the network plot column
                         forceNetworkOutput("network_plot", width = "100%")  # Adjust height as necessary
                       ),
                       column(
                         width = 3, 
                         plotOutput("network_legend", height = "600px")  # Set height for the legend plot
                       )
                     )
                   )
                 ))
             )),
    tabPanel("Module profile",
             fluidRow(
               column(
                 width = 12,
                 bsCollapse(
                   id = "instructions", open = "NULL",
                   bsCollapsePanel(
                     "Instructions",
                     style = "primary",
                     # style = "height: 500px, width = 1100px",
                     htmlOutput("module_instructions")
                   )
                 )
               ),
             ),
             fluidRow(
               column(
                 width = 3,
                 style = "position: sticky; top: 20px; height: 100vh; overflow-y: auto;",  # Sidebar styling for stickiness
                 bsCollapse(
                   id = "module_selection", open = "Module selection",
                   bsCollapsePanel(
                     "Module selection",
                     selectizeInput("module", "Choose Module:", choices = NULL, multiple = FALSE,
                                    options = list(placeholder = 'Type in your module', maxOptions = 10)),
                     # actionButton("module_heatmap", "Generate info module", status = "info"),
                     div(style = "margin-bottom: 40px;"),
                     # textOutput("error_message"),
                     # div(style = "margin-bottom: 40px;"), 
                     # textOutput("error_message_2"),
                     # div(style = "margin-bottom: 40px;"), 
                     # textOutput("threshold_message"),
                     # div(style = "margin-bottom: 40px;"), 
                     # sliderInput("weight_thr", "Edge Weight Threshold", min = 0.3, max = 1, value = 0.3, step = 0.05),
                     actionButton("submit_GO", "Get genes", status = "info"),
                     tags$script(HTML("
    // JavaScript function to open a specific bsCollapsePanel
    $(document).on('click', '#submit_GO', function() {
      $('#module_genes .panel-collapse').collapse('show');  
    });
  "))))), 
                     # div(style = "margin-bottom: 40px;"), 
                     # actionButton("plot_network", "Plot Network", status = "info")))),
               column(
                 width = 9,
                 bsCollapse(
                   id = "expression_module", open = NULL,
                   bsCollapsePanel(
                     "Expression profile",
                     style = "primary",
                     uiOutput("module_heatmap"),
                     div(style = "margin-bottom: 40px;"), 
                     div(style = "display: flex; justify-content: flex-start; margin-top: 20px;",
                         # downloadButton("download_plot_png", "Expression profile.png", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"),
                         # downloadButton("download_plot_svg", "Expression profile.svg", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"), 
                         downloadButton("download_module_pdf", "Module expression profile.pdf", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px; "))
                         # downloadButton("download_expression_data", "Expression data.csv", class = "btn-info", style = "width: 200px; padding: 5px 10px;")))
                 )
               )),
               column(
                 width = 9,
                 bsCollapse(
                   id = "GO", open = NULL,
                   bsCollapsePanel(
                     "GO terms",
                     style = "primary",
                     div(
                       style = "display: flex; flex-direction: column; align-items: flex-start; overflow-x: auto; width: 100%;",
                       htmlOutput("GO_module_footnote"),
                       DT::dataTableOutput("GO_table"),
                      
                       div(style = "display: flex; justify-content: flex-start; margin-top: 20px;", 
                           downloadButton("download_GO", "GO terms", class = "btn-info")) 
                                          # style = "width: 100px; margin-right: 10px; padding: 5px 10px;"),
                           # downloadButton("download_edges", "Edges", class = "btn-info", style = "width: 100px; padding: 5px 10px;"))
                     ),
                     tags$style(HTML("
                              #download_btn_wrapper {
                              margin-top: 20px;  /* Minimum space */
                              }
                              #my_table table {
                              margin-bottom: 20px;  /* Minimum space below table */
                              }"))
                   ))),
               column(
                 width = 9,
                 bsCollapse(
                   id = "module_genes", open = NULL,
                   bsCollapsePanel(
                     "Module genes",
                     style = "primary",
                     div(
                       style = "display: flex; flex-direction: column; align-items: flex-start; overflow-x: auto; width: 100%;",
                       htmlOutput("genes_module_footnote"),
                       checkboxInput("show_GO_term", "Show GO terms", value = FALSE),
                       DT::dataTableOutput("module_genes_table"),
                       div(style = "display: flex; justify-content: flex-start; margin-top: 20px;", 
                           downloadButton("download_genes", "GO genes", class = "btn-info")) 
                       # style = "width: 100px; margin-right: 10px; padding: 5px 10px;"),
                       # downloadButton("download_edges", "Edges", class = "btn-info", style = "width: 100px; padding: 5px 10px;"))
                     ),
                     tags$style(HTML("
                              #download_btn_wrapper {
                              margin-top: 20px;  /* Minimum space */
                              }
                              #my_table table {
                              margin-bottom: 20px;  /* Minimum space below table */
                              }"))
                   )))
             )),
    tabPanel("References and citations",
             fluidRow(
               column(
                 width = 6,
                 bsCollapse(
                   id = "datasets_tools", open = "Datasets and tools",
                   bsCollapsePanel(
                     "Datasets and tools",
                     htmlOutput("tools")
                   )
                 )),
               column(
                 width = 6,
                 bsCollapse(
                   id = "cite", open = "How to cite POPUL-R",
                   bsCollapsePanel(
                     "How to cite POPUL-R",
                     htmlOutput("citation")
                   )
                 )
               )
             )
    )
  )
)




# Server
server <- function(input, output, session) {
  
  output$Welcome <-renderText({
    paste("<p style='text-align:justify'>", "Welcome to a Shiny R application, POPUL-R, that allows the user to visualize the expression profile of the",  em("Populus"), "transcriptome throughout a yearly growth cycle. This app uses the RNA-Seq dataset we have recently published, obtained from different poplar tissues in various conditions indoors and outdoors. In POPUL-R, the
          user is able to obtain expression profiles of one gene or interactive heatmaps of several genes. In addition, they can also identify their genes of interest and their first neighbours in a co-expression interactive network that encompasses the whole dataset. The users can set a threshold parameter for the edge weight of the selected network, thereby filtering the connections between their genes of interest and their
          first neighbours by co-expression strength. The users can export the expression profiles and the filtered networks for further analysis if necessary.","</p>","<p style='text-align:justify'>",
          'In each tab, there are two sections called "Instructions" and "Material and methods". The first section explains how to use all the displayed features, and the second one provides methodological details.', "</p>","<p style='text-align:justify'>",
          'The next section "How it works" provides a quick description of all the tabs in POPUL-R.', "</p>", "<p style='text-align:justify'>",
          "For any questions or comments, please contact:", "<br>",
          "laura.garcia@slu.se, ove.nilsson@slu.se")
  })

  output$Summary <- renderText({
    paste("<p style='text-align:justify'>", 'This section summarizes the main functionalities of POPUL-R. There is a more detailed description under the section "Instructions" in each tab.', "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Gene profile"</b>: By choosing a transcription factor family from our dataset or providing one or several <i>Populus tremula</i> gene ids from <a href="https://plantgenie.org" target="_blank">PlantGenIE</a>, visualize the expression profile of the genes of interest indoors and outdoors, identify and select the first neighbours from a thresholded co-expression network and visualize the selected network.', "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Module profile"</b>: By choosing a module from our data set,  visualize the expression profile of the selected module, identify the GO terms associated to the selected module and obtain a list of genes found in the specific module.', "</p>","<p style='text-align:justify'>")
  })

  # output$single_genes_instructions <-  renderText({
  #   paste("<p style='text-align:justify'>", "This tab generates plots and tables from the dataset published in Marcon <i> et.al</i> (2024).", "</p>","<p style='text-align:justify'>",
  #         '- <b>Gene Selection and Expression Profile</b>: You need to enter a <i>Populus tremula </i>gene id from <a href="https://plantgenie.org" target="_blank">PlantGenIE</a> in the Gene Selection box such as Potra2n8c17315 (<i>FT1</i>) and click on "Plot Expression". You will obtain the expression profile of your gene of interest that you can download as .svg or .pdf. The raw data can be downloaded...', "</p>","<p style='text-align:justify'>",
  #         '- <b>First neighbours</b>: A message will appear in the Gene Selection panel with the highest edge weight of the gene: this value represents the highest connection between the selected genes and their first neighbours in the co-expression network. You can then choose an edge weight and click on "Get First Neighbours" to filter and obtain the first neighbours of your gene of interest. If the chosen value is higher than the maximum value provided, no connections are found from your gene of interest and the table appears empty. The table of first neighbours is ordered by edge weight in a descending order. If the gene name on the table is clicked,
  #         you will be directed to the <a href="https://plantgenie.org" target="_blank">PlantGenIE</a> entry of that gene. The centrality value displayed in the table is the centrality of the gene in the whole co-expression network. The module number corresponds to the module where the gene is found within the co-expression network. Once the table is displayed, you can download it as two separate .txt files, nodes.txt and edges.txt, that can be used in other network visualization softwares for further analyses.', "</p>","<p style='text-align:justify'>",
  #         '- <b>Co-expression network</b>: If you want to visualize the connections between your gene of interest and some first neighbours, you can select the genes directly by clicking on the table. The button "Plot Network" allows a visualization of your selected genes, displaying the nodes with colors that correspond to the Module number. The node size depends on the centrality of the gene.')
  # })
  # 

  output$genes_instructions <-  renderText({
    paste("<p style='text-align:justify'>", "This tab generates plots and tables from the dataset published in Marcon <i>et.al</i> (2024).", "</p>","<p style='text-align:justify'>",
          '- <b>Gene Selection and Combined expression profile</b>: You need to select a transcription factor family from the dataset or enter several <i>Populus tremula </i>gene ids from <a href "https://plantgenie.org" target="_blank">PlantGenIE</a> in the Gene Selection box, such as Potra2n8c17315 (<i>FT1</i>). The genes should be entered one per line, with no separators. By clicking on "Generate Heatmap", an interactive heatmap of your selected genes will be displayed. You can zoom-in sections of the heatmap by clicking on them, and they will be displayed in the panel "Selected sub-heatmap". The output section just displays
          detailed information about the sub-heatmap. You can download this heatmap and the sub-heatmap as .png, .svg or .pdf. The raw data from the subheatmap can be downloaded by clicking on "Export to table". ', "</p>","<p style='text-align:justify'>",
          '- <b>First neighbours</b>: A message will appear in the Gene Selection panel with an interval of the highest edge weights of the chosen genes: since multiple genes have been selected and each has a highest edge weight value, the lowest and highest maximum edge weight value from all selected genes are displayed. A maximum edge weight value represents the highest connection between the selected genes and their first neighbours in the co-expression network. You can then choose an edge weight and click on on "Get First Neighbours" to filter and obtain the first neighbours of your genes of interest. This is interval of maximum edge weight values is shown so that, if you choose a higher edge value than the lowest value of the interval, you are aware
          that you will start missing some of your selected genes and first neighbours because the chosen edge value is higher than their maximum edge weight value. If the chosen value is higher than the maximum value from the interval, no connections are found from your gene of interest and the table appears empty. 
          Once the table is displayed, you can download it as two separate .txt files, nodes.txt and edges.txt, that can be used in other network visualization softwares for further analyses.', "</p>","<p style='text-align:justify'>",
          '- <b>Co-expression network</b>: If you want to visualize the connections between your genes of interest and some first neighbours, you can select the genes directly by clicking on the table. The button "Plot Network" allows a quick visualization of your selected genes, displaying the nodes with colors that correspond to the Module number. The node size depends on the centrality of the gene.')
  })
  output$genes_methods <-  renderText({
    paste("The data was obtained from samples collected outdoors throughout a whole year and indoors throughout a growth cycle.", "</p>","<p style='text-align:justify'>",
          "The samples from outdoors were collected from ca. 1-year-old and 35-year-old local (Umeå, Sweden) aspen trees once a month around midday (buds from September to May, except November, and leaves from June to August).", "</p>","<p style='text-align:justify'>",
          "Our growth cycle conditions that simulate the yearly growth cycle of <i>Populus</i> trees have been previously described in <a href='https://pub.epsilon.slu.se/24748/1/andre_d_210629.pdf' target='_blank'>André (2021)</a>. Briefly, <i>in vitro </i>cultures of the trees were grown in jars with MS medium (Murashige and Skoog, 1962) until they were potted in soil for experiments. The trees were grown in long day (LD, 18h light/6h dark) for four weeks at ~20°C during the day and 18°C at night to simulate spring and summer. During this time, the plants were fertilized weekly. After four weeks, the trees were subjected to short day treatment
          (SD, 14h light/10h dark) to simulate autumn and induce dormancy. During this period, fertilization was stopped but the temperature remained constant at ~20/18 °C. After 15 weeeks, the trees were subjected to cold treatment (CT, 8h light/16h dark at 4°C) during 10 weeks to release dormancy. To induce bud flush, the trees were returned to the same LD conditions but fertilization only started after bud flush occured. Samples from buds and/or leaves were taken in each condition.", "</p>","<p style='text-align:justify'>",
          'RNA was extracted and pre-processed as described in <a href="https://www.sciencedirect.com/science/article/pii/S0960982222007825?via%3Dihub" target="_blank">André (2022)</a>. The expression heatmaps are generated with the packages <a href="https://academic.oup.com/bioinformatics/article/32/18/2847/1743594" target="blank">"ComplexHeatmap"</a> and <a href="https://academic.oup.com/bioinformatics/article/38/5/1460/6448211" target="blank">"InteractiveComplexHeatmap"</a>. All heatmaps are scaled by rows, and hierarchical clustering is used to cluster genes. Co-expression analysis was performed using the R package
          <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559" target="blank">"WGCNA"</a>. The visualization of the co-expression network is generated with the R package <a href= "https://cran.r-project.org/web/packages/networkD3/index.html"target="blank"> "networkD3"</a>. For more details, please see the Materials and Methods section in Marcon <i> et.al</i> (2024).')
  })
  
  output$module_instructions <- renderText({
    paste("<p style='text-align:justify'>", "This tab generates heatmaps and tables from the dataset published in Marcon <i>et.al</i> (2024).", "</p>","<p style='text-align:justify'>",
          '- <b>Module selection and Expression profile</b>: By selecting a module from the dataset, you will obtain the expression profile of the selected module. You can download the heatmap as a pdf.', "</p>","<p style='text-align:justify'>",
          '- <b>GO terms</b>: When selecting a module form the dataset, you will obtain a table with the GO terms associated to the selected module. The table can be donwloaded as a .txt file.', "</p>","<p style='text-align:justify'>",
          '- <b>Module genes</b>: If you want to see which genes in the module are associated with specific GO terms, you can select the GO terms directly by clicking on the table. The button "Get genes" will generate another table with the genes associated with the selected GO terms. The table can be downloaded as a .txt file')
  })

  output$tools <- renderText({#DT package!! And Github
    paste("The source code is available on GitHub...","</p>",
          "This app was inspired by CAST-R.", "<br>",
          "<b>Bonnot T, Gillard MB, Nagel DH </b> (2022). CAST-R: A shiny application to visualize circadian and heat stress-responsive genes in plants. <i>Plant Physiol</i> 190(2): 994-1004. <a href = 'https://academic.oup.com/plphys/article/190/2/994/6549534?login=false' target='_blank' > doi: 10.1093/plphys/kiac121 </a>", "</p>",
          "<b>Data </b>", "<br>",
          "<b>Marcon A, García Romañach L, André D, Delhomme N, Hvidsten T, Nilsson O </b> (2024). A transcriptional roadmap of the yearly growth cycle in Populus trees.", "</p>","<p style='text-align:justify'>",
          "<b>Tools </b>", "<br>",
          "<b>R Core Team </b>(2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <a href = 'https://www.r-project.org'target='_blank' >r-project.org</a> ","<br>",
          "<b>Chang W, Cheng J, Allaire J, Xie Y, McPherson J</b> (2020). shiny: Web Application Framework for R.","<br>",
          "<b>Wickham H</b> (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4. <a href = 'https://ggplot2.tidyverse.org' target='_blank'> https://ggplot2.tidyverse.org</a>" , "<br>",
          "<b>Gu, Z.</b> (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. <a href = 'https://academic.oup.com/bioinformatics/article/32/18/2847/1743594' target='_blank'> doi:10.1093/bioinformatics/btw313 </a>","<br>",
          "<b>Gu, Z. </b> (2022). Complex Heatmap Visualization. iMeta. <a href = 'https://onlinelibrary.wiley.com/doi/10.1002/imt2.43' target='_blank'> doi:10.1002/imt2.43 </a>", "<br>",
          "<b>Langfelder, P & Horvath, S</b> (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9, 559. <a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559' target='_blank'> doi: 10.1186/1471-2105-9-559 </a>", "<br>",
          "<b>Allaire J, Gandrud C, Russell K, Yetman C</b> (2017). networkD3: D3 JavaScript Network Graphs from R. <a href = 'https://cran.r-project.org/web/packages/networkD3/index.html' target='_blank'> https://CRAN.R-project.org/package=networkD3 </a>"
    )


  })

  output$citation <- renderText({
    paste("<b>Marcon A, García Romañach L, André D, Delhomme N, Hvidsten T, Nilsson O </b> (2024). A transcriptional roadmap of the yearly growth cycle in <i>Populus</i> trees.")
  })
  
  output$download_gois <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_GOIs_POPUL_R_mini.txt", sep = "")
    },
    content = function(file) {
      # write.csv(expression_data_plot(), file, row.names = FALSE)
      write.table(GOIs, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    })


  # updateSelectizeInput(session, "gene", choices = c("", unique(expression_data$Gene), server = TRUE))

  values <- reactiveValues(gene = NULL, 
                           selected_genes = NULL, 
                           max_weight = NULL,
                           min_weight = NULL,
                           threshold_message = NULL, 
                           # error_message = NULL, 
                           thr = NULL , 
                           list_genes = NULL,
                           selected_GO = NULL) 

  #Tab Gene profile
  ##From CAST-R app

  # Function to handle gene list creation and validation
  update_gene_list <- function(values, input) {
    if (input$List_gen == "fam") {
      ID <- input$TFfamil
      # ID <- "PHD"
      ID_list <- nodes %>% 
        filter(Family %in% ID) %>% 
        filter(`Gene name`%in% data$Gene) #Only take as input the genes that we have selected as GOIs (not its first neighbours too), which are the genes that have expression
      values$list_genes <- ID_list$`Gene name`
    } else if (input$List_gen == "paste") {
      values$list_genes <- unlist(strsplit(input$Genes_list_Tab2, "\n"))
    }
    
    valid_genes <- intersect(values$list_genes, nodes$`Gene name`)
    
    if (length(valid_genes) == 0) {
      print("Error found")
      # No valid genes found
      values$error_message <- c(
        "No data found for these genes. Please check the gene name"
      )
      values$list_genes <- NULL  
    } else {
      # Some valid genes found
      values$error_message <- NULL  
      values$list_genes <- valid_genes  
    }
    
    print(values$list_genes)
    return(values$list_genes)
  }
 
  # Function to process gene weights and thresholds
  
  process_genes <- function(values, edges_filtered) {
    values$threshold_message <- NULL
    values$max_weight <- NULL
    values$thr <- NULL
    
    if(is.null(values$error_message)) {
      if (nrow(edges_filtered) == 0) {
        values$max_weight <- NULL  
        values$threshold_message <- "This gene is below the minimum network threshold, but expression data is available for plotting."
      } #Gene correctly typed but it's below the minimum network threshold
    else {
      max_weight_goi <- edges_filtered %>%
        gather(key = "type", value = "Gene", fromNode, toNode) %>%
        filter(Gene %in% values$list_genes) %>%
        group_by(Gene) %>%
        summarise(max_weight = max(weight, na.rm = TRUE)) %>%
        ungroup()
      
      values$max_weight <- max(max_weight_goi$max_weight, na.rm = TRUE)
      values$max_min_weight <- min(max_weight_goi$max_weight, na.rm = TRUE)
      print(values$max_weight)
      print(values$max_min_weight)
      
      values$threshold_message <- paste("The interval of maximum edge weights of your genes of interest is",
                                                 trunc(values$max_min_weight * 10^2) / 10^2, "-",
                                                 trunc(values$max_weight * 10^2) / 10^2,
                                                 ". If the chosen network threshold is higher than the lowest value of the interval, some genes of interest will not be shown in the table. If the chosen network threshold is higher than the highest value of the interval, the table will appear empty.")
      values$error_message <- NULL
    }
    }
    
    values$thr <- input$weight_thr
    
  }
  
  
  filtered_edges <- eventReactive(input$show_heatmap, {
    df <- edges %>% 
      filter(fromNode %in% values$list_genes | toNode %in% values$list_genes)
    print("filtered edges done")
    return(df)
  })

  observeEvent(input$show_heatmap, {
    update_gene_list(values, input)  # Update gene list
    print("gene list updated")
    process_genes(values, filtered_edges()) # Process gene weights and thresholds
    print("gene list processed")
  })

  # Debugging observer
  observe({
    print("Observer triggered")  
    print(isolate(values$list_genes))  
  })

  # Display error message if any
  output$error_message <- renderText({
    values$error_message
  })

  observeEvent(input$submit_table, {
    values$thr <- input$weight_thr
    print(values$thr)
  })

  # Render max weight message
  output$threshold_message <- renderText({
    values$threshold_message
  })
  
  # list_genes <- "Potra2n6c14928"
  expression_data <- eventReactive(input$show_heatmap, {
    req(values$list_genes)
    print("Making expr data")
      df <- data %>% 
        filter(Gene %in% values$list_genes)
    print(head(df))
    return(df)
  })
  

  expression_data_plot <- eventReactive(input$show_heatmap, {
    req(values$list_genes)
    print("Making expr data plot")
    print(length(values$list_genes))
    
    
    if (length(values$list_genes) == 1) {
      
      df <- data %>% 
        filter(Gene %in% values$list_genes)
      print(df)
      return(df)
    } else {
      df <- data %>% 
        filter(Gene %in% values$list_genes) %>% 
        group_by(`Month/Treatment`, Location, Tissue) %>% 
        summarize(mean_all = mean(mean),
                  se = std.error(mean))
      colnames(df)[4] <- "mean"
    }
    print(head(df))
    return(df)
  })
  
  data_out <- eventReactive(input$show_heatmap, {
    req(values$list_genes)  
    
    # If there's an error, return NULL to prevent plotting
    if (nrow(expression_data()) == 0) {
      return(NULL)
    }
    
    df <- expression_data_plot() %>%
      filter(Location == "Outdoor")
    print("Outdoor data generated")
    return(df)
  })
  
  data_gh <- eventReactive(input$show_heatmap, {
    req(values$list_genes)  
    
    # If there's an error, return NULL to prevent plotting
    if (nrow(expression_data_plot()) == 0) {
      return(NULL)
    }
    
    df <- expression_data_plot() %>%
      filter(Location == "Indoor")
    print("Indoor data generated")
    return(df)
  })
  
  # Render Expression Plots
  output$expression_plots <- renderPlot({
    req(input$show_heatmap)
    
    if (nrow(expression_data_plot()) == 0) {
      ggplot() +
        theme_void() +
        labs(title = "")
    } else {
      
      p1 <- plot_expression(data_out(),labels_out, expression_data_plot()$mean) 
      
      p2 <- plot_expression(data_gh(), labels_gh, expression_data_plot()$mean)
      print("plots created")
      
      plot_grid(p1, p2, nrow = 1, align = "h", rel_widths = c(0.92, 1))
    }
  })
  
  
  
  output$plot_single_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. LD: long day. The number shown is the number of weeks in the treatment. The illustrations represent buds and leaves. Data are means +/- SE for n = 6 biological replicates.", "</p>"
    )
  })
  
  output$download_plot_svg <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_expression_plot.svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "svg", width = 16, height = 5)
    }
  )
  
  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_expression_plot.pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "pdf", width = 16, height = 5)
    }
  )
  
  output$download_expression_data <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_combined_expression_data.txt", sep = "")
    },
    content = function(file) {
      # write.csv(expression_data_plot(), file, row.names = FALSE)
      write.table(expression_data_plot(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  


  # Reactive expression to generate the heatmap
  heatmap_data <- eventReactive(input$show_heatmap, {
    # if (!is.null(values$error_message_multiple2)) {
    #   return(NULL)
    # }

    req(values$list_genes)
    print("Making heatmap data")

    mult <- expression_data() %>%
      dplyr::select(!c(se, Treatment2, Location, Tissue)) %>% 
      arrange(`Month/Treatment`) %>%
      pivot_wider(names_from =  `Month/Treatment`, values_from = mean) %>% 
      column_to_rownames("Gene") 

    return(mult)
  })

  heatmap_plot <- eventReactive(input$show_heatmap, {
    req(heatmap_data())
    
    if (nrow(heatmap_data()) == 1) {
      # Single row case: no clustering
      print("Rendering heatmap plot")
      hp <- pheatmap::pheatmap(
        mat = as.matrix(heatmap_data()),
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        scale = "row",
        legend = TRUE,
        border_color = NA,
        color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
        fontsize = 12,
        fontsize_row = 14,
        fontsize_col = 8,
        show_rownames = TRUE,
        show_colnames = FALSE,
        annotation_col = annot_col,
        annotation_colors = annot_colors)
    } else {
      # Multiple rows: perform clustering
      dist.obs <- as.dist(1 - cor(t(heatmap_data())))  # Distance metric
      dist.obs.tree <- hclust(dist.obs, method = "ward.D")
      print("Rendering heatmap plot")
      hp <- pheatmap::pheatmap(
        mat = as.matrix(heatmap_data()),
        cluster_rows = dist.obs.tree,
        cluster_cols = FALSE,
        scale = "row",
        legend = TRUE,
        border_color = NA,
        color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
        fontsize = 12,
        fontsize_row = 14,
        fontsize_col = 8,
        show_rownames = TRUE,
        show_colnames = FALSE,
        annotation_col = annot_col,
        annotation_colors = annot_colors)
    }
  })
  
  #For mini version, make static heatmap
  observeEvent(input$show_heatmap, {
    # n_rows <- nrow(heatmap_data())
    # height_dim <- max(600, n_rows *62.5)
    output$expression_heatmap <- renderPlot({
      # Show progress message
      withProgress(message = 'Generating heatmap...', value = 0, {
        # mult <- heatmap_data()  # Get heatmap data
        incProgress(0.5)  # Update progress
        grid.draw(heatmap_plot())
        # # Check the dimensions of the data
        # # print(dim(mult))
        # # print(head(mult))
        # if (nrow(heatmap_data()) == 1) {
        #   # Single row case: no clustering
        #   pheatmap::pheatmap(
        #     mat = as.matrix(heatmap_data()),
        #     cluster_rows = FALSE,
        #     cluster_cols = FALSE,
        #     scale = "row",
        #     legend = TRUE,
        #     border_color = NA,
        #     color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
        #     fontsize = 12,
        #     fontsize_row = 14,
        #     fontsize_col = 8,
        #     show_rownames = TRUE,
        #     show_colnames = FALSE,
        #     annotation_col = annot_col,
        #     annotation_colors = annot_colors)
        # } else {
        #   # Multiple rows: perform clustering
        #   dist.obs <- as.dist(1 - cor(t(heatmap_data())))  # Distance metric
        #   dist.obs.tree <- hclust(dist.obs, method = "ward.D")
        #   
        #   pheatmap::pheatmap(
        #     mat = as.matrix(heatmap_data()),
        #     cluster_rows = dist.obs.tree,
        #     cluster_cols = FALSE,
        #     scale = "row",
        #     legend = TRUE,
        #     border_color = NA,
        #     color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
        #     fontsize = 12,
        #     fontsize_row = 14,
        #     fontsize_col = 8,
        #     show_rownames = TRUE,
        #     show_colnames = FALSE,
        #     annotation_col = annot_col,
        #     annotation_colors = annot_colors)
        #   
        # 
        # 
        # if (nrow(heatmap_data()) == 1) {
        #   # Single row case: no clustering
        #   pheatmap::pheatmap(
        #     mat = as.matrix(heatmap_data()),
        #     cluster_rows = FALSE,
        #     cluster_cols = FALSE,
        #     scale = "row",
        #     legend = TRUE,
        #     border_color = NA,
        #     color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
        #     fontsize = 12,
        #     fontsize_row = 14,
        #     fontsize_col = 8,
        #     show_rownames = TRUE,
        #     show_colnames = FALSE,
        #     annotation_col = annot_col,
        #     annotation_colors = annot_colors)
        # } else {
        #   # Multiple rows: perform clustering
        #   dist.obs <- as.dist(1 - cor(t(heatmap_data())))  # Distance metric
        #   dist.obs.tree <- hclust(dist.obs, method = "ward.D")
        #   
        #   pheatmap::pheatmap(
        #     mat = as.matrix(heatmap_data()),
        #     cluster_rows = dist.obs.tree,
        #     cluster_cols = FALSE,
        #     scale = "row",
        #     legend = TRUE,
        #     border_color = NA,
        #     color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
        #     fontsize = 12,
        #     fontsize_row = 14,
        #     fontsize_col = 8,
        #     show_rownames = TRUE,
        #     show_colnames = FALSE,
        #     annotation_col = annot_col,
        #     annotation_colors = annot_colors)
        #   
          
        incProgress(0.5)  # Complete progress
      })
    },
    height = 600, width = 1100)
  })
  
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste0("heatmap_", Sys.Date(), ".pdf")  
    },
    content = function(file) {
      # Open a PDF graphics device
      pdf(file, width = 9, height = 10)

      grid.draw(heatmap_plot())

      # Close the PDF device
      dev.off()
    }
  )
  
  output$download_heatmap_svg <- downloadHandler(
    filename = function() {
      paste0("heatmap_", Sys.Date(), ".svg")  
    },
    content = function(file) {
      # Open a PDF graphics device
      svglite(file, width = 9, height = 10)
      
      grid.draw(heatmap_plot())
      
      # Close the PDF device
      dev.off()
    }
  )
  
  output$download_data_heatmap <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_heatmap_expression_data.txt", sep = "")
    },
    content = function(file) {
      # write.csv(heatmap_data(), file)
      write.table(heatmap_data(), file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    })

  # observeEvent(input$show_heatmap, {
  # 
  # 
  # 
  #   withProgress(message = 'Generating heatmap...', value = 0, {
  #     mult <- heatmap_data()  # Get heatmap data
  #     incProgress(0.5)  # Update progress
  # 
  # 
  #     print(dim(mult))
  #     print(head(mult))
  # 
  #     if(nrow(mult) == 1) {
  #       heatmap <- ComplexHeatmap::pheatmap(mat = as.matrix(mult),
  #                                           cluster_rows = FALSE,
  #                                           cluster_cols = FALSE, #dist.var.tree.GLOBAL,
  #                                           scale = "row",
  #                                           legend = TRUE,
  #                                           border_color = NA,
  #                                           color = colorRampPalette(c("dodgerblue","white","firebrick"))(10),
  #                                           fontsize = 30,
  #                                           fontsize_row = 14,
  #                                           fontsize_col = 12,
  #                                           # srtCol = 45,
  #                                           show_rownames = TRUE,
  #                                           show_colnames = FALSE,
  #                                           #labels_col = names,
  #                                           annotation_legend = TRUE,
  #                                           annotation_col = annot_col,
  #                                           #annotation_row = annot_row,
  #                                           annotation_colors = annot_colors,
  #                                           name = "Expression")
  #       ht = draw(heatmap)
  #       incProgress(0.5)
  #       makeInteractiveComplexHeatmap(input, output, session, ht_list = ht)
  #     } else {
  # 
  # 
  #       dist.obs <- as.dist(1-cor(t(mult)))
  #       dist.obs.tree <- hclust(dist.obs, method = "ward.D")
  # 
  #       heatmap <- ComplexHeatmap::pheatmap(mat = as.matrix(mult),
  #                                           cluster_rows = dist.obs.tree,
  #                                           cluster_cols = FALSE, #dist.var.tree.GLOBAL,
  #                                           scale = "row",
  #                                           legend = TRUE,
  #                                           border_color = NA,
  #                                           color = colorRampPalette(c("dodgerblue","white","firebrick"))(10),
  #                                           fontsize = 30,
  #                                           fontsize_row = 10,
  #                                           fontsize_col = 12,
  #                                           # srtCol = 45,
  #                                           show_rownames = TRUE,
  #                                           show_colnames = FALSE,
  #                                           #labels_col = names,
  #                                           annotation_legend = TRUE,
  #                                           annotation_col = annot_col,
  #                                           #annotation_row = annot_row,
  #                                           annotation_colors = annot_colors,
  #                                           name = "Expression")
  #       ht = draw(heatmap)
  #       # Draw the heatmap and make it interactive
  #       incProgress(0.5)  # Update progress
  #       makeInteractiveComplexHeatmap(input, output, session, ht_list = ht)
  #     }
  #   })
  # })
  
  

  output$plot_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. LD: long day. The number shown is the number of weeks in the treatment. Data are means for n = 6 biological replicates.", "</p>"
    )
  })

  # observeEvent(input$show_heatmap, {
  #   print("submit_table_observed")
  #   output$node_table <- DT::renderDT({
  #     # Create an empty dataframe with the same structure as filtered_nodes_ann
  #     empty_df <- data.frame("Gene name" = character(),
  #                            "Module" = character(),
  #                            "Centrality" = numeric(),
  #                            "ATG symbol" = character(),
  #                            "ATG full name" = character(),
  #                            stringsAsFactors = FALSE)
  #     DT::datatable(empty_df)
  #   })
  # })
  # list_genes <- c("Potra2n1c937",
  #   #                 "Potra2n1c935",
  #   #                 "Potra2n10c21797")
  
  filtered_nodes_thr <- eventReactive(input$submit_table, {
    
    filtered_edges_thr <- filtered_edges() %>%
      filter(weight >= values$thr)
    
    if (nrow(filtered_edges_thr) > 0) {
      
      goi_edge <- data.frame(fromNode = values$list_genes, 
                             toNode = values$list_genes, 
                             weight = 1)
      
      filtered_edges_thr <- rbind(goi_edge, filtered_edges_thr)
    }
   
    
    filtered_nodes <- data.frame(id = unique(c(filtered_edges_thr$fromNode, filtered_edges_thr$toNode))) 
    
    print("Rendered filtered_nodes_thr")
    print(head(filtered_nodes))
    return(filtered_nodes)
  })

  filtered_nodes_ann <- eventReactive(input$submit_table, {
    print(paste("submit_table_multiple triggered2"))
    req(values$list_genes)
    
    # ann_mult1 <- filtered_nodes_mult %>%
      #       filter(id %in% values$list_genes) %>%
      #       left_join(subannot, by = c("id" = "Gene name")) %>%
      #       separate(Module, c('Module_color', 'Module'), sep = "/") %>%
      #       arrange(desc(Centrality))
      # 
      #     ann_mult2 <- filtered_nodes_mult %>%
      #       filter(!(id %in% values$list_genes)) %>%
      #       left_join(subannot, by = c("id" = "Gene name")) %>%
      #       separate(Module, c('Module_color', 'Module'), sep = "/") %>%
      #       arrange(desc(Centrality))
      # 
      #     ann_mult <- rbind(ann_mult1, ann_mult2) %>%
      #       dplyr::select(-nodeID, -GOI, -Module_color)
      
    ann1 <- filtered_nodes_thr() %>% 
      filter(id %in% values$list_genes) %>% 
      left_join(nodes, by = c("id" = "Gene name")) %>% 
      select(!c(`GO id`, `GO term`)) %>% 
      arrange(desc(Centrality))
    
    ann2 <- filtered_nodes_thr() %>% 
      filter(!(id %in% values$list_genes)) %>% 
      left_join(nodes, by = c("id" = "Gene name")) %>% 
      select(!c(`GO id`, `GO term`)) %>% 
      arrange(desc(Centrality))
    
    ann <- rbind(ann1, ann2) 
    
    colnames(ann)[1] <- c("Gene name")

    head(ann)
    return(ann)

  })
  
  subnetwork <- eventReactive(input$submit_table, {
    df <- edges %>%
      filter((fromNode %in% filtered_nodes_ann()$`Gene name` & toNode %in% filtered_nodes_ann()$`Gene name`)) 
    print("subnetwork created")
    return(df)
  })

  # Render Node Table with selectable rows
  observeEvent(input$submit_table, {
    annotated_nodes <- filtered_nodes_ann()
    annotated_nodes[[1]] <- paste0('<a href="https://plantgenie.org/gene?id=', annotated_nodes[[1]], '" target="_blank">', annotated_nodes[[1]], '</a>')
    
    output$node_table <- DT::renderDT({
      req(filtered_nodes_ann())
      print("Rendering table with filtered nodes")# Ensure data is available
      DT::datatable(
        annotated_nodes,
        options = list(
          pageLength = 20,  
          lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
          searching = TRUE,  
          ordering = TRUE,  
          dom = 'lfrtip'),  
          # columnDefs = list(list(targets = ncol(filtered_nodes_ann()) - 1, orderable = FALSE))  
          # language = list(emptyTable = "Your genes of interest are below the chosen network threshold")  
        # ),
        escape = FALSE,
        selection = 'multiple',
        rownames = FALSE
      )
    })
  })
  
  output$node_table_footnote <-renderText({
    paste("<p style='text-align:justify'>","First neighbours of the selected genes according to the chosen edge weight threshold. The table is ordered by centrality of the gene in the co-expression network in descending order. There is no edge weight displayed on this table because several selected genes can be connected to the same first neigbours with different edge weight values. If the gene name on the table is clicked,
          you will be directed to the PlantGenIE entry of that gene. The module number corresponds to the module where the gene is found within the co-expression network.", "</p>"
    )
  })


  # Add download handler for node table
  output$download_table <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_nodes","_thr", values$thr, ".txt", sep = "")
    },
    content = function(file) {
      write.table(filtered_nodes_ann(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })

  # Add download for edges
  output$download_edges <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_edges", "_thr", values$thr, ".txt", sep = "")
    },
    content = function(file) {
      write.table(subnetwork(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })

  # Update selected genes when plot button is clicked
  observeEvent(input$plot_network, {
    req(input$submit_table, values$thr, input$List_gen)

    selected_rows <- input$node_table_rows_selected #???
    if (!is.null(selected_rows)) {
      # Extract the gene names and apply gsub to each element individually
      gene_names <- filtered_nodes_ann()[selected_rows, "Gene name"]

      # Use sapply to apply gsub to each gene name in the vector
      values$selected_genes <- sapply(gene_names, function(x) gsub('.*>(.*)<.*', '\\1', x))
    } else {
      values$selected_genes <- NULL
    }
    print(paste("Selected genes: ", paste(values$selected_genes, collapse = ", ")))
  })
  
  network_data <- eventReactive(input$plot_network, {
    req(values$list_genes, values$selected_genes, subnetwork())
    
    print(head(subnetwork()))
    if (length(values$selected_genes) == 0) {
      filtered_edges_network <- subnetwork()
    } else {
      filtered_edges_network <- subnetwork() %>%
        filter(fromNode %in% values$selected_genes & toNode %in% values$selected_genes)
    }
    
    if (nrow(filtered_edges_network) == 0) {
      showModal(modalDialog(
        title = "No edges found",
        "The selected genes are not connected to each other in the network.",
        easyClose = TRUE,
        footer = NULL
      ))
      
      # If there are no edges, return empty data (edges and nodes)
      return(list(
        edges = data.frame(),
        nodes = data.frame(`Gene name` = character(), Centrality = numeric(), Module = character(), nodeID = numeric())
      ))
    }
    
    filtered_nodes_network <- filtered_nodes_ann() %>% 
      filter(`Gene name` %in% values$selected_genes) %>% 
      select(`Gene name`, Centrality, Module) %>% 
      mutate(nodeID = 1:nrow(.) -1)
    
    filtered_edges_network$fromNodeID <- match(filtered_edges_network$fromNode, filtered_nodes_network$`Gene name`) - 1
    filtered_edges_network$toNodeID <- match(filtered_edges_network$toNode, filtered_nodes_network$`Gene name`) - 1
    
    print(head(filtered_edges_network))
    
    list(edges = filtered_edges_network, nodes = filtered_nodes_network)
  })
  
  output$network_plot <- renderForceNetwork({
    req(network_data())  
    
    if (nrow(network_data()$nodes) > 0) {
    
    network_plot <- forceNetwork(
      Links = network_data()$edges,
      Nodes = network_data()$nodes,
      Source = "fromNodeID",
      Target = "toNodeID",
      NodeID = "Gene name",
      Group = "Module",
      opacity = 0.8,
      zoom = TRUE,
      linkColour = "grey",
      # linkWidth =  "weight",
      # linkWidth = 2,
      fontSize = 20,
      linkDistance = 300,
      colourScale = networkD3::JS(color_scale_js),
      Nodesize = "Centrality",
      legend = FALSE
    )
    network_plot <- htmlwidgets::onRender(
      network_plot,
      '
  function(el, x) {
    d3.selectAll(".node circle")
      .style("stroke", "grey")
      .style("stroke-width", function(d) { return (d.weight); })

    d3.selectAll(".node text")
      .style("fill", "#000000");
  }
  '
    )
    
    return(network_plot)
    } else {
      # If no nodes, return NULL, effectively rendering nothing
      return(NULL)
    }
  })
  
  # Render legend plot
  observeEvent(input$plot_network, {
    req(values$list_genes, values$selected_genes, network_data())
    
    if (nrow(network_data()$nodes) > 0) {
    df <-  network_data()$nodes %>%
      select(Module) %>%
      left_join(legend_colors, by = c("Module" = "Number")) %>%
      unique()
    
    df <- df %>%
      mutate(x = 1,
             y = 1:nrow(df))
 
    plot_leg <- ggplot(df, aes(x, y, color = Module)) +
      geom_point() +
      scale_color_manual(values = unique(df$Color), name = "Module") +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme_classic() +
      theme(
        legend.key = element_rect(fill = "#FFFFFF"),
        legend.background = element_rect(fill = "#FFFFFF"),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.key.size = unit(2.5, 'cm'))
    
    values$legend_data <- cowplot::get_legend(plot_leg)
    } else {
      # If there are no nodes, clear the legend data
      values$legend_data <- NULL
    }
    
  })
  
  output$network_legend <- renderPlot({
    req(values$legend_data, input$plot_network)
    grid.newpage()
    pushViewport(viewport(x = 0.6, y = 0.5, just = "center"))  # Move the legend to the right
    
    grid.draw(values$legend_data)
  })
  
  ##Tab module profile
  updateSelectizeInput(session, "module", choices = c("", sort(as.numeric(unique(GO_modules_mini$Module)))), server = TRUE)
  
  # runjs("$('#GO .panel-collapse').collapse('hide');")
  # observeEvent(input$module, {
  #   runjs("$('#expression_module .panel-collapse').collapse('show');
  #         $('#GO .panel-collapse').collapse('show');")
  # })
  
  observeEvent(input$module, {
    if (!is.null(input$module) && input$module != "") {
      runjs("
      $('#expression_module .panel-collapse').collapse('show');
      $('#GO .panel-collapse').collapse('show');
    ")
    }
  })

  heatmap_path <- reactive({
    path <- file.path("heatmaps", paste0("heatmap_module_", input$module, ".pdf"))
    print(path)
    return(path)
  })
  
  
  # Render PDF
  output$module_heatmap <- renderUI({
    req(input$module)  # Ensure a module is selected
    
    print("module selected")
    # Verify the file exists
    if (file.exists(file.path("www", heatmap_path()))) {
      tags$iframe(
        src = heatmap_path(),
        width = "100%",
        height = "700px"
      )
    } else {
      tags$p("Selected heatmap file does not exist.")
    }
  })
  
  output$download_module_pdf <- downloadHandler(
    filename = function() {
      paste0(basename(heatmap_path()), ".pdf")
    },
    content = function(file) {
      file.copy(file.path("www", heatmap_path()), file)
    }
  )
  
  
  GO_module <- eventReactive(input$module,{
    print("GO_module table started")
    df <- GO_modules_mini %>% 
      filter(Module %in% input$module) %>% 
      arrange(`P-value`, desc(x))
    df$`P-value` <- sprintf("%.2e", df$`P-value`)
    print("GO module table created")
    print(df)
    return(df)
  })
  
  observeEvent(input$module, {
    
    output$GO_table <- DT::renderDT({
      req(input$module, GO_module())
      print("Rendering table GO")# Ensure data is available
      DT::datatable(
        GO_module(),
        options = list(
          pageLength = 20,  
          lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
          searching = TRUE,  
          ordering = TRUE,  
          dom = 'lfrtip'),  
          # columnDefs = list(list(targets = ncol(GO_module()) - 1, orderable = FALSE))  
          # language = list(emptyTable = "Your genes of interest are below the chosen network threshold")  
        # ),
        escape = FALSE,
        selection = 'multiple',
        rownames = FALSE
      )
    })
  })
  
  output$GO_module_footnote <-renderText({
    paste("<p style='text-align:justify'>","Gene Ontology enrichment analysis for the selected module.", "<br>",
    "n: number of expressed genes associated to the GO term in all modules.", "<br>",
    "x: number of expressed genes associated to the GO term in the selected module.", "<br>",
    "The table is ordered first by P-value and then by x in descending order.", "</p>"
    )
  })
  
  output$download_GO <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_module_", input$module, "_GO_terms.txt", sep = "")
    },
    content = function(file) {
      write.table(GO_module(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  
  observeEvent(input$submit_GO, {
    req(input$submit_GO, GO_module())
    
    selected_rows_go <- input$GO_table_rows_selected #???
    if (!is.null(selected_rows_go)) {
      # Extract the gene names and apply gsub to each element individually
      # gene_names <- filtered_nodes_ann()[selected_rows, "Gene name"]
      
      # Use sapply to apply gsub to each gene name in the vector
      values$selected_GO <- GO_module()[selected_rows_go, "GO id"]
    } else {
      values$selected_GO <- NULL
    }
    print(paste("Selected GO:", paste(values$selected_GO, collapse = ", ")))
  })
  
  # selected_GO <- c("GO:0015979")
  
  GO_genes <- eventReactive(input$submit_GO, {
    req(values$selected_GO, GO_module())
    print("rendering GO genes")
    if (length(values$selected_GO) == 0) {
        filtered_GO <- NULL
      
    } else {
      filtered_GO <- nodes %>%
        filter(Module %in% input$module) %>% 
          filter(
            str_detect(`GO id`, str_c(values$selected_GO, collapse = "|"))) %>% 
        select(!`GO id`) %>% 
        arrange(desc(Centrality))

    }

    print(head(filtered_GO))
    return(filtered_GO)
}) 
  
  observeEvent(input$submit_GO, {
    
    output$module_genes_table <- DT::renderDT({
      req(GO_genes())
      print("Rendering table genes GO")
      
      # go_column_visibility <- list(targets = ncol(GO_genes()) - 1)  # Target the last column
      
      # Conditionally set the visibility
      if (input$show_GO_term) {
        
        DT::datatable(
          GO_genes() ,
          options = list(
            pageLength = 20,  
            lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
            searching = TRUE,  
            ordering = TRUE,  
            dom = 'lfrtip'), 
            # columnDefs = list(
              # list(targets = ncol(GO_genes()) - 1, orderable = FALSE),
              # list(go_column_visibility)) 
            # language = list(emptyTable = "Your genes of interest are below the chosen network threshold")  
          # ),
          escape = FALSE,
          selection = 'multiple',
          rownames = FALSE
        )
        
        # go_column_visibility <- c(go_column_visibility, list(visible = TRUE))  # Hide the column if checkbox is unchecked
      } else {
        # 
        DT::datatable(
          GO_genes() %>%  select(!`GO term`),
          options = list(
            pageLength = 20,  
            lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
            searching = TRUE,  
            ordering = TRUE,  
            dom = 'lfrtip'), 
          # columnDefs = list(
          # list(targets = ncol(GO_genes()) - 1, orderable = FALSE),
          # list(go_column_visibility)) 
          # language = list(emptyTable = "Your genes of interest are below the chosen network threshold")  
          # ),
          escape = FALSE,
          selection = 'multiple',
          rownames = FALSE
        )
        # go_column_visibility <- c(go_column_visibility, list(visible = FALSE))
      } 
      
      # print("GO term clicked")
      
      # DT::datatable(
      #   GO_genes(),
      #   options = list(
      #     pageLength = 20,  
      #     lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
      #     searching = TRUE,  
      #     ordering = TRUE,  
      #     dom = 'lfrtip',  
      #     columnDefs = list(
      #       # list(targets = ncol(GO_genes()) - 1, orderable = FALSE),
      #       list(go_column_visibility)) 
      #     # language = list(emptyTable = "Your genes of interest are below the chosen network threshold")  
      #   ),
      #   escape = FALSE,
      #   selection = 'multiple',
      #   rownames = FALSE
      # )
    })
  })
  
  output$genes_module_footnote <-renderText({
    paste("<p style='text-align:justify'>","Genes expressed in the module associated with the selected GO terms.The table is ordered by centrality of the gene in the co-expression network in descending order.", "</p>"
    )
  })
  
  output$download_genes <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_module_", input$module, "_GO_genes.txt", sep = "")
    },
    content = function(file) {
      write.table(GO_genes(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  
}



# Run the application
shinyApp(ui = ui, server = server)