pacman::p_load(
  shiny,
  shinyBS,
  bslib,
  shinythemes,
  ggplot2,
  plotrix,
  networkD3,
  DT,
  dplyr,
  tidyverse,
  cowplot,
  bs4Dash,
  shinyWidgets,
  htmltools,
  grid,
  gridExtra,
  InteractiveComplexHeatmap,
  ComplexHeatmap,
  RColorBrewer,
  htmlwidgets,
  shinyjs,
  shinydashboard,
  svglite,
  data.table
)


#Timeout: https://stackoverflow.com/questions/33839543/shiny-server-session-time-out-doesnt-work

timeoutMinutes <- 30
inactivity <- sprintf("function idleTimer() {
  var t = setTimeout(logout, %s);
  window.onmousemove = resetTimer; // catches mouse movements
  window.onmousedown = resetTimer; // catches mouse movements
  window.onclick = resetTimer;     // catches mouse clicks
  window.onscroll = resetTimer;    // catches scrolling
  window.onkeypress = resetTimer;  // catches keyboard actions

  function logout() {
    Shiny.setInputValue('timeOut', '%s minutes')
  }

  function resetTimer() {
    clearTimeout(t);
    t = setTimeout(logout, %s);  // time is in milliseconds
  }
}
idleTimer();", 
timeoutMinutes * 60 * 1000,    
timeoutMinutes,                
timeoutMinutes * 60 * 1000)    

# Functions
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
    scale_x_discrete(name = NULL,
                     labels = labels)+
    ylim(-0.3, round(max(data_expr) + 1.25*max(data_expr))) +
    theme_classic() +
    theme(plot.title = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 0, b = 20, l = 0)),
          axis.title = element_text(size = 15,color = "black"),
          axis.title.y = element_text(margin = margin(r = 20)),
          axis.text.y = element_text(size = 13, color = "black"),
          axis.line = element_line(color = "black"),
          legend.position = "none",
          axis.text.x = ggtext::element_markdown(size = 13, color = "black"))
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

greenhouse_images <- c(rep("www/AB_0.png", 7), rep("www/Leaf.png", 5))
greenhouse_widths <- c(rep("16.5", 7), rep("33.5", 5))
greenhouse_names <- c("SD15", "CT2", "CT8", "CT10", "BB1", "BB2", "BB3", 
                      "LD", "SD1", "SD2", "SD3", "SD10")

labels_gh <- setNames(
  mapply(create_label, greenhouse_names, greenhouse_images, greenhouse_widths),
  greenhouse_names
)


outdoors_images <- c(rep("www/AB_0.png", 8), rep("www/Leaf.png", 3))
outdoors_widths <- c(rep("16.5", 8), rep("33.5", 3))
outdoors_names <- c("SEP", "OCT", "DEC", "JAN", "FEB", "MAR", "APR", 
                    "MAY", "JUN", "JUL", "AUG")

labels_out <- setNames(
  mapply(create_label, outdoors_names, outdoors_images, outdoors_widths),
  outdoors_names
)

data <- readRDS("expression_data/expression_data.rds") 
nodes <- readRDS("network/nodes.rds")
edges <- readRDS("network/edges.rds")
annot_col <- read_rds("annot_col.rds")
legend_colors <- readRDS("legend_colors.rds")
GO_modules <- readRDS("GO_modules.rds")

TF_list <- sort(discard(as.vector(unique(data$Family)), is.na))


custom_theme <- bslib::bs_theme(
  version = 3,          
  bootswatch = "flatly", 
  fg = "black",          
  primary = "#4A8B94",
  info = "#006F8E",
  link = "#006F8E", 
  bg = "#f0f8f7"         
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
    BB1 = treatment_palette[5],
    BB2 = treatment_palette[6],
    BB3 = treatment_palette[7],
    LD = treatment_palette[8],
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

# Parts of the UI have been adapted from the CAST-R app: https://github.com/Nagel-lab/CAST-R/blob/main/app.R

ui <- fluidPage(
  useShinyjs(),
  tags$script(inactivity),    
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
        font-size: 18px; /* Custom text size for panel body */
    color: #333333; /* Custom text color for panel body */
  }

  .panel-title {
    font-size: 22px; /* Custom font size for panel title */
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

  /* Change the text color in the dropdown options to black */
  .selectize-dropdown-content .option {
    color:black !important; /* Ensure default text color for dropdown options */
  }

    .title-container {
  display: flex;
  justify-content: space-between;
  align-items: center;
  flex-wrap: wrap; /* Allows items to wrap on smaller screens */
  padding: 10px;
}

.logo {
  height: 90px; /* Reduce size slightly for better adaptability */
  max-width: 100%; /* Ensure responsiveness */
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

.btn {
      width: 100% !important;  /* Makes button full width on smaller screens */
      max-width: 250px !important; /* Ensures it doesn’t get too large */
      font-size: 16px !important; /* Ensures text is readable */
      padding: 10px; /* Adds spacing for better UX */
      text-align: center;
    }

  
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
                     id = "concepts", open = "What is a co-expression network?",
                     bsCollapsePanel(
                       "What is a co-expression network?",
                       style = "primary",
                       htmlOutput("concepts"),
                       div(style = "text-align: center;", 
                           tags$img(src = "network_example.png", width = "70%", height = "60%", style = "margin-top: 20px;")),
                       div(style = "margin-bottom: 40px;"), 
                       htmlOutput("example_footnote", style = "font-size: 15px;")
                     )
                   )
                 )
                 ),
               fluidRow(
               column(
                 width = 6,
                 bsCollapse(
                   id = "info", open = "How it works",  
                   bsCollapsePanel(
                     "How it works",
                     htmlOutput("summary"),
                     div(style = "text-align: center;", 
                         tags$img(src = "how_it_works.png", width = "680px", height = "290px", style = "margin-top: 10px;")),
                     style = "primary"
                   )
                 )
               ),
               column(
                 width = 6,
                 bsCollapse(
                   id = "tutorial", open = "Tutorial video",  
                   bsCollapsePanel(
                     "Tutorial video",
                     div(class = "video-container",
                         tags$video(src = "Tutorial.mp4", type = "video/mp4", controls = TRUE, width = "650", style = "max-width: 100%;"), 
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
                     htmlOutput("genes_methods"),
                     div(style = "text-align: center;", 
                         tags$img(src = "Out_gh.png", width = "80%", style = "margin-top: 20px;")),
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("methods_footnote", style = "font-size: 15px;")
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
                     div(style = "text-align: center;", 
                         tags$img(src = "instructions_genes.png", width = "1150px", height = "1000px"))
                   )
                 )
               )
             ),
             fluidRow(
               column(
                 width = 3,
                 style = "position: sticky; top: 20px; height: 100vh; overflow-y: auto;",  
                 bsCollapse(
                   id = "gene_selection", open = "Gene selection",
                   bsCollapsePanel(
                     "Gene selection",
                     radioButtons(inputId="List_gen",label="1. Choose or paste a list of genes", selected=character(0),
                                  choices=c("Choose a gene family" = "fam",
                                            "Paste a list of genes" = "paste")),
                     conditionalPanel('input.List_gen === "fam"', selectizeInput(inputId = "TFfamil", label = "2. Pick a gene family",
                                                                                 choices = TF_list, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1))), 
                     conditionalPanel('input.List_gen === "paste"', textareaInput(id="Genes_list_Tab2","2. Paste a list of PlantGenIE gene ids (e.g. Potra2n1c28...)","",rows=20),
                                      div(HTML("<i>One gene per line, no separator</i>"),style = "margin-bottom:15px")),
                     
                     actionButton("show_heatmap", "Generate expression profiles", class = "btn btn-info"),
                     tags$script(HTML("
    $(document).on('click', '#show_heatmap', function() {
      $('#expression_single .panel-collapse').collapse('show');  
      $('#expression_mult .panel-collapse').collapse('show'); 
      $('#neigh .panel-collapse').collapse('hide'); 
      $('#network .panel-collapse').collapse('hide'); 
    });
  ")),
                     div(style = "margin-bottom: 40px;"), 
                     uiOutput("error_message"),
                     div(style = "margin-bottom: 40px;"), 
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("threshold_message"),
                     div(style = "margin-bottom: 40px;"), 
                     sliderInput("weight_thr", "Edge weight threshold", min = 0.1, max = 1, value = 0.3, step = 0.05),
                     actionButton("submit_table", "Get first neighbors", class = "btn btn-info"),
                     tags$script(HTML("
    $(document).on('click', '#submit_table', function() {
      $('#neigh .panel-collapse').collapse('show'); 
    });
  ")),
                     div(style = "margin-bottom: 40px;"), 
                     actionButton("plot_network", "Plot network", class = "btn btn-info"),
                     tags$script(HTML("
    $(document).on('click', '#plot_network', function() {
      $('#network .panel-collapse').collapse('show');  
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
                     InteractiveComplexHeatmapOutput("heatmap_output", height1 = 880, height2 = 880, width1 = 500, width2 = 500),
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("plot_footnote"),
                     div(style = "display: flex; justify-content: flex-start; margin-top: 20px"),
                         downloadButton("download_data_heatmap", "Expression data.txt", class = "btn-info", style = "width: 200px; padding: 5px 10px;"))
                 )
               ),
               column(
                 width = 9,
                 bsCollapse(
                   id = "neigh", open = NULL,
                   bsCollapsePanel(
                     "First neighbors",
                     style = "primary",
                     div(
                       style = "display: flex; flex-direction: column; align-items: flex-start; overflow-x: auto; width: 100%;",
                       htmlOutput("node_table_footnote"),
                       DT::dataTableOutput("node_table"),
                       div(style = "display: flex; justify-content: flex-start; margin-top: 20px;", 
                           downloadButton("download_table", "Nodes.txt", class = "btn-info", style = "width: 100px; margin-right: 10px; padding: 5px 10px;"),
                           downloadButton("download_edges", "Edges.txt", class = "btn-info", style = "width: 100px; padding: 5px 10px;"))
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
                         width = 8,  
                         
                         forceNetworkOutput("network_plot", width = "100%"),
                         htmlOutput("network_footnote")
                       ),
                       column(
                         width = 3, 
                         uiOutput("dynamic_legend")
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
                     div(style = "text-align: center;", 
                         tags$img(src = "instructions_modules.png", width = "1000px", height = "600px"))
                   )
                 )
               )
             ),
             fluidRow(
               column(
                 width = 3,
                 style = "position: sticky; top: 20px; height: 100vh; overflow-y: auto;",  
                 bsCollapse(
                   id = "module_selection", open = "Module selection",
                   bsCollapsePanel(
                     "Module selection",
                     selectizeInput("module", "Choose Module:", choices = NULL, multiple = FALSE,
                                    options = list(placeholder = 'Type in your module', maxOptions = 50)),
                     div(style = "margin-bottom: 40px;"),
                     actionButton("submit_GO", "Get genes", class = "btn btn-info"),
                     tags$script(HTML("
    // JavaScript function to open a specific bsCollapsePanel
    $(document).on('click', '#submit_GO', function() {
      $('#module_genes .panel-collapse').collapse('show');  
    });
  "))))),
               column(
                 width = 9,
                 bsCollapse(
                   id = "expression_module", open = NULL,
                   bsCollapsePanel(
                     "Expression profile",
                     style = "primary",
                     uiOutput("module_heatmap"),
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("heatmap_module_footnote"),
                     div(style = "display: flex; justify-content: flex-start; margin-top: 20px;",
                         downloadButton("download_module_pdf", "Expression profile.pdf", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px; "))
                   ))
               ),
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
                           downloadButton("download_GO", "GO terms.txt", class = "btn-info")) 
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
                           downloadButton("download_genes", "GO genes.txt", class = "btn-info")) 
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
  observeEvent(input$timeOut, { 
    print(paste0("Session (", session$token, ") timed out at: ", Sys.time()))
    showModal(modalDialog(
      title = "Timeout",
      paste("Session timeout due to", input$timeOut, "of inactivity -", Sys.time()),
      footer = NULL
    ))
    session$close()
  })
  
  #Tab Introduction
  output$Welcome <-renderText({
    paste("<p style='text-align:justify'>",
    "Welcome to POPUL-R, a Shiny app designed to explore the ",  em("Populus"), "transcriptome throughout a yearly growth cycle. This app uses a RNA-Seq dataset we have obtained from different poplar tissues in outdoor and indoor conditions.",
    "</p>",
    "<p style='text-align:justify'>",
    "In POPUL-R, you can perform exploratory analyses on your genes of interest, such as visualizing their expression profiles, identifying other genes with similar expression patterns, and exploring the biological processes they may be involved in. All data generated in the app can be exported in multiple formats for further analysis, if necessary.",
    "</p>",
    "<p style='text-align:justify'>",
          'Each tab includes a section called "Instructions", which explains how to use all the displayed features. In this tab, you will also find:', 
    "</p>",
    "<p style='text-align:justify'>",
          '- <b>What is a co-expression network?</b>: a summary of key concepts needed to use POPUL-R.',
    "</p>", 
    "<p style='text-align:justify'>",
          '- <b>How it works</b>: a quick overview of the main tabs in the app.', 
    "</p>", 
    "<p style='text-align:justify'>",
          '- <b>Tutorial video</b>: a step-by-step guide on how to use POPUL-R.', 
    "</p>", 
    "<p style='text-align:justify'>",
          '<i>Note: Your session will timeout after being 30 minutes of inactivity.</i>',
    "</p>", 
    "<p style='text-align:justify'>",
          "If you would like to incorporate your data to POPUL-R, please contact:", "<br>",
          "populr@proton.me")
  })
  
  output$concepts <-  renderText({
    paste("<p style='text-align:justify'>
            A co-expression network is a graphical representation where genes are connected based on the correlation or similarity of their expression patterns (Figure 1).
          </p>
            <p style='text-align:justify'><b>Key concepts:</b></p>
            
            <ul style='text-align:justify; list-style-type: disc; padding-left: 20px;'>
            <li><b>Node</b>: an individual gene.</li>
            <li><b>Edge</b>: a connection between nodes.</li>
            <li><b>Edge weight</b>: quantification of connection strength between two nodes. Higher edge weights indicate stronger co-expression between nodes.</li>
            <li><b>First neighbors</b>: nodes that are directly connected by an edge.</li>
            <li><b>Module</b>: group of co-expressed genes often involved in similar processes. They are identified by clustering genes with similar expression patterns.</li>
            <li><b>Centrality</b>: number of edges connected to a node. It indicates how important a gene is within the network.</li>
            </ul>"
    )
    
  })
  
  output$example_footnote <-  renderText({
    paste("<b>Figure 1. Example of a co-expression network illustrating some key concepts.</b>", "</br>","<p style='text-align:justify'>")
  })
  
  output$summary <- renderText({
    paste('This section summarizes the main tabs of POPUL-R. There is a more detailed description under the section "Instructions" in each tab.')
  })
  
  #Tab Materials and Methods
  output$genes_methods <-  renderText({
    paste("The data was obtained from samples collected outdoors throughout a whole year and indoors throughout a growth cycle.", "</p>","<p style='text-align:justify'>",
          "The samples from outdoors were collected from ca. 1-year-old and 35-year-old local (Umeå, Sweden) aspen trees once a month around midday (Figure 2).", "</p>","<p style='text-align:justify'>",
          "The indoor conditions that simulate the yearly growth cycle of <i>Populus</i> trees have been previously described in <a href='https://pub.epsilon.slu.se/24748/1/andre_d_210629.pdf' target='_blank'>André et al. (2021)</a> (Figure 2). Briefly, the plants were transferred to soil and grown in long day (LD, 18h light/6h dark, 20°C/18°C) for four weeks. The trees were then subjected to short day treatment
          (SD, 14h light/10h dark, 20°C/18°C) for 15 weeks to induce growth cessation and dormancy. The plants were moved to cold treatment (CT, 8h light/16h dark, 6°C/6°C) for 10 weeks to release dormancy. Lastly, plants were moved back to LD conditions to induce bud break. Samples from buds or leaves were taken in each treatment.", "</p>","<p style='text-align:justify'>",
          'RNA was extracted and pre-processed as described in <a href="https://www.sciencedirect.com/science/article/pii/S0960982222007825?via%3Dihub" target="_blank">André et al. (2022)</a>. The expression heatmaps are generated with the packages <a href="https://academic.oup.com/bioinformatics/article/32/18/2847/1743594" target="blank">"ComplexHeatmap"</a> and <a href="https://academic.oup.com/bioinformatics/article/38/5/1460/6448211" target="blank">"InteractiveComplexHeatmap"</a>, and co-expression analysis was performed using the R package
          <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559" target="blank">"WGCNA"</a>. The visualization of the co-expression network is generated with the R package <a href= "https://cran.r-project.org/web/packages/networkD3/index.html"target="blank"> "networkD3"</a>. For more details, please see the Materials and Methods section in Marcon <i> et.al</i> (2025).')
  })
  
  output$methods_footnote <-  renderText({
    paste("<b>Figure 2. Temperature and day length of outdoor and indoor conditions.</b>", "</br>","<p style='text-align:justify'>",
          "Outdoor temperature represents the average daily temperature from January to December 2013 (excluding November). Meteorological data were obtained from the Swedish Meteorological and Hydrological Institute (<a href='https://www.smhi.se' target='blank'>SMHI</a>).", "<br>",  
          "SD: short day. CT: cold treatment. BB: bud break. LD: long day. The number shown is the number of weeks in the treatment when the samples were taken (see Figure 1).")
  })
  

  #Tab Gene profile
  ## Inspired by CAST-R
  values <- reactiveValues(previous_gene = NULL,
                           list_genes = NULL,
                           selected_genes = NULL, 
                           max_weight = NULL,
                           min_weight = NULL,
                           threshold_message = NULL,  
                           error_message = NULL,
                           expression_plot = NULL,
                           thr = NULL , 
                           selected_GO = NULL) 
  
  observeEvent(input$show_heatmap, {
    updateSliderInput(session, "weight_thr", value = 0.3)
  })


  update_gene_list <- function(values, input) {
    req(input$List_gen)
    if (input$List_gen == "fam") {
      ID <- input$TFfamil
      values$list_genes <- data %>%
        filter(Family %in% ID) %>% 
        pull(Gene) %>% 
        unique()
    } else if (input$List_gen == "paste") {
      values$list_genes <- unlist(strsplit(input$Genes_list_Tab2, "\n"))
    }
    
    if (length(values$list_genes) == 0 | nrow(data[Gene %in% values$list_genes]) == 0) {
      
      values$error_message <- c(
        "No data found for these genes."
      )
      values$list_genes <- character(0)  
    } else {

      values$error_message <- NULL  
    }
    return(values$list_genes)
  }
  
  ## Function to process gene weights and give thresholds
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
        melted <- melt(edges_filtered,
                       measure.vars = c("fromNode", "toNode"),
                       variable.name = "type",
                       value.name = "Gene")
        
        melted_filtered <- melted[Gene %in% values$list_genes]
        max_weight_goi <- melted_filtered[, .(max_weight = max(weight, na.rm = TRUE)), by = Gene]
        
        values$max_weight <- max(max_weight_goi$max_weight, na.rm = TRUE)
        values$max_min_weight <- min(max_weight_goi$max_weight, na.rm = TRUE)
        
        values$threshold_message <- paste0("<b>Edge weight summary</b><br> The range of maximum edge weights of your genes of interest is ",
                                           trunc(values$max_min_weight * 10^2) / 10^2, "-",
                                           trunc(values$max_weight * 10^2) / 10^2,
                                           ". If the chosen edge weight threshold is higher than the lowest value of this range, edges from some genes of interest will not be shown. If the chosen edge weight threshold exceeds the highest value of this range, the table will appear empty.")
        values$error_message <- NULL
      }
    }
    
    values$thr <- input$weight_thr
    
  }
  
  observeEvent(
    list(input$TFfamil, input$Genes_list_Tab2,input$List_gen), 
    {
      update_gene_list(values, input)
    }
  )
  
    filtered_edges <- reactiveVal(data.table())
    
    observeEvent(input$show_heatmap, {
      # update_gene_list(values, input)

      filtered_from <- edges[J(values$list_genes), on = "fromNode", nomatch = 0]
      filtered_to <- edges[J(values$list_genes), on = "toNode", nomatch = 0]
      all_edges <- unique(rbind(filtered_from, filtered_to))
      
      filtered_edges(all_edges)
      
      process_genes(values, all_edges)
    })
    
   error_display <- eventReactive(input$show_heatmap, {
      values$error_message
    })
    
  output$error_message <- renderUI({
    # req(input$show_heatmap)
      tags$div(style = "color: red; font-weight: bold;", error_display())
  })
  
  observeEvent(input$submit_table, {
    values$thr <- input$weight_thr
  })
  
  output$threshold_message <- renderText({
    HTML(values$threshold_message)
  })
  
  expression_data <- eventReactive(input$show_heatmap, {
    req(values$list_genes)
    print("Making expr data")
    data[Gene %in% values$list_genes]
  })
  
  
  expression_data_plot <- eventReactive(input$show_heatmap, {
    req(values$list_genes)
    print("Making expr data plot")
    print(length(values$list_genes))
    
    if (length(values$list_genes) == 1) {
      
      df <- expression_data()
    } else {
      df <- expression_data() %>% 
        group_by(`Month/Treatment`, Location, Tissue) %>% 
        summarize(mean_all = mean(mean),
                  se = std.error(mean)) %>% 
        ungroup()
      colnames(df)[4] <- "mean"
    }
    return(df)
  })
  
  data_out <- eventReactive(input$show_heatmap, {
    req(values$list_genes)  
    
    if (nrow(expression_data()) == 0) {
      return(NULL)
    }
    
    expression_data_plot() %>%
      filter(Location == "Outdoor")
  })
  
  data_gh <- eventReactive(input$show_heatmap, {
    req(values$list_genes)  
    if (nrow(expression_data_plot()) == 0) {
      return(NULL)
    }
    
    expression_data_plot() %>%
      filter(Location == "Indoor")
  })
  
  #Expression plots
  output$expression_plots <- renderPlot({
    req(input$show_heatmap)
    
    df_plot <- expression_data_plot()
    
    if (nrow(df_plot) == 0) {
      ggplot() +
        theme_void() +
        labs(title = "")
    } else {
      
      p1 <- plot_expression(data_out(),labels_out, df_plot$mean) 
      
      p2 <- plot_expression(data_gh(), labels_gh, df_plot$mean)
      print("plots created")
      
      combined_plot <- plot_grid(p1, p2, nrow = 1, align = "h", rel_widths = c(0.92, 1))
      values$expression_plot <- combined_plot
    }
    values$expression_plot 
  })
  
  
  output$plot_single_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. BB: bud break. LD: long day. The number shown is the number of weeks in the treatment. The illustrations represent buds and leaves. If only one gene has been selected, data are means +/- SE for n = 6 biological replicates. If several genes have been selected, data are means +/- SE of all selected genes", "</p>"
    )
  })
  
  output$download_plot_svg <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_expression_plot.svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = values$expression_plot, device = "svg", width = 16, height = 5)
    }
  )
  
  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_expression_plot.pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = values$expression_plot, device = "pdf", width = 16, height = 5)
    }
  )
  
  output$download_expression_data <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_combined_expression_data.txt", sep = "")
    },
    content = function(file) {
      write.table(expression_data_plot(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  
  heatmap_data <- eventReactive(input$show_heatmap, {
    
    req(values$list_genes)
    print("Making heatmap data")

    df <- expression_data() %>%
      ungroup() %>% 
      dplyr::select(!c(se, Treatment2, Location, Tissue, Family)) 
    
    mult <- dcast(df, Gene ~ `Month/Treatment`, value.var = "mean")
    
    mult <- mult %>% 
      as.data.frame %>% 
      column_to_rownames("Gene")
    
    return(mult)
  })
  
  output$download_data_heatmap <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_heatmap_expression_data.txt", sep = "")
    },
    content = function(file) {
      write.table(expression_data() %>% 
                    dplyr::select(!(Treatment2)), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })

  observeEvent(input$show_heatmap, {
     
    withProgress(message = 'Generating heatmap...', value = 0, {
      incProgress(0.2)  

      
      if (is.null(heatmap_data()) || nrow(heatmap_data()) == 0) {
        
        # Case 1: No data — return blank heatmap or empty plot
        heatmap <- grid::grid.newpage()  # This makes the output blank
        
      } else if (nrow(heatmap_data()) == 1) {
      
      # if(nrow(heatmap_data()) == 1) {
        heatmap <- ComplexHeatmap::pheatmap(mat = as.matrix(heatmap_data()),
                                            cluster_rows = FALSE,
                                            cluster_cols = FALSE, 
                                            scale = "row",
                                            legend = TRUE,
                                            border_color = NA,
                                            color = colorRampPalette(c("dodgerblue","white","firebrick"))(10),
                                            fontsize = 12,
                                            fontsize_row = 14,
                                            fontsize_col = 8,
                                            show_rownames = TRUE,
                                            show_colnames = FALSE,
                                            annotation_legend = TRUE,
                                            annotation_col = annot_col,
                                            annotation_colors = annot_colors,
                                            name = "Expression")
      } else {
        
        dist.obs <- as.dist(1-cor(t(heatmap_data())))
        dist.obs.tree <- hclust(dist.obs, method = "ward.D")

        heatmap <- ComplexHeatmap::pheatmap(mat = as.matrix(heatmap_data()),
                                            cluster_rows = dist.obs.tree,
                                            cluster_cols = FALSE, 
                                            scale = "row",
                                            legend = TRUE,
                                            border_color = NA,
                                            color = colorRampPalette(c("dodgerblue","white","firebrick"))(10),
                                            fontsize = 30,
                                            fontsize_row = 10,
                                            fontsize_col = 12,
                                            show_rownames = TRUE,
                                            show_colnames = FALSE,
                                            annotation_legend = TRUE,
                                            annotation_col = annot_col,
                                            annotation_colors = annot_colors,
                                            name = "Expression"
        )
      }
      
      ht = draw(heatmap)
      incProgress(0.5)
      makeInteractiveComplexHeatmap(input, output, session, ht_list = ht)
      incProgress(0.3)
    })
  })
  
  output$plot_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. LD: long day. The number shown is the number of weeks in the treatment. Data are means for n = 6 biological replicates.", "</p>"
    )
  })
  
  filtered_nodes_thr <- eventReactive(input$submit_table, {
    
    filtered_edges_thr <- filtered_edges() %>%
      filter(weight >= values$thr)
    
    if (nrow(filtered_edges_thr) > 0) {
      
      valid_genes <- intersect(values$list_genes, nodes$`Gene name`)
      
      goi_edge <- data.frame(fromNode = valid_genes, 
                             toNode = valid_genes, 
                             weight = 1)
      
      filtered_edges_thr <- rbind(goi_edge, filtered_edges_thr)
    }
    
    filtered_nodes <- data.frame(id = unique(c(filtered_edges_thr$fromNode, filtered_edges_thr$toNode))) 
    return(filtered_nodes)
  })
  
  filtered_nodes_ann <- eventReactive(input$submit_table, {
    print(paste("submit_table_multiple triggered2"))
    req(values$list_genes)
    
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
    edges %>%
      filter((fromNode %in% filtered_nodes_ann()$`Gene name` & toNode %in% filtered_nodes_ann()$`Gene name`)) %>% 
      filter(weight >= values$thr)
  })
  
 #First neighbors table
  observeEvent(input$submit_table, {
    if(values$thr > values$max_weight) {
      output$node_table <- DT::renderDT({
        DT::datatable(
          data.frame("The chosen edge weight threshold exceeds the highest value of the edge weights range (see Edge weight summary)."),  # Message as cell content
          options = list(
            dom = 't',                  
            searching = FALSE,          
            paging = FALSE,             
            info = FALSE,               
            columnDefs = list(
              list(                     
                targets = 0,
                title = "",            
                className = "dt-center" 
              )
            )
          ),
          rownames = FALSE,
          colnames = ""                 
        )
      })

    } else {
    annotated_nodes <- filtered_nodes_ann()
    annotated_nodes[[1]] <- paste0('<a href="https://plantgenie.org/gene?id=', annotated_nodes[[1]], '" target="_blank">', annotated_nodes[[1]], '</a>')
    
    output$node_table <- DT::renderDT({
      req(filtered_nodes_ann())
      print("Rendering table with filtered nodes")
      DT::datatable(
        annotated_nodes,
        options = list(
          pageLength = 20,  
          lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
          searching = TRUE,  
          ordering = TRUE,  
          dom = 'lfrtip'),  
        escape = FALSE,
        selection = 'multiple',
        rownames = FALSE
      )
    })
    }
  })
  
  output$node_table_footnote <-renderText({
    paste("<p style='text-align:justify'>","First neighbors of the input genes according to the chosen edge weight threshold. The table is ordered by centrality of the gene in the co-expression network in descending order. The module number corresponds to the module where the gene is found within the co-expression network. <br>", "Note: some genes will have expression data but they will not be found in the co-expression network because they do not have any strong correlations with any gene in the dataset.", "</p>"
    )
  })
  
  output$download_table <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_nodes","_thr", values$thr, ".txt", sep = "")
    },
    content = function(file) {
      write.table(filtered_nodes_ann(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  
  output$download_edges <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_edges", "_thr", values$thr, ".txt", sep = "")
    },
    content = function(file) {
      write.table(subnetwork(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })
  
  # Update selected genes from the first neigbors table when plot button is clicked
  observeEvent(input$plot_network, {
    req(input$submit_table, values$thr, input$List_gen)
    
    selected_rows <- input$node_table_rows_selected 
    if (!is.null(selected_rows)) {
      gene_names <- filtered_nodes_ann()[selected_rows, "Gene name"]
      
      values$selected_genes <- sapply(gene_names, function(x) gsub('.*>(.*)<.*', '\\1', x))
    } else {
      values$selected_genes <- NULL
    }
  })
  
  #Sub co-expression network
  network_data <- eventReactive(input$plot_network, {
    req(values$list_genes, values$selected_genes, subnetwork())
    
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
    d3.selectAll(".node text")
      .style("fill", "#000000");
  }
  '
      )
      
      return(network_plot)
    } else {
      return(NULL)
    }
  })
  
  # Network Legend
  observeEvent(input$plot_network, {
    req(values$list_genes, values$selected_genes, network_data())

    
    df <- network_data()$nodes %>%
      select(Module) %>%
      left_join(legend_colors, by = c("Module" = "Number")) %>%
      distinct() %>% 
      mutate(
        Module_number = as.numeric(Module)) %>% 
      arrange(Module_number) %>% 
      mutate(y = -seq_along(Module_number))
    
    legend_plot <- 
      ggplot(df, aes(x = 1, y = y)) +
      geom_point(aes(fill = Module), shape = 21, size = 6, color = "grey") +
      geom_text(aes(label = Module), hjust = 0, nudge_x = 0.1, size = 6) +
      scale_fill_manual(values = setNames(df$Color, df$Module)) +
      coord_cartesian(xlim = c(0.9, 1.5)) +  
      labs(title = "Module") +
      theme_void() +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.2),
        title = element_text(face = "bold", size = 14)
      )
    
    
    values$legend_data <- legend_plot
    
    
  })
  
  output$network_legend <- renderPlot({
    req(values$legend_data)
    if (!is.null(values$legend_data)) {
      values$legend_data
    }
  }, height = function() {
    if (is.null(values$legend_data)) return(1)
    else return(400) 
  })
  
  output$dynamic_legend <- renderUI({
    tags$div(
      style = "height: 600px; overflow-y: scroll; border: 1px solid #ddd; padding: 5px;",
      plotOutput("network_legend", height = paste0(nrow(unique(network_data()$nodes["Module"])) * 50, "px"))
    )
  })
  
  
  output$network_footnote <-renderText({
    paste("<p style='text-align:justify'>","Subnetwork of selected genes. The color and size of the nodes correspond to the module number and the centrality of the gene, respectively. The gene id will appear when you hover over the nodes.", "</p>"
    )
  })
  
  
  # Tab module profile
  updateSelectizeInput(session, "module", choices = c("", sort(as.numeric(unique(GO_modules$Module)))), server = TRUE)
  
  observeEvent(input$module, {
    if (!is.null(input$module) && input$module != "") {
      runjs("
      $('#module_genes .panel-collapse').collapse('hide'); 
      $('#expression_module .panel-collapse').collapse('show');
      $('#GO .panel-collapse').collapse('show');
    ")
    }
  })
  
  ##Module heatmaps
  heatmap_path <- reactive({
    file.path("heatmaps", paste0("heatmap_module_", input$module, ".pdf"))
  })
  

  output$module_heatmap <- renderUI({
    req(input$module)  
    
    print("module selected")

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
  
  output$heatmap_module_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. LD: long day. The number shown is the number of weeks in the treatment. Data are means for n = 6 biological replicates.", "</p>"
    )
  })
  
  output$download_module_pdf <- downloadHandler(
    filename = function() {
      paste0(Sys.Date(), basename(heatmap_path()), ".pdf")
    },
    content = function(file) {
      file.copy(file.path("www", heatmap_path()), file)
    }
  )
  
  ## GO terms
  GO_module <- eventReactive(input$module,{
    print("GO_module table started")
    df <- GO_modules[Module %in% input$module][order(-`P-value`)]
    if(nrow(df) > 1) {
      df$`P-value` <- sprintf("%.2e", df$`P-value`)
    }
    print("GO module table created")
    return(df)
  })
  
  observeEvent(input$module, {
    
    output$GO_table <- DT::renderDT({
      req(input$module, GO_module())
      print("Rendering table GO")
      DT::datatable(
        GO_module(),
        options = list(
          pageLength = 20,  
          lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
          searching = TRUE,  
          ordering = TRUE,  
          dom = 'lfrtip'),  
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
    
    ##Update selected GO terms
    selected_rows_go <- input$GO_table_rows_selected 
    if (!is.null(selected_rows_go)) {
  
      values$selected_GO <- GO_module()[selected_rows_go, "GO id"][[1]]
    } else {
      values$selected_GO <- NULL
    }
    
    print(values$selected_GO)
  })
  
  ## Genes associated with selected GO terms
  GO_genes <- eventReactive(input$submit_GO, {
    req(values$selected_GO, GO_module())
    print("rendering GO genes")
    if (length(values$selected_GO) == 0) {
      filtered_GO <- NULL
      
    } else {
      nodes[Module %in% input$module &
                   str_detect(`GO id`, str_c(values$selected_GO, collapse = "|"))
      ][, !"GO id", with = FALSE
      ][order(-Centrality)]
    }
  }) 
  
  observeEvent(input$submit_GO, {
    
    output$module_genes_table <- DT::renderDT({
      req(GO_genes())
      print("Rendering table genes GO")
      
      GO_link <- GO_genes()
      print(GO_link)
      GO_link[[1]] <- paste0('<a href="https://plantgenie.org/gene?id=', GO_link[[1]], '" target="_blank">', GO_link[[1]], '</a>')
      
      if (input$show_GO_term) {
        
        DT::datatable(
          GO_link ,
          options = list(
            pageLength = 20,  
            lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
            searching = TRUE,  
            ordering = TRUE,  
            dom = 'lfrtip'), 
          escape = FALSE,
          selection = 'multiple',
          rownames = FALSE
        )
      } else {
        # 
        DT::datatable(
          GO_link %>%  select(!`GO term`),
          options = list(
            pageLength = 20,  
            lengthMenu = list(c(20, 30, 40, -1), c('20', '30', '40', 'All')),  
            searching = TRUE,  
            ordering = TRUE,  
            dom = 'lfrtip'), 
          escape = FALSE,
          selection = 'multiple',
          rownames = FALSE
        )
      } 
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
  
  # Tab References and citations
  output$tools <- renderText({
    paste("The source code is available on <a href = 'https://github.com/lauragarciaromanach/POPUL_R' target='_blank'>  GitHub</a>." ,"</p>",
          "This app was inspired by CAST-R:", "<br>",
          "<b>Bonnot T, Gillard MB, Nagel DH </b> (2022). CAST-R: A shiny application to visualize circadian and heat stress-responsive genes in plants. <i>Plant Physiol</i> 190(2): 994-1004. <a href = 'https://academic.oup.com/plphys/article/190/2/994/6549534?login=false' target='_blank' > doi: 10.1093/plphys/kiac121 </a>", "</p>",
          "<b>Data </b>", "<br>",
          "<b>Marcon A, García Romañach L, André D, Delhomme N, Hvidsten T, Nilsson O </b> (2025). A transcriptional roadmap of the yearly growth cycle in Populus trees.", "</p>",
          "<b>Tools </b>", "<br>",
          "<b>R Core Team </b>(2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. <a href = 'https://www.r-project.org'target='_blank' >r-project.org</a> ","<br>",
          "<b>Chang W, Cheng J, Allaire J, Xie Y, McPherson J</b> (2020). shiny: Web Application Framework for R.","<br>",
          "<b>Wickham H</b> (2016) ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4. <a href = 'https://ggplot2.tidyverse.org' target='_blank'> https://ggplot2.tidyverse.org</a>" , "<br>",
          "<b>Gu, Z.</b> (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. Bioinformatics. <a href = 'https://academic.oup.com/bioinformatics/article/32/18/2847/1743594' target='_blank'> doi:10.1093/bioinformatics/btw313 </a>","<br>",
          "<b>Gu, Z. </b> (2022). Complex Heatmap Visualization. iMeta. <a href = 'https://onlinelibrary.wiley.com/doi/10.1002/imt2.43' target='_blank'> doi:10.1002/imt2.43 </a>", "<br>",
          "<b>Xie Y., Cheng J., Tan X. </b> (2025). DT: A Wrapper of the JavaScript Library 'DataTables'. <a href = 'https://cran.r-project.org/web/packages/DT/index.html' target = '_blank'> https://CRAN.R-project.org/package=DT </a>", "<br>",
          "<b>Langfelder, P & Horvath, S</b> (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9, 559. <a href='https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559' target='_blank'> doi: 10.1186/1471-2105-9-559 </a>", "<br>",
          "<b>Allaire J, Gandrud C, Russell K, Yetman C</b> (2017). networkD3: D3 JavaScript Network Graphs from R. <a href = 'https://cran.r-project.org/web/packages/networkD3/index.html' target='_blank'> https://CRAN.R-project.org/package=networkD3 </a>"
    )
    
    
  })
  
  output$citation <- renderText({
    paste("<b>Marcon A, García Romañach L, André D, Ding J, Zhang B, Hvidsten T R, Nilsson O </b> (2025). A transcriptional roadmap of the yearly growth cycle in <i>Populus</i> trees.", "</p>",
          "I would like to thank the <a href = 'https://github.com/UPSCb' target='_blank' >Umeå Plant Science Bioinformatics Facility </a> and <a href = 'https://github.com/MaxiEstravis' target='_blank' >Maximiliano Estravis-Barcala </a> for all their support in developing POPUL-R.")
  })
  
  #Restart values and collect garbage at the end of the session
  session$onSessionEnded(function() {
    
    values$list_genes <- NULL
    values$selected_genes <- NULL
    values$max_weight <- NULL
    values$min_weight <- NULL
    values$threshold_message <- NULL
    values$error_message <- NULL   
    values$expression_plot <- NULL
    values$thr <- NULL 
    values$selected_GO <- NULL
    
    gc()
    message("Session ended and garbage collected.")
  })
  
}

shinyApp(ui = ui, server = server)