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

#paquete pacman: pacman load

#Functions
create_label <- function(name, image, width) {
  paste0("<img src = '", image, "' width = '", width, "' /><br>", name)
}

#Load data
load("RData/data_1_3_2_year_around_v2_app.RData")
load("RData/genetable_1_3_2_year_around_v2.RData")

samples.sep <- samples.sep %>% 
  unite("Treatment_Week", Treatment, Week, na.rm = TRUE, sep = "") 

samples.sep$Treatment_Week <- if_else(samples.sep$Location == "Greenhouse" , samples.sep$Treatment_Week, "")

samples.sep <- samples.sep %>% 
  unite("Month/Treatment", Month, Treatment_Week, na.rm = TRUE, sep = "")

treatment_mapping <- c(
  "SEP", "OCT", "DEC", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
  "AUG", "SD15", "CT2", "CT8", "CT10", "LD1", "LD2", "LD3", "LD4", "SD1", 
  "SD2", "SD3", "SD10"
)

samples.sep$Treatment2 <- match(samples.sep$`Month/Treatment`, treatment_mapping)

samples.sep$Location <- factor(samples.sep$Location,
                               levels = c("Outdoor", 
                                          "Greenhouse"),
                               labels = c("Outdoor",
                                          "Indoor"))

samples.sep <- samples.sep %>%
  arrange(desc(Location), Tissue) %>% 
  arrange(Treatment2)

expression_data <- data %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Gene") %>% 
  pivot_longer(!Gene, names_to = "Samples", values_to = "Expression") %>% 
  mutate(Sample = factor(Samples, levels = colnames(data))) %>% 
  left_join(samples.sep %>% dplyr::select(-Expression), by = "Samples") %>% 
  arrange(Treatment2)


TF_list <-sort(discard(as.vector(unique(subannot$Family)), is.na))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

nodes <- data.frame(read_delim("network/nodes-all_1_3_2_year_around_v2_0.1.txt"))
edges <- data.frame(read_delim("network/edges-all_1_3_2_year_around_v2_0.1.txt"))

custom_theme <- bslib::bs_theme(
  version = 3,           # Using Bootstrap 3
  bootswatch = "flatly", # Using the Flatly theme from Bootswatch
  fg = "black",          # Customize foreground color
  primary = "#4A8B94",
  info = "#006F8E",
  link = "#006F8E", # Customize primary color
  bg = "#f0f8f7"         # Light background color
)

font_add_google(name = "Lato", family = "lato")
showtext_auto()

#https://stackoverflow.com/questions/14452465/how-to-create-textarea-as-input-in-a-shiny-webapp-in-r
textareaInput <- function(id, label, value, rows=20, cols=35, class="form-control"){
  tags$div(
    class="form-group shiny-input-container",
    tags$label('for'=id,label),
    tags$textarea(id=id,class=class,rows=rows,cols=cols,value))
}

legend <- {
  df <- subannot %>% 
    dplyr::select(Module) %>% 
    separate(Module, c('Color', 'Module'), sep = "/") %>% 
    unique() %>% 
    mutate(y = 1:47,
           x = 1)
  df$Module <- factor(df$Module, 
                      levels = as.character(0:47))
  suppressWarnings({
    
    plot_leg <- ggplot(df, aes(x, y, color = Module)) +
      geom_point(show.legend = TRUE) +
      scale_color_manual(name = "Modules",
                         values = c("grey", "turquoise", "blue", "brown", "yellow","green","red",  "black", "pink", "magenta", "purple", "greenyellow", "tan", "salmon", "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen",
                                    "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise","darkgrey", "orange", "darkorange", "white","skyblue", "saddlebrown", "paleturquoise", "steelblue", "violet", "darkolivegreen", "darkmagenta", "sienna3", "yellowgreen", "skyblue3",        
                                    "plum1", "orangered4", "mediumpurple3", "lightsteelblue1", "lightcyan1", "ivory", "floralwhite", "darkorange2", "brown4")) +
      guides(color = guide_legend(override.aes = list(size = 4))) +
      theme_classic() +
      theme(
        legend.key = element_rect(fill = "#C1D6E0"),
        legend.background = element_rect(fill = "#C1D6E0"),
        legend.title = element_text(size = 18, family = "lato", face = "bold"),
        legend.text = element_text(size = 14, family = "lato"),
        legend.key.size = unit(0.8, 'cm'))
    
    cowplot::get_legend(plot_leg)
  })
  }
    




ui <- fluidPage(
  # theme = shinytheme("flatly"),
  theme = custom_theme, #Recheck this, title container might be too much
  tags$style(
    ".selectize-dropdown {position: static}",
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

  /* Custom styles for selectizeInput placeholder */
    .selectize-input {
      font-size: 16px !important; /* Increase the font size of the input field */
    }
    .selectize-dropdown { font-size: 16px; line-height: 16px; }
    
 /* Change the text color of the selected option in the input box */
  .selectize-input.single.has-items {
    color: black !important; /* Set the text color to black when an item is selected */
  }

  /* Change the text color in the dropdown options to black */
  .selectize-dropdown-content .option {
    color:#4A8B94 !important; /* Ensure default text color for dropdown options */
  }

  /* Change the text color on hover */
  .selectize-dropdown-content .option:hover {
    color: black !important; /* Text color on hover */
  }

  /* Change the text color when the input is focused */
  .selectize-input:focus {
    color: black !important; /* Focused text color */
  }
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
                   id = "welcome", open = "Welcome",  # Initial open panel
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
                   id = "info", open = NULL,  # Initial open panel
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
                   id = "tutorial", open = "Tutorial video",  # Initial open panel
                   bsCollapsePanel(
                     "Tutorial video",
                     div(class = "video-container",
                         tags$video(src = "Tutorial.mp4", type = "video/mp4", controls = TRUE, width = "650", style = "max-width: 100%;"), 
                         style = "margin: 0 auto;"  # Center the video
                     ),
                     style = "primary"
                   )
                 )
               )
    )
    ),
    tabPanel("Single genes",
             fluidRow(
               column(
                 width = 6,
                 bsCollapse(
                   id = "instructions_single", open = "NULL",
                   bsCollapsePanel(
                     "Instructions",
                     style = "primary",
                     htmlOutput("single_genes_instructions")
                   )
                 )
               ),
               column(
                 width = 6,
                 bsCollapse(
                   id = "methods_single", open = "NULL",
                   bsCollapsePanel(
                     "Materials and methods",
                     style = "primary",
                     htmlOutput("single_genes_methods")
                   )
                 )
               )
             ),
             
             fluidRow(
               column(
                 width = 3,
                 style = "position: sticky; top: 20px; height: 100vh; overflow-y: auto;",  # Sidebar styling for stickiness
                 bsCollapse(
                   id = "gene_selection", open = "Gene Selection",
                   bsCollapsePanel(
                     "Gene Selection",
                     selectizeInput("gene", "Choose Gene:", choices = NULL, multiple = FALSE, 
                                    options = list(placeholder = 'Type in your gene: ie. Potra2n...', maxOptions = 50)),
                     actionButton("submit", "Plot Expression", status = "info"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("error_message"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("threshold_message"),
                     div(style = "margin-bottom: 40px;"), 
                     sliderInput("weight_thr", "Edge Weight Threshold", min = 0.1, max = 1, value = 0.1, step = 0.05),
                     actionButton("submit_table", "Get first neighbours", status = "info"),
                     div(style = "margin-bottom: 20px;"), 
                     actionButton("plot_button", "Plot Network", status = "info")))),
               column(
                 width = 9,
                 bsCollapse(
                   id = "expression_single", open = NULL,
                   bsCollapsePanel(
                     "Expression Profile",
                     style = "primary",
                     plotOutput("expression_plots"),
                     div(style = "margin-bottom: 40px;"), 
                     htmlOutput("plot_single_footnote"),
                     div(style = "display: flex; justify-content: flex-start; margin-top: 20px;",
                         # downloadButton("download_plot_png", "Expression profile.png", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"),
                         downloadButton("download_plot_svg", "Expression profile.svg", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px;"), 
                         downloadButton("download_plot_pdf", "Expression profile.pdf", class = "btn-info", style = "width: 200px; margin-right: 10px; padding: 5px 10px; "),
                         downloadButton("download_expression_data", "Expression data.csv", class = "btn-info", style = "width: 200px; padding: 5px 10px;")))
                 )
               )
               ,
               column(
                 width = 9,
                 bsCollapse(
                   id = "neigh_single", open = NULL,
                   bsCollapsePanel(
                     "First neighbours",
                     style = "primary",
                     # style = "height: 1400px;",  # Set the height for the card containing the table
                     div(
                       style = "display: flex; flex-direction: column; align-items: flex-start; overflow-x: auto; width: 100%;",
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
                   id = "network_single", open = NULL,
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
    tabPanel("Multiple genes",
             fluidRow(
               column(
                 width = 6,
                 bsCollapse(
                   id = "instructions_mult", open = "NULL",
                   bsCollapsePanel(
                     "Instructions",
                     style = "primary",
                     # style = "height: 500px, width = 1100px",
                     htmlOutput("mult_genes_instructions")
                   )
                 )
               ),
               column(
                 width = 6,
                 bsCollapse(
                   id = "methods_mult", open = "NULL",
                   bsCollapsePanel(
                     "Materials and methods",
                     style = "primary",
                     # style = "height: 500px, width = 1100px",
                     htmlOutput("mult_genes_methods")
                   )
                 )
               )
             ),
             fluidRow(
               column(
                 width = 3,
                 style = "position: sticky; top: 20px; height: 100vh; overflow-y: auto;",  # Sidebar styling for stickiness
                 bsCollapse(
                   id = "gene_selection_mult", open = "Gene Selection",
                   bsCollapsePanel(
                     "Gene Selection",
                     radioButtons(inputId="List_gen",label="1. Choose or paste a list of genes", selected=character(0),
                                  choices=c("Choose a transcription factor family" = "fam",
                                            "Paste a list of genes" = "paste")),
                     conditionalPanel('input.List_gen === "fam"', selectizeInput(inputId = "TFfamil", label = "2. Pick a TF family",
                                                                                 choices = TF_list, selected = character(0), multiple =T, options=list(placeholder = '', maxItems=1))),
                     conditionalPanel('input.List_gen === "paste"', textareaInput(id="Genes_list_Tab2","2. Paste a list of genes","",rows=20),
                                      div(HTML("<i>One gene per line, no separator</i>"),style = "margin-bottom:15px")),
                     actionButton("show_heatmap", "Generate heatmap", status = "info"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("error_message_multiple"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("error_message_multiple2"),
                     div(style = "margin-bottom: 40px;"), 
                     textOutput("threshold_message_multiple"),
                     div(style = "margin-bottom: 40px;"), 
                     sliderInput("weight_thr_multiple", "Edge Weight Threshold", min = 0.1, max = 1, value = 0.1, step = 0.05),
                     actionButton("submit_table_multiple", "Get first neighbours", status = "info"),
                     div(style = "margin-bottom: 40px;"), 
                     actionButton("plot_button_multiple", "Plot Network", status = "info")))),
               
               column(
                 width = 9,
                 bsCollapse(
                   id = "expression_mult", open = NULL,
                   bsCollapsePanel(
                     "Expression Profile",
                     style = "primary",
                     InteractiveComplexHeatmapOutput("heatmap_output", height1 = 880, height2 = 880, width1 = 500, width2 = 500),
                     htmlOutput("plot_multiple_footnote"))
                 )
               ),
               column(
                 width = 9,
                 bsCollapse(
                   id = "neigh_mult", open = NULL,
                   bsCollapsePanel(
                     "First neighbours",
                     style = "primary",
                     div(
                       style = "display: flex; flex-direction: column; align-items: flex-start; overflow-x: auto; width: 100%;",
                       DT::dataTableOutput("node_table_multiple"),
                       div(style = "display: flex; justify-content: flex-start; margin-top: 20px;", # Add some margin to create space between the table and button
                           downloadButton("download_table_multiple", "Nodes", class = "btn-info", style = "width: 100px; margin-right: 10px; padding: 5px 10px;"),
                           downloadButton("download_edges_multiple", "Edges", class = "btn-info", style = "width: 100px; padding: 5px 10px;"))
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
                   id = "network_mult", open = NULL,
                   bsCollapsePanel(
                     "Co-expression network",
                     style = "primary",
                     fluidRow(
                       column(
                         width = 8,  # Adjust the width of the network plot column
                         forceNetworkOutput("network_plot_multiple", width = "100%")  # Adjust height as necessary
                       ),
                       column(
                         width = 3, 
                         plotOutput("network_legend_multiple", height = "600px")  # Set height for the legend plot
                       )
                     )
                   )
                 ))
             )),
    tabPanel("References and citations",
             fluidRow(
               column(
                 width = 6,
                 bsCollapse(
                   id = "datasets_tools", open = NULL,
                   bsCollapsePanel(
                     "Datasets and tools",
                     # style = "height: 500px, width = 1100px",
                     htmlOutput("tools")
                   )
                 )),
               column(
                 width = 6,
                 bsCollapse(
                   id = "cite", open = NULL,
                   bsCollapsePanel(
                     "How to cite POPUL-R",
                     # style = "height: 500px, width = 1100px",
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
    paste("<p style='text-align:justify'>", 'This section summarizes the main functionalities of POPUL-R. There is a more detailed description under the section "Instructions and methods" in each tab.', "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Single genes"</b>: By providing a <i>Populus tremula</i> gene id from <a href="https://plantgenie.org" target="_blank">PlantGenIE</a>, visualize the expression profile of the gene of interest indoors and outdoors, identify and select the first neighbours from a thresholded co-expression network and visualize the selected network.', "</p>","<p style='text-align:justify'>",
          '- <b>Tab "Multiple genes"</b>: By choosing a transcription factor family from our dataset or providing several <i>Populus tremula</i> gene ids from <a href="https://plantgenie.org" target="_blank">PlantGenIE</a>, visualize the expression profile of the selected genes indoors and outdoors as an interactive heatmap, identify and select the first neighbours from a thresholded co-expression network and visualize the selected network.', "</p>","<p style='text-align:justify'>")
  })

  output$single_genes_instructions <-  renderText({
    paste("<p style='text-align:justify'>", "This tab generates plots and tables from the dataset published in Marcon <i> et.al</i> (2024).", "</p>","<p style='text-align:justify'>",
          '- <b>Gene Selection and Expression Profile</b>: You need to enter a <i>Populus tremula </i>gene id from <a href="https://plantgenie.org" target="_blank">PlantGenIE</a> in the Gene Selection box such as Potra2n8c17315 (<i>FT1</i>) and click on "Plot Expression". You will obtain the expression profile of your gene of interest that you can download as .svg or .pdf. The raw data can be downloaded...', "</p>","<p style='text-align:justify'>",
          '- <b>First neighbours</b>: A message will appear in the Gene Selection panel with the highest edge weight of the gene: this value represents the highest connection between the selected genes and their first neighbours in the co-expression network. You can then choose an edge weight and click on "Get First Neighbours" to filter and obtain the first neighbours of your gene of interest. If the chosen value is higher than the maximum value provided, no connections are found from your gene of interest and the table appears empty. The table of first neighbours is ordered by edge weight in a descending order. If the gene name on the table is clicked,
          you will be directed to the <a href="https://plantgenie.org" target="_blank">PlantGenIE</a> entry of that gene. The centrality value displayed in the table is the centrality of the gene in the whole co-expression network. The module number corresponds to the module where the gene is found within the co-expression network. Once the table is displayed, you can download it as two separate .txt files, nodes.txt and edges.txt, that can be used in other network visualization softwares for further analyses.', "</p>","<p style='text-align:justify'>",
          '- <b>Co-expression network</b>: If you want to visualize the connections between your gene of interest and some first neighbours, you can select the genes directly by clicking on the table. The button "Plot Network" allows a visualization of your selected genes, displaying the nodes with colors that correspond to the Module number. The node size depends on the centrality of the gene.')
  })

  output$single_genes_methods <-  renderText({
    paste("The data was obtained from samples collected outdoors throughout a whole year and indoors throughout a growth cycle.", "</p>","<p style='text-align:justify'>",
          "The samples from outdoors were collected from ca. 1-year-old and 35-year-old local (Umeå, Sweden) aspen trees once a month around midday (buds from September to May, except November, and leaves from June to August).", "</p>","<p style='text-align:justify'>",
          "Our growth cycle conditions that simulate the yearly growth cycle of <i>Populus</i> trees have been previously described in <a href='https://pub.epsilon.slu.se/24748/1/andre_d_210629.pdf' target='_blank'>André (2021)</a>. Briefly, <i>in vitro </i>cultures of the trees were grown in jars with MS medium (Murashige and Skoog, 1962) until they were potted in soil for experiments. The trees were grown in long day (LD, 18h light/6h dark) for four weeks at ~20°C during the day and 18°C at night to simulate spring and summer. During this time, the plants were fertilized weekly. After four weeks, the trees were subjected to short day treatment
          (SD, 14h light/10h dark) to simulate autumn and induce dormancy. During this period, fertilization was stopped but the temperature remained constant at ~20/18 °C. After 15 weeeks, the trees were subjected to cold treatment (CT, 8h light/16h dark at 4°C) during 10 weeks to release dormancy. To induce bud flush, the trees were returned to the same LD conditions but fertilization only started after bud flush occured. Samples from buds and/or leaves were taken in each condition.", "</p>","<p style='text-align:justify'>",
          'RNA was extracted and pre-processed as described in <a href="https://www.sciencedirect.com/science/article/pii/S0960982222007825?via%3Dihub" target="_blank">André (2022)</a>. The expression heatmaps are generated with the packages <a href="https://academic.oup.com/bioinformatics/article/32/18/2847/1743594" target="blank">"ComplexHeatmap"</a> and <a href="https://academic.oup.com/bioinformatics/article/38/5/1460/6448211" target="blank">"InteractiveComplexHeatmap"</a>. All heatmaps are scaled by rows, and hierarchical clustering is used to cluster genes. Co-expression analysis was performed using the R package 
          <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559" target="blank">"WGCNA"</a>. The visualization of the co-expression network is generated with the R package <a href= "https://cran.r-project.org/web/packages/networkD3/index.html"target="blank"> "networkD3"</a>. For more details, please see the Materials and Methods section in Marcon <i> et.al</i> (2024).')
  })

  output$mult_genes_instructions <-  renderText({
    paste("<p style='text-align:justify'>", "This tab generates plots and tables from the dataset published in Marcon <i>et.al</i> (2024).", "</p>","<p style='text-align:justify'>",
          '- <b>Gene Selection and Expression Profile</b>: You need to select a transcription factor family from the dataset or enter several <i>Populus tremula </i>gene ids from <a href "https://plantgenie.org" target="_blank">PlantGenIE</a> in the Gene Selection box, such as Potra2n8c17315 (<i>FT1</i>). The genes should be entered one per line, with no separators. By clicking on "Generate Heatmap", an interactive heatmap of your sleected genes will be displayed. You can zoom-in sections of the heatmap by clicking on them, and they will be displayed in the panel "Selected sub-heatmap". The output section just displays
          detailed information about the sub-heatmap. You can download this heatmap and the sub-heatmap as .png, .svg or .pdf. The raw data fromt eh subheatmap can be downloaded by clicking on "Export to table". ', "</p>","<p style='text-align:justify'>",
          '- <b>First neighbours</b>: A message will appear in the Gene Selection panel with an interval of the highest edge weights of the chosen genes: since multiple genes have been selected and each has a highest edge weight value, the lowest and highest maximum edge weight value from all selected genes are displayed. A maximum edge weight value represents the highest connection between the selected genes and their first neighbours in the co-expression network. You can then choose an edge weight and click on on "Get First Neighbours" to filter and obtain the first neighbours of your genes of interest. This is interval of maximum edge weight values is shown so that, if you choose a higher edge value than the lowest value of the interval, you are aware
          that you will start missing some of your selected genes and first neighbours because the chosen edge value is higher than their maximum edge weight value. If the chosen value is higher than the maximum value from the interval, no connections are found from your gene of interest and the table appears empty. The table of first neighbours is ordered by centrality of the gene in the whole co-expression network in a descending order. There is no edge weight displayed on this table because several selected genes can be connected to the same first neigbours with different edge weight values. If the gene name on the table is clicked,
          you will be directed to the PlantGenIE entry of that gene. The module number corresponds to the module where the gene is found within the o-expression network. Once the table is displayed, you can download it as two separate .txt files, nodes.txt and edges.txt, that can be used in other network visualization softwares for further analyses.', "</p>","<p style='text-align:justify'>",
          '- <b>Co-expression network</b>: If you want to visualize the connections between your genes of interest and some first neighbours, you can select the genes directly by clicking on the table. The button "Plot Network" allows a quick visualization of your selected genes, displaying the nodes with colors that correspond to the Module number. The node size depends on the centrality of the gene.')
  })
    output$mult_genes_methods <-  renderText({
      paste("The data was obtained from samples collected outdoors throughout a whole year and indoors throughout a growth cycle.", "</p>","<p style='text-align:justify'>",
            "The samples from outdoors were collected from ca. 1-year-old and 35-year-old local (Umeå, Sweden) aspen trees once a month around midday (buds from September to May, except November, and leaves from June to August).", "</p>","<p style='text-align:justify'>",
            "Our growth cycle conditions that simulate the yearly growth cycle of <i>Populus</i> trees have been previously described in <a href='https://pub.epsilon.slu.se/24748/1/andre_d_210629.pdf' target='_blank'>André (2021)</a>. Briefly, <i>in vitro </i>cultures of the trees were grown in jars with MS medium (Murashige and Skoog, 1962) until they were potted in soil for experiments. The trees were grown in long day (LD, 18h light/6h dark) for four weeks at ~20°C during the day and 18°C at night to simulate spring and summer. During this time, the plants were fertilized weekly. After four weeks, the trees were subjected to short day treatment
          (SD, 14h light/10h dark) to simulate autumn and induce dormancy. During this period, fertilization was stopped but the temperature remained constant at ~20/18 °C. After 15 weeeks, the trees were subjected to cold treatment (CT, 8h light/16h dark at 4°C) during 10 weeks to release dormancy. To induce bud flush, the trees were returned to the same LD conditions but fertilization only started after bud flush occured. Samples from buds and/or leaves were taken in each condition.", "</p>","<p style='text-align:justify'>",
            'RNA was extracted and pre-processed as described in <a href="https://www.sciencedirect.com/science/article/pii/S0960982222007825?via%3Dihub" target="_blank">André (2022)</a>. The expression heatmaps are generated with the packages <a href="https://academic.oup.com/bioinformatics/article/32/18/2847/1743594" target="blank">"ComplexHeatmap"</a> and <a href="https://academic.oup.com/bioinformatics/article/38/5/1460/6448211" target="blank">"InteractiveComplexHeatmap"</a>. All heatmaps are scaled by rows, and hierarchical clustering is used to cluster genes. Co-expression analysis was performed using the R package 
          <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-559" target="blank">"WGCNA"</a>. The visualization of the co-expression network is generated with the R package <a href= "https://cran.r-project.org/web/packages/networkD3/index.html"target="blank"> "networkD3"</a>. For more details, please see the Materials and Methods section in Marcon <i> et.al</i> (2024).')
  })

  output$tools <- renderText({#DT package!! And Github
    paste("The source code is available on GitHub...","</p>","<p style='text-align:justify'>",
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


  updateSelectizeInput(session, "gene", choices = c("", unique(expression_data$Gene), server = TRUE))
  
  values <- reactiveValues(gene = NULL, selected_genes = NULL, max_weight = NULL, threshold_message = NULL, error_message = NULL, thr = NULL , list_genes = NULL, max_weight_mult = NULL, max_min_weight_mult = NULL, threshold_message_multiple = NULL, error_message_multiple = NULL, thr_mult = NULL, selected_genes_mult = NULL)

  process_gene <- function(values, input_gene, input_weight_thr, edges) {

    values$error_message <- NULL
    values$threshold_message <- NULL
    values$max_weight <- NULL
    values$thr <- NULL

    # Validate the gene input
    if (is.null(input_gene) || input_gene == "") {
      values$error_message <- "This gene is not in the dataset"
    } else {
      # Gene is valid; proceed with processing
      values$error_message <- NULL

      # Filter edges related to the gene
      filtered_edges <- edges %>%
        filter(fromNode == input_gene | toNode == input_gene)

      # Compute the maximum weight
      max_wt <- max(filtered_edges$weight, na.rm = TRUE)
      if (is.infinite(max_wt)) {
        values$max_weight <- NULL
        values$threshold_message <- "This gene is below the minimum network threshold."
      } else {
        values$max_weight <- max_wt
        print(values$max_weight)
        values$threshold_message <- paste(
          "The highest edge weight of your gene of interest is",
          trunc(values$max_weight * 10^2) / 10^2,
          ". If the chosen network threshold is higher than this value, the table and the network will appear empty."
        )
      }
    }

    # Update the threshold value
    values$thr <- input_weight_thr
  }

  observeEvent(input$submit, {
    values$gene <- input$gene  # Update values$gene
    print("Gene captured in submit")
    print(values$gene)

    process_gene(values, values$gene, input$weight_thr, edges)
  })

  # Update values$gene when submit_table button is clicked
  observeEvent(input$submit_table, {
    values$gene <- input$gene  # Update values$gene
    print("Gene captured in submit_table")
    print(values$gene)

    process_gene(values, values$gene, input$weight_thr, edges)
  })

  observe({
    print("Observer triggered")  # Debugging: see if observer is triggered
    # req(values$gene)
    print(isolate(values$gene))  # Safely print the gene value
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

  #Create reactive plot data
  data_single_out <- eventReactive(input$submit, {
    req(values$gene)  # Ensure that gene is set

    # If there's an error, return NULL to prevent plotting
    if (!is.null(values$error_message)) {
      return(NULL)
    }

    df <- expression_data %>%
      filter(Gene == values$gene) %>%
      group_by(`Month/Treatment`, Location, Treatment2) %>%
      summarize(mean = mean(Expression),
                se = std.error(Expression)) %>%
      arrange(Treatment2) %>%
      filter(Location == "Outdoor")
    df$`Month/Treatment` <- factor(df$`Month/Treatment`, levels = c("SEP",
                                                              "OCT",
                                                              "DEC",
                                                              "JAN",
                                                              "FEB",
                                                              "MAR",
                                                              "APR",
                                                              "MAY",
                                                              "JUN",
                                                              "JUL",
                                                              "AUG"))
    df
  })

  data_single_gh <- eventReactive(input$submit, {
    req(values$gene)  # Ensure that gene is set

    # If there's an error, return NULL to prevent plotting
    if (!is.null(values$error_message)) {
      return(NULL)
    }

    df <- expression_data %>%
      filter(Gene == values$gene) %>%
      group_by(`Month/Treatment`, Location, Treatment2) %>%
      summarize(mean = mean(Expression),
                se = std.error(Expression)) %>%
      arrange(Treatment2) %>%
      filter(Location == "Indoor")
    df$`Month/Treatment` <- factor(df$`Month/Treatment`, levels = c("SD15","CT2","CT8","CT10",
                                                              "LD1","LD2","LD3", "LD4",
                                                              "SD1","SD2","SD3","SD10"))
    df

  })
  
  # Render Expression Plots
  output$expression_plots <- renderPlot({
    req(input$submit)
    if (!is.null(values$error_message)) {
      ggplot() +
        theme_void() +
        labs(title = "")
    } else {
      
   
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
      

    filtered_data <- expression_data %>%
      dplyr::filter(Gene == values$gene)
    
    plot_expression <- function(data, labels) {
      ggplot(data, aes(x = `Month/Treatment`, y = mean, ymin = mean - se,
                       ymax = mean + se)) +
        geom_line(aes(group = 1), color = "#009E73",
                  linewidth = 1.7) +
        geom_ribbon(aes(group = 1, y = mean, ymin = mean - se, ymax = mean + se), fill = "#009E73", alpha = 0.2) +
        labs(x = NULL,
             y = "Expression",
             title = values$gene) +
        scale_x_discrete(name = NULL,
                         labels = labels)+
        ylim(-0.3, round(max(filtered_data$Expression) + 2)) +
        theme_classic() +
        theme(plot.title = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 0, b = 20, l = 0)),
              axis.title = element_text(size = 15,color = "black"),
              axis.title.y = element_text(margin = margin(r = 20)),
              axis.text.y = element_text(size = 13, color = "black"),
              axis.line = element_line(color = "black"),
              legend.position = "none",
              axis.text.x = ggtext::element_markdown(size = 13, color = "black"))
    }
    

    # p1 <- ggplot(data_single_out(), aes(x = `Month/Treatment`, y = mean, ymin = mean - se,
    #                                     ymax = mean + se)) +
    #   geom_line(aes(group = 1), colour = "#67C77B",
    #             linewidth = 1.7) +
    #   geom_ribbon(aes(group = 1, y = mean, ymin = mean - se, ymax = mean + se), fill = "#009E73", alpha = 0.2) +
    #   labs(x = NULL,
    #        y = "Expression",
    #        title = values$gene) +
    #   scale_x_discrete(name = NULL,
    #                    labels = labels_outdoors)+
    #   ylim(-0.3, round(max(filtered_data$Expression)) + 2) +
    #   theme_classic() +
    #   theme(plot.title = element_text(size = 24, face = "bold", margin = margin(t = 0, r = 0, b = 20, l = 0)),
    #         axis.title = element_text(size = 15,color = "black"),
    #         axis.title.y = element_text(margin = margin(r = 20)),
    #         axis.text.y = element_text(size = 13, color = "black"),
    #         axis.line = element_line(color = "black"),
    #         legend.position = "none",
    #         axis.text.x = ggtext::element_markdown(size = 13, color = "black"))
    
  p1 <- plot_expression(data_single_out(),labels_out) 

    # p2 <- ggplot(data_single_gh(), aes(x = `Month/Treatment`, y = mean, ymin = mean - se,
    #                                    ymax = mean + se)) +
    #   geom_line(aes(group = 1), colour = "#009E73",
    #             linewidth = 1.7) +
    #   geom_ribbon(aes(group = 1, y = mean, ymin = mean - se, ymax = mean + se), fill = "#009E73", alpha = 0.2) +
    #   labs(x = NULL,
    #        y = "Expression",
    #        title = values$gene) +
    #   scale_x_discrete(name = NULL,
    #                    labels = labels_greenhouse)+
    #   ylim(-0.3, round(max(filtered_data$Expression) + 2)) +
    #   theme_classic() +
    #   theme(plot.title = element_blank(),
    #         axis.title = element_text(size = 15,color = "black"),
    #         axis.title.y = element_text(margin = margin(r = 20)),
    #         axis.text.y = element_text(size = 13, color = "black"),
    #         axis.line = element_line(color = "black"),
    #         legend.position = "none",
    #         axis.text.x = ggtext::element_markdown(size = 13, color = "black"))
  
  p2 <- plot_expression(data_single_gh(), labels_gh) +
    theme(plot.title = element_blank())

    plot_grid(p1, p2, nrow = 1, align = "h", rel_widths = c(0.92, 1))
    }
  })



  output$plot_single_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. LD: long day. The number shown is the number of weeks in the treatment. The illustrations represent buds and leaves. Data are means +/- SE for n = 6 biological replicates.", "</p>"
    )
  })

  output$download_plot_svg <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_expression_plot_", values$gene, ".svg", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "svg", width = 16, height = 5)
    }
  )

  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste(Sys.Date(),"_expression_plot_", values$gene, ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "pdf", width = 16, height = 5)
    }
  )
  
  output$download_expression_data <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_expression_data_", values$gene,".csv", sep = "")
    },
    content = function(file) {
      write.csv(rbind(data_single_out(), data_single_gh()), file, row.names = FALSE)
    })

  observeEvent(input$submit, {
    output$node_table <- DT::renderDT({
      # Create an empty dataframe with the same structure as filtered_nodes_ann
      empty_df <- data.frame("Gene name" = character(),
                             "Module" = character(),
                             "Centrality" = numeric(),
                             "ATG symbol" = character(),
                             "ATG full name" = character(),
                             stringsAsFactors = FALSE)
      DT::datatable(empty_df)
    })
  })

  # #Make table with selected gene and 1st neighbours
  filtered_nodes_ann <- eventReactive(input$submit_table, {
    req(values$gene)  # Ensure gene is not NULL
    gene_to_check <- input$gene
    # gene_to_check <- "Potra2n1c3477"
    filtered_edges <- edges %>%
      filter((fromNode == gene_to_check | toNode == gene_to_check) & weight >= values$thr) # Combine filtering conditions

      if (nrow(filtered_edges) > 0) {

        goi_edge <- data.frame(fromNode = gene_to_check,
                               toNode = gene_to_check,
                               weight = 1,
                               direction = "undirected",
                               fromAltName = gene_to_check,
                               toAltName = gene_to_check)

        filtered_edges <- rbind(goi_edge, filtered_edges)
      }

    filtered_nodes <- data.frame(id = unique(c(filtered_edges$fromNode, filtered_edges$toNode)),
                                 # Module = nodes$Module[nodes$nodeName %in% unique(c(filtered_edges$fromNode, filtered_edges$toNode))],
                                 # Degree = nodes$Degree[nodes$nodeName %in% unique(c(filtered_edges$fromNode, filtered_edges$toNode))],
                                 Edge_weight = filtered_edges$weight) #Nodes are in the same order in node and edge file

    filtered_nodes$nodeID <- match(filtered_nodes$id, unique(filtered_nodes$id)) - 1

    ann <- filtered_nodes %>%
      left_join(subannot, by = c("id" = "Gene name")) %>%
      separate(Module, c('Module_color', 'Module'), sep = "/") %>%
      arrange(desc(Edge_weight)) %>%
      dplyr::select(-nodeID, -GOI, -Module_color)

    colnames(ann)[1] <- c("Gene name")

    if (nrow(filtered_edges) > 0) {

    ann[[1]] <- paste0('<a href="https://plantgenie.org/gene?id=', ann[[1]], '" target="_blank">', ann[[1]], '</a>')
    }
    return(ann)
})

  # Render Node Table with selectable rows
  observeEvent(input$submit_table, {
  output$node_table <- DT::renderDT({
    req(filtered_nodes_ann())
    print("Rendering table with filtered nodes")# Ensure data is available
    DT::datatable(
      filtered_nodes_ann(),
        options = list(
          pageLength = 10,  # Show 10 rows per page by default
          lengthMenu = list(c(10, 15, 20, -1), c('10', '15', '20', 'All')),  # Options to show 10, 15, 20 or all rows
          order = list(list(1, 'desc')),  # Order by the first column initially
          searching = TRUE,  # Enable search box
          ordering = TRUE,  # Enable column ordering
          dom = 'lfrtip',  # Show length menu, filtering input, and pagination controls
          columnDefs = list(list(targets = ncol(filtered_nodes_ann()) - 1, orderable = FALSE)),  # Disable ordering on the checkbox column
          language = list(emptyTable = "Your gene of interest is below the chosen network threshold")  # Custom error message
        ),
        escape = FALSE,
        selection = 'multiple',
        rownames = FALSE
      )
  })
  })


  #Create nodes without links
  filtered_nodes_download <- eventReactive(input$submit_table, {
    req(filtered_nodes_ann())
    print("Triggered filtered_nodes_download")

    `Gene name` <- gsub('.*>(.*)<.*', '\\1', filtered_nodes_ann()$`Gene name`)
    print("Gene names for nodes extracted")
    cbind(`Gene name`,
          filtered_nodes_ann() %>%
            select(!`Gene name`))
  })

 #Create edges from table
  filtered_edges_download <- eventReactive(input$submit_table, {
    req(filtered_nodes_ann())
    print("Triggered filtered_edges_download")

    ann_gene_names <- gsub('.*>(.*)<.*', '\\1', filtered_nodes_ann()$`Gene name`)
    print("Gene names for edges extracted")
    edges %>%
      filter(fromNode %in% ann_gene_names & toNode %in% ann_gene_names)
  })

  # Add download handler for node table
  output$download_table <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_node_table_", values$gene,"_thr", values$thr, ".txt", sep = "")
    },
    content = function(file) {
      write.table(filtered_nodes_download(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })

  # Add download for edges
  output$download_edges <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_edges_table_", values$gene, "_thr", values$thr, ".txt", sep = "")
    },
    content = function(file) {
      write.table(filtered_edges_download(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })



  # # Update selected genes when plot button is clicked
  observeEvent(input$plot_button, {
    req(input$submit_table, values$gene, values$thr)

        selected_rows <- input$node_table_rows_selected #???
        if (!is.null(selected_rows)) {
          # Extract the gene names and apply gsub to each element individually
          gene_names <- filtered_nodes_ann()[selected_rows, "Gene name"]

          # Use sapply to apply gsub to each gene name in the vector
          values$selected_genes <- sapply(gene_names, function(x) gsub('.*>(.*)<.*', '\\1', x))
          # values$selected_genes <- gsub('.*>(.*)<.*', '\\1', filtered_nodes_ann()[selected_rows, "Gene name"])
        } else {
          values$selected_genes <- NULL
        }
        print(paste("Selected genes: ", paste(values$selected_genes, collapse = ", ")))
      })

  # Create a reactive expression to hold the network data
  network_data <- reactiveVal(NULL)

  # Observe button click to update network data
  observeEvent(input$plot_button, {
    # Ensure that inputs are available and valid
    req(values$gene, values$selected_genes, input$plot_button)

    selected_genes <- values$selected_genes
    gene_to_check <- input$gene

    filtered_edges <- edges %>%
      filter((fromNode == gene_to_check | toNode == gene_to_check) & weight >= values$thr)

    if (length(selected_genes) == 0) {
      filtered_edges_network <- filtered_edges
    } else {
      filtered_edges_network <- edges %>%
        filter(fromNode %in% selected_genes & toNode %in% selected_genes)
    }

    if (nrow(filtered_edges_network) == 0) {
      network_data(NULL)  # Set to NULL if no data
      return(NULL)
    }

    filtered_nodes <- data.frame(id = unique(c(filtered_edges$fromNode, filtered_edges$toNode)))
    filtered_nodes$nodeID <- match(filtered_nodes$id, unique(filtered_nodes$id)) - 1

    filtered_nodes_network <- nodes %>%
      dplyr::select(-c(Degree, Module)) %>%
      filter(nodeName %in% unique(c(filtered_edges_network$fromNode, filtered_edges_network$toNode))) %>%
      left_join(filtered_nodes, by = c("nodeName" = "id")) %>%
      left_join(subannot %>% select(`Gene name`, Centrality, Module), by = c("nodeName" = "Gene name")) %>%
      separate(Module, c('Module_color', 'Module'), sep = "/")

    filtered_edges_network$fromNodeID <- match(filtered_edges_network$fromNode, filtered_nodes_network$nodeName) - 1
    filtered_edges_network$toNodeID <- match(filtered_edges_network$toNode, filtered_nodes_network$nodeName) - 1

    network_data(list(edges = filtered_edges_network, nodes = filtered_nodes_network))
  })

  # Render network plot based on the reactive network_data
  output$network_plot <- renderForceNetwork({
    network_info <- network_data()
    req(network_info)  # Ensure network data is available

    network_plot <- forceNetwork(
      Links = network_info$edges,
      Nodes = network_info$nodes,
      Source = "fromNodeID",
      Target = "toNodeID",
      NodeID = "nodeName",
      Group = "Module",
      opacity = 0.8,
      zoom = TRUE,
      linkColour = "grey",
      linkWidth = 2,
      fontSize = 20,
      linkDistance = 300,
      colourScale = networkD3::JS('d3.scaleOrdinal().domain(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46"]).range(["#40E0D0", "#0000FF", "#A52A2A", "#FFFF00", "#00FF00", "#FF0000", "#000000", "#FFC0CB", "#FF00FF", "#A020F0", "#ADFF2F", "#D2B48C", "#FA8072", "#00FFFF", "#191970", "#E0FFFF", "#999999", "#90EE90", "#FFFFE0", "#4169E1", "#8B0000", "#006400", "#00CED1", "#A9A9A9", "#FFA500", "#FF8C00", "#EEE9E9", "#87CEEB", "#8B4513", "#AFEEEE", "#4682B4", "#EE82EE", "#556B2F", "#8B008B", "#CD6839", "#9ACD32", "#6CA6CD", "#FFBBFF", "#8B2500", "#8968CD", "#CAE1FF", "#E0FFFF", "#FFFFF0", "#FFFAF0", "#EE7600", "#8B2323"])'),
      Nodesize = "Centrality",
      legend = FALSE
    )
    network_plot <- htmlwidgets::onRender(
      network_plot,
      '
  function(el, x) {
    d3.selectAll(".node circle")
      .style("stroke", "grey")
      .style("stroke-width", "2px")

    d3.selectAll(".node text")
      .style("fill", "#B3BAB9");
  }
  '
    )

    network_plot
  })


  # Render legend plot
  output$network_legend <- renderPlot({
    req(values$gene, values$selected_genes, input$plot_button)
    suppressWarnings({
    grid.newpage()

    grid.roundrect(
      x = 0.6, y = 0.5,
      width = unit(6, "cm"), height = unit(14.5, "cm"),  # Adjust width and height of the legend background
      r = unit(0.5, "cm"),  # Set the corner radius, tweak as needed
      gp = gpar(fill = "#C1D6E0", col = NA)  # Fill color for the legend background
    )
    pushViewport(viewport(x = 0.6, y = 0.5, just = "center"))  # Move the legend to the right
    grid.draw(legend)
    })

  })


#Tab Multiple genes
  ##From CAST-R app

  # Function to handle gene list creation and validation
  update_gene_list <- function(values, input, subannot) {
    if (input$List_gen == "fam") {
      # validate(
      #   need(input$TFfamil != "", "Select a TF family before clicking on any button")
      # )
      ID <- input$TFfamil
      ID_list <- filter(subannot, Family %in% ID)
      values$list_genes <- ID_list$`Gene name`
    } else if (input$List_gen == "paste") {
      # validate(
      #   need(input$Genes_list_Tab2 != "", "Provide a list of genes before clicking on any button")
      # )
      values$list_genes <- unlist(strsplit(input$Genes_list_Tab2, "\n"))
    }
    print(values$list_genes)
    return(values$list_genes)
  }

  # Function to process gene weights and thresholds
  process_gene_weights <- function(values, edges) {
    # Initialize error and threshold messages every time the function is called
    values$error_message_multiple <- NULL
    values$error_message_multiple2 <- NULL
    values$threshold_message_multiple <- NULL
    values$max_weight_mult <- NULL
    values$max_min_weight_mult <- NULL

    if (is.null(values$list_genes) || any(values$list_genes == "")) {
      values$error_message_multiple <- "These genes are not in the dataset."
    } else if (length(values$list_genes) == 1) {
      values$error_message_multiple <- "Please provide at least two genes."
    } else {
      values$error_message_multiple <- NULL
    }

    # Proceed with the rest of the logic only if there's no error
    if (is.null(values$error_message_multiple)) {
      filtered_edges_mult <- edges %>%
        filter(fromNode %in% values$list_genes | toNode %in% values$list_genes)

      if (is.infinite(max(filtered_edges_mult$weight, na.rm = TRUE))) {
        values$max_weight_mult <- NULL
        values$threshold_message_multiple <- "These genes are below the minimum network threshold."
      } else {
        max_weight_goi <- filtered_edges_mult %>%
          gather(key = "type", value = "Gene", fromNode, toNode) %>%
          filter(Gene %in% values$list_genes) %>%
          group_by(Gene) %>%
          summarise(max_weight = max(weight, na.rm = TRUE)) %>%
          ungroup()

        values$max_weight_mult <- max(max_weight_goi$max_weight, na.rm = TRUE)
        values$max_min_weight_mult <- min(max_weight_goi$max_weight, na.rm = TRUE)
        print(values$max_weight_mult)
        print(values$max_min_weight_mult)

        values$threshold_message_multiple <- paste("The interval of maximum edge weights of your genes of interest is",
                                                   trunc(values$max_min_weight_mult * 10^2) / 10^2, "-",
                                                   trunc(values$max_weight_mult * 10^2) / 10^2,
                                                   ". If the chosen network threshold is higher than the lowest value of the interval, some genes of interest will not be shown in the table or the network. If the chosen network threshold is higher than the highest value of the interval, the table and the network will appear empty.")
      }
    }
    values$thr_mult <- input$weight_thr_multiple
  }

  # Observe heatmap button
  observeEvent(input$show_heatmap, {
    update_gene_list(values, input, subannot)  # Update gene list
    process_gene_weights(values, edges)          # Process gene weights and thresholds
  })

  # Observe submit_table_multiple button
  observeEvent(input$submit_table_multiple, {
    update_gene_list(values, input, subannot)  # Update gene list
    process_gene_weights(values, edges)          # Process gene weights and thresholds
  })

  # Debugging observer
  observe({
    print("Observer triggered")  # Debugging: see if observer is triggered
    print(isolate(values$list_genes))  # Safely print the gene value
  })

  # Display error message if any
  output$error_message_multiple <- renderText({
    values$error_message_multiple
  })

  observeEvent(input$submit_table_multiple, {
    values$thr_mult <- input$weight_thr_multiple
    print(values$thr_mult)
  })



  # Render max weight message
  output$threshold_message_multiple <- renderText({
    values$threshold_message_multiple
  })



  # Reactive expression to generate the heatmap
  heatmap_data <- reactive({
    # if (!is.null(values$error_message_multiple2)) {
    #   return(NULL)
    # }

    req(input$List_gen)
    print("Making heatmap data")

    mult <- filter(expression_data, Gene %in% values$list_genes)

    mult <- mult %>%
      group_by(Gene, `Month/Treatment`) %>%
      summarize(mean = mean(Expression))

    mult$`Month/Treatment` <- factor(mult$`Month/Treatment`, levels = c("SEP", "OCT", "DEC", "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SD15","CT2","CT8","CT10", "LD1","LD2","LD3", "LD4", "SD1","SD2","SD3","SD10"))

    mult <- mult %>%
      arrange(`Month/Treatment`) %>%
      pivot_wider(names_from =  `Month/Treatment`, values_from = mean) %>%
      column_to_rownames("Gene") %>%
      mutate(Total = rowSums(.)) %>%
      filter(Total != 0) %>%
      select(-Total)

    return(mult)
  })

  # Observe heatmap generation when button is clicked
  observeEvent(input$show_heatmap, {


    # values$error_message_multiple2 <- NULL

    withProgress(message = 'Generating heatmap...', value = 0, {
    mult <- heatmap_data()  # Get heatmap data
    incProgress(0.5)  # Update progress

    # Check the length of input$List_gen
    if (input$List_gen == "paste" && length(values$list_genes) == 1) {
      mult <- NULL  # Assign NULL to mult for a single gene
    }

    if (nrow(mult) == 0 || is.null(mult)) {
      output$error_message_multiple2 <- renderText({"No data available for the selected genes."})

      # Create a blank heatmap with a message
      blank_matrix <- matrix(0, nrow = 1, ncol = 1)  # 1x1 matrix to create a blank heatmap
      colnames(blank_matrix) <- c("")  # Optional: Give a name to the column
      rownames(blank_matrix) <- c("")  # No row name for blank
      heatmap <- ComplexHeatmap::pheatmap(mat = blank_matrix,
                                          border_color = NA,
                                          color = colorRampPalette(c("dodgerblue", "white", "firebrick"))(10),
                                          breaks = seq(-1, 1, length.out = 10),
                                          name = "Expression")

      ht = draw(heatmap)
      incProgress(0.5)
      makeInteractiveComplexHeatmap(input, output, session, ht_list = ht)  # Assuming you have a plot output for the heatmap
      return(NULL)
    }

    missing_genes <- setdiff(values$list_genes, rownames(mult))
    if (length(missing_genes) > 0) {
      output$error_message_multiple2 <- renderText({
        paste("Some of these genes are not in the dataset:",
              paste(missing_genes, collapse = ", "),
              ". The heatmap will be rendered with only the genes that appear in the dataset.")
      })
    } else {
      output$error_message_multiple2 <- renderText(NULL)
    }
    # Debugging: print data dimensions and first few rows
    print(dim(mult))  # Check dimensions of the matrix
    print(head(mult))  # Check first few rows of the matrix

    annot_col <- data.frame(
      Tissue   = samples.sep$Tissue,
      `Month/Treatment` = samples.sep$`Month/Treatment`,
      Location = samples.sep$Location,
      check.names = FALSE) %>%
      unique()
    rownames(annot_col) <- colnames(mult)
    annot_col$`Month/Treatment` <- factor(annot_col$`Month/Treatment`, levels = c("SEP",
                                                                                  "OCT",
                                                                                  "DEC",
                                                                                  "JAN",
                                                                                  "FEB",
                                                                                  "MAR",
                                                                                  "APR",
                                                                                  "MAY",
                                                                                  "JUN",
                                                                                  "JUL",
                                                                                  "AUG",
                                                                                  "SD15","CT2","CT8","CT10",
                                                                                  "LD1","LD2","LD3", "LD4",
                                                                                  "SD1","SD2","SD3","SD10"))
    annot_col$Location <- factor(annot_col$Location, levels = c("Outdoor",
                                                                "Indoor"))


    vals <- unique(samples.sep$Location)
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
    print("Annotation colors saved")

    if(nrow(mult) == 1) {
      heatmap <- ComplexHeatmap::pheatmap(mat = as.matrix(mult),
                                          cluster_rows = FALSE,
                                          cluster_cols = FALSE, #dist.var.tree.GLOBAL,
                                          scale = "row",
                                          legend = TRUE,
                                          border_color = NA,
                                          color = colorRampPalette(c("dodgerblue","white","firebrick"))(10),
                                          fontsize = 30,
                                          fontsize_row = 14,
                                          fontsize_col = 12,
                                          # srtCol = 45,
                                          show_rownames = TRUE,
                                          show_colnames = FALSE,
                                          #labels_col = names,
                                          annotation_legend = TRUE,
                                          annotation_col = annot_col,
                                          #annotation_row = annot_row,
                                          annotation_colors = annot_colors,
                                          name = "Expression")
      ht = draw(heatmap)
      incProgress(0.5)
      makeInteractiveComplexHeatmap(input, output, session, ht_list = ht)
    } else {


    dist.obs <- as.dist(1-cor(t(mult)))
    dist.obs.tree <- hclust(dist.obs, method = "ward.D")

    heatmap <- ComplexHeatmap::pheatmap(mat = as.matrix(mult),
                                        cluster_rows = dist.obs.tree,
                                        cluster_cols = FALSE, #dist.var.tree.GLOBAL,
                                        scale = "row",
                                        legend = TRUE,
                                        border_color = NA,
                                        color = colorRampPalette(c("dodgerblue","white","firebrick"))(10),
                                        fontsize = 30,
                                        fontsize_row = 10,
                                        fontsize_col = 12,
                                        # srtCol = 45,
                                        show_rownames = TRUE,
                                        show_colnames = FALSE,
                                        #labels_col = names,
                                        annotation_legend = TRUE,
                                        annotation_col = annot_col,
                                        #annotation_row = annot_row,
                                        annotation_colors = annot_colors,
                                        name = "Expression")
    ht = draw(heatmap)
    # Draw the heatmap and make it interactive
    incProgress(0.5)  # Update progress
    makeInteractiveComplexHeatmap(input, output, session, ht_list = ht)
    }
  })
})

  output$plot_multiple_footnote <-renderText({
    paste("<p style='text-align:justify'>","Expression data in VST counts. SD: short day. CT: cold treatment. LD: long day. The number shown is the number of weeks in the treatment. Data are means for n = 6 biological replicates.", "</p>"
    )
  })

  observeEvent(input$show_heatmap, {
    print("submit_table_multiple observed")
    output$node_table_multiple <- DT::renderDT({
      # Create an empty dataframe with the same structure as filtered_nodes_ann
      empty_df <- data.frame("Gene name" = character(),
                             "Module" = character(),
                             "Centrality" = numeric(),
                             "ATG symbol" = character(),
                             "ATG full name" = character(),
                             stringsAsFactors = FALSE)
      DT::datatable(empty_df)
    })
  })

  #Make table with selected gene and 1st neighbours
  filtered_edges_mult <- eventReactive(input$submit_table_multiple, {
    print(paste("submit_table_multiple triggered"))
    print("Rendering data for edges multiple")
    print(values$list_genes)
    print(values$thr_mult)

    req(input$List_gen)
    filtered_edges_mult <- edges %>%
           filter((fromNode %in% values$list_genes | toNode %in% values$list_genes) & weight >= values$thr_mult)

    if (nrow(filtered_edges_mult) > 0) {

      goi_edge_mult <- data.frame(fromNode =  values$list_genes[values$list_genes %in% unique(c(filtered_edges_mult$fromNode, filtered_edges_mult$toNode))],
                                  toNode = values$list_genes[values$list_genes %in% unique(c(filtered_edges_mult$fromNode, filtered_edges_mult$toNode))],
                                  weight = 1,
                                  direction = "undirected",
                                  fromAltName = values$list_genes[values$list_genes %in% unique(c(filtered_edges_mult$fromNode, filtered_edges_mult$toNode))],
                                  toAltName = values$list_genes[values$list_genes %in% unique(c(filtered_edges_mult$fromNode, filtered_edges_mult$toNode))])

      filtered_edges_mult <- rbind(goi_edge_mult, filtered_edges_mult)
    }
    head(filtered_edges_mult)
    return(filtered_edges_mult)

  })

  filtered_nodes_ann_mult <- eventReactive(input$submit_table_multiple, {
    print(paste("submit_table_multiple triggered2"))
    print("Rendering data for table multiple")
    req(input$List_gen)

    filtered_edges_data <- filtered_edges_mult()

    filtered_nodes_mult <- data.frame(id = unique(c(filtered_edges_data$fromNode, filtered_edges_data$toNode)))

    filtered_nodes_mult$nodeID <- match(filtered_nodes_mult$id, unique(filtered_nodes_mult$id)) - 1

    ann_mult1 <- filtered_nodes_mult %>%
      filter(id %in% values$list_genes) %>%
      left_join(subannot, by = c("id" = "Gene name")) %>%
      separate(Module, c('Module_color', 'Module'), sep = "/") %>%
      arrange(desc(Centrality))

    ann_mult2 <- filtered_nodes_mult %>%
      filter(!(id %in% values$list_genes)) %>%
      left_join(subannot, by = c("id" = "Gene name")) %>%
      separate(Module, c('Module_color', 'Module'), sep = "/") %>%
      arrange(desc(Centrality))

    ann_mult <- rbind(ann_mult1, ann_mult2) %>%
      dplyr::select(-nodeID, -GOI, -Module_color)

    colnames(ann_mult)[1] <- c("Gene name")

    if (nrow(filtered_edges_mult()) > 0) {

      ann_mult[[1]] <- paste0('<a href="https://plantgenie.org/gene?id=', ann_mult[[1]], '" target="_blank">', ann_mult[[1]], '</a>')
    }
    head(ann_mult)
    return(ann_mult)

  })

  # Render Node Table with selectable rows
  observeEvent(input$submit_table_multiple, {
    output$node_table_multiple <- DT::renderDT({
      req(filtered_nodes_ann_mult())
      print("Rendering table with filtered nodes")# Ensure data is available
      DT::datatable(
        filtered_nodes_ann_mult(),
        options = list(
          pageLength = 10,  # Show 10 rows per page by default
          lengthMenu = list(c(10, 15, 20, -1), c('10', '15', '20', 'All')),  # Options to show 10, 15, 20 or all rows
          #order = list(list(1, 'desc')),  # Order by the first column initially
          searching = TRUE,  # Enable search box
          ordering = TRUE,  # Enable column ordering
          dom = 'lfrtip',  # Show length menu, filtering input, and pagination controls
          columnDefs = list(list(targets = ncol(filtered_nodes_ann_mult()) - 1, orderable = FALSE)),  # Disable ordering on the checkbox column
          language = list(emptyTable = "Your genes of interest are below the chosen network threshold")  # Custom error message
        ),
        escape = FALSE,
        selection = 'multiple',
        rownames = FALSE
      )
    })
  })

  #Create nodes without links
  filtered_nodes_mult_download <- eventReactive(input$submit_table_multiple, {
    req(filtered_nodes_ann_mult())

    `Gene name` <- gsub('.*>(.*)<.*', '\\1', filtered_nodes_ann_mult()$`Gene name`)
    cbind(`Gene name`,
          filtered_nodes_ann_mult() %>%
            select(!`Gene name`))
  })

  filtered_edges_mult_download <- eventReactive(input$submit_table_multiple, {
    req(filtered_nodes_ann_mult())

    ann_gene_names_mult <- gsub('.*>(.*)<.*', '\\1', filtered_nodes_ann_mult()$`Gene name`)
    edges %>%
      filter(fromNode %in% ann_gene_names_mult & toNode %in% ann_gene_names_mult)
  })

  # Add download handler for node table
  output$download_table_multiple <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_node_table_multiple_genes","_thr", values$thr_mult, ".txt", sep = "")
    },
    content = function(file) {
      write.table(filtered_nodes_mult_download(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })

  # Add download for edges
  output$download_edges_multiple <- downloadHandler(
    filename = function() {
      paste(Sys.Date(), "_edges_table_multiple_genes", "_thr", values$thr_mult, ".txt", sep = "")
    },
    content = function(file) {
      write.table(filtered_edges_mult_download(), file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    })

# Update selected genes when plot button is clicked
  observeEvent(input$plot_button_multiple, {
    req(input$submit_table_multiple, values$thr_mult, input$List_gen)

    selected_rows <- input$node_table_multiple_rows_selected #???
    if (!is.null(selected_rows)) {
      # Extract the gene names and apply gsub to each element individually
      gene_names <- filtered_nodes_ann_mult()[selected_rows, "Gene name"]

      # Use sapply to apply gsub to each gene name in the vector
      values$selected_genes_mult <- sapply(gene_names, function(x) gsub('.*>(.*)<.*', '\\1', x))
    } else {
      values$selected_genes_mult <- NULL
    }
    print(paste("Selected genes: ", paste(values$selected_genes_mult, collapse = ", ")))
  })

  # Create a reactive expression to hold the network data
  network_data_mult <- reactiveVal(NULL)

  # Observe button click to update network data
    observeEvent(input$plot_button_multiple, {
    # Ensure that inputs are available and valid
    req(values$selected_genes_mult, input$plot_button_multiple, input$List_gen)

    filtered_edges_data <- filtered_edges_mult()

    selected_genes <- values$selected_genes_mult
    # print(selected_genes)
    # selected_genes <-c("Potra2n8c17315", "Potra2n10c20839")

    if (length(selected_genes) == 0) {
      filtered_edges_network <- filtered_edges_data
    } else {
      filtered_edges_network <- edges %>%
        filter(fromNode %in% selected_genes & toNode %in% selected_genes)
    }

    if (nrow(filtered_edges_network) == 0) {
      network_data_mult(NULL)  # Set to NULL if no data
      return(NULL)
    }


    filtered_nodes <- data.frame(id = unique(c(filtered_edges_data$fromNode, filtered_edges_data$toNode)))
    filtered_nodes$nodeID <- match(filtered_nodes$id, unique(filtered_nodes$id)) - 1

    filtered_nodes_network <- nodes %>%
      dplyr::select(-c(Degree, Module)) %>%
      filter(nodeName %in% unique(c(filtered_edges_network$fromNode, filtered_edges_network$toNode))) %>%
      left_join(filtered_nodes, by = c("nodeName" = "id")) %>%
      left_join(subannot %>% select(`Gene name`, Centrality, Module), by = c("nodeName" = "Gene name")) %>%
      separate(Module, c('Module_color', 'Module'), sep = "/")

    filtered_edges_network$fromNodeID <- match(filtered_edges_network$fromNode, filtered_nodes_network$nodeName) - 1
    filtered_edges_network$toNodeID <- match(filtered_edges_network$toNode, filtered_nodes_network$nodeName) - 1

    network_data_mult(list(edges = filtered_edges_network, nodes = filtered_nodes_network))
  })

  # Render network plot based on the reactive network_data
  output$network_plot_multiple <- renderForceNetwork({
    network_info <- network_data_mult()
    req(network_info)  # Ensure network data is available

    network_plot <- forceNetwork(
      Links = network_info$edges,
      Nodes = network_info$nodes,
      Source = "fromNodeID",
      Target = "toNodeID",
      NodeID = "nodeName",
      Group = "Module",
      opacity = 0.8,
      zoom = TRUE,
      linkColour = "grey",
      linkWidth = 2,
      fontSize = 20,
      linkDistance = 400,
      colourScale = networkD3::JS('d3.scaleOrdinal().domain(["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46"]).range(["#40E0D0", "#0000FF", "#A52A2A", "#FFFF00", "#00FF00", "#FF0000", "#000000", "#FFC0CB", "#FF00FF", "#A020F0", "#ADFF2F", "#D2B48C", "#FA8072", "#00FFFF", "#191970", "#E0FFFF", "#999999", "#90EE90", "#FFFFE0", "#4169E1", "#8B0000", "#006400", "#00CED1", "#A9A9A9", "#FFA500", "#FF8C00", "#EEE9E9", "#87CEEB", "#8B4513", "#AFEEEE", "#4682B4", "#EE82EE", "#556B2F", "#8B008B", "#CD6839", "#9ACD32", "#6CA6CD", "#FFBBFF", "#8B2500", "#8968CD", "#CAE1FF", "#E0FFFF", "#FFFFF0", "#FFFAF0", "#EE7600", "#8B2323"])'),
      Nodesize = "Centrality",
      legend = FALSE
    )
    network_plot <- htmlwidgets::onRender(
      network_plot,
      '
  function(el, x) {
    d3.selectAll(".node circle")
      .style("stroke", "grey")
      .style("stroke-width", "2px")

    d3.selectAll(".node text")
      .style("fill", "#B3BAB9");
  }
  '
    )

    network_plot
  })

  # Render legend plot
  output$network_legend_multiple <- renderPlot({
    req(values$selected_genes_mult, input$plot_button_multiple)
    suppressWarnings({
    grid.newpage()
    grid.roundrect(
      x = 0.6, y = 0.5,
      width = unit(6, "cm"), height = unit(14.5, "cm"),  # Adjust width and height of the legend background
      r = unit(0.5, "cm"),  # Set the corner radius, tweak as needed
      gp = gpar(fill = "#C1D6E0", col = NA)  # Fill color for the legend background
    )
    pushViewport(viewport(x = 0.6, y = 0.5, just = "center"))  # Move the legend to the right
    grid.draw(legend)
    })
  })
 }  



# Run the application
shinyApp(ui = ui, server = server)