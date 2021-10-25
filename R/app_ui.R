#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinydashboard
#' @import shinyjs
#' @import shinyBS
#' @import shinycssloaders
#' @import shinyFiles
#' @import DT
#' @noRd
app_ui <- function(request) {
  fluidPage(
    
    dashboardPage(
      dashboardHeader(title = "scWizard"),
      dashboardSidebar(
        sidebarMenu(
          id = "tabs",
          menuItem("User Guide", tabName = "introTab", icon = icon("info-circle")),
          menuItem("Input Data", tabName = "datainput", icon = icon("upload")),
          menuItem("Batch Processing", tabName = "batchTab", icon = icon("th")),
          menuItem("Cell Annotation", tabName = "annoTab", icon = icon("th")),
          menuItem("GSVA", tabName = "gsvaTab", icon = icon("th")),
          menuItem("Find Markers", tabName = "findMarkersTab", icon = icon("th")),
          menuItem("Monocle", tabName = "monocleTab", icon = icon("th")),
          menuItem("TF-SCENIC", tabName = "tfScenicTab", icon = icon("th")),
          menuItem("Correlation Analysis", tabName = "corTab", icon = icon("th")),
          menuItem("CellphoneDB", tabName = "cellphoneDBTab", icon = icon("th")),
          menuItem("Plot", tabName = "plotTab", icon = icon("bar-chart"))
          
          
        )
      ),
      dashboardBody(
        shinyjs::useShinyjs(),
        extendShinyjs(script = "app/www/custom.js", functions = c("winprint")),
        tags$head(
          tags$style(HTML(
            " .shiny-output-error-validation {color: darkred; } "
          )),
          tags$link(rel = "stylesheet", type = "text/css", href = "app/www/custom.css"),
          tags$link(rel = "stylesheet", type = "text/css", href = "app/www/buttons.css")
        ),
        tabItems(
          tabItem(tabName = "introTab",
                  fluidRow(
                    
                    box(title = "User Guide", width = 11, solidHeader = T, status = "primary",
                        column(12,
                               includeMarkdown(system.file("app/www/intro.Rmd", package='scWizard'))
                        )
                    )
                  )
          ),
          
          tabItem(tabName = "datainput",
                  hr(),
                  fluidRow(column(3,
                                  box(title = "Upload Data", solidHeader = T, status = "primary", width = 12, collapsible = T,id = "uploadbox",
                                      
                                      radioButtons('data_file_type','Use example file or upload your own data',
                                                   c(
                                                     'Data-rds'="data_rds",
                                                     'Data-example'="data_example"
                                                   ),selected = "data_example"),
                                      
                                      conditionalPanel(condition="input.data_file_type=='data_rds'",
                                                       p("please input .rds file"),
                                                       fileInput('datafile', 'Choose File(s) Containing Data', multiple = TRUE)
                                      )
                                  )
                                  
                  ),#column
                  
                  column(9,
                         bsCollapse(id="input_collapse_panel",open="data_panel",multiple = FALSE,
                                    bsCollapsePanel(title="Data Contents Table:",value="data_panel",
                                                    p("Note: if there are more than 20 columns, only the first 20 will show here"),
                                                    textOutput("inputInfo"),
                                                    withSpinner(dataTableOutput('countdataDT'))
                                    )
                         )#bscollapse
                  )#column
                  )#fluidrow
          ),#tabpanel
          tabItem(tabName = "batchTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Batch Processing")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("Different factors such as the collection time of sample data set, collection organization and sequencing platform may automatically form different batches, thus affecting the real data. Therefore, the batch effect should be checked and removed before subsequent analysis. Otherwise, all subsequent analysis results will not be used."),
                                           p("Introduction parameters(Seurat):"),
                                           p("-dims: Which dimensions to use from the CCA to specify the neighbor search space"),
                                           p("-normalization.method: Name of normalization method used"),
                                           p("-reduction: Dimensional reduction to perform when finding anchors"),
                                           p("-k.weight: Number of neighbors to consider when weighting anchors"),
                                           p("Introduction parameters(Harmony):"),
                                           p("-group.by.vars: Which variable(s) to remove"),
                                           p("-max.iter.harmony: Maximum number of rounds to run Harmony"),
                                           p("-dims.use: Which PCA dimensions to use for Harmony"),
                                           p("-lambda: Ridge regression penalty parameter"),
                                           p("-theta: Diversity clustering penalty parameter"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Seurat CCA",
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of Seurat parameter:"),
                                                                         column(4, numericInput('dims', 'Dims',value = 30)),
                                                                         column(4, numericInput('kweight', 'k.weight',value = 100)),
                                                                         column(4,selectInput("normalizationmethod", "normalization",
                                                                                              choices = c("LogNormalize", "SCT")
                                                                                              , selected = "LogNormalize")),
                                                                         column(4,selectInput("reductionmethod", "reduction",
                                                                                              choices = c("cca", "rpca")
                                                                                              , selected = "cca")),
                                                                         column(4, uiOutput("myselectboxbatch")),
                                                                         
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startCCA","start CCA",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.CCAAvailable",
                                                                                        downloadButton('downloadCCAPlot','Save Results as Plot File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(plotOutput(outputId = "CCAPlot", height = 350))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       ),
                                                       tabPanel("Harmony",
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of Harmony parameter:"),
                                                                         column(4,numericInput("maxiterharmony", "max.iter.harmony", value = 10)),
                                                                         column(4,numericInput("lambda", "Lambda", value = 2)),
                                                                         column(4,numericInput("theta", "Theta", value = 1)),
                                                                         column(4,numericInput("dimsuse", "dims.use", value = 2000)),
                                                                         column(4,uiOutput("myselectboxbatch2")),
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startHarmony","start Harmony",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.HarmonyAvailable",
                                                                                        downloadButton('downloadHarmonyPlot','Save Results as Plot File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(plotOutput(outputId = "HarmonyPlot", height = 350))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          tabItem(tabName = "annoTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Cell Annotation")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("Automatically assign a cell type to each cell and return the result."),
                                           p("Introduction parameters(cell annotion):"),
                                           p("-learning rate: Adjusting the convergence speed of neural network"),
                                           p("-Number of hidden layer nodes: A reasonable number of nodes can effectively avoid over fitting"),
                                           p("-Regularization rate: Avoid over fitting"),
                                           p("-PCA.k: Dimension to which PCA is reduced"),
                                           p("Introduction parameters(Cell classification):"),
                                           p("-resolution: Use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities"),
                                           p("Introduction parameters(Cell subcluster annotation):"),
                                           p("Parameters are the same as cell annotation"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Cell Annotion",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of BP3 parameter:"),
                                                                         column(4,numericInput("learning", "learning rate", value = 0.001)),
                                                                         column(4,numericInput("layer1", "Number of hidden layer1 nodes", value = 100)),
                                                                         column(4,numericInput("regularization", "Regularization rate", value = 0.05)),
                                                                         column(4,numericInput("PCAk", "PCA.k", value = 200)),
                                                                         div(style = "clear:both;"),
                                                                         actionButton("installPython","install Python",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%"),
                                                                         actionButton("startAnnotion","start Annotion",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.AnnotionAvailable",
                                                                                        downloadButton('downloadAnnotionPlot','Save Results as Plot File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(plotOutput(outputId = "AnnotionPlot", height = 700))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       ),
                                                       tabPanel("Cell Classification",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of Classification parameter:"),
                                                                         column(4,numericInput("resolution", "Resolution", value = 0.5)),
                                                                         column(4, uiOutput("myselectboxanno1")),
                                                                         
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startClassification","start Classification",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.ClassificationAvailable",
                                                                                        downloadButton('downloadClassificationPlot','Save Results as Plot File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(plotOutput(outputId = "ClassificationPlot", height = 700))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       ),
                                                       tabPanel("Cell subcluster annotation",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of BP5 parameter:"),
                                                                         column(4,numericInput("learning2", "learning rate", value = 0.001)),
                                                                         column(4,numericInput("layer2_1", "Number of hidden layer1 nodes", value = 100)),
                                                                         column(4,numericInput("layer2_2", "Number of hidden layer2 nodes", value = 50)),
                                                                         column(4,numericInput("layer2_3", "Number of hidden layer3 nodes", value = 25)),
                                                                         column(4,numericInput("regularization2", "Regularization rate", value = 0.05)),
                                                                         column(4,numericInput("PCAk2", "PCA.k", value = 2000)),
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startSubannotion","start Subannotion",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.SubannotionAvailable",
                                                                                        downloadButton('downloadSubannotionPlot','Save Results as Plot File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(plotOutput(outputId = "SubannotionPlot", height = 700))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          tabItem(tabName = "gsvaTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Gene Set Variation Analysis")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("GSVA is a non parametric and unsupervised analysis method, which is mainly used to evaluate the gene set enrichment results of transcriptome. In order to evaluate whether different metabolic pathways are enriched in different products, the expression matrix of genes among different products is transformed into the expression matrix of gene sets among different samples."),
                                           p("Introduction parameters:"),
                                           p("-kcdf: Choose method to estimate"),
                                           p("-gmtfile: Geneset file"),
                                           p("-celltype: Choose the cell type you want to calculate"),
                                           p("-mix.diff: Offers two approaches to calculate the enrichment statistic"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("GSVA",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of gsva parameter:"),
                                                                         
                                                                         #column(4,numericInput("threshAll", "Logfc Thresh", value = 0.25)),
                                                                         #column(4,numericInput("minsz", "min % (min.sz)", value = 1)),
                                                                         column(4,selectInput("kcdfmethod", "kcdf",
                                                                                              choices = c("Gaussian", "Poisson", "none")
                                                                                              , selected = "Gaussian")),
                                                                         column(4, fileInput('gmtfile', 'Choose gmtFile Containing Data', multiple = TRUE)),
                                                                         column(4, uiOutput("myselectbox2")),
                                                                         column(4, checkboxInput("mxdiff","mx.diff"), value = TRUE),
                                                                         
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startGSVA","start GSVA",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.GSVAAvailable",
                                                                                        downloadButton('downloadGSVACSV','Save Results as CSV File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(dataTableOutput('GSVA'))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          
          tabItem(tabName = "findMarkersTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Finding differentially expressed genes (cluster biomarkers)")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells."),
                                           p("Introduction parameters:"),
                                           p("-Logfc Thresh: Limit testing to genes"),
                                           p("-min.pct: Only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations"),
                                           p("-Test.use: Denotes which test to use"),
                                           p("-min.diff.pct(-Inf~0): only test genes that show a minimum difference in the fraction of detection between the two groups"),
                                           p("-max.cells.per.ident(0~Inf): Down sample each identity class to a max number"),
                                           p("-only.pos: Only return positive markers (FALSE by default)"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Find markers by cluster",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Select clusters to find markers:"),
                                                                         
                                                                         column(4,numericInput("threshAll", "Logfc Thresh", value = 0.25)),
                                                                         column(4,numericInput("minPct", "min.pct", value = 0.1)),
                                                                         column(4,selectInput("testuse", "Test.use",
                                                                                              choices = c("wilcox","bimod","roc","t","tobit","poisson","negbinom","MAST","DESeq2")
                                                                                              , selected = "wilcox")),
                                                                         column(4,numericInput("mindiffpct", "min.diff.pct", value = -1/0)),
                                                                         column(4,numericInput("maxcellsperident", "max.cells.per.ident", value = 1/0)),
                                                                         
                                                                         column(4, checkboxInput("onlypos","only.pos"), value = FALSE),
                                                                         
                                                                         div(style = "clear:both;"),
                                                                         actionButton("findClusterMarkers","Find Cluster Markers",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.clusterMarkersAvailable",
                                                                                        downloadButton('downloadClusterMarkersCSV','Save Results as CSV File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(dataTableOutput('clusterMarkers'))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          
          tabItem(tabName = "monocleTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Pseudotime analysis by monocle")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("Single cell pseudo time analysis software, monocle, is an installation package based on R language. Its function is based on the expression matrix of single cell transcriptome, which can simulate the biological process of cell population by placing cells on different branches of development trajectory through unsupervised learning. That is, we often say quasi time series analysis, also known as cell trajectory analysis."),
                                           p("Introduction parameters:"),
                                           p("-lowerDetectionLimit: The minimum expression level that consistitutes true expression"),
                                           p("-max_components: The dimensionality of the reduced space"),
                                           p("-mean_expression: Select data by expression quantity"),
                                           p("-method: The algorithm to use for dimensionality reduction"),
                                           p("-celltype: Choose the cell type you want to calculate"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Pseudotime analysis",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of Monocle parameter:"),
                                                                         
                                                                         column(4,numericInput("lowerdetectionlimit", "lowerDetectionLimit", value = 0.5)),
                                                                         column(4,numericInput("maxcomponents", "max_components", value = 2)),
                                                                         column(4,numericInput("meanexpression", "mean_expression", value = 0.1)),
                                                                         column(4,selectInput("rmethod", "method",
                                                                                              choices = c("DDRTree",
                                                                                                          "ICA", "tSNE", "SimplePPT", "L1-graph", "SGL-tree"), selected = "DDRTree")),
                                                                         column(4,uiOutput("myselectbox3")),
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startMonocle","start Monocle",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.MonocleAvailable",
                                                                                        downloadButton('downloadMonocleRDS','Save Results as RDS File', class = "btn btn-primary"),
                                                                                        p("Calculation completed, please download data...")
                                                                       ),
                                                                       br()
                                                                       #withSpinner(dataTableOutput('Monocle'))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          
          tabItem(tabName = "tfScenicTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Single-Cell Regulatory Network Inference And Clustering")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("Scenic analysis is to study the transcription factors (TFs) in the data of scrna SEQ, and finally screen the TFs with significant regulatory intensity and core role. The results are usually displayed in the form of heat map. Especially in oncology, Scenic analysis can help to find the key 'driver' related to the occurrence and development of tumor, so as to lay the foundation for exploring its pathogenesis."),
                                           p("introduction parameters"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("TF-Scenic",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of Scenic parameter:"),
                                                                         # column(4,selectInput("clusterNum", "Cluster Num (ident.1)",
                                                                         #                      choices = c(1,2,3,4), selected = 1)),
                                                                         column(4,numericInput("minCountsPerGene1", "minCountsPerGene_1", value = 3)),
                                                                         column(4,numericInput("minCountsPerGene2", "minCountsPerGene_2", value = 0.1)),
                                                                         column(4,numericInput("minSamples", "choose minSamples", value = 0.1)),
                                                                         column(4,selectInput("org", "choose org",
                                                                                              choices = c("mgi", "hgnc", "dmel")
                                                                                              , selected = "hgnc")),
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startCoexpression","First start Coexpression",class = "button button-block button-pill button-primary button-large", style = "width: 100%"),
                                                                         br(),
                                                                         actionButton("startGRN","Second start GRN",class = "button button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       # conditionalPanel("output.TFscenicAvailable",
                                                                       #                  actionButton("startGRN","Second start GRN",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       # ),
                                                                       br(),
                                                                       
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          
          tabItem(tabName = "corTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Correlation Analysis")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("Correlation analysis refers to the analysis of two or more variable elements with correlation, so as to measure the degree of correlation between two variable factors. Correlation analysis can only be carried out when there is a certain connection or probability between the elements of correlation. Relevance is not equal to causality, nor is it simple personalization. The scope and field of relevance almost cover all aspects we have seen. The definition of relevance in different disciplines is also very different."),
                                           p("Introduction parameters:"),
                                           p("-use: Giving a method for computing covariances in the presence of missing values"),
                                           p("-method: Indicating which correlation coefficient is to be computed"),
                                           p("-celltype: Choose the cell type you want to calculate"),
                                           p("-gene: Choose the genes you want to calculate"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("Correlation analysis",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of Cor parameter:"),
                                                                         
                                                                         column(4,selectInput("use", "choose use",
                                                                                              choices = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")
                                                                                              , selected = "everything")),
                                                                         column(4,selectInput("method", "choose method",
                                                                                              choices = c("pearson", "kendall", "spearman")
                                                                                              , selected = "pearson")),
                                                                         
                                                                         column(4,uiOutput("myselectbox5")),
                                                                         column(4,uiOutput("myselectgroupbox5")),
                                                                         
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startCor","start Cor",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       conditionalPanel("output.CorAvailable",
                                                                                        downloadButton('downloadCorCSV','Save Results as CSV File', class = "btn btn-primary")
                                                                       ),
                                                                       br(),
                                                                       withSpinner(dataTableOutput('Cor'))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          
          tabItem(tabName = "cellphoneDBTab",
                  
                  fluidRow(
                    
                    column(12,
                           h3(strong("Intercellular communication based on ligand receptor pairs")),
                           hr(),
                           
                           column(12,
                                  tags$div(class = "BoxArea2",
                                           p("CellPhoneDB is a publicly available repository of curated receptors, ligands and their interactions. Subunit architecture is included for both ligands and receptors, representing heteromeric complexes accurately. This is crucial, as cell-cell communication relies on multi-subunit protein complexes that go beyond the binary representation used in most databases and studies."),
                                           p("Introduction parameters:"),
                                           p("-countsdata: Type of gene identfiers in counts data"),
                                           p("-cellphonedbin: Select inputfile location"),
                                           p("-cellphonedbin: Select outputfile location"),
                                           
                                           tabsetPanel(type = "tabs",
                                                       tabPanel("CellphoneDB",
                                                                
                                                                column(12,
                                                                       wellPanel(
                                                                         h4("Set the value of cellphonedb parameter:"),
                                                                         
                                                                         column(4,selectInput("countsdata", "counts-data",
                                                                                              choices = c("ensembl","gene_name", "hgnc_symbol")
                                                                                              , selected = "gene_name")),
                                                                         
                                                                         column(4,shinyDirButton("cellphonedbin", "Choose a cellphonedbin_folder" ,
                                                                                                 title = "",
                                                                                                 buttonType = "default", class = NULL,
                                                                                                 icon = icon("folder", lib = "font-awesome"), multiple = TRUE)),
                                                                         column(4,shinyDirButton("cellphonedbout", "Choose a cellphonedbout_folder" ,
                                                                                                 title = "",
                                                                                                 buttonType = "default", class = NULL,
                                                                                                 icon = icon("folder", lib = "font-awesome"), multiple = TRUE)),
                                                                         
                                                                         div(style = "clear:both;"),
                                                                         actionButton("startcellphonedb","Start CellphoneDB",class = "button button-3d button-block button-pill button-primary button-large", style = "width: 100%")
                                                                       ),
                                                                       # conditionalPanel("output.cellphonedbAvailable",
                                                                       #                  downloadButton('downloadcellphonedb','Save Results as CSV File', class = "btn btn-primary")
                                                                       # ),
                                                                       br(),
                                                                       #withSpinner(dataTableOutput('cellphonedb'))
                                                                ),
                                                                tags$div(class = "clearBoth")
                                                       )
                                                       
                                           )
                                  )
                                  
                           ),
                           tags$div(class = "clearBoth")
                           
                    )
                  )
          ),
          
          tabItem(tabName = "plotTab",
                  fluidRow(
                    column(12,h3(strong("Plot Output:")),
                           column(3,
                                  wellPanel(
                                    h4("Choose graph to plot:"),
                                    selectInput("methodsUsed", "Methods used",
                                                choices = c("GSVA", "Monocle", "Scenic", "Cor", "CellphoneDB")
                                                , selected = "GSVA"),
                                    uiOutput("sel_box"),
                                    fileInput('plotfile1', 'Choose plotFile1 Containing Data', multiple = TRUE),
                                    fileInput('plotfile2', 'Choose plotFile2 Containing Data', multiple = TRUE)
                                    
                                  ),
                                  actionButton("startPlot","Start Plot",class = "btn btn-success btn-lg", style = "width: 100%"),
                                  
                                  conditionalPanel("output.saveAvailable",
                                                   downloadButton('downloadplot','Save Plot', class = "btn btn-success btn-lg", style = "width: 100%")
                                  )
                           ),
                           column(9,
                                  column(12,
                                         withSpinner(plotOutput(outputId = "plot", height = 700))
                                  )
                           ),
                           
                    )
                  )
          )
          
        )
        
      )
      
    ),
    tags$footer(
      wellPanel(
        HTML(
          '
      <p align="center" width="4">Developed and maintained by: School of bioengineering and science, South China University of Technology</p>
      <p align="center" width="4">Copyright (C) 2021  </p>
      '
        )
      ),
      tags$script(src = "imgModal.js")
    )
    
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'scWizard'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

