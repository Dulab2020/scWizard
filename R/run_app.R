#' Run the Shiny Application
#'
#' @param ... A series of options to be used inside the app.
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(options = list()) {
  shiny::shinyApp(ui = app_ui,
                  server = app_server,
                  options = options) 
}

#golem::use_recommended_deps(recommended = c("shinydashboard","shinyjs","shinyBS","shinycssloaders","shinyFiles","DT","Seurat","dplyr","Matrix","V8","sodium","GSVA","GSVAdata","future","ggplot2","patchwork","monocle","SCENIC","reticulate","doParallel","pheatmap"))
