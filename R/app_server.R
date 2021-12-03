#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @import shinyjs
#' @import shinyBS
#' @import shinycssloaders
#' @import shinyFiles
#' @import shinydashboard
#' @import DT
#' @import Seurat
#' @import dplyr
#' @import Matrix
#' @import V8
#' @import harmony
#' @import sodium
#' @import GSVA
#' @import GSVAdata
#' @import future
#' @import ggplot2
#' @import patchwork
#' @import monocle
#' @import SCENIC
#' @import GSEABase
#' @import infercnv
#' @import reticulate
#' @import doParallel
#' @import BiocGenerics
#' @import pheatmap
#' @noRd
app_server <- function( input, output, session ) {
  options(shiny.maxRequestSize=-1) # Remove limit of upload   
  options(shiny.deprecation.messages=F)   
  options(warn =-1)
  
  ### input data
  observe({
    inputDataReactive()
  })
  inputDataReactive <- reactive({
    print("inputting data")
    if (input$data_file_type == "data_example") {
      inFile=system.file("app/www/example/data_300cell.RDS", package='scWizard')
    }
    else{
      inFile = input$datafile
    }
    if (!is.null(inFile)) {
      if(input$data_file_type == "data_rds"){
        seqdata <- readRDS(inFile$datapath)
      }
      else{
        seqdata <- readRDS(inFile)
      }
      print('uploaded seqdata')
      shiny::validate(need(ncol(seqdata)>1,
                           message="File appears to be one column. Check that it is a comma or tab delimited file."))
      return(list('data'=seqdata))
    }
    return(NULL)
  })
  output$countdataDT <- renderDataTable({
    tmp <- inputDataReactive()
    if(!is.null(tmp))
    {
      df = tmp$data@assays$RNA@counts
      if(ncol(df) > 20)
        return(as.matrix(df[,1:20]))
      return(as.matrix(df))
    }
    return(matrix(nrow=0,ncol=0))
  },
  options = list(scrollX = TRUE))
  output$inputInfo <- renderText({
    tmp <- inputDataReactive()$data
    if(!is.null(tmp))
    {
      df = tmp@assays$RNA@counts
      outStr = paste0(
        paste("dense size: ", object.size(x = as.matrix(x = df))),
        '\n',
        paste("sparse size: ", object.size(x = df)))
    }
    else{
      paste("dense size: NULL")
    }
  })
  
  ### Quality Control
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  observe({
    QCReactive()
  })
  QCReactive <- eventReactive(input$startQC, {
    withProgress(message = "Processing,please wait",{
      data_rds = inputDataReactive()$data
      tryCatch({
        data_rds[["percent.mt"]] <- PercentageFeatureSet(data_rds, pattern = "^MT-")
        data_rds <- subset(data_rds, subset = nFeature_RNA>input$featurelow & nFeature_RNA<input$featurehigh & nCount_RNA>input$countlow & nCount_RNA<input$counthigh & percent.mt<input$percent)
        p1 = FeatureScatter(data_rds, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'orig.ident')
        p2 = FeatureScatter(data_rds, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'orig.ident')
        shiny::setProgress(value = 0.8, detail = "Done.")
        res_plot = p1+p2
        return(list("plot" = res_plot,"data" = data_rds))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
      return(NULL)
    })
  })
  observe({
    viewReactive()
  })
  viewReactive <- eventReactive(input$viewQC, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        if(input$startQC > 0)
          data_rds = QCReactive()$data
        else
          data_rds = inputDataReactive()$data
        data_rds[["percent.mt"]] <- PercentageFeatureSet(data_rds, pattern = "^MT-")
        p = VlnPlot(data_rds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident')
        shiny::setProgress(value = 0.8, detail = "Done.")
        res_plot = p
        return(list("plot" = res_plot))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
      return(NULL)
    })
  })
  output$QCPlot <- renderPlot({
    tmp <- QCReactive()
    if(!is.null(tmp)){
      tmp$plot
    }
  })
  output$viewPlot <- renderPlot({
    tmp <- viewReactive()
    if(!is.null(tmp)){
      tmp$plot
    }
  })
   
  ### GSVA
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  #choose celltype
  output$myselectbox2 <-
    renderUI({
      data_rds = inputDataReactive()$data
      label=names(summary(data_rds@active.ident))
      selectInput("celltype", "choose celltype",
                  choices =c('all' ,label), selected = 'all')
    })
  observe({
    GSVAReactive()
  })
  GSVAReactive <- eventReactive(input$startGSVA, {
    withProgress(message = "Processing,please wait",{
      data_rds = inputDataReactive()$data
      tryCatch({  
        if(input$celltype=='all')
        {
          new_data_rds = data_rds
        }
        else
        {
          new_data_rds = subset(data_rds, idents = input$celltype)
        }
        data_counts = as.matrix(new_data_rds@assays$RNA@data)
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        gmt<-getGmt(input$gmtfile$datapath)
        res_gsva <- GSVA::gsva(data_counts, gmt, kcdf=input$kcdfmethod, mx.diff=input$mxdiff)
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = t(res_gsva)))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NA)
      })
    })
  })
  output$GSVA <- renderDataTable({
    tmp <- GSVAReactive()
    if(!is.null(tmp)){
      tmp$data
    }
  })
  output$GSVAAvailable <-
    reactive({
      return(!is.null(GSVAReactive()$data))
    })
  outputOptions(output, 'GSVAAvailable', suspendWhenHidden=FALSE)
  output$downloadGSVACSV <- downloadHandler(
    filename = function()  {paste0(GSVAReactive()$data,".csv")},
    content = function(file) {
      write.csv(GSVAReactive()$data, file, row.names=TRUE)}
  )
  
  ### batch processing
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  # return selectbox
  output$myselectboxbatch <-
    renderUI({
      if(input$startQC > 0)
        data_rds = QCReactive()$data
      else
        data_rds = inputDataReactive()$data
      label=colnames(data_rds@meta.data)
      selectInput("batch", "choose batch",
                  choices =c(' ' ,label), selected = ' ')
    })
  # return selectbox
  output$myselectboxbatch2 <-
    renderUI({
      if(input$startQC > 0)
        data_rds = QCReactive()$data
      else
        data_rds = inputDataReactive()$data
      label=colnames(data_rds@meta.data)
      selectInput("batch2", "choose batch",
                  choices =c(' ' ,label), selected = ' ')
    })
  # CCA
  observe({
    CCAReactive()
  })
  CCAReactive <- eventReactive(input$startCCA, {
    withProgress(message = "Processing,please wait",{
      data_rds = inputDataReactive()$data
      tryCatch({
        p1 = DimPlot(data_rds, reduction = "tsne",  pt.size = .1,group.by = input$batch)
        ifnb.list <- SplitObject(data_rds, split.by = input$batch)
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
          x <- NormalizeData(x)
          x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
        })
        features <- SelectIntegrationFeatures(object.list = ifnb.list)
        immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)
        immune.combined <- IntegrateData(anchorset = immune.anchors)
        DefaultAssay(immune.combined) <- "integrated"
        immune.combined <- ScaleData(immune.combined, verbose = FALSE)
        immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
        immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:30)
        immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
        immune.combined <- FindClusters(immune.combined, resolution = 0.5)
        p2 = DimPlot(immune.combined, reduction = "tsne",  pt.size = .1,group.by = input$batch)
        shiny::setProgress(value = 0.8, detail = "Done.")
        res_plot = p1+p2
        return(list("plot" = res_plot,"data" = immune.combined))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NA)
      })
    })
  })
  output$CCAPlot <- renderPlot({
    tmp <- CCAReactive()
    if(!is.null(tmp)){
      tmp$plot
    }
  })
  output$CCAAvailable <-
    reactive({
      return(!is.null(CCAReactive()$plot))
    })
  outputOptions(output, 'CCAAvailable', suspendWhenHidden=FALSE)
  output$downloadCCAPlot <- downloadHandler(
    filename = function()  {'batchplot.pdf'},
    content = function(file) {
      ggsave(file,CCAReactive()$plot)
    }
  )
  # harmony
  observe({
    HarmonyReactive()
  })
  HarmonyReactive <- eventReactive(input$startHarmony, {
    withProgress(message = "Processing,please wait",{
      data_rds = inputDataReactive()$data
      tryCatch({
        p1 = DimPlot(data_rds, reduction = "tsne",  pt.size = .1,group.by = input$batch2)
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        data_rds <- data_rds %>% RunHarmony(input$batch2, plot_convergence = TRUE)
        data_rds <- RunTSNE(data_rds,reduction = "harmony", dims = 1:30)
        p2 = DimPlot(data_rds, reduction = "tsne",  pt.size = .1,group.by = input$batch2)
        shiny::setProgress(value = 0.8, detail = "Done.")
        res_plot = p1+p2
        return(list("plot" = res_plot, "data" = data_rds))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NA)
      })
    })
  })
  output$HarmonyPlot <- renderPlot({
    tmp <- HarmonyReactive()
    if(!is.null(tmp)){
      tmp$plot
    }
  })
  output$HarmonyAvailable <-
    reactive({
      return(!is.null(HarmonyReactive()$plot))
    })
  outputOptions(output, 'HarmonyAvailable', suspendWhenHidden=FALSE)
  output$downloadHarmonyPlot <- downloadHandler(
    filename = function()  {paste0(HarmonyReactive()$data,".pdf")},
    content = function(file) {
      pdf(file)
      print(HarmonyReactive()$plot)
      dev.off()}
  )
  
  ### cell annotion
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  observe({
    installpythonReactive()
  })
  installpythonReactive <- eventReactive(input$installPython, {
    withProgress(message = "Processing,please wait",{
      shiny::setProgress(value = 0.4, detail = "insatlling ...")
      tryCatch({
        if(!file.exists(system.file("miniconda", package='scWizard')))
        {
          dir.create(paste0(system.file("", package='scWizard'),'/miniconda'))
          conda_path = system.file("miniconda", package='scWizard')
          install_miniconda(path = conda_path)
        }
        if(file.exists(system.file("miniconda/envs/r-reticulate", package='scWizard')))
        {
          cellphonedb_path = system.file("app/www/CellPhoneDB-2.1.4.tar.gz", package='scWizard')
          envs = system.file("miniconda/envs/r-reticulate", package='scWizard')
          conda_install(envname = envs, packages = 'rpy2==3.4.2', pip = T)
          conda_install(envname = envs, packages = cellphonedb_path, pip = T)
          conda_install(envname = envs, packages = 'scikit-learn==0.22', pip = T)
          conda_install(envname = envs, packages = 'tensorflow-gpu==2.4.1', pip = T)
          conda_install(envname = envs, packages = 'tables', pip = T)
        }
        reticulate::use_python(system.file("miniconda/envs/r-reticulate", package='scWizard'), required = F)
        py_config()
        shiny::setProgress(value = 0.8, detail = "Done.")
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NA)
      })
    })
  })
  # anno celltype
  observe({
    AnnotionReactive()
  })
  AnnotionReactive <- eventReactive(input$startAnnotion, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        source_python(system.file("app/www/python/BP3_new.py", package='scWizard'))
        if(input$startCCA > 0)
          data_rds = CCAReactive()$data
        else if(input$startHarmony > 0)
          data_rds = HarmonyReactive()$data
        else
          data_rds = inputDataReactive()$data
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        X_total_path=system.file('app/www/python/trainset/trainx_all.h5', package='scWizard')
        
        Y_total=read.csv(system.file('app/www/python/trainset/trainy_all.csv', package='scWizard'))
        Y_total$X=NULL
        colnames(Y_total)=c('celltype')
        num_classes = length(unique(Y_total$celltype))
        X_verify=t(as.data.frame(data_rds@assays$RNA@data))
        X_verify=as.data.frame(X_verify)
        res_celltype=get_BP3_res(X_total_path, Y_total$celltype, X_verify, num_classes, input$PCAk, input$layer1, input$regularization, input$learning)
        data_rds@meta.data$pred_cell = res_celltype
        res_plot = DimPlot(data_rds,reduction = "tsne",pt.size = .1,group.by = 'pred_cell')
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = res_plot, "data_rds" = data_rds))
      },
      error=function(cond) {
        message("Here's the original error.")
        #message(cond)
        return(NA)
      })
    })
  })
  output$AnnotionPlot <- renderPlot({
    tmp <- AnnotionReactive()
    if(!is.null(tmp)){
      tmp$data
    }
  })
  output$AnnotionAvailable <-
    reactive({
      return(!is.null(AnnotionReactive()$data))
    })
  outputOptions(output, 'AnnotionAvailable', suspendWhenHidden=FALSE)
  output$downloadAnnotionPlot <- downloadHandler(
    filename = function()  {'plot_pred_cell.pdf'},
    content = function(file) {
      pdf(file)
      print(AnnotionReactive()$data)
      dev.off()}
  )
  
  # return selectbox
  output$myselectboxanno1 <-
    renderUI({
      data_rds = AnnotionReactive()$data_rds
      label = unique(data_rds@meta.data$pred_cell)
      selectInput("celltype1", "choose celltype",
                  choices =c(NULL ,label), selected = NULL)
    })
  observe({
    ClassificationReactive()
  })
  ClassificationReactive <- eventReactive(input$startClassification, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        #data_rds = inputDataReactive()$data
        if(input$startAnnotion > 0)
          data_rds = AnnotionReactive()$data_rds
        else
          data_rds = inputDataReactive()$data
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        tmp_data = subset(data_rds, subset = pred_cell==input$celltype1)
        tmp_counts = tmp_data@assays$RNA@counts
        tmp_meta = tmp_data@meta.data
        tmp_data = CreateSeuratObject(counts = tmp_counts, project = input$celltype1, min.cells = 3, min.features = 200)
        tmp_data@meta.data = tmp_meta
        tmp_data <- NormalizeData(tmp_data, normalization.method = "LogNormalize", scale.factor = 10000)
        tmp_data <- FindVariableFeatures(tmp_data, selection.method = "vst", nfeatures = 2000)
        all.genes <- rownames(tmp_data)
        tmp_data <- ScaleData(tmp_data, features = all.genes)
        tmp_data <- RunPCA(tmp_data, features = VariableFeatures(object = tmp_data), npcs=30)
        tmp_data <- FindNeighbors(tmp_data, dims = 1:30)
        tmp_data <- FindClusters(tmp_data, resolution = input$resolution)
        tmp_data <- RunTSNE(tmp_data, dims = 1:30, perplexity=5)
        res_plot = DimPlot(tmp_data,reduction = "tsne",pt.size = .1)
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = res_plot, "tmp_data" = tmp_data))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
  output$ClassificationPlot <- renderPlot({
    tmp <- ClassificationReactive()
    if(!is.null(tmp)){
      tmp$data
    }
  })
  output$ClassificationAvailable <-
    reactive({
      return(!is.null(ClassificationReactive()$data))
    })
  outputOptions(output, 'ClassificationAvailable', suspendWhenHidden=FALSE)
  output$downloadClassificationPlot <- downloadHandler(
    filename = function()  {'plot_pred_subclusters.pdf'},
    content = function(file) {
      pdf(file)
      print(ClassificationReactive()$data)
      dev.off()}
  )
  # anno_cellsubtype
  observe({
    SubannotionReactive()
  })
  SubannotionReactive <- eventReactive(input$startSubannotion, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = ClassificationReactive()$tmp_data
        #cur_celltype = data_rds@meta.data$pred_cell[1]
        cur_celltype = input$celltype1
        source_python(system.file("app/www/python/BP5_new.py", package='scWizard'))
        if(cur_celltype == 'End')
        {
          X_total_path=system.file('app/www/python/trainset/trainx_End.h5', package='scWizard')
          Y_total=read.csv(system.file('app/www/python/trainset/trainy_End.csv', package='scWizard'))
        }
        else if(cur_celltype == 'Fib')
        {
          X_total_path=system.file('app/www/python/trainset/trainx_Fib.h5', package='scWizard')
          Y_total=read.csv(system.file('app/www/python/trainset/trainy_Fib.csv', package='scWizard'))
        }
        else if(cur_celltype == 'Mye')
        {
          X_total_path=system.file('app/www/python/trainset/trainx_Mye.h5', package='scWizard')
          Y_total=read.csv(system.file('app/www/python/trainset/trainy_Mye.csv', package='scWizard'))
        }
        else if(cur_celltype == 'T&NK')
        {
          X_total_path=system.file('app/www/python/trainset/trainx_T.h5', package='scWizard')
          Y_total=read.csv(system.file('app/www/python/trainset/trainy_T.csv', package='scWizard'))
        }
        else if(cur_celltype == 'CD4T')
        {
          X_total_path=system.file('app/www/python/trainset/trainx_CD4T.h5', package='scWizard')
          Y_total=read.csv(system.file('app/www/python/trainset/trainy_CD4T.csv', package='scWizard'))
        }
        else
        {
          X_total_path=system.file('app/www/python/trainset/trainx_CD8T.h5', package='scWizard')
          Y_total=read.csv(system.file('app/www/python/trainset/trainy_CD8T.csv', package='scWizard'))
        }
        Y_total$X=NULL
        colnames(Y_total)=c('celltype')
        num_classes = length(unique(Y_total$celltype))
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        X_verify=t(as.data.frame(data_rds@assays$RNA@data))
        X_verify=as.data.frame(X_verify)
        subclusters=as.vector(data_rds@active.ident)
        subclusters=as.data.frame(subclusters)
        res_celltype=get_BP5_res(X_total, Y_total$celltype, X_verify, subclusters, num_classes, input$PCAk2, input$layer2_1, input$layer2_2, input$layer2_3, input$regularization2, input$regularization2, input$learning2)
        names(res_celltype) <- levels(data_rds)
        data_rds <- RenameIdents(data_rds, res_celltype)
        data_rds@meta.data$pred_sub_cell = as.vector(data_rds@active.ident)
        res_plot = DimPlot(data_rds,reduction = "tsne",pt.size = .1,group.by = 'pred_sub_cell')
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = res_plot, "data_rds" = data_rds))
      },
      error=function(cond){
        message("Here's the original error.")
        return(NULL)
      })
      
    })
  })
  output$SubannotionPlot <- renderPlot({
    tmp <- SubannotionReactive()
    if(!is.null(tmp)){
      tmp$data
    }
  })
  output$SubannotionAvailable <-
    reactive({
      return(!is.null(SubannotionReactive()$data))
    })
  outputOptions(output, 'SubannotionAvailable', suspendWhenHidden=FALSE)
  output$downloadSubannotionPlot <- downloadHandler(
    filename = function()  {'plot_pred_sub_cell.pdf'},
    content = function(file) {
      pdf(file)
      print(SubannotionReactive()$data)
      dev.off()}
  )
  
  ### Find markers
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  observe({
    findClusterMarkersReactive()
  })
  findClusterMarkersReactive <- eventReactive(input$findClusterMarkers, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = inputDataReactive()$data
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        if(is.na(input$mindiffpct) && is.na(input$maxcellsperident))
        {
          cluster.markers <- FindAllMarkers(object = data_rds, logfc.threshold = input$threshAll, min.pct = input$minPct, 
                                            test.use = input$testuse, only.pos = input$onlypos)
        }
        else if(is.na(input$mindiffpct) && (!is.na(input$maxcellsperident)))
        {
          cluster.markers <- FindAllMarkers(object = data_rds, logfc.threshold = input$threshAll, min.pct = input$minPct, 
                                            test.use = input$testuse, only.pos = input$onlypos,max.cells.per.ident = input$maxcellsperident)
        }
        else if((!is.na(input$mindiffpct)) && is.na(input$maxcellsperident))
        {
          cluster.markers <- FindAllMarkers(object = data_rds, logfc.threshold = input$threshAll, min.pct = input$minPct, 
                                            test.use = input$testuse, only.pos = input$onlypos,min.diff.pct = input$mindiffpct)
        }
        else
        {
          cluster.markers <- FindAllMarkers(object = data_rds, logfc.threshold = input$threshAll, min.pct = input$minPct, 
                                            test.use = input$testuse, only.pos = input$onlypos,min.diff.pct = input$mindiffpct,max.cells.per.ident = input$maxcellsperident)
        
        }
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = cluster.markers))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
  output$clusterMarkers <- renderDataTable({
    tmp <- findClusterMarkersReactive()
    if(!is.null(tmp)){
      tmp$data
    }
  })
  output$clusterMarkersAvailable <-
    reactive({
      return(!is.null(findClusterMarkersReactive()$data))
    })
  outputOptions(output, 'clusterMarkersAvailable', suspendWhenHidden=FALSE)
  
  output$downloadClusterMarkersCSV <- downloadHandler(
    filename = function()  {paste0(findClusterMarkersReactive()$data,".csv")},
    content = function(file) {
      write.csv(findClusterMarkersReactive()$data, file, row.names=TRUE)}
  )
  
  ### inferCNV
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  #choose celltype
  output$myselectinfercnvbox <-
    renderUI({
      if(input$startAnnotion > 0)
        data_rds = AnnotionReactive()$data
      else
        data_rds = inputDataReactive()$data
      label=names(summary(data_rds@active.ident))
      selectInput("refgroupnames", "ref_group_names",
                  choices =c(NULL ,label), selected = NULL, multiple = TRUE)
    })
  observe({
    infercnvReactive()
  })
  infercnvReactive <- eventReactive(input$startinfercnv, {
    withProgress(message = "Processing,please wait",{
      if(input$startAnnotion > 0)
        data_rds = AnnotionReactive()$data
      else
        data_rds = inputDataReactive()$data
      tryCatch({  
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        cellinfo = data_rds@meta.data
        cellinfo$cellname = row.names(cellinfo)
        annofile = cellinfo[, c('cellname', 'celltype')]
        row.names(annofile) = NULL
        write.table(annofile, './annofile.txt', quote = F, row.names = F, col.names = F, sep = '\t')
        infercnv_obj = CreateInfercnvObject(raw_counts_matrix=as.matrix(data_rds@assays$RNA@counts),
                                            annotations_file='./annofile.txt',
                                            delim="\t",
                                            gene_order_file=input$geneorderfile$datapath,
                                            ref_group_names=input$refgroupnames) 
        
        infercnv_obj = infercnv::run(infercnv_obj,
                                     cutoff=input$cutoff, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                                     out_dir='./infercnv_res', 
                                     cluster_by_groups=TRUE, 
                                     denoise=TRUE,
                                     HMM=TRUE)
        file.remove('./annofile.txt')
        shiny::setProgress(value = 0.8, detail = "Done.")
      },
      error=function(cond) {
        message("Here's the original error.")
        #message(cond)
        return(NA)
      })
    })
  })
  
  ### pseudotime
  observe({
    
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
    
  })
  
  # return selectbox
  output$myselectbox3 <-
    renderUI({
      data_rds = inputDataReactive()$data
      label=names(summary(data_rds@active.ident))
      selectInput("celltype", "choose celltype",
                  choices =c('all' ,label), selected = 'all')
    })
  observe({
    MonocleReactive()
  })
  MonocleReactive <- eventReactive(input$startMonocle, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = inputDataReactive()$data
        if(input$celltype=='all')
        {
          new_data_rds = data_rds
        }
        else
        {
          new_data_rds = subset(data_rds, idents = input$celltype)
        }
      
        cellinfo<-new_data_rds@meta.data
      
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
      
        gene_annotation<-data.frame(gene_short_name=row.names(new_data_rds))#“gene_short_name”
        row.names(gene_annotation)<-gene_annotation[,1]
        pd <- new("AnnotatedDataFrame", data = cellinfo)#cell information
        fd <- new("AnnotatedDataFrame", data = gene_annotation)#gene information
        HSMM <- newCellDataSet(new_data_rds@assays$RNA@counts,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = input$lowerdetectionlimit,
                              expressionFamily = negbinomial.size())
        # pre-process
        HSMM <- estimateSizeFactors(HSMM)
        HSMM <- estimateDispersions(HSMM)
        # choose special gene
        disp_table <- dispersionTable(HSMM)
        unsup_clustering_genes <- subset(disp_table, mean_expression >= input$meanexpression)
        HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
        # Dimensionality reduction by DDRtree
        HSMM_myo <- reduceDimension(HSMM, max_components = input$maxcomponents, method = input$rmethod)
        # calculation pseudotime
        HSMM_myo <- orderCells(HSMM_myo)
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = HSMM_myo))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
  output$MonocleAvailable <-
    reactive({
      return(!is.null(MonocleReactive()$data))
    })
  outputOptions(output, 'MonocleAvailable', suspendWhenHidden=FALSE)
  
  output$downloadMonocleRDS <- downloadHandler(
    filename = function()  {"res_moncle.rds"},
    content = function(file) {
      saveRDS(MonocleReactive()$data, file)}
  )
  
  
  ### TF-scenic
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  observe({
    startCoexpressionReactive()
    startGRNReactive()
  })
  startCoexpressionReactive <- eventReactive(input$startCoexpression, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = inputDataReactive()$data
        cellalltype = as.vector(data_rds@active.ident)
        cellalltype = as.data.frame(cellalltype)
        meta_data <- cbind(data_rds@meta.data, cellalltype)
        cellInfo = meta_data[,c(7,2,3)]
        colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
        saveRDS(cellInfo,"./int/cellInfo.RDS")
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        # init
        exprMat <- as.matrix(data_rds@assays$RNA@counts)
        exprMat <- exprMat[which(rowSums(exprMat)>0),]
        mydbDIR <- "./cisTarget"
        mydbs <- c("hg19-500bp-upstream-7species.mc9nr.feather",
                  "hg19-tss-centered-10kb-7species.mc9nr.feather")
        names(mydbs) <- c("500bp", "10kb")
        scenicOptions <- initializeScenic(org=input$org, dbDir=mydbDIR, dbs = mydbs, nCores=1)
        saveRDS(scenicOptions, "./int/scenicOptions.rds")
        # build co-expression net
        genesKept <- geneFiltering(exprMat, scenicOptions, 
                                  minCountsPerGene = input$minCountsPerGene1 * input$minCountsPerGene2 * ncol(exprMat), 
                                  minSamples = ncol(exprMat) * input$minSamples)
        exprMat_filtered <- exprMat[genesKept, ]
        runCorrelation(exprMat_filtered, scenicOptions)
        exprMat_filtered_log <- log2(exprMat_filtered+1)
        runGenie3(exprMat_filtered_log, scenicOptions)
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
  startGRNReactive <- eventReactive(input$startGRN, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = inputDataReactive()$data
        exprMat <- as.matrix(data_rds@assays$RNA@counts)
        exprMat <- exprMat[which(rowSums(exprMat)>0),]
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        # Build and score the GRN
        scenicOptions <- readRDS("./int/scenicOptions.rds")
        exprMat_log <- log2(exprMat+1)
        scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
        scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
        scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                                    coexMethod=c("top5perTarget")) # Toy run settings
        library(doParallel)
        scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
        scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
        tsneAUC(scenicOptions, aucType="AUC") # choose settings
        #motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
 
  ### Correlation
  observe({
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  # return selectbox
  output$myselectbox5 <-
    renderUI({
      data_rds = inputDataReactive()$data
      label=names(summary(data_rds@active.ident))
      
      selectInput("celltype", "choose celltype",
                  choices =c('all' ,label), selected = 'all')
    })
  # return selectgroupbox
  output$myselectgroupbox5 <-
    renderUI({
      data_rds = inputDataReactive()$data
      dt_mat = as.matrix(data_rds@assays$RNA@counts)
      new_dt_mat = dt_mat[which(rowSums(dt_mat)>0),]
      rowName = row.names(new_dt_mat)
      selectInput("geneName", "choose gene",
                  choices =c('all' ,rowName), selected = 'all', multiple = TRUE)
    })
  observe({
    startCorReactive()
  })
  startCorReactive <- eventReactive(input$startCor, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = inputDataReactive()$data
        if(input$celltype=='all')
        {
          new_data_rds = data_rds
        }
        else
        {
          new_data_rds = subset(data_rds, idents = input$celltype)
        }
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        dt_mat = as.matrix(new_data_rds@assays$RNA@counts)
        if(input$geneName=='all')
        {
          new_dt_mat = dt_mat
        }
        else
        {
          new_dt_mat = dt_mat[input$geneName,]
        }
        new_dt_mat = t(new_dt_mat)
        resCor <- cor(new_dt_mat, use = input$use,
                      method = input$method)
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = resCor))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
  output$Cor <- renderDataTable({
    tmp <- startCorReactive()
    if(!is.null(tmp)){
      tmp$data
    }
  })
  output$CorAvailable <-
    reactive({
      return(!is.null(startCorReactive()$data))
    })
  outputOptions(output, 'CorAvailable', suspendWhenHidden=FALSE)
  output$downloadCorCSV <- downloadHandler(
    filename = function()  {"res_cor.csv"},
    content = function(file) {
      write.csv(startCorReactive()$data, file, row.names=TRUE)}
  )
  
  
  ### cell-cell communication
  volumes = getVolumes()()
  observe({
    shinyDirChoose(input, "cellphonedbin", roots = volumes)
    shinyDirChoose(input, "cellphonedbout", roots = volumes)
    if(!is.null(inputDataReactive()))
    {
      data_rds = inputDataReactive()$data
    }
  })
  observe({
    cellphonedbReactive()
  })
  cellphonedbReactive <- eventReactive(input$startcellphonedb, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        data_rds = inputDataReactive()$data
        shiny::setProgress(value = 0.4, detail = "Calculating ...")
        data_rds = inputDataReactive()$data
        counts_file = paste0(as.character(parseDirPath(volumes, input$cellphonedbin)), "/count.txt")
        meta_file = paste0(as.character(parseDirPath(volumes, input$cellphonedbin)), "/meta.txt")
        dir.create("cellphonedb_in")
        count_data = as.matrix(data_rds@assays$RNA@data)
        write.table(as.matrix(count_data), counts_file, sep='\t', quote=F)
        cellalltype = as.vector(data_rds@active.ident)
        cellalltype = as.data.frame(cellalltype)
        meta_data <- cbind(rownames(data_rds@meta.data), cellalltype)
        meta_data <- as.matrix(meta_data)
        meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
        write.table(meta_data, meta_file, sep='\t', quote=F, row.names=F)
        out_file = as.character(parseDirPath(volumes, input$cellphonedbout))
        counts_data = input$countsdata
        command = paste("cellphonedb method statistical_analysis",meta_file)
        command = paste(command,counts_file)
        command = paste(command,"--output-path=")
        command = paste0(command,out_file)
        command = paste(command,"--counts-data=")
        command = paste0(command,counts_data)
        command = paste(command,"--threads=")
        command = paste0(command,10)
        print(command)
        os = import("os")
        os$system(command)
        shiny::setProgress(value = 0.8, detail = "Done.")
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      }) 
    })
  })
  
  ### Plot
  # return selectbox
  output$sel_box <-
    renderUI({
      if(input$methodsUsed=='GSVA')
      {
        selectInput("Gragh", "choose gragh",
                    choices = c('heatmap'), selected = 'heatmap')
      }
      else if(input$methodsUsed=='CellphoneDB')
      {
        selectInput("Gragh", "choose gragh",
                    choices = c('heatmap1', 'heatmap2', 'bubble'), selected = 'heatmap1')
      }
      else if(input$methodsUsed=='Monocle')
      {
        selectInput("Gragh", "choose gragh",
                    choices = c('cellTrajectory_pse','cellTrajectory_sta','cellTrajectory_seu'), selected = 'cellTrajectory_pse')
      }
      else if(input$methodsUsed=='Scenic')
      {
        selectInput("Gragh", "choose gragh",
                    choices = c('heatmap_step3', 'heatmap_step4', ), selected = 'heatmap')
      }
      else if(input$methodsUsed=='Cor')
      {
        selectInput("Gragh", "choose gragh",
                    choices = c('heatmap'), selected = 'heatmap')
      }
    })
  observe({
    plotReactive()
  })
  plotReactive <- eventReactive(input$startPlot, {
    withProgress(message = "Processing,please wait",{
      tryCatch({
        shiny::setProgress(value = 0.4, detail = "Ploting ...")
        if(input$methodsUsed=='GSVA')
        {
          gsva_dt = read.csv(input$plotfile1$datapath)
          gsva_dt = cor_dt[-1,-1]
          gsva_dt_matrix = data.matrix(cor_dt)
          res_plot = heatmap(gsva_dt_matrix, Rowv=NA, Colv=NA, col=cm.colors(256), revC=TRUE, scale='column')
        }
        else if(input$methodsUsed=='CellphoneDB')
        {
          source("plot_heatmap.R", local = TRUE)
          source("plot_dot.R", local = TRUE)
          pval = input$plotfile1$datapath
          if(input$Gragh=='heatmap1')
          {
            meta = input$plotfile2$datapath
            res_plot = heatmaps_plot(meta, pval,count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf", count_network_filename = "count_network.txt", interaction_count_filename = "interaction_count.txt",
                                    count_network_separator = "\t", interaction_count_separator = "\t")
            res_plot = res_plot$p1
          }
          else if(input$Gragh=='heatmap2')
          {
            meta = input$plotfile2$datapath
            res_plot = heatmaps_plot(meta, pval,count_filename = "heatmap_count.pdf", log_filename = "heatmap_log_count.pdf", count_network_filename = "count_network.txt", interaction_count_filename = "interaction_count.txt",
                                    count_network_separator = "\t", interaction_count_separator = "\t")
            res_plot = res_plot$p2
          }
          else if(input$Gragh=='bubble')
          {
            mean = input$plotfile2$datapath
          
            res_plot = dot_plot(means_path = mean, pvalues_path = pval)
          }
        }
        else if(input$methodsUsed=='Monocle')
        {
          monocle_dt = readRDS(input$plotfile1$datapath)
          if(input$Gragh=='cellTrajectory_pse')
          {
            res_plot = plot_cell_trajectory(monocle_dt, color_by = "Pseudotime")
          }
          else if(input$Gragh=='cellTrajectory_sta')
          {
            res_plot = plot_cell_trajectory(monocle_dt, color_by = "State")
          }
          else if(input$Gragh=='cellTrajectory_seu')
          {
            res_plot = plot_cell_trajectory(monocle_dt, color_by = "seurat_clusters")
          }
        }
        else if(input$methodsUsed=='Scenic')
        {
        
        }
        else if(input$methodsUsed=='Cor')
        {
          cor_dt = read.csv(input$plotfile$datapath)
          cor_dt = cor_dt[-1,-1]
          cor_dt_matrix = data.matrix(cor_dt)
          res_plot = heatmap(cor_dt_matrix, Rowv=NA, Colv=NA, col=cm.colors(256), revC=TRUE, scale='column')
        }
        shiny::setProgress(value = 0.8, detail = "Done.")
        return(list("data" = res_plot))
      },
      error=function(cond) {
        message("Here's the original error.")
        return(NULL)
      })
    })
  })
  output$plot <- renderPlot({
    tmp = plotReactive()
    if(!is.null(tmp))
    {
      tmp$data
    }
  })
  output$saveAvailable <-
    reactive({
      return(!is.null(plotReactive()$data))
    })
  outputOptions(output, 'saveAvailable', suspendWhenHidden=FALSE)
  output$downloadplot <- downloadHandler(
    filename = function()  {"plot.pdf"},
    content = function(file) {
      ggsave(file, plotReactive()$data )}
  )
}
