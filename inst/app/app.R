library(shiny)
library(FunCZIDM)
library(rhandsontable)
library(reticulate)
library(shinyMatrix)
library(shinyjqui)

# Allow up to 100 MB uploads
options(shiny.maxRequestSize = 400*1024^2)

useCondaEnv <- shiny::getShinyOption("useCondaEnv", FALSE)
if (useCondaEnv) {
  reticulate::use_condaenv("shinyPy", required = TRUE)
} else {
  print("Configured Python:")
  py_config()
}

pyPackages <- c("numpy", "pandas", "matplotlib", "seaborn")
missing <- pyPackages[!vapply(pyPackages, reticulate::py_module_available, 
                              logical(1))]

if (useCondaEnv) {
  if (length(missing)) {
    stop("Missing Python modules in env shinyPy: ",
         paste(missing, collapse = ", "))
  }
} else {
  py_require(missing)
}

source("helpers.R")
source_python("plotFunctions.py")

# Load a default output
defaultData <- readRDS(system.file("extdata", "caseStudySamples.rds",
                                    package = "FunCZIDM"))

ui <- fluidPage(
  titlePanel("FunC-ZIDM Plots"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        inputId = "file", 
        label = "Upload Output rds (optional)", 
        accept = c(".rds")
      ),
      selectInput(
        inputId = "plotType", 
        label = "Choose a plot:",
        choices = c("Relative Abundance", 
                    "Multiplicative Change in Relative Abundance",
                    "Diversity",
                    "Multiplicative Change in Diversity"),
        selected = "Relative Abundance"
      ),
      conditionalPanel(
        condition = "['Multiplicative Change in Relative Abundance', 
                      'Relative Abundance'].includes(input.plotType)",
        selectInput(
          inputId = "category", 
          label = "Choose a category:",
          choices = defaultData$catNames,
          selected = defaultData$catNames[1],
          multiple = TRUE
        )
      ),
      conditionalPanel(
        condition = "['Multiplicative Change in Relative Abundance',
                'Multiplicative Change in Diversity'].includes(input.plotType)",
        selectInput(
          inputId = "covariate", 
          label = "Choose a covariate:",
          choices = unlist(defaultData$colMapping)[-1],
          selected = unlist(defaultData$colMapping)[2]
        ),
        uiOutput("changeUI")
      ),
      uiOutput("plotBetas_ui"),
      conditionalPanel(
        condition = "['Multiplicative Change in Diversity',
                      'Diversity'].includes(input.plotType)",
        tags$div(
          style = "margin-bottom: 4px; display: inline-block;",
          tags$span(withMathJax("\\(l\\):"),
                    style = "font-weight: bold; margin-right: 4px;"),
          tags$span(icon("question-circle"),
                    title = 
"This parameter controls the weight given towards smaller categories. 1 gives equal weight to all categories, while 0 weights larger categories more.",
                    style = "cursor: help; color: #555;")
        ),
        numericInput(inputId = "l", label = NULL, value = 0, min = 0, max = 1)
      ),
      tags$div(
        style = "display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
        actionButton("createPlot", "Create Plot"),
        downloadButton("downloadPlot", "Download plot"),
        downloadButton("downloadPlotData", "Download plot data")
      ),
    ),
    mainPanel(
      jqui_resizable(
        tags$div(
          id    = "imgContainer",
          style = "width:800px; height:500px; overflow:auto; 
                   border:1px solid #ccc;",
          imageOutput("pyPlot")
        ),
        options = list(handles = "n, e, s, w, ne, se, sw, nw")
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      wellPanel(
        h4("Covariate Profile"),
        rHandsontableOutput("covProfileTable")
      )
    )
  ),
  fluidRow(
    column(
      width = 12,
      wellPanel(
        h4("Plot Settings"),
        uiOutput("plotOptions")
      )
    )
  )
)


server <- function(input, output, session) {
  # Reactive that holds either the uploaded data or defaultData
  data_reactive <- reactive({
    req <- input$file
    if (is.null(req)) {
      data <- defaultData
    } else {
      data <- readRDS(req$datapath)
    }
    data
  })
  
  # Whenever data_reactive() changes, update the covariate dropdown choices
  observeEvent(data_reactive(), {
    data <-  data_reactive()
    updateSelectInput(
      session,
      inputId = "covariate",
      choices = unlist(data$colMapping)[-1],
      selected = unlist(data$colMapping)[2]
    )
    updateSelectInput(
      session,
      inputId = "category",
      choices = data$catNames,
      selected = data$catNames[1]
    )
  })

  # Build a data.frame of covariate names and default values
  covProfileDF <- reactive({
    data <- data_reactive()
    df <- data.frame(matrix(ncol=0, nrow = 1))
    startingVals <- data$centerScaleList
    covLabels <- unlist(data$colMapping)
    for (i in 2:length(covLabels)) {
      if (is.character(startingVals[[i]]) || is.null(startingVals[[i]])) {
        df[covLabels[i]] = FALSE
      } else {
        df[covLabels[i]] = startingVals[[i]]$center
      } 
    }
    df
  })

  # Render the editable table
  output$covProfileTable <- renderRHandsontable({
    df <- covProfileDF()
    rhandsontable(df, rowHeaders = NULL)
  })

  output$changeUI <- renderUI({
    req(input$covariate)
    req(input$covProfileTable)
    data <- data_reactive()
    df <- if (is.null(input$covProfileTable)) covProfileDF() else 
          hot_to_r(input$covProfileTable)
    
    covIdx <- which(colnames(df) == input$covariate) + 1 
    # grab that one cell (it might be numeric or logical)
    val <- df[[input$covariate]]

    # if it’s NOT logical, show the “Amount Change” widget
    if (!is.character(data$centerScaleList[[covIdx]])) {
      min <- data$centerScaleList[[covIdx]]$center - 
              data$centerScaleList[[covIdx]]$scale*4
      max <- data$centerScaleList[[covIdx]]$center + 
              data$centerScaleList[[covIdx]]$scale*4
      minFromVal <- round(min - val, 4)
      maxFromVal <- round(max - val, 4)
      tagList(
        checkboxInput(
          inputId = "changeRange",
          label =  withMathJax("Range of Changes"), 
          value = FALSE
        ),
        conditionalPanel(
          condition = "input.changeRange == 1",
          checkboxInput(inputId="inCI", label="1 if mean inside C.I.", 
                        value = FALSE),
          tags$div(
            style = "margin-bottom: 4px; display: inline-block;",
            tags$span("Range of Change:",
                      style = "font-weight: bold; margin-right: 4px;"),
            tags$span(icon("question-circle"),
                      title = 
"The amount change in the covariate selected from its starting point in covariate profile.",
                      style = "cursor: help; color: #555;")
          ),
          sliderInput(inputId="changeSlide", label=NULL, min=minFromVal, 
                      max=maxFromVal, value=c(minFromVal, maxFromVal)),
          tags$div(
            style = "margin-bottom: 4px; display: inline-block;",
            tags$span("Number of changes:",
                      style = "font-weight: bold; margin-right: 4px;"),
            tags$span(icon("question-circle"),
                      title = 
"The amount of points to estimate the change. 2 would just be the end points selected.",
                      style = "cursor: help; color: #555;")
          ),
          numericInput(inputId="numChangeSlide", label=NULL, value=2, min=2,
                       max=50)
        ),
        conditionalPanel(
          condition = "input.changeRange == 0",
          tags$div(
            style = "margin-bottom: 4px; display: inline-block;",
            tags$span("Amount Change:",
                      style = "font-weight: bold; margin-right: 4px;"),
            tags$span(icon("question-circle"),
                      title = 
"The amount change in the covariate selected from its starting point in covariate profile.",
                      style = "cursor: help; color: #555;")
          ),
          numericInput(inputId="changeNum", label=NULL, value=1)
        )
      )
    } else {
      NULL
    }
  })

  output$plotBetas_ui <- renderUI({
    req(input$covariate)
    req(input$covProfileTable)
    data <- data_reactive()
    df <- if (is.null(input$covProfileTable)) covProfileDF() else 
          hot_to_r(input$covProfileTable)
    
    covIdx <- which(colnames(df) == input$covariate) + 1 

    # When the checkbox should be visible
    raShow  <- input$plotType == "Relative Abundance"
    deltaRAShow <- input$plotType == 
                                   "Multiplicative Change in Relative Abundance" 

    # Preserve current state if it exists; default TRUE on first render
    current <- isolate(input$plotBetas)
    val <- if (is.null(current)) TRUE else isTRUE(current)

    if (isTRUE(deltaRAShow)) {
      if (!is.character(data$centerScaleList[[covIdx]])) {
        conditionalPanel(
          condition = "input.changeRange == 0",
          checkboxInput(
            inputId = "plotBetas",
            label =  withMathJax("Plot \\(\\beta_{jp}(t)\\)"), 
            value = val
          )
        )
      } else {
        checkboxInput(
          inputId = "plotBetas",
          label =  withMathJax("Plot \\(\\beta_{jp}(t)\\)"), 
          value = val
        )
      }
    } else if (isTRUE(raShow)) {
      checkboxInput(
        inputId = "plotBetas",
        label = withMathJax("Plot \\(\\beta_{j0}(t)\\)"),
        value = val
      )
    } else {
      NULL
    }
  })

  output$plotOptions <- renderUI({
    # Preserve current state if it exists
    suptitle <- isolate(input$suptitle)
    if (is.null(suptitle)) suptitle <- "FunC-ZIDM Plots"
    titleSize <- isolate(input$titleSize)
    if (is.null(titleSize)) titleSize <- 16
    xLabel <- isolate(input$xLabel)
    if (is.null(xLabel)) xLabel <- "Days After Birth"
    fontSize <- isolate(input$fontSize)
    if (is.null(fontSize)) fontSize <- 12
    plotWidth <- isolate(input$plotWidth)
    if (is.null(plotWidth)) plotWidth <- 10
    plotHeight <- isolate(input$plotHeight)
    if (is.null(plotHeight)) plotHeight <- 6

    # --- static controls (always shown) ---
    base <- list(
      textInput("suptitle", "Plot Title:", suptitle),
      numericInput("titleSize", "Title Font Size:", value=titleSize, min=1),
      textInput("xLabel", "X-axis Label:", xLabel),
      numericInput("fontSize",   "Font Size:", value = fontSize, min = 1),
      numericInput("plotWidth",  "Plot Width (inches):", 
                   value = plotWidth, min = 1),
      numericInput("plotHeight", "Plot Height (inches):", 
                   value = plotHeight,  min = 1)
    )

    # --- dynamic controls (conditionally shown) ---
    dyn <- list()

    changeRA <- input$plotType == "Multiplicative Change in Relative Abundance"
    divPlot <- input$plotType == "Multiplicative Change in Diversity" ||
                input$plotType == "Diversity"
    changeRanges <- !is.null(input$changeRange) && isTRUE(input$changeRange) && 
                    input$plotType != "Relative Abundance" &&
                    input$plotType != "Diversity"
    plotBetas <- !is.null(input$plotBetas) && isTRUE(input$plotBetas)

    # Preserve current state if it exists
    lineWidth <- isolate(input$lineWidth)
    if (is.null(lineWidth)) lineWidth <- 1
    xAxisRange <- isolate(input$xAxisRange)
    if (is.null(xAxisRange)) xAxisRange <- c(0, 50)
    suptitleLoc <- isolate(input$suptitleLoc)
    if (is.null(suptitleLoc)) suptitleLoc <- c(0.55, 1.025)
    gridDims <- isolate(input$gridDims)
    if (is.null(gridDims)) gridDims <- c(length(input$category), 1)
    sameYaxis <- input$sameYaxis
    if (is.null(sameYaxis)) sameYaxis <- FALSE

    if (!changeRanges) {
      dyn <- c(dyn, list(
        numericInput("lineWidth", "Line Width:", value=lineWidth,  min=0.1),
        checkboxInput("sameYaxis", "Same Y-axis for all categories", 
                      value=sameYaxis)
      ))
      
      if (isTRUE(sameYaxis)) {
        dyn <- c(dyn, list(
          numericInput(
            inputId = "by",
            label = "Number of Y-axis Ticks:",
            value = 5,
            min = 1,
            step = 1
          )))
      } else {
        val <- ifelse((input$plotType == "Relative Abundance"), 0.05, 0.5)
        dyn <- c(dyn, list(
          numericInput(
            inputId = "by",
            label = "Y-axis Tick Spacing:",
            value = val,
            min = .001
          )))
        if (plotBetas && !divPlot) {
          dyn <- c(dyn, list(
            numericInput("byBeta", "Beta's Y-axis Tick Spacing:", value = 0.5, 
                         min = 0.001)
          ))
        }
      } 

      if (!plotBetas && !divPlot) {
        dyn <- c(dyn, list(
          matrixInput(
            "gridDims", "Grid dimensions:",
            value = matrix(gridDims, nrow = 1,
                           dimnames = list(NULL, c("Rows","Cols"))),
            class = "numeric", rows = list(names = FALSE), 
            cols = list(names = TRUE)
          )
        ))
      }

      dyn <- c(dyn, list(
        matrixInput("legendLoc", "Legend Location:",
          value = matrix(c(.55, -0.05), nrow = 1, 
          dimnames = list(NULL, c("x","y"))),
          class = "numeric", rows = list(names = FALSE), 
          cols = list(names = TRUE)
        )))
    } else if (changeRA || 
               input$plotType == "Multiplicative Change in Diversity") {
      dyn <- c(dyn, list(
        textInput("yLabel", "Y-axis Label:", 
                  "Weeks from 27 Weeks of Gestational Age"),
        numericInput("by", "Number of Y-axis Ticks:",
                     value = input$numChangeSlide, min = 1),
        numericInput("byX", "Number of X-axis Ticks:",
                     value = 50, min = 1)
      ))
    }

    # You can also append your other matrixInputs here:
    dyn <- c(dyn, list(
      matrixInput("xAxisRange", "X-axis Range:",
        value = matrix(xAxisRange, nrow = 1,
                       dimnames = list(NULL, c("min","max"))),
        class = "numeric", rows = list(names = FALSE), 
        cols = list(names = TRUE)
      ),
      matrixInput("suptitleLoc", "Title Location:",
        value = matrix(suptitleLoc, nrow = 1, 
                       dimnames = list(NULL, c("x","y"))),
        class = "numeric", rows = list(names = FALSE), cols = list(names = TRUE)
      )
    ))

    # Return one full flowLayout
    do.call(flowLayout, c(base, dyn))
  })

  plotRequest <- eventReactive(input$createPlot, {
    output   <- data_reactive()
    covTable  <- if (is.null(input$covProfileTable)) covProfileDF() else
                 hot_to_r(input$covProfileTable)
    covVals   <- c(1, unlist(covTable[1, ], use.names = TRUE))
    cov    <- input$covariate
    cats    <- input$category
    pType      <- input$plotType
    covIdx <- which(colnames(covTable) == cov) + 1 
    if (!is.null(input$changeRange) && input$changeRange && 
        !is.character(output$centerScaleList[[covIdx]])) {
      changeEnds <- input$changeSlide
      numChanges <- input$numChangeSlide
      inCI <- input$inCI
      covChange <- seq(changeEnds[1], changeEnds[2], length.out=numChanges) 
    } else {
      covChange <- input$changeNum
    }
    lVal      <- input$l
    plotBetas <- input$plotBetas

    # Pull plot settings
    suptitle <- input$suptitle
    titleSize <- input$titleSize
    xLabel <- input$xLabel
    figsize <- c(input$plotWidth, input$plotHeight)
    fontSize <- input$fontSize
    lineWidth <- input$lineWidth
    by <- input$by
    byBeta <- input$byBeta
    byX <- input$byX
    if (plotBetas) {
      figRowCol <- NULL
    } else {
      figRowCol <- c(input$gridDims[1,]["Rows"], input$gridDims[1,]["Cols"])
    }
    legendLoc <- c(input$legendLoc[1,]["x"], input$legendLoc[1,]["y"])
    suptitleLoc <- c(input$suptitleLoc[1,]["x"], input$suptitleLoc[1,]["y"])
    xAxisRange <- c(input$xAxisRange[1,]["min"], input$xAxisRange[1,]["max"])
    sameYaxis <- input$sameYaxis
    
    # ---- Processing the Data ----
    covVals <- getCenterScaledCovProfile(covVals, output$centerScaleList)
    if (pType == "Multiplicative Change in Relative Abundance" ||
        pType == "Multiplicative Change in Diversity") {
      heatmapIndex <- covChange
      covChange <- getCenterScaledChange(covChange, output$centerScaleList, 
                                        which(unlist(output$colMapping) == cov),
                                         covVals)
    }

    if (pType == "Relative Abundance") {
      plotData <- RAData(output, covVals, xAxisRange, interceptBetas=plotBetas)
      dataDict <- list()

      for (cat in cats) {
        colIdx <- which(output$catNames == cat)
        confInt <- t(apply(plotData$RA[,colIdx,], 1,
                         function(x) quantile(x, c(0.025, 0.975))))
        dataDict[[cat]] = list()
        dataDict[[cat]][["RA"]] <- cbind(confInt,
                                     rowMeans(plotData$RA[,colIdx,]),
                                     plotData$points)
        colnames(dataDict[[cat]][["RA"]]) <- c("lowerCI", "upperCI", "mean", 
                                               "time")
        dataDict[[cat]][["RA"]] <- data.frame(dataDict[[cat]][["RA"]])
        if (plotBetas) {
          confInt <- t(apply(plotData$intercepts[,colIdx,], 1,
                         function(x) quantile(x, c(0.025, 0.975))))
          dataDict[[cat]][["beta"]] <- cbind(confInt,
                                        rowMeans(plotData$intercepts[,colIdx,]),
                                        plotData$points)
          colnames(dataDict[[cat]][["beta"]]) <- c("lowerCI", "upperCI", "mean",
                                                   "time")
          dataDict[[cat]][["beta"]] <- data.frame(dataDict[[cat]][["beta"]])
        }
      }
      raTitlePrefix <- if (plotBetas) "$\\mathrm{RA}_{cp}(t)$ of " else ""
      betaTitlePrefix <- "$\\beta_{cp}(t)$ of "
      raYlabel <- "$\\mathrm{RA}_{cp}(t)$"
      betaYlabel <- "$\\beta_{cp}(t)$"
      changeIn <- FALSE
      fig <- curvePlots(dataDict, suptitle, raTitlePrefix, betaTitlePrefix,
                        raYlabel, betaYlabel, xLabel, changeIn, plotBetas, 
                        lineWidth=lineWidth, byRA=by, byBeta=byBeta, 
                        fontsize=fontSize, titleSize=titleSize, figsize=figsize,
                        xAxisRange=xAxisRange, figRowCol=figRowCol, 
                        legendLoc=legendLoc, suptitleLoc=suptitleLoc, 
                        sameYaxis=sameYaxis)
    } else if (pType == "Multiplicative Change in Relative Abundance") {
      dataDict <- deltaRAData(output, cov, cats, covVals, covChange, xAxisRange,
                              betas=plotBetas)
      raTitlePrefix <- if (plotBetas) 
                          "$\\Delta_{v}\\mathrm{RA}_{cp}(t)$ of " else ""
      betaTitlePrefix <- "$\\beta_{cp}(t)$ of "
      raYlabel <- "$\\Delta_{v}\\mathrm{RA}_{cp}(t)$"
      betaYlabel <- "$\\beta_{cp}(t)$"
      changeIn <- TRUE
      if (input$changeRange && 
          !is.character(output$centerScaleList[[covIdx]])) {
        cbarLabel <- "$\\Delta_{v} \\mathrm{RA}_{cp}(t)$"
        ylabel <- input$yLabel
        fig <- heatmapPlots(dataDict, suptitle, heatmapIndex, inCI, 
                            figsize=figsize, fontsize=fontSize, xlabel=xLabel, 
                            ylabel=ylabel, cmap='PuOr', vmin=0, vcenter=1, 
                            vmax=2, cbarLabel=cbarLabel, 
                            suptitleLoc=suptitleLoc, titleSize=titleSize, 
                            numXticks=byX, numYticks=by)
      } else {
        fig <- curvePlots(dataDict, suptitle, raTitlePrefix, betaTitlePrefix, 
                          raYlabel, betaYlabel, xLabel, changeIn, plotBetas, 
                          lineWidth=lineWidth, byRA=by,  byBeta=byBeta, 
                          fontsize=fontSize, titleSize=titleSize, 
                          figsize=figsize, xAxisRange=xAxisRange, 
                          figRowCol=figRowCol,  legendLoc=legendLoc, 
                          suptitleLoc=suptitleLoc, sameYaxis=sameYaxis)
      }
    } else if (pType == "Diversity") {
      plotData <- diversityData(output, covVals, xAxisRange, l = lVal)
      confInt <- t(apply(plotData$div, 1,
                         function(x) quantile(x, c(0.025, 0.975))))
      dataDict = list()
      dataDict[["div"]] <- cbind(confInt, rowMeans(plotData$div), 
                                plotData$testpoints)
      dataDict[["div"]] <- data.frame(dataDict[["div"]])
      
      ylabel <- "$\\alpha(t)$"
      changeIn <- FALSE
      fig <- divCurvePlots(dataDict, suptitle, ylabel, xLabel, changeIn, 
                           lineWidth=lineWidth, byY=by, fontsize=fontSize, 
                           titleSize=titleSize, figsize=figsize, 
                           xAxisRange=xAxisRange)
    } else if (pType == "Multiplicative Change in Diversity") {
      dataDict <- deltaDivData(output, cov, covVals, covChange, lVal, 
                               xAxisRange, betas = plotBetas)
      ylabel <- "$\\Delta_{v}\\alpha_{p}(t)$"
      changeIn <- TRUE
      if (input$changeRange && 
          !is.character(output$centerScaleList[[covIdx]])) {
        cbarLabel <- "$\\Delta_{v}\\alpha_{p}(t)$"
        ylabel <- input$yLabel
        fig <- divHeatmapPlots(dataDict, suptitle, heatmapIndex, inCI, 
                               figsize=figsize, fontsize=fontSize, 
                               xlabel=xLabel, ylabel=ylabel, cmap='PuOr', 
                               vmin=0, vcenter=1, vmax=2, cbarLabel=cbarLabel, 
                               suptitleLoc=suptitleLoc, titleSize=titleSize,  
                               numXticks=byX, numYticks=by)
      } else {
        fig <- divCurvePlots(dataDict, suptitle, ylabel, xLabel, changeIn, 
                             multiChange=TRUE, lineWidth=lineWidth, byY=by, 
                             fontsize=fontSize, titleSize=titleSize,  
                             figsize=figsize, xAxisRange=xAxisRange, 
                             legendLoc=legendLoc, suptitleLoc=suptitleLoc)
      }
    }

    tmpfile <- tempfile(fileext = ".png")
    fig$savefig(tmpfile, dpi = 96, bbox_inches = "tight")
    list(path = tmpfile, data = dataDict) 
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(gsub("\\s+", "_", input$plotType), "_", Sys.Date(), ".png")
    },
    contentType = "image/png",
    content = function(file) {
      req(res  <- plotRequest()) # make sure a plot exists (Create Plot clicked)
      file.copy(res$path, file, overwrite = TRUE)
    }
  )

  output$downloadPlotData <- downloadHandler(
    filename = function() {
      paste0(gsub("\\s+", "_", input$plotType), "_dataDict_", Sys.Date(),
             ".rds")
    },
    contentType = "application/octet-stream",
    content = function(file) {
      req(res <- plotRequest())
      saveRDS(res$data, file)
    }
  )

  output$pyPlot <- renderImage({
    # only triggered by plotRequest()
    req(tmpfile <- plotRequest()$path)
    px_w <- input$plotWidth*96  
    px_h <- input$plotHeight*96
    list(src = tmpfile, contentType = "image/png", width = px_w, height = px_h, 
         alt = "FunC-ZIDM plot")
    }, deleteFile = FALSE)
}

shinyApp(ui, server)
