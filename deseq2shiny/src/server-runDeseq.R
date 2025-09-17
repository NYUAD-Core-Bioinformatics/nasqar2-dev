output$ddsColData <- renderTable(
    {
        # print(myValues)
        # print(myValues$dds)
        print("ddsColData")
        print(colData(myValues$dds))
        colDataTable <- colData(myValues$dds)
        colsToSelect <- labels(terms(design(myValues$dds)))
        print(colsToSelect)
        print(length(colsToSelect))

        if (length(colsToSelect) == 0) {
            return(head(colDataTable, n = 5))
        }
        if (length(colsToSelect) == 1) {
            print(colsToSelect[[1]])
            h <- head(colDataTable[, colsToSelect], n = 5)
            h <- as.data.frame(h)
            colnames(h) <- c(colsToSelect[[1]])
            return(h)
        }

        head(colDataTable[, colsToSelect], n = 5)
    },
    bordered = TRUE,
    spacing = "xs",
    rownames = T
)

output$ddsDesignFormula <- renderText({
    dds <- myValues$dds
    paste(as.character(design(dds)))
})

observe({
    ddsReactive()
})

ddsReactive <- eventReactive(input$run_deseq2, {
    # print('ddsReactive')
    dds <- myValues$dds # ()

    withProgress(message = "Running DESeq2 , please wait", {
        js$addStatusIcon("deseqTab", "loading")
        shiny::setProgress(value = 0.4, detail = " ...")

        removeNotification("errorNotify")
        removeNotification("errorNotify1")
        # nasqar use 3 cores
        BiocParallel::register(MulticoreParam(3))

        validate(need(
            tryCatch(
                {
                    dds <- DESeq(dds, parallel = T)
                },
                error = function(e) {
                    myValues$status <- paste("DESeq2 Error: ", e$message)

                    showNotification(id = "errorNotify", myValues$status, type = "error", duration = NULL)
                    showNotification(id = "errorNotify1", "If this is intended, please select 'No Replicates' in Input Data step. OR use ~ 1 as the design formula", type = "error", duration = NULL)


                    js$addStatusIcon("deseqTab", "fail")

                    return(NULL)
                }
            ),
            "Error"
        ))

        BiocParallel::register(SerialParam())


        if (input$computeRlog) {
            shiny::setProgress(value = 0.5, detail = "Calculating RLog transformation ...")
            rld <- rlog(dds)
            myValues$rld <- rld
            myValues$rlogMat <- assay(rld)

            myValues$rldColNames <- colnames(rld)
        }


        shiny::setProgress(value = 0.7, detail = "Computing Variance Stabilizing Transformation ...")

        vsd <- varianceStabilizingTransformation(dds)
        myValues$vsd <- vsd
        myValues$vstMat <- assay(vsd)
        myValues$vsdColNames <- colnames(vsd)

        shiny::setProgress(value = 0.8, detail = "Formatting data ...")

        shiny::setProgress(value = 1, detail = "...")

        js$addStatusIcon("deseqTab", "done")

        myValues$dds <- dds

        shinyjs::show(selector = "a[data-value=\"boxplotTab\"]")
        shinyjs::show(selector = "a[data-value=\"resultsTab\"]")
        shinyjs::show(selector = "a[data-value=\"heatmapTab\"]")
        shinyjs::show(selector = "a[data-value=\"vstTab\"]")

        if (input$computeRlog) {
            shinyjs::show(selector = "a[data-value=\"rlogTab\"]")
        } else {
            shinyjs::hide(selector = "a[data-value=\"rlogTab\"]")
        }


        factorChoices <- colnames(colData(dds))

        # print('factorChoices')
        # print(factorChoices)
        factorChoices <- factorChoices[!grepl("^SV[::digit::]*", factorChoices)]

        updateSelectInput(session, "rlogIntGroupsInput", choices = factorChoices, selected = factorChoices[1])
        updateSelectInput(session, "vsdIntGroupsInput", choices = factorChoices, selected = factorChoices[1])

        factorChoices <- factorChoices[!(factorChoices %in% c("sizeFactor", "replaceable"))]

        # Filter factorChoices to only include categorical factors that are in the design formula
        design_terms <- labels(terms(design(dds)))
        print("Design formula terms:")
        print(design_terms)
        print("Available factor choices:")
        print(factorChoices)
        
        validFactorChoices <- getValidCategoricalFactors(factorChoices, design_terms, myValues$DF)
        
        print("Valid categorical factor choices:")
        print(validFactorChoices)
        
        # If no valid factors found, keep original list but add warning
        if (length(validFactorChoices) == 0) {
            print("Warning: No categorical factors from design formula found for differential analysis")
            validFactorChoices <- factorChoices
        }

        updateSelectizeInput(session, "resultNamesInput", choices = resultsNames(dds), selected = NULL)
        updateSelectizeInput(session, "factorNameInput", choices = validFactorChoices, selected = if(length(validFactorChoices) > 0) validFactorChoices[1] else NULL)

        # Update condition dropdowns now that DESeq2 analysis is complete
        if (!is.null(input$factorNameInput) && input$factorNameInput != "" && 
            input$factorNameInput %in% validFactorChoices) {
            factor_levels <- levels(myValues$DF[, input$factorNameInput])
            if (length(factor_levels) >= 2) {
                updateSelectInput(session, "condition1", choices = factor_levels)
                updateSelectInput(session, "condition2", choices = factor_levels)
            }
        }

        disable("data_file_type")
        disable("no_replicates")

        js$addStatusIcon("conditionsTab", "done")
    })

    # }
})

observeEvent(input$factorNameInput,
    {
        print("factorNameInput")
        print(input$factorNameInput)
        print(myValues$DF)
        
        # Validate that the selected factor still exists in myValues$DF
        if (is.null(myValues$DF) || is.null(input$factorNameInput) || input$factorNameInput == "") {
            return()
        }
        
        # Check if the selected column exists in the dataframe
        if (!(input$factorNameInput %in% colnames(myValues$DF))) {
            print(paste("Warning: Column", input$factorNameInput, "not found in myValues$DF"))
            return()
        }
        
        # Only populate conditions if DESeq2 analysis has been completed and results exist
        if (is.null(myValues$dds)) {
            print("DESeq2 analysis not completed - conditions will be populated after analysis")
            updateSelectInput(session, "condition1", choices = NULL)
            updateSelectInput(session, "condition2", choices = NULL)
            return()
        }
        
        myValues$DF[] <- lapply(myValues$DF, function(x) {
            if (is.character(x) || is.factor(x)) {
                as.factor(x)
            } else {
                x
            }
        })
        
        # Verify that the factor is present in the design formula and results
        design_terms <- labels(terms(design(myValues$dds)))
        if (!(input$factorNameInput %in% design_terms)) {
            print(paste("Factor", input$factorNameInput, "not found in design formula"))
            updateSelectInput(session, "condition1", choices = NULL)
            updateSelectInput(session, "condition2", choices = NULL)
            return()
        }
        
        # Verify that the factor is actually categorical (not continuous)
        factor_data <- myValues$DF[, input$factorNameInput]
        if (!isCategoricalFactor(factor_data, nrow(myValues$DF))) {
            print(paste("Factor", input$factorNameInput, "appears to be continuous, not suitable for categorical comparison"))
            updateSelectInput(session, "condition1", choices = NULL)
            updateSelectInput(session, "condition2", choices = NULL)
            return()
        }
        
        # print(levels(factor(myValues$DF[,input$factorNameInput])))
        print("factorNameInput changed")
        
        # Get factor levels properly, converting to factor if necessary
        factor_data <- myValues$DF[, input$factorNameInput]
        
        # Convert to factor if it's not already one
        if (!is.factor(factor_data)) {
            factor_data <- as.factor(factor_data)
        }
        
        factor_levels <- levels(factor_data)
        print(paste("Factor levels:", paste(factor_levels, collapse = ", ")))
        
        # Only populate if there are at least 2 levels for comparison
        if (length(factor_levels) < 2) {
            print(paste("Factor", input$factorNameInput, "must have at least 2 levels for differential analysis"))
            updateSelectInput(session, "condition1", choices = NULL)
            updateSelectInput(session, "condition2", choices = NULL)
            return()
        }
        
        updateSelectInput(session, "condition1", choices = factor_levels)
        updateSelectInput(session, "condition2", choices = factor_levels)
    },
    ignoreInit = T
)


observe({
    if (input$goto_svaTab > 0) {
        GotoTab("svaseqTab")
    }
})

output$rlogData <- renderDataTable(
    {
        if (!is.null(myValues$rlogMat)) {
            myValues$rlogMat
        }
    },
    options = list(scrollX = TRUE, pageLength = 5)
)

output$rlogPlot <- renderPlotly({
    if (!is.null(myValues$rlogMat)) {
        sampleDists <- dist(t(myValues$rlogMat))
        print("rlogPlot")
        print(myValues$rldColNames)
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- paste(myValues$rldColNames)
        colnames(sampleDistMatrix) <- NULL
        colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
        # rlogHeat = pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
        # browser()
        # rlogHeat
        heatmaply::heatmaply(sampleDistMatrix, showticklabels = c(F, T))
    }
})

output$rlogPcaPlot <- renderPlotly({
    if (!is.null(myValues$rld)) {
        intgroups <- input$rlogIntGroupsInput
        if (is.null(intgroups) | intgroups == "") {
            intgroups <- names(colData(myValues$dds))[1]
        }



        # Use DESeq2's plotPCA and returnData=TRUE to get the data
        rld <- myValues$rld
        
        # Get PCA data from DESeq2
        pca_data <- DESeq2::plotPCA(rld, intgroup = intgroups, returnData = TRUE)
        percentVar <- attr(pca_data, "percentVar")
        
        # Add sample names for hover
        pca_data$SampleName <- rownames(pca_data)
        
        # Create group column for coloring and display
        if (length(intgroups) == 1) {
            pca_data$group <- pca_data[[intgroups]]
        } else {
            # If multiple groups, combine them
            pca_data$group <- pca_data$group  # DESeq2 already creates this
        }
        
        # Prepare hover text with all available information
        hover_text <- paste(
            "Sample:", pca_data$SampleName,
            "<br>Group of Interest:", pca_data$group
        )
        
        # Add additional columns to hover if they exist
        sample_data <- as.data.frame(colData(rld))
        additional_cols <- setdiff(names(sample_data), c(intgroups, "sizeFactor", "replaceable"))
        if (length(additional_cols) > 0) {
            for (col in additional_cols) {
                if (!is.null(sample_data[[col]])) {
                    hover_text <- paste(hover_text, 
                                      paste0("<br>", col, ":"), sample_data[[col]])
                }
            }
        }
        
        # Create plotly visualization
        p <- plot_ly(
            data = pca_data,
            x = ~PC1, 
            y = ~PC2,
            color = ~group,
            text = hover_text,
            hoverinfo = "text",
            type = "scatter",
            mode = "markers",
            marker = list(size = 8, opacity = 0.8)
        ) %>%
        layout(
            title = "PCA Plot",
            xaxis = list(title = paste0("PC1: ", round(percentVar[1] * 100, 1), "% variance")),
            yaxis = list(title = paste0("PC2: ", round(percentVar[2] * 100, 1), "% variance")),
            showlegend = TRUE,
            legend = list(title = list(text = "Group of Interest"))
        )
        
        return(p)
    }
})


output$vsdPlot <- renderPlotly({
    if (!is.null(myValues$vstMat)) {
        print("vsdPlot")
        print(myValues$rldColNames)
        print(myValues$vsdColNames)
        sampleDists <- dist(t(myValues$vstMat))
        sampleDistMatrix <- as.matrix(sampleDists)
        rownames(sampleDistMatrix) <- paste(myValues$vsdColNames)
        colnames(sampleDistMatrix) <- paste(myValues$vsdColNames)
        colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
        # vstHeat <- pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors)
        # vstHeat
        # heatmaply::heatmaply(sampleDistMatrix, showticklabels = c(F,T))
        # label both axes
        heatmaply::heatmaply(sampleDistMatrix, showticklabels = c(F, T))
    }
})

output$vsdPcaPlot <- renderPlotly({
    if (!is.null(myValues$vsd)) {
        intgroups <- input$vsdIntGroupsInput
        if (is.null(intgroups) | intgroups == "") {
            intgroups <- names(colData(myValues$dds))[1]
        }

        print("colData")
        print(colData(myValues$dds))


        vsd <- myValues$vsd
        
        # Get PCA data from DESeq2  (default)
        pca_data <- DESeq2::plotPCA(vsd, intgroup = intgroups, returnData = TRUE)
        percentVar <- attr(pca_data, "percentVar")
        
        # Add sample names for hover
        pca_data$SampleName <- rownames(pca_data)
        
        # Create group column for coloring and display
        if (length(intgroups) == 1) {
            pca_data$group <- pca_data[[intgroups]]
        } else {
            # If multiple groups, combine them
            pca_data$group <- pca_data$group  # DESeq2 already creates this
        }
        
        # Prepare hover text with all available information
        hover_text <- paste(
            "Sample:", pca_data$SampleName,
            "<br>Group of Interest:", pca_data$group
        )
        
        # Add additional columns to hover if they exist
        sample_data <- as.data.frame(colData(vsd))
        additional_cols <- setdiff(names(sample_data), c(intgroups, "sizeFactor", "replaceable"))
        if (length(additional_cols) > 0) {
            for (col in additional_cols) {
                if (!is.null(sample_data[[col]])) {
                    hover_text <- paste(hover_text, 
                                      paste0("<br>", col, ":"), sample_data[[col]])
                }
            }
        }
        
        # Create plotly visualization
        p <- plot_ly(
            data = pca_data,
            x = ~PC1, 
            y = ~PC2,
            color = ~group,
            text = hover_text,
            hoverinfo = "text",
            type = "scatter",
            mode = "markers",
            marker = list(size = 8, opacity = 0.8)
        ) %>%
        layout(
            title = "PCA Plot",
            xaxis = list(title = paste0("PC1: ", round(percentVar[1] * 100, 1), "% variance")),
            yaxis = list(title = paste0("PC2: ", round(percentVar[2] * 100, 1), "% variance")),
            showlegend = TRUE,
            legend = list(title = list(text = "Group of Interest"))
        )
        
        return(p)
    }
})

output$vsdData <- renderDataTable(
    {
        if (!is.null(myValues$vstMat)) {
            myValues$vstMat
        }
    },
    options = list(scrollX = TRUE, pageLength = 5)
)

output$downloadRlogCsv <- downloadHandler(
    filename = function() {
        paste0("rlog", ".csv")
    },
    content = function(file) {
        csv <- myValues$rlogMat

        write.csv(csv, file, row.names = T)
    }
)

output$downloadVsdCsv <- downloadHandler(
    filename = function() {
        paste0("vsd", ".csv")
    },
    content = function(file) {
        csv <- myValues$vstMat

        write.csv(csv, file, row.names = T)
    }
)

output$ddsComputed <- reactive({
    return(!is.null(myValues$dds))
})
outputOptions(output, "ddsComputed", suspendWhenHidden = FALSE)

output$deseqError <- reactive({
    return(!is.null(myValues$status))
})
outputOptions(output, "deseqError", suspendWhenHidden = FALSE)
