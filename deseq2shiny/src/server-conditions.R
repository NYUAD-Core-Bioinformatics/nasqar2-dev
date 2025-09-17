observe({
    colnamesChoices <- colnames(myValues$DF[!(names(myValues$DF) %in% c("Samples", "Groups"))])
    updateSelectInput(session, "colToRemove", choices = colnamesChoices, selected = NULL)
    
    # Update sample exclusion choices when DF changes
    if (!is.null(myValues$DF)) {
        sampleChoices <- rownames(myValues$DF)
        updateSelectizeInput(session, "samplesToExclude", choices = sampleChoices, selected = input$samplesToExclude)
    }
})

observeEvent(input$removeCol, {
    validate(
        need(input$colToRemove != "", message = "need to select column to remove")
    )

    if (length(colnames(myValues$DF)) == 1) {
        showNotification(id = "errorNotify1", "At least two column required in the table, to remove one column", type = "error", duration = NULL)
    }
    validate(
        need(length(colnames(myValues$DF)) > 1, message = "At least two column required in the table, to remove one column")
    )
    removeNotification("errorNotify1")
    myValues$DF <- hot_to_r(input$table)
    myValues$DF[, input$colToRemove] <- NULL
    updateDesignFormula()
    
    # Update factorNameInput choices only if DESeq2 analysis is completed
    if (exists("input") && !is.null(input) && "factorNameInput" %in% names(input) && !is.null(myValues$dds)) {
        factorChoices <- colnames(myValues$DF)
        factorChoices <- factorChoices[!(factorChoices %in% c("sizeFactor", "replaceable"))]
        
        # Filter to only include categorical factors that are in the design formula
        design_terms <- labels(terms(design(myValues$dds)))
        validFactorChoices <- getValidCategoricalFactors(factorChoices, design_terms, myValues$DF)
        
        if (length(validFactorChoices) > 0) {
            updateSelectizeInput(session, "factorNameInput", choices = validFactorChoices, selected = validFactorChoices[1])
        }
    }
})

observe({
    tableEditReactive()
})

tableEditReactive <- reactive({
    if (input$addConditions > 0) {
        isolate({
            print("input$table start")
            # print(input$table)
            myValues$DF <- hot_to_r(input$table)
            DF <- myValues$DF

            print("input$table stop")
            print(myValues$DF)


            validate(
                need(!(input$conditionName %in% colnames(DF)), message = "Condition name already exists"),
                need(trimws(input$conditionName) != "", message = "Condition name empty"),
                need(trimws(input$conditions) != "", message = "Conditions empty")
            )

            newDF <- data.frame(newCol = character(dim(DF)[1]))
            names(newDF) <- c(input$conditionName)

            DF <- cbind(DF, newDF)

            myValues$DF <- DF

            print("myValues$conditions")
            print(DF)

            myValues$conditions[[dim(DF)[2] - 1]] <- input$conditions

            print(myValues$conditions)

            updateTextInput(session, "conditionName", value = "")
            updateTextInput(session, "conditions", value = "")
            updateDesignFormula()
        })
    }
})

# output$tableCreated <-
#   reactive({
#     return(!is.null(myValues$DF))
#   })
# outputOptions(output, 'tableCreated', suspendWhenHidden=FALSE)

output$table <- renderRHandsontable({
    DF1 <- myValues$DF
    if (is.null(DF1)) {
        return()
    }


    table <- rhandsontable(DF1, rowHeaderWidth = 100) %>%
        hot_cols(colWidths = 150)

    table <- table %>% hot_table(highlightCol = TRUE, highlightRow = TRUE, colHeaders = NULL)

    for (i in 2:dim(DF1)[2])
    {
        if (dim(DF1)[2] < 2) {
            break()
        }
        if (!is.null(myValues$conditions[[i - 1]])) {
            table <- table %>% hot_col(col = colnames(DF1)[i], type = "dropdown", source = getConditionsListFromStr(myValues$conditions[[i - 1]]))
        }
    }

    table
})

output$downloadCSV <- downloadHandler(
    filename = paste0("metadatatable_", format(Sys.time(), "%y-%m-%d_%H-%M-%S"), ".csv"),
    content = function(file) {
        write.csv(hot_to_r(input$table), file, row.names = T)
    }
)

observe({
    metadataFileReactive()
})
metadataFileReactive <- reactive({
    # Check if example selected, or if not then ask to upload a file.
    shiny::validate(
        need((!is.null(input$metadatafile)),
            message = "Please select a file"
        )
    )

    inFile <- input$metadatafile
    if (is.null(inFile)) {
        return(NULL)
    }

    inFile <- inFile$datapath

    # print('metadataFileReactive')
    # print(inFile)
    # select file separator, either tab or comma
    sep <- "\t"
    if (length(inFile) > 0) {
        testSep <- read.csv(inFile[1], header = TRUE, sep = "\t")
        if (ncol(testSep) < 2) {
            sep <- ","
        }
    } else {
        return(NULL)
    }

    fileContent <- read.csv(inFile[1], header = TRUE, sep = sep)

    sampleN <- colnames(fileContent)[-1]
    print("sampleN")

    # metaData <- fileContent[,sampleN]
    # if meta dat has only one entry returns as vector.
    metaData <- data.frame(fileContent[, colnames(fileContent)[-1]])
    colnames(metaData) <- colnames(fileContent)[-1]

    metaData[] <- lapply(metaData, function(x) {
        if (is.character(x) || is.factor(x)) {
            as.factor(x)
        } else {
            x
        }
    })
    print("metaData")
    print(fileContent[, 1])
    print("metaData")
    # metaData <- data.frame(sapply( metaData, as.factor ))
    # print(levels(metaData$Stage))
    row.names(metaData) <- fileContent[, 1]

    myValues$DF <- metaData

    updateDesignFormula()

    return(metaData)
})


observe({
    ddsInitReactive()
})

# Handle clearing excluded samples
observeEvent(input$clearExcludedSamples, {
    updateSelectizeInput(session, "samplesToExclude", selected = character(0))
})

# Initialize excluded samples reactive value
if (!exists("myValues$excludedSamples")) {
    myValues$excludedSamples <- character(0)
}

# Update excluded samples when selection changes
observeEvent(input$samplesToExclude, {
    myValues$excludedSamples <- input$samplesToExclude
}, ignoreNULL = FALSE)


ddsInitReactive <- eventReactive(input$init_deseq2, {
    withProgress(message = "Initializing DESeq2 ...", {
        # print("Initializing DESeq2 ...")
        removeNotification("errorNotify")
        removeNotification("errorNotify1")
        # shinyjs::hide(selector = "a[data-value=\"deseqTab\"]")
        shinyjs::hide(selector = "a[data-value=\"rlogTab\"]")
        shinyjs::hide(selector = "a[data-value=\"vstTab\"]")

        shinyjs::hide(selector = "a[data-value=\"resultsTab\"]")
        shinyjs::hide(selector = "a[data-value=\"boxplotTab\"]")
        shinyjs::hide(selector = "a[data-value=\"heatmapTab\"]")

        myValues$status <- NULL
        myValues$dds <- NULL

        js$addStatusIcon("conditionsTab", "loading")
        shiny::setProgress(value = 0.2, detail = "...")


        myValues$DF <- hot_to_r(input$table)

        samples <- myValues$DF
        dataCounts <- myValues$dataCounts

        # print('dataCounts...')

        # print(rownames(myValues$dataCounts))

        # Set up sample names properly first
        rownames(samples) <- rownames(myValues$DF)

        # Exclude selected samples after sample names are properly set
        excludedSamples <- input$samplesToExclude
        if (!is.null(excludedSamples) && length(excludedSamples) > 0) {
            # Filter out excluded samples from both samples and dataCounts
            samplesToKeep <- !rownames(samples) %in% excludedSamples
            samples <- samples[samplesToKeep, , drop = FALSE]
            
            # Filter dataCounts by column names (sample names)
            columnsToKeep <- !colnames(dataCounts) %in% excludedSamples
            dataCounts <- dataCounts[, columnsToKeep, drop = FALSE]
            
            # Store excluded samples info for display
            myValues$excludedSamplesInfo <- excludedSamples
            
            print(paste("Excluded samples:", paste(excludedSamples, collapse = ", ")))
            print(paste("Remaining samples:", paste(rownames(samples), collapse = ", ")))
        } else {
            myValues$excludedSamplesInfo <- character(0)
        }

        # convert factors to unordered
        # factor(samples, ordered = F)

        for (i in 1:ncol(samples)) {
            if (all(class(samples[, i]) %in% c("ordered", "factor"))) {
                samples[, i] <- factor(samples[, i], ordered = F)
            }
        }

        print("samples")
        print(samples)

        isolate({
            validate(need(
                tryCatch(
                    {
                        dds <- DESeqDataSetFromMatrix(dataCounts, colData = samples, design = as.formula(input$designFormula))
                    },
                    error = function(e) {
                        myValues$status <- paste("DESeq2 Error: ", e$message)

                        showNotification(id = "errorNotify", myValues$status, type = "error", duration = NULL)
                        showNotification(id = "errorNotify1", "Fix design formula OR Factors/Conditions", type = "error", duration = NULL)


                        js$addStatusIcon("conditionsTab", "fail")
                    }
                ),
                "Error"
            ))
        })

        updateTextInput(session, "designFormulaSva", value = paste(as.character(design(dds)), collapse = " "))

        js$addStatusIcon("conditionsTab", "done")
        shinyjs::show(selector = "a[data-value=\"deseqTab\"]")
        updateTabItems(session, "tabs", "deseqTab")

        shiny::setProgress(value = 0.8, detail = "...")
        myValues$dds <- dds
        return(dds)
    })
})

output$ddsInitAvailable <- reactive({
    return(!is.null(ddsInitReactive()))
})
outputOptions(output, "ddsInitAvailable", suspendWhenHidden = FALSE)

# Display excluded samples
output$excludedSamplesList <- renderText({
    if (!is.null(myValues$excludedSamplesInfo) && length(myValues$excludedSamplesInfo) > 0) {
        paste(myValues$excludedSamplesInfo, collapse = ", ")
    } else {
        "None"
    }
})

# For the DESeq tab display
output$excludedSamplesList2 <- renderText({
    if (!is.null(myValues$excludedSamplesInfo) && length(myValues$excludedSamplesInfo) > 0) {
        paste(myValues$excludedSamplesInfo, collapse = ", ")
    } else {
        ""
    }
})

output$hasExcludedSamples <- reactive({
    return(!is.null(myValues$excludedSamplesInfo) && length(myValues$excludedSamplesInfo) > 0)
})
outputOptions(output, "hasExcludedSamples", suspendWhenHidden = FALSE)





updateDesignFormula <- function() {
    isolate({
        print("updateDesignFormula")
        
        # Check if myValues$DF exists and has columns
        if (is.null(myValues$DF) || ncol(myValues$DF) == 0) {
            return()
        }
        
        groupvars <- colnames(myValues$DF)

        if (length(groupvars) == 1) {
            # Validate that the column still exists
            if (groupvars[1] %in% colnames(myValues$DF) && 
                length(levels(myValues$DF[, groupvars])) == nrow(myValues$DF)) {
                designFormula <- "~ 1"
            } else {
                designFormula <- paste("~ ", groupvars)
            }
        } else {
            designFormula <- paste(" ~ ", paste(groupvars, collapse = " + "))
        }

        print(designFormula)
        updateTextInput(session, "designFormula", value = designFormula)
        updateTextInput(session, "designFormulaSva", value = designFormula)
    })
}


getConditionsListFromStr <- function(conditonsStr) {
    conditions <- isolate(unlist(strsplit(conditonsStr, ",")))
    conditions <- trimws(conditions)
    conditions <- conditions[conditions != ""]
    conditions <- unique(conditions)
    return(conditions)
}
