# Max upload size
options(shiny.maxRequestSize = 600 * 1024^2)
suppressPackageStartupMessages(library(kableExtra))

# Define server

# Function to generate n distinct random colors
generate_random_colors <- function(n) {
  colors <- grDevices::colors()[sample(1:657, n)]
  return(colors)
}

# Helper function to check if a variable is categorical (suitable for factor analysis)
isCategoricalFactor <- function(factor_data, sample_count) {
    # Remove NA values for analysis
    factor_data_clean <- factor_data[!is.na(factor_data)]
    
    if (length(factor_data_clean) == 0) {
        return(FALSE)  # All NA values
    }
    
    # Explicit factor or character variables are categorical
    if (is.factor(factor_data) || is.character(factor_data)) {
        unique_count <- length(unique(factor_data_clean))
        # Must have at least 2 levels and no more than half the sample size
        return(unique_count >= 2 && unique_count <= max(2, sample_count * 0.5))
    }
    
    # For numeric variables, be very strict
    if (is.numeric(factor_data)) {
        unique_count <- length(unique(factor_data_clean))
        # Only consider numeric as categorical if:
        # 1. Very few unique values (≤ 5)
        # 2. At least 2 levels for comparison
        # 3. No more than 20% of sample size
        return(unique_count >= 2 && unique_count <= min(5, sample_count * 0.2))
    }
    
    return(FALSE)
}

# Helper function to get valid categorical factors from design formula
getValidCategoricalFactors <- function(factorChoices, design_terms, metadata_df) {
    # First filter by design formula presence
    validFactorChoices <- factorChoices[factorChoices %in% design_terms]
    
    # Then check if they are actually categorical factors
    if (length(validFactorChoices) > 0) {
        actualFactors <- c()
        for (factor_name in validFactorChoices) {
            if (factor_name %in% colnames(metadata_df)) {
                factor_data <- metadata_df[, factor_name]
                factor_data_clean <- factor_data[!is.na(factor_data)]
                unique_count <- length(unique(factor_data_clean))
                data_type <- class(factor_data)[1]
                
                if (isCategoricalFactor(factor_data, nrow(metadata_df))) {
                    print(paste("✓ Factor", factor_name, "is categorical:", data_type, "with", unique_count, "levels"))
                    actualFactors <- c(actualFactors, factor_name)
                } else {
                    print(paste("✗ Skipping", factor_name, "-", data_type, "with", unique_count, "unique values (not suitable for categorical analysis)"))
                }
            }
        }
        validFactorChoices <- actualFactors
    }
    
    return(validFactorChoices)
}


server <- function(input, output, session) {
    source("server-inputdata.R", local = TRUE)

    source("server-conditions.R", local = TRUE)

    source("server-svaseq.R", local = TRUE)

    source("server-runDeseq.R", local = TRUE)

    source("server-analysisRes.R", local = TRUE)
    source("server-venndiagram.R", local = TRUE)
    source("server-volcanoplot.R", local = TRUE)

    source("server-boxplot.R", local = TRUE)

    source("server-heatmap.R", local = TRUE)

    GotoTab <- function(name) {
        shinyjs::show(selector = paste0("a[data-value=\"", name, "\"]"))

        shinyjs::runjs("window.scrollTo(0, 0)")
    }
}
