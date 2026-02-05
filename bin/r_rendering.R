#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(knitr)
})

# Parse command-line arguments manually
args <- commandArgs(trailingOnly = TRUE)

# Function to display help
show_help <- function() {
  cat("Process and format MultiQC general stats table\n\n")
  cat("Usage: r_rendering.R -i INPUT [-o OUTPUT] [-f FORMAT]\n\n")
  cat("Options:\n")
  cat("  -i, --input FILE    Input TSV file path (required)\n")
  cat("  -o, --output FILE   Output file path (optional, prints to stdout if not specified)\n")
  cat("  -f, --format FORMAT Output format: 'tsv' or 'markdown' [default: tsv]\n")
  cat("  -h, --help          Show this help message\n\n")
  quit(save = "no", status = 0)
}

# Initialize options with defaults
opt <- list(input = NULL, output = NULL, format = "tsv")

# Parse arguments
i <- 1
while (i <= length(args)) {
  arg <- args[i]
  if (arg %in% c("-h", "--help")) {
    show_help()
  } else if (arg %in% c("-i", "--input")) {
    opt$input <- args[i + 1]
    i <- i + 1
  } else if (arg %in% c("-o", "--output")) {
    opt$output <- args[i + 1]
    i <- i + 1
  } else if (arg %in% c("-f", "--format")) {
    opt$format <- args[i + 1]
    i <- i + 1
  }
  i <- i + 1
}

# Check if input file is provided and exists
if (is.null(opt$input)) {
  cat("Error: Input file is required\n\n")
  show_help()
}

if (!file.exists(opt$input)) {
  stop(sprintf("Error: Input file '%s' not found", opt$input))
}

# Read the TSV file
df <- read.delim(opt$input, check.names = FALSE)

# Clean column names: extract metric name after last dash
# samtools_stats_bbmap_stats-error_rate -> error_rate
colnames(df)[-1] <- colnames(df)[-1] |> 
  stringr::str_extract("[^-]+$")

# clean sample name to remove suffix _*_samtoolsstats
df$Sample <- df$Sample |> stringr::str_remove_all("_\\d+_samtoolsstats")

# sample name as row name
rownames(df) <- df$Sample

# remove Sample column and clean up the column names
tableout <- cbind(ID = rownames(df), stack(df[-1])) |> 
  transform(ind = as.character(ind) |> stringr::str_remove_all("\\.\\d+"))

# remove na values
tableout <- tableout[!is.na(tableout$values),]
# remove . values
tableout$values <- tableout$values |> stringr::str_remove_all("^\\.$")

# pivot data
tableout <- tableout |> pivot_wider(id_cols = ID , names_from = ind, values_from = values, 
              values_fn = \(x) paste(unique(x), collapse = ""))

# round each value to 4 decimals
tableout <- tableout |> mutate(across(-ID, ~round(as.numeric(.), 4)))

# Output results
if (tolower(opt$format) == "markdown") {
  # Markdown format
  if (!is.null(opt$output)) {
    output_table <- knitr::kable(tableout, format = "markdown", align = 'r')
    writeLines(output_table, con = opt$output)
    cat(sprintf("Output written to: %s\n", opt$output))
  } else {
    cat(knitr::kable(tableout, format = "markdown", align = 'r'), sep = "\n")
  }
} else {
  # TSV format (default)
  if (!is.null(opt$output)) {
    write.table(tableout, file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("Output written to: %s\n", opt$output))
  } else {
    write.table(tableout, file = stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}