/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Please provide your markdown file to be compiled and your output file", call = FALSE)
}

rmarkdown::render(input = args[1],
                  output_format = "html_notebook",
                  output_file = args[2],
                  clean = TRUE,
                  run_pandoc = TRUE)
