#!/usr/bin/env Rscript
#
# Command line tool to decode a RAPPOR data set.  It is a simple wrapper for
# Decode() in decode.R.

library(optparse)
library(RJSONIO)
setwd("./decode/")
source("read_input.R")
source("decode.R")
source("util.R")
source("alternative.R")


Decode_main <- function(params_file, map_file, counts_file, out_result) {
    params <- ReadParameterFile(params_file)
    counts <- ReadCountsFile(counts_file, params, adjust_counts = TRUE)
    map <- LoadMapFile(map_file, params)

    res <- Decode(counts, map$map, params, correction = "FDR", alpha = .05)
    
    if (nrow(res$fit) == 0) {
        Log("FATAL: Analysis returned no strings.")
        quit(status = 1)
    }
    
    # Write analysis results as CSV.
    write.csv(res$fit, file = out_result, row.names = FALSE)

}
