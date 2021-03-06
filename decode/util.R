#!/usr/bin/Rscript
#
# Common utility library for all R scripts.

# Log message with timing.  Example:
#
# _____ 1.301 My message
#
# The prefix makes it stand out (vs R's print()), and the number is the time so far.
#
# NOTE: The shell script log uses hyphens.

Log <- function(...) {
  cat(sprintf('_____ %.3f ', proc.time()[['elapsed']]))
  cat(sprintf(...))
  cat('\n')
}
