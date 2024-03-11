library(tidyverse)

#Process Data
source("original/00_PHOS_DataPreparation.R")

#Load Model Code and Functions
source("original/01_PHOS_functions_and_JAGSmodel.R")

## Fix the missing 4 states after capture in the trap (and removal from the system)
fix_trapped_state <- function(input) {
  z <- input$z
  K <- ncol(z)
  for (r in seq_len(nrow(z))) {
    merwin <- which(z[r, ] == 4)
    if (length(merwin) > 0) {
      first_merwin <- min(merwin)
      z[r, first_merwin:K] <- 4
    }
  }
  input$z <- z
  input
}
pHOS.Input <- fix_trapped_state(pHOS.Input)

y <- pHOS.Input$y

tibble(
  Date = seq(as.Date("2016-03-08"), as.Date("2016-05-10"), by = "1 week"),
  Captures = apply(pHOS.Input$y, 2, \(.) sum(. %in% c(2, 3), na.rm = TRUE))[2:11]
)
