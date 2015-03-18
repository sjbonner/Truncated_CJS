library(devtools)

packdir <- "TruncatedCJS"

## Create new package (First time only)
## devtools::setup(packdir)

## Build documentation
devtools::document(packdir)

## Build package
devtools::build(packdir)

## Install package
devtools::install(packdir)
