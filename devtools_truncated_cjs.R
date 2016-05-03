library(devtools)

packdir <- "TruncatedCJS"

## Create new package (First time only)
## devtools::setup(packdir)

## Reload files from package directory
packdir <- "TruncatedCJS"
devtools::reload(packdir)

## Build documentation
devtools::document(packdir)

## Build package
devtools::build(packdir)

## Install package
devtools::install(packdir)
