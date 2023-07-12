# pkgLoad function adapted from rosscova on Stack Overflow
# https://stackoverflow.com/a/38928678
pkgLoad <- function( packages = "favourites" ) {

    if( length( packages ) == 1L && packages == "favourites" ) {
        packages <- c( "data.table", "chron", "plyr", "dplyr", "shiny",
                       "shinyjs", "parallel", "devtools", "doMC", "utils",
                       "stats", "microbenchmark", "ggplot2", "readxl",
                       "feather", "googlesheets", "readr", "DT", "knitr",
                       "rmarkdown", "Rcpp"
        )
    }

    packagecheck <- match( packages, utils::installed.packages()[,1] )
    packagestoinstall <- packages[ is.na( packagecheck ) ]

    if( length( packagestoinstall ) > 0L ) {
        utils::install.packages(
            packagestoinstall,
            repos = "https://mirror.las.iastate.edu/CRAN/"
        )
    } else {
        print( "All requested packages already installed" )
    }

    for( package in packages ) {
        suppressPackageStartupMessages(
            library( package, character.only = TRUE, quietly = TRUE )
        )
    }
}

# Choice of functions specific to Proximal ID code
pkgLoad(c("R.utils", "jsonlite", "ipw", "cubature", "nnet", "numDeriv", "data.table"))
