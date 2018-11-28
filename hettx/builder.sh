#!/bin/bash

# Build documentation
R -e 'devtools::document()'
cd ..

# Build and run CRAN checks
R CMD BUILD hettx --resave-data 
R CMD CHECK hettx_*.tar.gz --as-cran
R CMD INSTALL hettx_*.tar.gz
