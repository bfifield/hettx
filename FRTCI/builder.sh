#!/bin/bash

# Build documentation
R -e 'devtools::document()'
cd ..

# Build and run CRAN checks
R CMD BUILD FRTCI --resave-data 
R CMD CHECK FRTCI_*.tar.gz --as-cran
R CMD INSTALL FRTCI_*.tar.gz
