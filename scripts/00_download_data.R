# ----------------------------------------------------------
# Download and prepare microarray data for E-MTAB-1620
# Source: https://www.ebi.ac.uk/biostudies/ArrayExpress/studies/E-MTAB-1620
# ----------------------------------------------------------

# The dataset corresponds to:
# “Morphological and Molecular Characterization of Adult Midgut
#  Compartmentalization in Drosophila” (Buchon et al., 2013)
#  Platform: Affymetrix Drosophila Genome 2.0 Array

# 1. Download the ZIP file manually from the ArrayExpress link above.
#    File name: E-MTAB-1620.zip
# 2. Unzip and place it in the `raw_data/` folder of this repository.

file_dir <- "raw_data/E-MTAB-1620/"

cel_files <- list.files(file_dir, pattern = "\\.CEL$")
cel_files
length(cel_files)
