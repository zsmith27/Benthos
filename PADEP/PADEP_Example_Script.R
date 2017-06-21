#==============================================================================
#==============================================================================
# Title: PADEP 
# Author: Zachary M. Smith
# Organization: ICPRB
# Email: zsmith@icprb.org
# Date: 2/09/2017
# Purpose: Script written to help PADEP calculate metrics.
#==============================================================================
#==============================================================================
# Install the package devtools, which makes it easy to install the Benthos
# package.
install.packages(c("devtools","curl", "httr"))
# Install the Benthos package.
devtools::install_github("zsmith27/Benthos", build_vignettes = TRUE)
# Load the Benthos package.
library(Benthos)
# Load the Master Taxa List contained within the Benthos package.
data(master)
vignette(package = "Benthos")
#==============================================================================
# Set your working directory.
# This is where your files are stored on your computer.
#setwd("C:/Users/zsmith/Desktop/DEC_LAKE")
# Import csv containing taxonomic counts.
#dec <- read.csv("MasterSpeciesTable_R.csv", stringsAsFactors = FALSE)
#==============================================================================
# Import PADEP Master Taxa List.
setwd("C:/Users/zsmith/Desktop/Benthos_R/master/PADEP")
master.pa <- read.csv("PADEP_Original.csv", stringsAsFactors = FALSE)
names(master.pa)[names(master.pa) %in% "STR_TAXA_NAME"] <- "AGENCY_ID"
master.pa$AGENCY_ID <- toupper(master.pa$AGENCY_ID)

test <- merge(master.pa, master, by = "AGENCY_ID", all.x = TRUE)




