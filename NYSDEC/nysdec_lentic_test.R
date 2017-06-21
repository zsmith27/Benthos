#==============================================================================
#==============================================================================
# Title: NYSDEC Lentic Metrics
# Author: Zachary M. Smith
# Organization: ICPRB
# Email: zsmith@icprb.org
# Date: 2/09/2017
# Purpose: Script written to help E. Stieber calculate benthic
# macroinvertebrate metrics for NYS lentic data.
#==============================================================================
#==============================================================================
# Install the package devtools, which makes it easy to install the Benthos
# package.
install.packages(c("devtools","curl", "httr"))
# Install the Benthos package.
devtools::install_github("zsmith27/Benthos",  force = TRUE)
# Load the Benthos package.
library(Benthos)
# Load the Master Taxa List contained within the Benthos package.
data(master)
#vignette(package = "Benthos")
#==============================================================================
# Set your working directory.
# This is where your files are stored on your computer.
#setwd("C:/Users/zsmith/Desktop/DEC_LAKE")
setwd("C:/Users/Owner/Desktop/NYSDEC/NYSDEC_LAKE")
# Import csv containing taxonomic counts.
dec <- read.csv("MasterSpeciesTable.csv", stringsAsFactors = FALSE)
#==============================================================================
# Format data to meet the formating requirments of the Benthos package.
long <- nysdec_data_prep_old(dec, master)
# We need to updated the master taxa list to include the missing taxa indicated
# by the ISSUE column.

# For now remove the rows that are causing issues to test calculating metrics.
long.sub <- long[long$ISSUE %in% "NO", ]
# If you look at long you will see NAs in the taxonomic hierarchy.
# These values are missing because the taxon was not identified beyond
# a specific rank. For example, if the taxa was identified to the family-level
# then it is not possible to assign a genus or species name.
# The names might also be missing because they are not specified on ITIS.gov.
# The fill_taxa function fills in the NAs with proceeding taxonomic level.
long.fill <- fill_taxa(long.sub)
#==============================================================================
# Calculate all of the available metrics.
metrics.df <- all_metrics(long.fill, master, "SPECIES", 
                          tv.col = "NYSDEC_TV",
                          ffg.col = "NYSDEC_FFG")
#==============================================================================
# Change your working directory if you want the calculated metrics to be
# stored in a different file location then where the taxonomic count file
# was stored.
setwd("C:/Users/Owner/Desktop/NYSDEC/NYSDEC_LAKE")
# Export as a csv file.
write.csv(metrics.df, "my_metrics.csv", row.names = FALSE)
#==============================================================================
#==============================================================================
# Use the data frame below to compare to your master taxa list.
master.check <- unique(master[c("FINAL_ID", "NYSDEC_TV", "NYSDEC_FFG", 
                                "NYSDEC_NBI.P", "NYSDEC_NBI.N")])

long.sub <- data.frame(unique(long[long$ISSUE %in% "YES", "AGENCY_ID"]))

test <- seq_taxa_rich(long.fill, rank = "SPECIES", master)
names(test)[grepl("\\.", names(test))]


