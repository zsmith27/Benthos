#==============================================================================
#==============================================================================
# Title: NYSDEC Lentic Metrics
# Author: Zachary M. Smith
# Organization: ICPRB
# Email: zsmith@icprb.org
# Date: 2/06/2017
# Purpose: Script written to help E. Stieber calculate benthic
# macroinvertebrate metrics for NYS lentic data.
#==============================================================================
#==============================================================================
# Install the package devtools, which makes it easy to install the Benthos
# package.
install.packages("devtools")
devtools::install_github("zsmith27/Benthos")
# Load the Benthos package.
library(Benthos)
# Load the Master Taxa List contained within the Benthos package.
data(master)
#==============================================================================
# Set your working directory.
# This is where your files are stored on your computer.
setwd("C:/Users/Owner/Desktop/NYSDEC/NYSDEC_LAKE")
# Import csv containing taxonomic counts.
dec <- read.csv("Taxa Counts_Metric Dev_Lakes.csv", stringsAsFactors = FALSE)
#==============================================================================
# Format data to meet the formating requirments of the Benthos package.
long <- nysdec_data_prep(dec)
# If you look at long you will see NAs in the taxonomic hierarchy.
# These values are missing because the taxon was not identified beyond
# a specific rank. For example, if the taxa was identified to the family-level
# then it is not possible to assign a genus or species name.
# The names might also be missing because they are not specified on ITIS.gov.
# The fill_taxa function fills in the NAs with proceeding taxonomic level.
long.fill <- fill_taxa(long)
#==============================================================================
# Calculate all of the available metrics.
metrics.df <- all_metrics(long.fill, master, "FAMILY")
#==============================================================================
# Change your working directory if you want the calculated metrics to be
# stored in a different file location then where the taxonomic count file
# was stored.
setwd("C:/Users/Owner/Desktop/NYSDEC/NYSDEC_LAKE")
# Export as a csv file.
write.csv(metrics.df, "my_metrics.csv", row.names = FALSE)
#==============================================================================
#==============================================================================

