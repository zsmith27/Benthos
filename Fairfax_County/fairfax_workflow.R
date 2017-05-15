#==============================================================================
#==============================================================================
# Title: Fairfax County Metrics
# Author: Zachary M. Smith
# Organization: ICPRB
# Email: zsmith@icprb.org
# Created: 5/14/2017
# Update: 5/15/2017
# Purpose: Script written for J. Witt to calculate metrics using the 
# Benthos package.
#==============================================================================
#==============================================================================
# Install the packages devtools, curl, and httr which enable you to install the
# Benthos package from GitHub.
# It is not necesarry to install these packages with each use.
install.packages(c("devtools","curl", "httr"))
# Install the Benthos package from GitHub.
# Best practice would to be install the latest version of the Benthos package
# with each use.
devtools::install_github("zsmith27/Benthos",  force = TRUE)
# Load the Benthos package.
library(Benthos)
# Load the Master Taxa List contained within the Benthos package.
data(master)
#==============================================================================
# Set your working directory.
# This is where your files are stored on your computer.
# You will have to update the file pathway.
# Notice that the folders are seperate by a forward slash (/) and not a 
# back slash (\). Back slashes will cause issues.
setwd("C:/Users/Owner/Desktop/Benthos/Benthos/Fairfax_County")
# Import csv containing taxonomic counts.
counts.df <- read.csv("fairfax_counts.csv", stringsAsFactors = FALSE)
#==============================================================================
# Format data to meet the formating requirments of the Benthos package.
long.df <- fairfax_data_prep(counts.df)

# If you look at long.df you will see NAs in the taxonomic hierarchy.
# These values are missing because the taxon does not have a valid taxonomic 
# name at a specified level (e.g., Genus) or the taxon was not identified
# beyond a specific rank. For example, if the taxa was identified to the
# family-level then it is not possible to assign a genus name.
# The fill_taxa function fills in the NAs with proceeding taxonomic level.
long.fill <- fill_taxa(long.df)
#==============================================================================
# Calculate all of the available metrics.
metrics.df <- all_metrics(long.fill, master, "FAMILY", 
                          tv.col = "BIBI_TV",
                          ffg.col = "BIBI_FFG",
                          hab.col = "BIBI_HABIT",
                          beck.col = "BECK_CLASS")
#==============================================================================
# Change your working directory if you want the calculated metrics to be
# stored in a different file location then where the taxonomic count file
# was stored.
setwd("C:/Users/Owner/Desktop/Benthos/Benthos/Fairfax_County")
# Create a file name including the current date.
file.name <- paste0(paste("metrics", format(Sys.Date(), "%m_%d_%Y"),
                          sep = "_"), ".csv")
# Export as a csv file.
write.csv(metrics.df, file.name, row.names = FALSE)
#==============================================================================
#==============================================================================
env.col <- names(metrics.df[, 1:7])
metric.col <- c("RICH", "RICH_EPT", "PCT_EPT", "PCT_EPHEMEROPTERA",
                "PCT_CLING", "PCT_COLEOPTERA", "PCT_DOM1")
sub.metrics <- metrics.df[, c(env.col, metric.col)]


