#==============================================================================
#==============================================================================
# Title: DRBC Metrics
# Author: Zachary M. Smith
# Organization: ICPRB
# Email: zsmith@icprb.org
# Date: 4/23/2017
# Purpose: Script written to help R. Limbeck calculate metrics using the 
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
setwd("C:/Users/Owner/Desktop/Benthos/DRBC")
# Import csv containing taxonomic counts.
drbc.df <- read.csv("Results.csv", stringsAsFactors = FALSE)
#==============================================================================
# Format data to meet the formating requirments of the Benthos package.
long.df <- drbc_data_prep(drbc.df)

# If you look at long.df you will see NAs in the taxonomic hierarchy.
# These values are missing because the taxon does not have a valid taxonomic 
# name at a specified level (e.g., Infraorder) or the taxon was not identified
# beyond a specific rank. For example, if the taxa was identified to the
# family-level then it is not possible to assign a genus or species name.
# The fill_taxa function fills in the NAs with proceeding taxonomic level.
long.fill <- fill_taxa(long.df)
#==============================================================================
# Benthos has a built in table (i.e., master) that has a pretty extensive list
# of taxonomic attributes. If you want to test out the Benthos "master" table
# let me know. I will need to add some additional script.

# For now I am assuming you want to use the taxa attributes provided in your
# data table. The function below will extract the applicable taxa attributes 
# columns from the supplied table.

drbc.master <- drbc_master_prep(drbc.df)
#==============================================================================
# Calculate all of the available metrics.
metrics.df <- all_metrics(long.fill, drbc.master, "GENUS", 
                          tv.col = "TV",
                          ffg.col = "FFG",
                          hab.col = "HABIT")
system.time(
metrics.df <- all_metrics(long.fill, master, "GENUS", 
                          tv.col = "BIBI_TV",
                          ffg.col = "BIBI_FFG",
                          hab.col = "BIBI_HABIT"))

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


