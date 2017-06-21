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
# Uncomment to view the Benthos attributes table.
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
# Import csv of taxonomic attributes.
rbp.df <- read.csv("taxa_values_rbp.csv", stringsAsFactors = FALSE)
#==============================================================================
# Prepare Attributes Table
#==============================================================================
attributes.df <- fairfax_attributes_prep(rbp.df)
#==============================================================================
# Prepare Taxonomic Counts
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
#------------------------------------------------------------------------------
# Change the long data format to a wide data format necessary for metric 
# calculations below.
ord.df <- wide(long.fill, "ORDER")
fam.df <- wide(long.fill, "FAMILY")
gen.df <- wide(long.fill, "GENUS")
#==============================================================================
# Calculate Metrics
#==============================================================================
# Create a data frame to store the metrics.
metrics.df <- gen.df[, 1:7]
#------------------------------------------------------------------------------
# Richness
metrics.df$RICH <- vegan::specnumber(gen.df[, 8:ncol(gen.df)])
#------------------------------------------------------------------------------
# EPT Richness
metrics.df$RICH_EPT <- rich_ept(long.fill, "GENUS")
#------------------------------------------------------------------------------
# Percent EPT minus the percentage of Hydrosychidae
metrics.df$PCT_EPT <- pct_ept(ord.df) - pct_taxon(fam.df, "HYDROPSYCHIDAE")
#------------------------------------------------------------------------------
# Percent Ephemeroptera
metrics.df$PCT_EPHEMEROPTERA <- pct_taxon(ord.df, "EPHEMEROPTERA")
#------------------------------------------------------------------------------
# %Trichoptera minus %Hydropsychidae.
metrics.df$PCT_TRICHOP_NO_HYDRO <- pct_trichoptera_no_hydro(ord.df, fam.df)
#------------------------------------------------------------------------------
# Percent Clinger
metrics.df$PCT_CLING <- pct_attribute(gen.df, attributes.df,
                                   "HABITAT", "CG", "GENUS")
#------------------------------------------------------------------------------
# Percent Coleoptera
metrics.df$PCT_COLEOPTERA <- pct_taxon(ord.df, "COLEOPTERA")
#------------------------------------------------------------------------------
# Percent Clinger
metrics.df$PCT_PREDATOR <- pct_attribute(gen.df, attributes.df,
                                      "FFG", "P", "GENUS")
#------------------------------------------------------------------------------
# Percent Clinger
metrics.df$PCT_SHRED <- pct_attribute(gen.df, attributes.df,
                                      "FFG", "SC", "GENUS")
#------------------------------------------------------------------------------
# Percent Clinger
metrics.df$PCT_SHRED <- pct_attribute(gen.df, attributes.df,
                                      "FFG", "SC", "GENUS")
#------------------------------------------------------------------------------
# Percent Dominant Taxon
metrics.df$PCT_DOM1 <- pct_dom(gen.df, 1)
#------------------------------------------------------------------------------
# Genus-level Biotic Index (HBI).
metrics.df$HBI <- tol_index(long.df = long.fill,
                            master.df = attributes.df,
                            tolerance.value = "MBSS.UTV",
                            taxa.rank = "GENUS",
                            remove_na = TRUE)
#------------------------------------------------------------------------------
# Family-level Biotic Index (FBI).
metrics.df$FBI <- tol_index(long.df = long.fill,
                             master.df = attributes.df,
                             tolerance.value = "MBSS.UTV",
                             taxa.rank = "FAMILY",
                             remove_na = TRUE)
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
# Calculate all of the available metrics.
#==============================================================================
all.df <- all_metrics(long.fill, attributes.df, "GENUS", 
                          tv.col = "MBSS.UTV",
                          ffg.col = "FFG",
                          hab.col = "HABITAT")
pct_taxon(ord.df, c("EPHEMEROPTERA", "TRICHOPTERA", "PLECOPTERA"))
blank_col(c("EPHEMEROPTERA", "TRICHOPTERA", "PLECOPTERA"))
wide.df <- ord.df
taxon <- c("EPHEMEROPTERA", "TRICHOPTERA", "PLECOPTERA")

identical(all.df$GOLD, all.df$GOLD2)
identical(all.df$PCT_COTE, all.df$PCT_COTE2)
identical(all.df$PCT_POTEC, all.df$PCT_POTEC2)
test <- all.df[, c("UNIQUE_ID", "GOLD", "GOLD2")]
test$DIFF <- test[, 2] - test[, 3]

tapply(fam.df$CHIRONOMIDAE, fam.df$UNIQUE_ID + fam.df$STATION_ID, sum)
taxon <- "TRICHOPTERA"

sub.df <- long.fill[long.fill$FAMILY %in% taxon, ]
sapply(unique(long.fill$UNIQUE_ID), function(x) {
  if (x %in% sub.df$UNIQUE_ID) {
    sub.df[]
  }
})
tapply(sub.df$REPORTING_VALUE, sub.df$UNIQUE_ID, sum) /
  tapply(long.fill$REPORTING_VALUE, long.fill$UNIQUE_ID,  sum) *100

sum(long.fill[long.fill$FAMILY %in% taxon, "REPORTING_VALUE"])
tapply(long.fill[long.fill$FAMILY %in% taxon, "REPORTING_VALUE"],
       long.fill[long.fill$FAMILY %in% taxon, "UNIQUE_ID"], sum, na.rm = TRUE) 
taxon <- "CHIRONOMIDAE"
taxon.rank <- "FAMILY"

    
    
    
    
  }
})