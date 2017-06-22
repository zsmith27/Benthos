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
#------------------------------------------------------------------------------
# Set your working directory.
# Example: setwd("C:/Users/Owner/Desktop/Benthos/DRBC")
# This is where your files are stored on your computer.
# You will have to update the file pathway.
# Notice that the folders are seperate by a forward slash (/) and not a 
# back slash (\). Back slashes will cause issues.
# Do not use the follow line of code, it has been streamlined me.
setwd("DRBC")
# Import csv containing taxonomic counts.
drbc.df <- read.csv("Results.csv", stringsAsFactors = FALSE)
#------------------------------------------------------------------------------
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
# Calculating All Metrics
#==============================================================================
# Benthos has a built in table (i.e., master) that has a pretty extensive list
# of taxonomic attributes. If you want to test out the Benthos "master" table
# let me know. You can test out the Chessie BIBI attributes like this:
all.metrics.bibi <- all_metrics(long.df = long.fill,
                          master.df = master,
                          taxa.rank = "GENUS", 
                          tv.col = "BIBI_TV",
                          ffg.col = "BIBI_FFG",
                          hab.col = "BIBI_HABIT",
                          beck.col = "BECK_CLASS")
#------------------------------------------------------------------------------
# However, I am assuming you want to use the taxa attributes provided in your
# data table. The function below will extract the applicable taxa attributes 
# columns from the supplied table.
drbc.master <- drbc_master_prep(drbc.df)
# Currently, there are duplicated final IDs when I extract the taxa attributes
# from the orginal table you provide. This causes issues with the Benthos
# functions becuase they need a unique set of attributes for each taxon. 
# Differing attributes for the same taxon will throw an error. 
# Run the script below and you can see duplicates.
library(dplyr)
duplicated.final.id <- drbc.master %>% 
  select(FINAL_ID, TV, FFG, HABIT) %>% 
  distinct() %>% 
  group_by(FINAL_ID) %>% 
  filter(n() > 1) %>% 
  arrange(FINAL_ID)
#------------------------------------------------------------------------------
# Many of the duplicates were the result of a taxon having assigned taxonomic
# attributes in one row and then the same taxon having no assigned taxonomic
# attributes in the another row (ie., TV, FFG, and HABIT are NA). 
# The following script removes the row that is filled with NAs.
# Additionally, Polypedilum has two TV values (5 and 6) and
# two FFG variables (SH and SC). There was only one instance where Polypedilum
# was assigned a TV value of 5 and a FFG variable of SC. The script removes,
# what I am assumming, is a typo.
# I recommend that you thuroughly review the output of the script below to
# verify I have not made any mistakes and that I have not made any changes
# that you disaggree with.
working.drbc.master <- drbc.master %>% 
  filter(case_when(
    FINAL_ID %in% duplicated.final.id$FINAL_ID & 
      is.na(TV) & is.na(FFG) & is.na(HABIT)  ~ FALSE,
    FINAL_ID %in% duplicated.final.id$FINAL_ID & 
      any(!is.na(TV), !is.na(FFG), !is.na(HABIT)) &
      !FINAL_ID %in% "POLYPEDILUM" ~ TRUE,
    !FINAL_ID %in% duplicated.final.id$FINAL_ID  ~ TRUE,
    FINAL_ID %in% "POLYPEDILUM" & FFG %in% "SC"  ~ FALSE
  ))
#------------------------------------------------------------------------------
# Calculate all of the available metrics using the drbc taxa attribute table.
all.metrics.drbc <- all_metrics(long.df = long.fill,
                          master.df = working.drbc.master,
                          taxa.rank = "GENUS", 
                          tv.col = "TV",
                          ffg.col = "FFG",
                          hab.col = "HABIT", 
                          beck.col = "BECK")
#------------------------------------------------------------------------------
# Change your working directory if you want the calculated metrics to be
# stored in a different file location then where the taxonomic count file
# was stored.
setwd("C:/Users/zsmith/Desktop/Benthos_R/Benthos")
# Create the name of the output file.
file.name <- paste0("metrics_", Sys.Date(), ".csv")
# Export as a csv file.
write.csv(all.metrics.drbc, file.name, row.names = FALSE)
#==============================================================================
# Calculating Specific Metrics
#==============================================================================
# If you are not interested in exploring new metrics or you want to calculate
# metrics individually, then you can use the script below. 
# I tried to calculate all of the metrics listed in your "Interim Methodology
# for Bioassessment of the Delaware River for the DRBC 2010 Integrated 
# Assessment" report.
#------------------------------------------------------------------------------
# Change the long data format to a wide data format necessary for metric 
# calculations below.
ord.df <- wide(long.fill, "ORDER")
fam.df <- wide(long.fill, "FAMILY")
gen.df <- wide(long.fill, "GENUS")
#------------------------------------------------------------------------------
# Create a data frame to store the metrics.
metrics.drbc <- gen.df[, 1:7]
#------------------------------------------------------------------------------
# Richness
metrics.drbc$RICH <- vegan::specnumber(gen.df[, 8:ncol(gen.df)])
#------------------------------------------------------------------------------
# EPT Richness
metrics.drbc$RICH_EPT <- rich_ept(long.fill, "FINAL_ID")
#------------------------------------------------------------------------------
# Ephemeroptera Richness
metrics.drbc$RICH_EPHEMEROPTERA <- taxon_richness(long = long.fill,
                                                taxon = "EPHEMEROPTERA",
                                                low.res.rank = "ORDER",
                                                high.res.rank = "FINAL_ID")
#------------------------------------------------------------------------------
# Plecoptera Richness
metrics.drbc$RICH_PLECOPTERA <- taxon_richness(long = long.fill,
                                                taxon = "PLECOPTERA",
                                                low.res.rank = "ORDER",
                                                high.res.rank = "FINAL_ID")
#------------------------------------------------------------------------------
# Trichoptera Richness
metrics.drbc$RICH_TRICHOPTERA <- taxon_richness(long = long.fill,
                                                taxon = "TRICHOPTERA",
                                                low.res.rank = "ORDER",
                                                high.res.rank = "FINAL_ID")
#------------------------------------------------------------------------------
# Invertebrate Richness
metrics.drbc$RICH_INVERT <- vegan::specnumber(gen.df[, 8:ncol(gen.df)]) - 
  taxon_richness(long = long.fill,
                 taxon = "INSECTA",
                 low.res.rank = "CLASS",
                 high.res.rank = "FINAL_ID")
  
#------------------------------------------------------------------------------
# Percent EPT
metrics.drbc$PCT_EPT <- pct_ept(ord.df)
#------------------------------------------------------------------------------
# Shannon-Wiener Diversity
metrics.drbc$SHANNON <- shannon(gen.df)
#------------------------------------------------------------------------------
# Percentage of 3 Most Dominant Taxon
metrics.drbc$PCT_DOM1 <- pct_dom(gen.df, 3)
#------------------------------------------------------------------------------
# Genus-level Biotic Index (HBI).
metrics.drbc$HBI <- tol_index(long.df = long.fill,
                            master.df = working.drbc.master,
                            tolerance.value = "TV",
                            taxa.rank = "FINAL_ID",
                            remove_na = TRUE)
#------------------------------------------------------------------------------
# Beck's Index
metrics.drbc$BECK <- becks(taxa.wide = gen.df,
                           taxa.rank = "GENUS",
                           master.df= working.drbc.master,
                           beck.column = "BECK",
                           beck.version = 1)
#------------------------------------------------------------------------------
# Intolerant Richness
metrics.drbc$RICH_INTOL <- rich_tolerance(taxa.wide = gen.df,
                                        master.df = working.drbc.master, 
                                        tolerance.value = "TV",
                                        lower.value = 0,
                                        upper.value = 2)
#------------------------------------------------------------------------------
# Intolerant Percent Richness
# I will develop a generic function for performing this task.
metrics.drbc$PCT_RICH_INTOL <- metrics.drbc$RICH_INTOL / metrics.drbc$RICH * 100
#------------------------------------------------------------------------------
# Intolerant Percent Abundance
metrics.drbc$PCT_INTOL <- pct_tol_val(taxa.wide = gen.df,
                                    master.df = working.drbc.master,
                                    tolerance.value = "TV",
                                    lower.value = 0,
                                    upper.value = 2)
#------------------------------------------------------------------------------
# Tolerant Richness
metrics.drbc$RICH_TOL <- rich_tolerance(taxa.wide = gen.df,
                                        master.df = working.drbc.master, 
                                        tolerance.value = "TV",
                                        lower.value = 8,
                                        upper.value = 10)
#------------------------------------------------------------------------------
# Tolerant Percent Richness
# I will develop a generic function for performing this task.
metrics.drbc$PCT_RICH_TOL <- metrics.drbc$RICH_TOL / metrics.drbc$RICH * 100
#------------------------------------------------------------------------------
# Tolerant Percent Abundance
metrics.drbc$PCT_TOL <- pct_tol_val(taxa.wide = gen.df,
                                    master.df = working.drbc.master,
                                    tolerance.value = "TV",
                                    lower.value = 8,
                                    upper.value = 10)
#------------------------------------------------------------------------------
# Scraper Richness
metrics.drbc$RICH_SHRED <- rich_attribute(taxa.wide = gen.df,
                                       master.df = working.drbc.master,
                                       attribute.column = "FFG",
                                       attribute.interest = "SC",
                                       rank = "FINAL_ID")
#------------------------------------------------------------------------------
# Scraper Percent Richness
# I will develop a generic function for performing this task.
metrics.drbc$PCT_RICH_SHRED <- metrics.drbc$RICH_SHRED / metrics.drbc$RICH * 100
#------------------------------------------------------------------------------


