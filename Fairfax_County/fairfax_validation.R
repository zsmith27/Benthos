#==============================================================================
#==============================================================================
# Title: Fairfax County Metric Validation
# Author: Zachary M. Smith
# Organization: ICPRB
# Email: zsmith@icprb.org
# Created: 5/14/2017
# Update: 5/15/2017
# Purpose: Script written for J. Witt to verify the metric are being calculated
# correctly.
#==============================================================================
#==============================================================================
# Set working directory.
setwd("C:/Users/Owner/Desktop/Benthos/Benthos/Fairfax_County")
# Import the previously calculated metrics (See fairfax_workflow.R).
metrics.df <- read.csv("metrics_05_16_2017.csv", stringsAsFactors = FALSE)
#==============================================================================
# Prepare Metric values
#==============================================================================
# Identify the columns containing sample information.
samp.col <- names(metrics.df[, 1:7])
# Identify the metric columns to keep.
metric.col <- c("RICH", "RICH_EPT", "PCT_EPT", "PCT_EPHEMEROPTERA",
                "PCT_TRICHOP_NO_HYDRO", "PCT_CLING",
                #"PCT_CLING_PLECOPTERA",
                "PCT_COLEOPTERA", "PCT_PREDATOR",
                "PCT_SHRED", "HBI", "FBI", "PCT_DOM1")
# Create a new dataframe containing only the columns of interest.
sub.metrics <- metrics.df[, c(samp.col, metric.col)]
#==============================================================================
# Check Metric Values
#==============================================================================
# Import metric values calculated by Fairfax County.
check.wide <- read.csv("metric_check.csv", stringsAsFactors = FALSE)
# Transform the check.wide dataframe from a wide data format to a long data format.
check.long <- tidyr::gather(check.wide, "METRIC", "CHECK", 2:ncol(check.wide))
#------------------------------------------------------------------------------
# Keep only the UNIQUE_ID and metric columns.
benthos.wide <- sub.metrics[, c(1, 8:ncol(sub.metrics))]
# Transform the check.wide dataframe from a wide data format to a long data format.
benthos.long <- tidyr::gather(benthos.wide, "METRIC", "BENTHOS", 2:ncol(benthos.wide))
#------------------------------------------------------------------------------
# Merge the metric calculated by Fairfax County with the metric calculate with
# the Benthos package.
final.df <- merge(check.long, benthos.long, by = c("UNIQUE_ID", "METRIC"))
#------------------------------------------------------------------------------
# Round the metric values to the nearest hundreth.
final.df[, c("CHECK", "BENTHOS")] <- sapply(final.df[, c("CHECK", "BENTHOS")],
                                          function(x) round(x, 2))
#------------------------------------------------------------------------------
# Find the difference between the two metric values.
# Values should be equal, and thus, DIFF should equal zero.
final.df$DIFF <- final.df$CHECK - final.df$BENTHOS
# Sort the dataframe by UNIQUE_ID and METRIC.
final.df <- final.df[order(final.df$UNIQUE_ID, final.df$METRIC), ]
#------------------------------------------------------------------------------
