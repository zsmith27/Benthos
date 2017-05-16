#==============================================================================
#'Prepare the NYSDEC SBU Data
#'
#'@param taxa.df Taxonomic counts in a long data format.
#'@param master.df A Master Taxa List containing the appropriate and applicable
#'taxonomic hierarchy and taxonomic traits. Specify 'master' to use the Master
#'Taxa List built in to the Benthos package.
#'@return Prepare taxonomic data from the New York State Department of 
#'Environmental Conservation (NYSDEC) Stream Biomonitoring Unit (SBU).
#'@export

nysdec_data_prep_old <- function(taxa.df, master.df){
  # Column names to uppercase.
  names(taxa.df) <- toupper(names(taxa.df))
  # Check for an AGENCY_CODE column and add one if missing.
  if(!"AGENCY_CODE" %in% names(taxa.df)) {
    taxa.df$AGENCY_CODE <- "NYSDEC_SBU"
  } 
  # Check for a REPLICATE column. If absent, create a SAMPLE_NUMBER column
  # and fill with ones. If present, change the column name to SAMPLE_NUMBER.
  if(!"REPLICATE" %in% names(taxa.df)){
    taxa.df$SAMPLE_NUMBER <- 1
  } else {
    names(taxa.df)[names(taxa.df) %in% "REPLICATE"] <- "SAMPLE_NUMBER"
    taxa.df$SAMPLE_NUMBER[is.na(taxa.df$SAMPLE_NUMBER)] <- 1
  }
  # Change column name COLL_DATE to DATE.
  names(taxa.df)[names(taxa.df) %in% "COLL_DATE"] <- "DATE"
  # Change column name SITE_ID to STATION_ID.
  names(taxa.df)[names(taxa.df) %in% "SITE_ID"] <- "STATION_ID"
  # Make sure the STATION_ID column is numeric.
  taxa.df$STATION_ID <- as.character(taxa.df$STATION_ID)
  # Column GENSPECIES to uppercase.
  taxa.df$GENSPECIES <- trimws(toupper(as.character(taxa.df$GENSPECIES)))
  # Change column name INDIV to REPORTING_VALUE.
  names(taxa.df)[names(taxa.df) %in% "INDIV"] <- "REPORTING_VALUE"
  taxa.df$REPORTING_VALUE <- as.numeric(as.character(taxa.df$REPORTING_VALUE))
  # Check for UNIQUE_ID column. If absent, create one.
  if(!"UNIQUE_ID" %in% names(taxa.df)) {
    taxa.df$UNIQUE_ID <- with(taxa.df, paste(STATION_ID, DATE, SAMPLE_NUMBER, METHOD, sep = "_"))
  }
  
  if(!"CONDITION" %in% names(taxa.df)) taxa.df$CONDITION <- "NONE"
  
  
  agency.df <- taxa.df[, c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                          "METHOD", "SAMPLE_NUMBER", "CONDITION", "GENSPECIES",
                          "REPORTING_VALUE")]
  agency.df$GENUS <- agency.df$GENSPECIES
  # Replace "Undet. Tubificidae w/ cap. setae" with Tubificidae.
  tubificidae.cols <- c("UNDET. TUBIFICIDAE W/O CAP. SETAE",
                        "UNDET. TUBIFICIDAE W/ CAP. SETAE")
  agency.df$GENUS <- ifelse(agency.df$GENUS %in% tubificidae.cols,
                            "TUBIFICIDAE", agency.df$GENUS)
  # Remove any UNDETERMINED
  agency.df$GENUS <- gsub("UNDETERMINED ","", agency.df$GENUS)
  # Remove any text contained within parentheses
  agency.df$GENUS <- gsub("\\([^\\)]+\\)","", agency.df$GENUS)
  # Remove NR.
  agency.df$GENUS <- gsub("\ NR\\.","", agency.df$GENUS)
  # Remove GR.
  agency.df$GENUS <- gsub("\ GR\\.","", agency.df$GENUS)
  # Remove ?
  agency.df$GENUS <- gsub("\\?","", agency.df$GENUS)
  agency.df$GENUS <- trimws(agency.df$GENUS)
  agency.df$GENUS <- ifelse(grepl(" ",  agency.df$GENUS),
                           gsub( " .*$", "", agency.df$GENUS),
                           agency.df$GENUS)
  # Complexes are not allowed. Pick the first taxon listed in the complex.
  agency.df$GENUS <- ifelse(grepl("/",  agency.df$GENUS),
                           gsub( "/.*$", "", agency.df$GENUS),
                           agency.df$GENUS)
  
  # Remove all genera, groups, undetermineds, complexes, and uncertainties
  agency.df$FINAL_ID <- sapply(agency.df$GENSPECIES, function(x){
    remove <- c("SP\\.", "SPP\\.", "CF\\.", "UNDET\\.", "UNDETERMINED", "/")
    ifelse(grepl(paste(remove, collapse="|"), x), "", paste(x))
  })
  # Remove any text contained within parentheses
  agency.df$FINAL_ID <- gsub("\\([^\\)]+\\)","", agency.df$FINAL_ID)
  # Remove NR.
  agency.df$FINAL_ID <- gsub("\ NR\\.","", agency.df$FINAL_ID)
  # Remove GR.
  agency.df$FINAL_ID <- gsub("\ GR\\.","", agency.df$FINAL_ID)
  # Remove ?
  agency.df$FINAL_ID <- gsub("\\?","", agency.df$FINAL_ID)
  # Replace the space between genus and FINAL_ID with "_"
  agency.df$FINAL_ID <- gsub(" ","_", agency.df$FINAL_ID)
  # Replace the blanks with "UNDETERMINED"
  #agency.df$FINAL_ID <- sub("^$", "UNDETERMINED", agency.df$FINAL_ID)
  agency.df$FINAL_ID <- sub("^$", NA, agency.df$FINAL_ID)
  #============================================================================
  agency.df[, c("GENUS", "FINAL_ID")] <- t(apply(agency.df[, c("GENUS", "FINAL_ID")], 1, zoo::na.locf))
  agency.df <- agency.df[c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                           "METHOD","SAMPLE_NUMBER", "CONDITION", "FINAL_ID",
                           "REPORTING_VALUE")]
  #============================================================================
  # Use the generic data_prep function to merge the taxonomic counts with the
  # Master Taxa List (master.df).
  final.df <- data_prep(agency.df, master.df)
  
  return(final.df)
}

#==============================================================================
#==============================================================================
#'Prepare the NYSDEC SBU Data
#'
#'@param taxa.df Taxonomic counts in a long data format.
#'@param master.df A Master Taxa List containing the appropriate and applicable
#'taxonomic hierarchy and taxonomic traits. Specify 'master' to use the Master
#'Taxa List built in to the Benthos package.
#'@return Prepare taxonomic data from the New York State Department of 
#'Environmental Conservation (NYSDEC) Stream Biomonitoring Unit (SBU).
#'@export

nysdec_data_prep <- function(taxa.df, master.df){
  # Column names to uppercase.
  names(taxa.df) <- toupper(names(taxa.df))
  # If GENSPECIES is a column name, change the name to FINAL_ID.
  names(taxa.df)[names(taxa.df) %in% "GENSPECIES"] <- "FINAL_ID"
  # Check for an AGENCY_CODE column and add one if missing.
  if(!"AGENCY_CODE" %in% names(taxa.df)) {
    taxa.df$AGENCY_CODE <- "NYSDEC_SBU"
  } 
  # Check for a REPLICATE column. If absent, create a SAMPLE_NUMBER column
  # and fill with ones. If present, change the column name to SAMPLE_NUMBER.
  if(!"REPLICATE" %in% names(taxa.df)){
    taxa.df$SAMPLE_NUMBER <- 1
  } else {
    names(taxa.df)[names(taxa.df) %in% "REPLICATE"] <- "SAMPLE_NUMBER"
    taxa.df$SAMPLE_NUMBER[is.na(taxa.df$SAMPLE_NUMBER)] <- 1
  }
  # Change column name COLL_DATE to DATE.
  names(taxa.df)[names(taxa.df) %in% "COLL_DATE"] <- "DATE"
  # Change column name SITE_ID to STATION_ID.
  names(taxa.df)[names(taxa.df) %in% "SITE_ID"] <- "STATION_ID"
  # Make sure the STATION_ID column is class character.
  taxa.df$STATION_ID <- as.character(taxa.df$STATION_ID)
  # Column FINAL_ID to uppercase.
  taxa.df$FINAL_ID <- trimws(toupper(as.character(taxa.df$FINAL_ID)))
  # Change column name INDIV to REPORTING_VALUE.
  names(taxa.df)[names(taxa.df) %in% "INDIV"] <- "REPORTING_VALUE"
  taxa.df$REPORTING_VALUE <- as.numeric(as.character(taxa.df$REPORTING_VALUE))
  # Check for UNIQUE_ID column. If absent, create one.
  if(!"UNIQUE_ID" %in% names(taxa.df)) {
    taxa.df$UNIQUE_ID <- with(taxa.df, paste(STATION_ID, DATE, SAMPLE_NUMBER, METHOD, sep = "_"))
  }
  
  if(!"CONDITION" %in% names(taxa.df)) taxa.df$CONDITION <- "NONE"

  #============================================================================
  taxa.df <- taxa.df[, c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                           "METHOD","SAMPLE_NUMBER", "CONDITION", "FINAL_ID",
                           "REPORTING_VALUE")]
  #============================================================================
  # Use the generic data_prep function to merge the taxonomic counts with the
  # Master Taxa List (master.df).
  #final.df <- data_prep(taxa.df, master.df)
  #============================================================================
  # Subset the master taxa list to only include taxonomic hierarchy.
  sub.master <- unique(master.df[, c("FINAL_ID", "PHYLUM", "CLASS",
                                     "ORDER", "FAMILY", "SUBFAMILY",
                                     "GENUS", "SPECIES")])
  
  #============================================================================
  # Merge the taxonomic counts with the master taxa list.
  final.df <- merge(taxa.df, sub.master, by = "FINAL_ID", all.x = TRUE)
  #============================================================================
  # If there are any NAs in the FINAL_ID column, a taxon or taxa are missing 
  # from the master taxa list. This issue must be fixed.
  final.df$ISSUE <- ifelse(is.na(final.df$PHYLUM), "YES", "NO")
  if (nrow(final.df[final.df$ISSUE %in% "YES", ]) > 0) {
    
    issue.rows <- paste("\n", "There was an issue matching ", 
                        nrow(final.df[final.df$ISSUE %in% "YES", ]),
                        " rows with the Master Taxa List (master.df).",  "\n",
                        "Please review the rows in the 'ISSUE' column equal to 'YES'.",  "\n",
                        "The missing taxon or taxa must be added to master.df before proceeding.",
                        "\n",
                        sep = "")
    warning(issue.rows)
  } else {
    print("No issues merging the Master Taxa List (master.df) with your taxonomic counts (long.df).")
  }
  
  return(final.df)
}

#==============================================================================

