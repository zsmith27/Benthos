#==============================================================================
#'Prepare the DRBC Data
#'
#'@param taxa.df = Taxonomic counts in a long data format.
#'@return Prepare taxonomic data from the Delaware River Basin Commission 
#'(DRBC).
#'@importFrom magrittr "%>%"
#'@export

drbc_data_prep <- function(taxa.df){
  # Column names to uppercase and remove any leading/trailing white space.
  names(taxa.df) <- toupper(names(taxa.df)) %>% trimws()
  #----------------------------------------------------------------------------
  # If "SUBJECT.TAXONOMIC.NAME" is a column name, change the name to FINAL_ID.
  names(taxa.df)[names(taxa.df) %in% "SUBJECT.TAXONOMIC.NAME"] <- "FINAL_ID"
  #----------------------------------------------------------------------------
  # Check for an AGENCY_CODE column and add one if missing.
  if(!"AGENCY_CODE" %in% names(taxa.df)) {
    taxa.df$AGENCY_CODE <- "DRBC"
  } 
  #----------------------------------------------------------------------------
  # Check for a REPLICATE column. If absent, create a SAMPLE_NUMBER column
  # and fill with ones. If present, change the column name to SAMPLE_NUMBER.
  # Any NAs will be replaced with a 1.
  if(!"REPLICATE" %in% names(taxa.df)){
    taxa.df$SAMPLE_NUMBER <- 1
  } else {
    names(taxa.df)[names(taxa.df) %in% "REPLICATE"] <- "SAMPLE_NUMBER"
    taxa.df$SAMPLE_NUMBER[is.na(taxa.df$SAMPLE_NUMBER)] <- 1
  }
  #----------------------------------------------------------------------------
  # Change column name "ACTIVITY.START.DATE" to DATE.
  names(taxa.df)[names(taxa.df) %in% "ACTIVITY.START.DATE"] <- "DATE"
  #----------------------------------------------------------------------------
  # Change column name "MONITORING.LOCATION.ID" to STATION_ID.
  names(taxa.df)[names(taxa.df) %in% "MONITORING.LOCATION.ID"] <- "STATION_ID"
  # Make sure the STATION_ID column is class character.
  taxa.df$STATION_ID <- as.character(taxa.df$STATION_ID)
  #----------------------------------------------------------------------------
  # Change column name "SAMPLE.COLLECTION.METHOD.ID" to METHOD.
  names(taxa.df)[names(taxa.df) %in% "SAMPLE.COLLECTION.METHOD.ID"] <- "METHOD"
  # Make sure the METHOD column is class character.
  taxa.df$METHOD <- as.character(taxa.df$METHOD)
  #----------------------------------------------------------------------------
  # Make sure the final IDs are class character, all uppercase, and
  # leading/trailing white space has been removed.
  taxa.df$FINAL_ID <- as.character(taxa.df$FINAL_ID) %>% toupper() %>% trimws()
  # Replace any remaining spaces int he FINAL_ID column with "_". This should
  # only be applicable to species.
  taxa.df$FINAL_ID <- gsub(" ", "_", taxa.df$FINAL_ID)
  #----------------------------------------------------------------------------
  # Change column name "RESULT.VALUE" to REPORTING_VALUE.
  names(taxa.df)[names(taxa.df) %in% "RESULT.VALUE"] <- "REPORTING_VALUE"
  # Make sure the REPORTING_VALUE is class numeric.
  taxa.df$REPORTING_VALUE <- as.character(taxa.df$REPORTING_VALUE) %>% as.numeric()
  #----------------------------------------------------------------------------
  # Change column name "UNIQUEID" to "UNIQUE_ID".
  names(taxa.df)[names(taxa.df) %in% "UNIQUEID"] <- "UNIQUE_ID"
  #----------------------------------------------------------------------------
  # Refers to assigned conditions, such as Reference or Degraded.
  if(!"CONDITION" %in% names(taxa.df)) taxa.df$CONDITION <- "NONE"
  #----------------------------------------------------------------------------
  # Specify taxonomic columns of interest.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "INFRACLASS",
                 "SUPERORDER", "ORDER", "SUBORDER", "INFRAORDER",
                 "SUPERFAMILY", "FAMILY", "SUBFAMILY", 
                 "TRIBE", "GENUS", "SPECIES")
  # For each taxonomic column convert the values to class character, make all 
  # of the characters uppercase, and remove any leading/trailing white space.
  taxa.df[, taxa.cols] <- lapply(taxa.df[, taxa.cols], function(x){
    as.character(x) %>% toupper() %>% trimws()
  })
  
  # Replace blanks ("") with NA.
  taxa.df[, taxa.cols] <- lapply(taxa.df[, taxa.cols], function(x){
    ifelse(x %in% "", NA, x)
  })

  # For the SPECIES column replace the "." sperating the genus and species 
  # names with "_".
  taxa.df$SPECIES <- gsub(" ", "_", taxa.df$SPECIES)
  #----------------------------------------------------------------------------
  # Keep only the necessary columns.
  final.df <- taxa.df[, c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                         "METHOD","SAMPLE_NUMBER", "CONDITION", "FINAL_ID",
                         "REPORTING_VALUE", taxa.cols)]
  #============================================================================
  # End drbc_data_prep function.
  return(final.df)
}

#==============================================================================
#'Prepare the DRBC Attributes Table
#'
#'@param taxa.df = Taxonomic counts in a long data format.
#'@return Prepare taxa attributes table from the data provided by the Delaware 
#'River Basin Commission (DRBC).
#'@importFrom magrittr "%>%"
#'@export

drbc_master_prep <- function(taxa.df){
  # Column names to uppercase and remove any leading/trailing white space.
  names(taxa.df) <- toupper(names(taxa.df)) %>% trimws()
  #----------------------------------------------------------------------------
  # Specify taxonomic columns of interest.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", "INFRACLASS",
                 "SUPERORDER", "ORDER", "SUBORDER", "INFRAORDER",
                 "SUPERFAMILY", "FAMILY", "SUBFAMILY", 
                 "TRIBE", "GENUS", "SPECIES")
  #----------------------------------------------------------------------------
  # Keep only the Final ID and taxonomic attribute columns.
  final.df <- taxa.df[, c("SUBJECT.TAXONOMIC.NAME",
                          "TAXON.POLLUTION.TOLERANCE",
                          "TAXON.FUNCTIONAL.FEEDING.GROUP", "TAXON.HABIT",
                          taxa.cols)]
  # Rename the columns.
  names(final.df) <- c("FINAL_ID", "TV", "FFG", "HABIT", taxa.cols)
  #----------------------------------------------------------------------------
  # Identify the columns of class character.
  char.cols <- c("FINAL_ID", "FFG", "HABIT", taxa.cols)
  # Make sure all values are class character, are all uppercase, and any
  # leading/trailing white space has been removed.
  final.df[, char.cols] <- lapply(final.df[, char.cols], function(x){
    as.character(x) %>% toupper() %>% trimws()
  })
  #----------------------------------------------------------------------------
  # Replace any remaining spaces int he FINAL_ID column with "_". This should
  # only be applicable to species.
  final.df$FINAL_ID <- gsub(" ", "_", final.df$FINAL_ID)
  #----------------------------------------------------------------------------
  # For FFG, HABIT, and taxa.cols, replace all blanks ("") with NA.
  sub.cols <- c("FFG", "HABIT", taxa.cols)
  final.df[, sub.cols] <- lapply(final.df[, sub.cols], function(x){
    ifelse(x %in% "", NA, x)
  })
  #----------------------------------------------------------------------------
  # For the SPECIES column replace the "." sperating the genus and species 
  # names with "_".
  taxa.df$SPECIES <- gsub(" ", "_", taxa.df$SPECIES)
  #----------------------------------------------------------------------------
  # Convert the FFG values to match the Benthos syntax.
  final.df$FFG <- ifelse(final.df$FFG %in% "C-F", "CF",
                          ifelse(final.df$FFG %in% "C-G", "CG",
                                 ifelse(final.df$FFG %in% "PRD", "PR",
                                        ifelse(final.df$FFG %in% "SCR", "SC",
                                               ifelse(final.df$FFG %in% "SHR", "SH",
                                                      ifelse(is.na(final.df$FFG), NA, "ERROR"))))))
  #----------------------------------------------------------------------------
  # Remove any duplicate rows.
  final.df <- unique(final.df)
  #----------------------------------------------------------------------------
  # End drbc_master_prep function.
  return(final.df)
}

