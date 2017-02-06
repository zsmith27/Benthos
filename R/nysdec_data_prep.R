#==============================================================================
#'Prepare the NYSDEC SBU Data
#'
#'@param Taxa.df = A data frame containing raw taxonomic counts.
#'@return Prepare taxonomic data from the New York State Department of 
#'Environmental Conservation (NYSDEC) Stream Biomonitoring Unit (SBU).
#'@export

nysdec_data_prep <- function(Taxa.df){
  # Column names to uppercase.
  names(Taxa.df) <- toupper(names(Taxa.df))
  # Check for an AGENCY_CODE column and add one if missing.
  if(!"AGENCY_CODE" %in% names(Taxa.df)) {
    Taxa.df$AGENCY_CODE <- "NYSDEC_SBU"
  } 
  # Check for a REPLICATE column. If absent, create a SAMPLE_NUMBER column
  # and fill with ones. If present, change the column name to SAMPLE_NUMBER.
  if(!"REPLICATE" %in% names(Taxa.df)){
    Taxa.df$SAMPLE_NUMBER <- 1
  } else {
    names(Taxa.df)[names(Taxa.df) %in% "REPLICATE"] <- "SAMPLE_NUMBER"
    Taxa.df$SAMPLE_NUMBER[is.na(Taxa.df$SAMPLE_NUMBER)] <- 1
  }
  # Change column name COLL_DATE to DATE.
  names(Taxa.df)[names(Taxa.df) %in% "COLL_DATE"] <- "DATE"
  # Change column name SITE_ID to STATION_ID.
  names(Taxa.df)[names(Taxa.df) %in% "SITE_ID"] <- "STATION_ID"
  # Make sure the STATION_ID column is numeric.
  Taxa.df$STATION_ID <- as.character(Taxa.df$STATION_ID)
  # Column GENSPECIES to uppercase.
  Taxa.df$GENSPECIES <- toupper(as.character(Taxa.df$GENSPECIES))
  # Change column name INDIV to REPORTING_VALUE.
  names(Taxa.df)[names(Taxa.df) %in% "INDIV"] <- "REPORTING_VALUE"
  Taxa.df$REPORTING_VALUE <- as.numeric(as.character(Taxa.df$REPORTING_VALUE))
  # Check for EVENT_ID column. If absent, create one.
  if(!"EVENT_ID" %in% names(Taxa.df)) {
    Taxa.df$EVENT_ID <- with(Taxa.df, paste(STATION_ID, DATE, SAMPLE_NUMBER, sep = "_"))
  }
  
  
  
  final.df <- Taxa.df[, c("EVENT_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                          "SAMPLE_NUMBER", "GENSPECIES", "REPORTING_VALUE")]
  final.df$GENUS <- final.df$GENSPECIES
  # Remove any UNDETERMINED
  final.df$GENUS <- gsub("UNDETERMINED ","", final.df$GENUS)
  # Remove any text contained within parentheses
  final.df$GENUS <- gsub("\\([^\\)]+\\)","", final.df$GENUS)
  # Remove NR.
  final.df$GENUS <- gsub("\ NR\\.","", final.df$GENUS)
  # Remove GR.
  final.df$GENUS <- gsub("\ GR\\.","", final.df$GENUS)
  # Remove ?
  final.df$GENUS <- gsub("\\?","", final.df$GENUS)
  final.df$GENUS <- trimws(final.df$GENUS)
  final.df$GENUS <- ifelse(grepl(" ",  final.df$GENUS),
                           gsub( " .*$", "", final.df$GENUS),
                           final.df$GENUS)
  # Complexes are not allowed. Pick the first taxon listed in the complex.
  final.df$GENUS <- ifelse(grepl("/",  final.df$GENUS),
                           gsub( "/.*$", "", final.df$GENUS),
                           final.df$GENUS)
  
  # Remove all genera, groups, undetermineds, complexes, and uncertainties
  final.df$FINAL_ID <- sapply(final.df$GENSPECIES, function(x){
    remove <- c("SP\\.", "SPP\\.", "CF\\.", "UNDET\\.", "UNDETERMINED", "/")
    ifelse(grepl(paste(remove, collapse="|"), x), "", paste(x))
  })
  # Remove any text contained within parentheses
  final.df$FINAL_ID <- gsub("\\([^\\)]+\\)","", final.df$FINAL_ID)
  # Remove NR.
  final.df$FINAL_ID <- gsub("\ NR\\.","", final.df$FINAL_ID)
  # Remove GR.
  final.df$FINAL_ID <- gsub("\ GR\\.","", final.df$FINAL_ID)
  # Remove ?
  final.df$FINAL_ID <- gsub("\\?","", final.df$FINAL_ID)
  # Replace the space between genus and FINAL_ID with "_"
  final.df$FINAL_ID <- gsub(" ","_", final.df$FINAL_ID)
  # Replace the blanks with "UNDETERMINED"
  #final.df$FINAL_ID <- sub("^$", "UNDETERMINED", final.df$FINAL_ID)
  final.df$FINAL_ID <- sub("^$", NA, final.df$FINAL_ID)
  #============================================================================
  
  final.df[, c("GENUS", "FINAL_ID")] <- t(apply(final.df[, c("GENUS", "FINAL_ID")], 1, zoo::na.locf))
  final.df <- final.df[c("EVENT_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                         "SAMPLE_NUMBER", "FINAL_ID", "REPORTING_VALUE")]
  
  return(final.df)
}

#==============================================================================

