#==============================================================================
#'Prepare the Fairfax County Data
#'
#'@param taxa.df = Taxonomic counts in a long data format.
#'@return Prepare taxonomic data from the Fairfax County, VA.
#'@importFrom magrittr "%>%"
#'@export

fairfax_data_prep <- function(taxa.df){
  # Column names to uppercase and remove any leading/trailing white space.
  names(taxa.df) <- toupper(names(taxa.df)) %>% trimws()
  #----------------------------------------------------------------------------
  # Check for an AGENCY_CODE column and add one if missing.
  if(!"AGENCY_CODE" %in% names(taxa.df)) {
    taxa.df$AGENCY_CODE <- "FAIRFAX_COUNTY"
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
  # Change column name "MONITORING.LOCATION.ID" to STATION_ID.
  names(taxa.df)[names(taxa.df) %in% "SITEID"] <- "STATION_ID"
  # Make sure the STATION_ID column is class character.
  taxa.df$STATION_ID <- as.character(taxa.df$STATION_ID)
  #----------------------------------------------------------------------------
  # Change column name "COUNT" to REPORTING_VALUE.
  names(taxa.df)[names(taxa.df) %in% "COUNT"] <- "REPORTING_VALUE"
  # Make sure the REPORTING_VALUE is class numeric.
  taxa.df$REPORTING_VALUE <- as.character(taxa.df$REPORTING_VALUE) %>% as.numeric()
  #----------------------------------------------------------------------------
  # Create a UNIQUE_ID column by concatenating STATION_ID and DATE.
  taxa.df$UNIQUE_ID <- paste(taxa.df$STATION_ID, taxa.df$DATE, sep = "_")
  #----------------------------------------------------------------------------
  # If the METHOD column is absent, create it and fill it with NAs.
  if(!"METHOD" %in% names(taxa.df)){
    taxa.df$METHOD <- "NONE"
  }
  # Make sure the METHOD column is class character.
  taxa.df$METHOD <- as.character(taxa.df$METHOD)
  #----------------------------------------------------------------------------
  # Change column name "REACH" to "CONDITION".
  names(taxa.df)[names(taxa.df) %in% "REACH"] <- "CONDITION"
  #----------------------------------------------------------------------------
  # Specify taxonomic columns of interest.
  taxa.cols <- c("PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS")
  # For each taxonomic column convert the values to class character, make all 
  # of the characters uppercase, and remove any leading/trailing white space.
  taxa.df[, taxa.cols] <- lapply(taxa.df[, taxa.cols], function(x){
    as.character(x) %>% toupper() %>% trimws()
  })
  
  # Replace blanks ("") with NA.
  taxa.df[, taxa.cols] <- lapply(taxa.df[, taxa.cols], function(x){
    ifelse(x %in% "", NA, x)
  })
  #----------------------------------------------------------------------------
  # Create a FINAL_ID column by using the fill_taxa function and extracting
  # the lowest taxonomic resolution (GENUS).
  taxa.df$FINAL_ID <- fill_taxa(taxa.df)$GENUS
  # Make sure the final IDs are class character, all uppercase, and
  # leading/trailing white space has been removed.
  taxa.df$FINAL_ID <- as.character(taxa.df$FINAL_ID) %>% toupper() %>% trimws()
  # Replace any remaining spaces int he FINAL_ID column with "_". This should
  # only be applicable to species.
  taxa.df$FINAL_ID <- gsub(" ", "_", taxa.df$FINAL_ID)
  #----------------------------------------------------------------------------
  # Keep only the necessary columns.
  final.df <- taxa.df[, c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                          "METHOD","SAMPLE_NUMBER", "CONDITION", "FINAL_ID",
                          "REPORTING_VALUE", taxa.cols)]
  #============================================================================
  # End fairfax_data_prep function.
  return(final.df)
}
