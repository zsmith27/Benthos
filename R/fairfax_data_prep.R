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

#==============================================================================
#'Prepare the Fairfax County Taxonomic Attributes Table
#'
#'@param attributes.df = Taxonomic attributes table.
#'@return Prepare Fairfax County's Taxonomic attribute table for use in the
#'Benthos package.
#'@importFrom magrittr "%>%"
#'@export

fairfax_attributes_prep <- function(attributes.df) {
  # All column names to uppercase and remove leading and trailing white space
  # from column names, character fields, and factor fields. 
  attributes.clean <- clean_up(attributes.df)
  #----------------------------------------------------------------------------
  # Replace strings with trailing "F" and "G"  with NA from the FAMILY and 
  # GENUS columns, respectively.
  attributes.clean$FAMILY <- ifelse(endsWith(attributes.clean$FAMILY, "F"), NA,
                                    attributes.clean$FAMILY)
  attributes.clean$GENUS <- ifelse(endsWith(attributes.clean$GENUS, "G"), NA,
                                   attributes.clean$GENUS)
  #----------------------------------------------------------------------------
  # Change clinger abreviation (CL) to the standard Benthos abreviation (CN).
  attributes.clean$HABITAT <- ifelse(attributes.clean$HABITAT %in% "CL", "CN",
                                   attributes.clean$HABITAT)
  #----------------------------------------------------------------------------
  # Identify the taxonomic rank columns present in the attributes table.
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS","SUBCLASS", "INFRACLASS",
                 "SUPERORDER", "ORDER", "SUBORDER", "INFRAORDER",
                 "SUPERFAMILY", "FAMILY", "SUBFAMILY", "TRIBE",
                 "GENUS", "SPECIES")
  taxa.cols <- taxa.cols[taxa.cols %in% names(attributes.clean)]
  # Make sure all taxonomic names are uppercase, have no leading/trailing
  # white space, and are of class character.
  attributes.clean[, taxa.cols] <- sapply(attributes.clean[, taxa.cols], function(x){
    toupper(x) %>% trimws() %>% as.character()
  })
  # Use the fill_taxa function to fill NAs in the taxonomic rank columns
  # with the previously identified taxonomic rank.
  attributes.clean$FINAL_ID <- fill_taxa(attributes.clean)$GENUS
  #----------------------------------------------------------------------------
  # Sort the columns so that the FINAL_ID column is first.
  attributes.final <- attributes.clean[, c(ncol(attributes.clean),
                                           1:(ncol(attributes.clean) - 1))]
  #----------------------------------------------------------------------------
  # End fairfax_attributes_prep function.
  return(attributes.final)
}
