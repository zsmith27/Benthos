#==============================================================================
#'Prepare the NYSDEC Master Taxa List
#'
#'@param dec.master A Master Taxa List containing the appropriate and applicable
#'taxonomic hierarchy and taxonomic traits. Specify 'master' to use the Master
#'Taxa List built in to the Benthos package.
#'@return Prepare master taxa list from the New York State Department of 
#'Environmental Conservation (NYSDEC) Stream Biomonitoring Unit (SBU).
#'@export
nysdec_master_prep <- function(dec.master){
  dec.master <- clean_up(dec.master)
  # Identify columns of class character and factor.
  char.cols <- sapply(dec.master, class) %in% c('character', 'factor')
  # Remove leading and trailing white space from characters and factors.
  dec.master[, char.cols] <- apply(dec.master[, char.cols], 2, toupper)
  #============================================================================
  names(dec.master)[names(dec.master) %in% "CLAS"] <- "CLASS"
  names(dec.master)[names(dec.master) %in% "ORDR"] <- "ORDER"
  names(dec.master)[names(dec.master) %in% "GENSPECIES"] <- "FINAL_ID"
  #============================================================================
  # Create a Genus column.
  dec.master$GENUS <- dec.master$FINAL_ID
  # Remove any UNDETERMINED
  rm.1 <- c("UNDETERMINED", "UNDET")
  dec.master$GENUS <- ifelse(grepl(paste(rm.1, collapse = "|"), dec.master$GENUS), "", dec.master$GENUS)
  
  rm.2 <- c("\\([^\\)]+\\)", "\ NR\\.", "\ GR\\.", "\\?")
  # Remove any text contained within parentheses
  dec.master$GENUS <- gsub(paste(rm.2, collapse = "|"),"", dec.master$GENUS)
  dec.master$GENUS <- trimws(dec.master$GENUS)
  dec.master$GENUS <- ifelse(grepl(" ",  dec.master$GENUS),
                            gsub( " .*$", "", dec.master$GENUS),
                            dec.master$GENUS)
  # Complexes are not allowed. Pick the first taxon listed in the complex.
  #dec.master$GENUS <- ifelse(grepl("/",  dec.master$GENUS),
  #                          gsub( "/.*$", "", dec.master$GENUS),
  #                          dec.master$GENUS)
  #============================================================================
  # Create a Species column.
  dec.master$SPECIES <- dec.master$FINAL_ID
  # Remove any UNDETERMINED
  rm.3 <- c("SP\\.", "SPP\\.", "CF\\.", "UNDET\\.", "UNDETERMINED", "/")
  dec.master$SPECIES <- ifelse(grepl(paste(rm.3, collapse = "|"), dec.master$SPECIES), "", dec.master$SPECIES)
  
  rm.4 <- c("\\([^\\)]+\\)", "\ NR\\.", "\ GR\\.", "\\?")
  # Remove any text contained within parentheses
  dec.master$SPECIES <- gsub(paste(rm.4, collapse = "|"),"", dec.master$SPECIES)
  dec.master$SPECIES <- trimws(dec.master$SPECIES)
  # Replace the space between genus and FINAL_ID with "_"
  dec.master$SPECIES <- gsub(" ","_", dec.master$SPECIES)
  #============================================================================
  dec.master[dec.master == ""] <- NA
  #============================================================================
  dec.master$FEEDINGHAB <- ifelse(dec.master$FEEDINGHAB %in% "PRD", "PR",
                                  ifelse(dec.master$FEEDINGHAB %in% "C-F", "CF",
                                         ifelse(dec.master$FEEDINGHAB %in% "C-G", "CG",
                                                ifelse(dec.master$FEEDINGHAB %in% "SCR", "SC",
                                                       ifelse(dec.master$FEEDINGHAB %in% "SHR", "SH",
                                                              ifelse(is.na(dec.master$FEEDINGHAB), NA, "ERROR"))))))
  


  return(dec.master)
}
