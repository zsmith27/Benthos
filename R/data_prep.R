#==============================================================================
# Title: Prepare Data for Analysis
#==============================================================================
clean_up <- function(x) {
  # Change all names to uppercase and remove leading and trailing white space.
  names(x) <- trimws(toupper(names(x)))
  # Identify columns of class character and factor.
  char.cols <- sapply(x, class) %in% c('character', 'factor')
  # Remove leading and trailing white space from characters and factors.
  x[, char.cols] <- apply(x[, char.cols], 2, trimws)
  return(x)
}



#==============================================================================
#'Fill empty rows with lowest level of
#'
#'@param Taxon_List = Master taxon list in wide format
#'@return Fills in blank spaces with the lowest level of taxonomic
#' identification.
#'@export
#requires package zoo
#applies lowest identification to each row
fill_taxa <- function(Taxon_List){
  
  taxa <- c("PHYLUM", "SUBPHYLUM", "CLASS",
            "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE",
            "GENUS", "SPECIES")
  
  Taxon_List[Taxon_List == ""]  <- NA
  final.df <- Taxon_List
  final.df[, taxa] <- data.frame(t(apply(final.df[, taxa], 1, zoo::na.locf)))
  names(final.df) <- toupper(colnames(final.df))
  final.df[, taxa] <- lapply(final.df[, taxa], toupper)
  #final.df[,"PHYLUM"] <- toupper(final.df[,"PHYLUM"])
  #final.df[,"SUBPHYLUM"] <- toupper(final.df[,"SUBPHYLUM"])
  #final.df[,"CLASS"] <- toupper(final.df[,"CLASS"])
  #final.df[,"SUBCLASS"] <- toupper(final.df[,"SUBCLASS"])
  #final.df[,"ORDER"] <- toupper(final.df[,"ORDER"])
  #final.df[,"SUBORDER"] <- toupper(final.df[,"SUBORDER"])
  #final.df[,"FAMILY"] <- toupper(final.df[,"FAMILY"])
  #final.df[,"SUBFAMILY"] <- toupper(final.df[,"SUBFAMILY"])
  #final.df[,"TRIBE"] <- toupper(final.df[,"TRIBE"])
  #final.df[,"GENUS"] <- toupper(final.df[,"GENUS"])
  #final.df[,"SPECIES"] <- toupper(final.df[,"SPECIES"])
  return(final.df)
}


#==============================================================================
#'Long Data Frame
#'Transform taxa count data from wide to a long format.
#'@param x = Taxa count data
#'@param Level = Taxonomic Level (PHYLUM, CLASS, ORDER, FAMILY, GENUS)
#'@return Taxa counts in a long data format.
#'@export
long <- function (wide.df, taxa.rank = "FAMILY") {
  
  wide.df <- clean_up(wide.df)
  
  long <- tidyr::gather_(wide.df, taxa.rank, "REPORTING_VALUE", names(wide.df[, 6:ncol(wide.df)]))
  

  return(wide)
}


#==============================================================================
#'Wide Data Frame
#'Transform taxa count data from long to wide format
#'@param long = Taxa count data
#'@param taxa.rank = Taxonomic Level (PHYLUM, CLASS, ORDER, FAMILY, GENUS)
#'@param pct.unid = A threshold can be established to exclude samples that do
#'have too many unidentified taxa in a sample. Enter the percentage of the
#'sample that you are comfortable with being unidentified.  If not specified
#'the function will not excluded any samples.
#'@return Taxa counts in wide format
#'@export

wide <- function (long.df, taxa.rank, pct.unid = NULL) {
  print("[1/2] Aggregating data for transformation")
  # Use data.table to speed up the aggregation process.
  long.df <- data.table::data.table(long.df)
  # List of columns to aggregate by.
  agg.list <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                "AGENCY_CODE", taxa.rank)
  # Aggregate the taxonomic counts.
  agg <- long.df[, sum(REPORTING_VALUE), by = agg.list]
  # Update column names.
  names(agg) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                     "AGENCY_CODE", taxa.rank, "REPORTING_VALUE")

  #============================================================================
  print("[2/2] Transforming from long.df data format to wide data format.")
  wide.df <- tidyr::spread_(agg, key = taxa.rank, "REPORTING_VALUE" )
  

  # Fill all NA's with zeros.
  wide.df[is.na(wide.df)] <- 0 #NA = zero
  # Sort the dataframe.
  wide.df <- wide.df[order(EVENT_ID, STATION_ID, DATE, SAMPLE_NUMBER, AGENCY_CODE), ]
  # All columns to uppercase. Easier for specification latter.
  names(wide.df) <- toupper(colnames(wide.df))
  #============================================================================
  # If pct.unid is specified the appropriate rows are removed and a message
  # is returned specifying the number of rows removed.
  if(!is.null(pct.unid) & "UNIDENTIFIED" %in% names(wide.df)){
    cat("Samples with >=", pct.unid, "% taxa unidentified at the specified 
        taxonomic level were excluded from the data set (N = ",
        nrow(wide.df) - sum((wide.df$UNIDENTIFIED / 
                               rowSums(wide.df[, 6:ncol(wide.df)]) * 100 >= pct.unid)),
        "). \n The number of samples with >= ", pct.unid, "% unidentified taxa: ",
        sum((wide.df$UNIDENTIFIED / rowSums(wide.df[, 6:ncol(wide.df)]) *
               100 >= pct.unid)), " (N = ", nrow(wide.df), "; ",
        round((sum((wide.df$UNIDENTIFIED / rowSums(wide.df[, 6:ncol(wide.df)]) * 
                      100 >= pct.unid)) / nrow(wide.df)) * 100, 2), "%)", sep ="")
    
    wide.df <- wide.df[!((wide.df$UNIDENTIFIED /
                            rowSums(wide.df[, 6:ncol(wide.df)])) * 100 >= pct.unid), ]
  }
  
  final.df <- data.frame(wide.df)
  return(final.df)
}

#==============================================================================
