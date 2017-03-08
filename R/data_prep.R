#==============================================================================
#==============================================================================
# Title: Prepare Data for Analysis
#==============================================================================
#==============================================================================
#'Data Prep
#'
#'@param long.df = Taxonomic counts in a long data format.
#'@param master.df = A Master Taxa List containing the appropriate and applicable
#'taxonomic hierarchy and taxonomic traits. Specify 'master' to use the Master
#'Taxa List built in to the Benthos package.
#'@return The provided taxonomic counts are merged with the provided Master 
#'Taxa List taxonomic hierarchy.
#'@export
#'


data_prep <- function(long.df, master.df){
  # Make sure all column names are uppercase.
  #names(long.df) <- toupper(names(long.df))
  # All column names to uppercase and remove leading and trailing white space
  # from column names, character fields, and factor fields.
  long.df <- clean_up(long.df)
  #============================================================================
  # Necessary column names.
  benthos.cols <- c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE", "METHOD",
                    "SAMPLE_NUMBER", "FINAL_ID", "REPORTING_VALUE")
  # Check if any of the necessary columns are missing.
  if(any((benthos.cols %in% names(long.df)) == FALSE)){
    missing.cols <- benthos.cols[!benthos.cols %in% names(long.df)]
    error.1 <- paste("You are missing the following column(s):",
                     paste(missing.cols, collapse = ", "))
    stop(error.1)
  }
  #============================================================================
  # If the condition column does not exist, create it and fill it with none.
  if(!"CONDITION" %in% names(long.df)) long.df$CONDITION <- "NONE"
  #============================================================================
  # Subset the master taxa list to only include taxonomic hierarchy.
  sub.master <- unique(master[, c("FINAL_ID", "AGENCY_ID",
                                  "PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                                  "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                                  "TRIBE", "GENUS", "SPECIES")])
  # Change the name of the FINAL_ID column to AGENCY_ID to make it easier to
  # merge with the master taxa list.
  names(long.df)[names(long.df) %in% "FINAL_ID"] <- "AGENCY_ID"
  long.df$AGENCY_ID <- toupper(long.df$AGENCY_ID)
  #============================================================================
  # Merge the taxonomic counts with the master taxa list.
  merged <- merge(long.df, sub.master, by = "AGENCY_ID", all.x = TRUE)
  #============================================================================
  # If there are any NAs in the FINAL_ID column, a taxon or taxa are missing 
  # from the master taxa list. This issue must be fixed.
  merged$ISSUE <- ifelse(is.na(merged$FINAL_ID), "YES", "NO")
  if (nrow(merged[merged$ISSUE %in% "YES", ]) > 0) {
    
    issue.rows <- paste("\n", "There was an issue matching ", 
                        nrow(merged[merged$ISSUE %in% "YES", ]),
                        " rows with the Master Taxa List (master.df).",  "\n",
                        "Please review the rows in the 'ISSUE' column equal to 'YES'.",  "\n",
                        "The missing taxon or taxa must be added to master.df before proceeding.",
                        "\n",
                        sep = "")
    warning(issue.rows)
  } else {
    print("No issues merging the Master Taxa List (master.df) with your taxonomic counts (long.df).")
  }
  #============================================================================
  # If the FINAL_ID does not match the AGENCY_ID, a warning is provided.
  # The taxon has been changed to the most up-to-date taxonomic name in ITIS
  # but the user may disagree. No action is required.
  merged$WARNING <- ifelse(merged$FINAL_ID %in% merged$AGENCY_ID, "NO", "YES")
  if (nrow(merged[merged$WARNING %in% "YES", ]) > 0) {
    
    warning.rows <- paste("\n", nrow(merged[merged$WARNING %in% "YES", ]),
                          " rows contain taxonomic final IDs that were changed to the current valid taxonomic name.",
                          "\n", 
                          "Please review the rows in the 'WARNING' column equal to 'YES'.",  "\n",
                          "If you disagree with the change(s) please email me at zsmith@icprb.org and I will address your concerns.",
                          "\n", "No action is necessary.", 
                          "\n",
                          sep = "")
    warning(warning.rows)
  } else {
    cat(paste("There were no differences between your final ID and the final ID found", 
              "\n", "in the Master Taxa List (master.df)."))
  }
  #============================================================================
  keep.cols <- c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE", "METHOD",
                 "SAMPLE_NUMBER", "CONDITION", "ISSUE", "WARNING",
                 "FINAL_ID", "AGENCY_ID",
                 "REPORTING_VALUE",
                 "PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS", 
                 "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES") 
  keep.cols[!keep.cols %in% names(merged)]
  final.df <- merged[, keep.cols]
  final.df[final.df == ""] <- NA
  # Sort columns by "UNIQUE_ID".
  final.df <- final.df[order(final.df$UNIQUE_ID), ]
  return(final.df)
}



#if((!class(long.df[, c("UNIQUE_ID", "STATION_ID", "FINAL_ID")]) %in% c("character", "factor")) == FALSE){}

#If(length(long.df[!long.df$FINAL_ID %in% master.df$FINAL_ID, ]) > 0) {
#  missing.taxa <- sort(unique(long.df[!long.df$FINAL_ID %in% master.df$FINAL_ID, "FINAL_ID"]))
#  error.2 <- paste0("The following taxon or taxa were not found in the Master Taxa List: ",
#                   paste(missing.taxa, collapse = ", "), ". Check the spelling of each name,
#                   an exact match is required.")
#  print(error.2)
#}
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
  
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS",
            "SUBCLASS", "ORDER", "SUBORDER",
            "FAMILY", "SUBFAMILY", "TRIBE",
            "GENUS", "SPECIES")
  taxa.cols <- taxa.cols[taxa.cols %in% names(Taxon_List)]
  
  
  Taxon_List[Taxon_List == ""]  <- NA
  final.df <- Taxon_List
  final.df[, taxa.cols] <- data.frame(t(apply(final.df[, taxa.cols], 1, zoo::na.locf)))
  names(final.df) <- toupper(colnames(final.df))
  final.df[, taxa.cols] <- lapply(final.df[, taxa.cols], toupper)
  
  return(final.df)
}


#==============================================================================
#'Long Data Frame
#'Transform taxa count data from wide to a long format.
#'@param x = Taxa count data
#'@param Level = Taxonomic Level (PHYLUM, CLASS, ORDER, FAMILY, GENUS)
#'@return Taxa counts in a long data format.
#'@export
#'
long <- function (wide.df, taxa.rank = "FAMILY") {
  
  wide.df <- clean_up(wide.df)
  
  long <- tidyr::gather_(wide.df, taxa.rank, "REPORTING_VALUE",
                         names(wide.df[, 8:ncol(wide.df)]))
  

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
  #print("[1/2] Aggregating data for transformation")
  # Use data.table to speed up the aggregation process.
  #long.dt <- data.table::data.table(long.df)
  # List of columns to aggregate by.
  agg.list <- c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE", "DATE",
                "METHOD", "SAMPLE_NUMBER", "CONDITION",
                 taxa.rank, "REPORTING_VALUE")
  # Aggregate the taxonomic counts.
  #agg <- long.dt[, sum(REPORTING_VALUE), by = agg.list]
  sub.df <- long.df[, agg.list]
  sub.df[, taxa.rank] <- ifelse(is.na(sub.df[, taxa.rank]) | sub.df[, taxa.rank] %in% "",
                             "UNIDENTIFIED", sub.df[, taxa.rank])
  agg <- aggregate(REPORTING_VALUE ~ .,
                   data = sub.df, FUN = sum)
  
  
  # Update column names.
  #names(agg) <- c("UNIQUE_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
  #                   "AGENCY_CODE", "CONDITION", taxa.rank, "REPORTING_VALUE")

  #============================================================================
  #print("[2/2] Transforming from long.df data format to wide data format.")
  wide.df <- tidyr::spread_(agg, key = taxa.rank, "REPORTING_VALUE" )
  if(nrow(long.df[!long.df$UNIQUE_ID %in% wide.df$UNIQUE_ID, ]) > 0){
    missing.long <- unique(long.df[!long.df$UNIQUE_ID %in% wide.df$UNIQUE_ID, 1:7])
    missing.wide <- cbind(missing.long, wide.df[1, 8:ncol(wide.df)])
    missing.wide[, 8:ncol(missing.wide)] <- NA
    wide.df <- rbind(wide.df, missing.wide)
  }

  # Fill all NA's with zeros.
  wide.df[is.na(wide.df)] <- 0 #NA = zero
  # Sort the dataframe.
  wide.df <- wide.df[order(wide.df$UNIQUE_ID, wide.df$STATION_ID,
                           wide.df$AGENCY_CODE,
                           wide.df$DATE, wide.df$METHOD,
                           wide.df$SAMPLE_NUMBER), ]
  # All columns to uppercase. Easier for specification latter.
  names(wide.df) <- toupper(colnames(wide.df))
  #============================================================================
  # If pct.unid is specified the appropriate rows are removed and a message
  # is returned specifying the number of rows removed.
  if(!is.null(pct.unid) & "UNIDENTIFIED" %in% names(wide.df)){
    cat("Samples with >=", pct.unid, "% taxa unidentified at the specified 
        taxonomic level were excluded from the data set (N = ",
        nrow(wide.df) - sum((wide.df$UNIDENTIFIED / 
                               rowSums(wide.df[, 8:ncol(wide.df)]) * 100 >= pct.unid)),
        "). \n The number of samples with >= ", pct.unid, "% unidentified taxa: ",
        sum((wide.df$UNIDENTIFIED / rowSums(wide.df[, 8:ncol(wide.df)]) *
               100 >= pct.unid)), " (N = ", nrow(wide.df), "; ",
        round((sum((wide.df$UNIDENTIFIED / rowSums(wide.df[, 8:ncol(wide.df)]) * 
                      100 >= pct.unid)) / nrow(wide.df)) * 100, 2), "%)", sep ="")
    
    wide.df <- wide.df[!((wide.df$UNIDENTIFIED /
                            rowSums(wide.df[, 8:ncol(wide.df)])) * 100 >= pct.unid), ]
  }
  
  if("UNIDENTIFIED" %in% names(wide.df)) {
    names(wide.df)[names(wide.df) %in% "UNIDENTIFIED"] <- paste("UNIDENTIFIED", taxa.rank, sep = "_")
  }
  
  final.df <- data.frame(wide.df)
  return(final.df)
}

#==============================================================================
#'Benthos Cheat
#'
#'@param long.df = Taxonomic counts in a long data format.
#'@return The Benthos package requires eight columns to exist before the metrics
#'can be calculated.  Although it is generally not recommended, five of these
#'columns (i.e., "STATION_ID", "AGENCY_CODE", "DATE", "METHOD", and
#' "SAMPLE_NUMBER") can be automatically generated and filled with a character
#'  sting (i.e., "blank"). Some of the required columns are not applicable to 
#'  all users and would have to be filled with a constant before proceeding;
#'  this function does the job for you.
#'@export
#'

benthos_cheat <- function(long.df){
  long.df <- clean_up(long.df)
  need.cols <- c("UNIQUE_ID", "FINAL_ID", "REPORTING_VALUE")
  if (!need.cols %in% names(long.df)){
    stop(paste("UNIQUE_ID, FINAL_ID, and REPORTING_VALUE must exist as column names."))
  } 
  benthos.cols <- c("STATION_ID", "AGENCY_CODE", "DATE", "METHOD", "SAMPLE_NUMBER")
  sub.benthos.cols <- benthos.cols[!benthos.cols %in% names(long.df)]
  if (length(sub.benthos.cols) > 0) {
    sub.benthos.cols <- benthos.cols[!benthos.cols %in% names(long.df)]
    fill.cols <- lapply(sub.benthos.cols, function(x){
      long.df[, x] <- rep("blank", nrow(long.df))
    })
    
    new.cols <- data.frame(do.call(cbind, fill.cols))
    names(new.cols) <- sub.benthos.cols
    final.df <- cbind(long.df, new.cols)
    
  } else {
    final.df <- long.df
  }
  
  return(final.df)
}