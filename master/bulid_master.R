#==============================================================================
# Build Master Taxa List
#==============================================================================
# Upload the current master table in the Benthos package.
library(Benthos)
data(master)
#==============================================================================
# Keep only the taxonomic hierarchy columns.
master.0 <- master[1:21, ]
#==============================================================================
setwd("C:/Users/zsmith/Desktop/Benthos_R/master")
#==============================================================================
# Function to remove leading and trailing white space of character columns.
rmws <- function(org.df){
  # All column names to uppercase.
  names(org.df) <- toupper(names(org.df))
  # Identify all factor columns.
  factor.cols <- names(org.df)[sapply(org.df, class) == 'factor']
  # Convert all factor columns to class character.
  org.df[, factor.cols] <- as.character(org.df[, factor.cols])
  # Identify all character columns.
  char.cols <- names(org.df)[sapply(org.df, class) == 'character']
  # All characters to uppercase.
  org.df[, char.cols] <- apply(org.df[char.cols], 2, toupper)
  if(length(char.cols) == 0) final.df <- org.df
  if(length(char.cols) == 1) {
    org.df[, char.cols] <- trimws(org.df[, char.cols])
    final.df <- org.df
  }
  if(length(char.cols) > 1) {
    org.df[, char.cols] <- apply(org.df[, char.cols], 2, trimws)
    final.df <- org.df
  }
  return(final.df)
}
#==============================================================================
# BIBI (Buchanan et al. 2011)
#==============================================================================
# Import 2011 BIBI master taxa list.
bibi.2011 <- read.csv("BIBI_2011_master_2_17_2017.csv", stringsAsFactors = FALSE)
# Remove the "TSN" column.
bibi.2011 <- bibi.2011[, !names(bibi.2011) %in% "TSN"]
# Remove whitespace from character columns.
bibi.2011 <- rmws(bibi.2011)
# Merge with master.
master.1 <- merge(master.0, bibi.2011, by = "FINAL_ID", all.x = TRUE)
test <- master.1[duplicated(master.1$TSN_FINAL), ]
#==============================================================================
# NYSDEC
#==============================================================================
# Import NYSDEC Stream Biomonitoring Unit master
nysdec <- read.csv("NYSDEC_2_17_2017.csv", stringsAsFactors = FALSE)
nysdec <- rmws(nysdec)
# Remove "CF" and "/" from column "GENSPECIES".
nysdec.rm1 <- c( "CF\\.", "UNDET\\.", "UNDETERMINED", "/")
nysdec2 <- nysdec[!grepl(paste(nysdec.rm1, collapse = "|"), nysdec$GENSPECIES), ]
# Split GENSPECIES column to just represent Genus.
nysdec$GENUS <- gsub(" .*$", "", nysdec$GENSPECIES)
nysdec$GENUS <- sapply(nysdec$GENUS, function(x){
  nysdec.rm2 <- c("SP.", "CF\\.", "UNDET.", "UNDETERMINED", "/", "\\?")
  ifelse(grepl(paste(nysdec.rm2, collapse="|"), x), "", paste(x))
})
# Remove all genera, groups, undeterminds, complexes, and uncertainties
nysdec$SPECIES <- sapply(nysdec$GENSPECIES, function(x){
  nysdec.rm3 <- c("SP\\.", "CF\\.")
  ifelse(grepl(paste(nysdec.rm3, collapse="|"), x), "", paste(x))
})
#nysdec.rm4 <- c("UNDET\\.", "UNDETERMINED", "/")
#nysdec <- nysdec[!(grepl(paste(nysdec.rm4, collapse="|"), nysdec$GENSPECIES)), ]
# Remove any text contained within parentheses
nysdec.rm4 <- c("\ \\([^\\)]+\\)", "\ NR\\.", "\ GR\\.", "\\?")
nysdec$SPECIES <- gsub("\ \\([^\\)]+\\)","", nysdec$SPECIES)
# Remove NR.
nysdec$SPECIES <- gsub("\ NR\\.","", nysdec$SPECIES)
# Remove GR.
nysdec$SPECIES <- gsub("\ GR\\.","", nysdec$SPECIES)
# Remove ?
nysdec$SPECIES <- gsub("\\?","", nysdec$SPECIES)
# Replace the space between genus and species with "_"
nysdec$SPECIES <- gsub(" ","_", nysdec$SPECIES)
# Reorder the dataframe
nysdec <- nysdec[, c("PHYLUM", "CLASS", "ORDER", "FAMILY",
                     "SUBFAMILY", "GENUS", "SPECIES", "GENSPECIES",
                     "TOLERANCE", "FEEDINGHAB", "NBI.P_TOLERANCE",
                     "NBI.N_TOLERANCE")]

nysdec.fill <- fill_taxa2(nysdec)


nysdec$FINAL_ID <- nysdec.fill$SPECIES

nysdec.final <- nysdec[, c("FINAL_ID", "TOLERANCE", "FEEDINGHAB",
                           "NBI.P_TOLERANCE","NBI.N_TOLERANCE")]

names(nysdec.final) <- c("FINAL_ID", "NYSDEC_TV", "NYSDEC_FFG",
                         "NYSDEC_NBI.P", "NYSDEC_NBI.N")

nysdec.final2 <- unique(nysdec.final)
ny <- nysdec.final2[, - which(names(nysdec.final2) %in% "NYSDEC_FFG")]
ny$NYSDEC_TV <- as.numeric(as.character(levels(ny$NYSDEC_TV)))[ny$NYSDEC_TV]
ny$NYSDEC_NBI.P <- as.numeric(as.character(levels(ny$NYSDEC_NBI.P)))[ny$NYSDEC_NBI.P]
ny$NYSDEC_NBI.N <- as.numeric(as.character(levels(ny$NYSDEC_NBI.N)))[ny$NYSDEC_NBI.N]
agg_ny <- aggregate(ny[,2:ncol(ny)] , by = list(ny$FINAL_ID), data = ny, FUN = mean, na.rm = T)
names(agg_ny)[1] <- "FINAL_ID"
agg_ny$NYSDEC_TV[is.nan(agg_ny$NYSDEC_TV)] <- ""
agg_ny$NYSDEC_NBI.P[is.nan(agg_ny$NYSDEC_NBI.P)] <- ""
agg_ny$NYSDEC_NBI.N[is.nan(agg_ny$NYSDEC_NBI.N)] <- ""
nysdec.final3 <- merge(agg_ny, unique(nysdec.final2[, c("FINAL_ID", "NYSDEC_FFG")]),
                       by = "FINAL_ID")
#colnames(nysdec.final3) <- paste("NYSDEC", colnames(nysdec.final3), sep = "_")

merged <- merge(master2, nysdec.final3, by = "FINAL_ID", all.x = T)


test <- nysdec.final3[duplicated(nysdec.final3$FINAL_ID),]

