#==============================================================================
#GENERA METRICS
#==============================================================================
#'Acid Tolerance Index (ATI)
#'
#'@param genus.wide = Taxonomic counts aggregated at the genus level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The number of Acid Tolerant Individuals (ATI) as described
#' by the NYSDEC (2014). The following genera are considered ATI taxa:
#' Epeorus (Ephemeroptera; Heptageniidae), Amphinemura (Plecoptera; Nemouridae),
#'  Leuctra (Plecoptera; Leuctridae), Isoperla (Plecoptera; Perlodidae),
#'   Rhyacophila (Trichoptera; Rhyacophilidae), Simulium (Diptera; Simuliidae),
#' Conchapelopia (Diptera; Chironomidae), Cricotopus (Diptera; Chironomidae),
#' Eukiefferiella (Diptera; Chironomidae), and Heterotrissocladius (Diptera; Chironomidae).
#'@export

pct_acid_tol <- function(genus.wide){
  acid_taxa <- c("EPEORUS", "AMPHINEMURA", "LEUCTRA", "ISOPERLA", "RHYACOPHILA",
                 "SIMULIUM", "CONCHAPELOPIA", "CRICOTOPUS", "EUKIEFFERIELLA",
                 "HETEROTRISSOCLADIUS")
  final.vec <- rowSums(genus.wide[, names(genus.wide) %in% acid_taxa])
  return(final.vec)
}

#==============================================================================
#'Percentage of Chironomids Identified as Chironomus and Cricoptopus
#'
#'@param genus.wide = Taxonomic counts aggregated at the genus level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param family.wide = Taxonomic counts aggregated at the family level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of chironomids (Family: Chironomidae) individuals
#'that were identified as Chironomus (Diptera; Chironomidae) and
#'Cricotopus (Diptera; Chironomidae). This metric was found in the
#'Macroinvertebrate Metric Calculation Guidelines provided by GADNR/EPD
#'Watershed Protection Branch.
#'@export

pct_cc_chironomidae <- function(family.wide, genus.wide){
  sum_cc <- blank_col("CHIRONOMUS", genus.wide) + blank_col("CRICOTOPUS", genus.wide)
  final.vec <- ifelse(family.wide$CHIRONOMIDAE > 0, 
                      (sum_cc / family.wide$CHIRONOMIDAE) * 100, 0)
  return(final.vec)
}

#==============================================================================
#'Ephemeropteran Richness minus Epeorus
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param genus.wide = Taxonomic counts aggregated at the genus level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The number of genera identified as ephemeropteran taxa
#' minus the genus Epeorus (Ephemeroptera; Heptageniidae),
#' as described by the NYSDEC.
#'@export

rich_ephem_epeorus <- function(long, genus.wide){
  if(!("EPEORUS" %in% names(genus.wide))){
    final.vec <- 0
  }
  
  if("EPEORUS" %in% names(genus.wide)){
    final.vec <- ifelse(genus.wide$EPEORUS > 0,
                        rich_ephemeroptera(long, "GENUS") - 1,
                        rich_ephemeroptera(long, "GENUS"))
  }
  
  return(final.vec)
}

#==============================================================================
#'Percent EPT Minus Cheumatopsyche
#'
#'@param order.wide = Taxonomic counts aggregated at the genus level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param genus.wide = Taxonomic counts aggregated at the genus level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as EPT (Orders: Ephemeroptera,
#'Plecoptera, and Trichoptera) minus Cheumatopsyche (Trichoptera; Hydropsychidae)
#'individuals.
#'@export

pct_ept_cheumatopsyche <- function(order.wide, genus.wide){
  pct_cheumatopsyche <- (blank_col("CHEUMATOPSYCHE", genus.wide) /
                           rowSums(genus.wide[, 7:ncol(genus.wide)])) * 100
  final.vec <- pct_ept(order.wide) - pct_cheumatopsyche
  return(final.vec)
}

#==============================================================================
#'Percent EPT minus Hydropsyche
#'
#'@param order.wide = Taxonomic counts aggregated at the order level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@param genus.wide = Taxonomic counts aggregated at the genus level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as EPT (Orders: Ephemeroptera,
#'Plecoptera, and Trichoptera) minus Hydropsyche (Trichoptera; Hydropsychidae)
#'individuals.
#'@export

pct_ept_hydropsyche <- function(order.wide, genus.wide){
  pct_hydropsyche <- (blank_col("HYDROPSYCHE", genus.wide) /
                        rowSums(genus.wide[, 7:ncol(genus.wide)])) * 100
  final.df <- pct_ept(order.wide) - pct_hydropsyche
  return(final.df)
}

#==============================================================================
#'Percentage of Tanytarsini
#'
#'@param tribe.wide = Taxonomic counts aggregated at the tribe level
#' classification in a wide data format. Use the wide function to
#' prepare the data.
#'@return The percentage of individuals identified as Tanytarsini
#'(Diptera; Chironomidae) individuals. Tanytarsini is a tribe from the family
#'Chironomidae.
#'@export

pct_tanytarsini <- function(tribe.wide){
  return((blank_col("TANYTARSINI", tribe.wide) /
            rowSums(tribe.wide[, 7:ncol(tribe.wide)])) * 100)
  
}

#==============================================================================
#'Percentage of Orthocladiinae
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@return The percentage of Orthocladiinae (Diptera; Chironomidae) individuals.
#'Orthocladiinae is a subfamily of the family Chironomidae.
#'@export

pct_orthocladiinae <- function(long){
  subfam.wide <- wide(long, "SUBFAMILY")
  final.vec <- (blank_col("ORTHOCLADIINAE", subfam.wide) /
                  rowSums(subfam.wide[, 7:ncol(subfam.wide)])) * 100
  return(final.vec)
}
