#==============================================================================
# Diversity Metrics
#============================================================================
#'Taxon Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. This
#'function requires a rank below the Order level taxonomic classification.
#'@return The number of taxa identified as ephemeropterans (Order: Ephemeroptera).
#'@export

taxon_richness <- function(long, taxon, low.res.rank, high.res.rank){
  taxa.split <- split(long[, high.res.rank], long[, low.res.rank])
  taxa.list <- unique(unlist(taxa.split[taxon]))
  taxa.wide <- wide(long, high.res.rank)
  final.vec <- group_rich(taxa.list, taxa.wide)
  return(final.vec)
}
#==============================================================================
#'Margalef's Index
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return Margalef's Index: (S - 1)/ln(N)
#'where:
#'S = taxa richness
#'N = the total number of individuals observed
#'
#'Taxa richness is calculated using the vegan package specnumber function.
#'@export

margalefs <- function(taxa.wide) {
  Rich <- vegan::specnumber(taxa.wide[, 8:ncol(taxa.wide)]) #Requires the Vegan Package
  Marg <- function(marg_rich.df, marg_fam.df){
    # (Richness - 1) / log(Total Count)
    return((marg_rich.df - 1) / log(rowSums(marg_fam.df)))
  }
  final.vec <- ifelse(vegan::specnumber(taxa.wide[, 8:ncol(taxa.wide)]) > 1,
                      Marg(Rich, taxa.wide[, 8:ncol(taxa.wide)]), 0)
  return(final.vec)
}

#==============================================================================
#'Menhinick's Index
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return Menhinick's Index: S/sqrt(N)
#'where:
#'S = taxa richness
#'N = the total number of individuals observed
#'
#'Taxa richness is calculated using the vegan package specnumber function.
#'@export

menhinicks <- function(taxa.wide) {
  #Requires the Vegan Package
  Rich <- vegan::specnumber(taxa.wide[, 8:ncol(taxa.wide)])
  # Richness / square root(Total Count)
  final.vec <- Rich / sqrt(rowSums(taxa.wide[, 8:ncol(taxa.wide)]))
  return(final.vec)
}

#==============================================================================
#'Percentage of the Most Dominant Taxon (Taxa)
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@param dom.level = A numeric value 1-5 indicating the number
#'@return Percent of individuals that represent the most abundant taxon or taxa.
#'dom.level can be used to specify 1st-5th most abundant taxa by specifying
#'the corresponding numeric value (1-5).  Values >1 include all of the previous
#' dominance levels. For example, dom.level = 3 is the percentage of the most
#' dominant taxon, the second most dominant taxon, and the third most
#' dominant taxon. This measure is related to taxa evenness. Typically
#' degradation is associated with elevated levels of the
#' most dominant taxon (taxa); therefore, this metric typically increases
#' with degradation.
#'@export

pct_dom <- function(taxa.wide, dom.level){
  final.vec <- apply(taxa.wide[, 8:ncol(taxa.wide)], 1, function(x){
    sum(order(x, decreasing = TRUE)[1:dom.level])
  })
  return(final.vec)
}

#==============================================================================
#'Simpson's Diversity Index
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return Simpson's Divieristy Index: 1 - sum(n/N)^2
#'where:
#'n = the abundance of a particular taxon
#'N = the total abundance of organisms
#'
#'Simpson's diversity is calculated using the vegan package function
#'diversity with the index set to "simpson".
#'@export

simpsons <- function(taxa.wide) {
  final.vec <- vegan::diversity(taxa.wide[, 8:ncol(taxa.wide)], "simpson")
  return(final.vec)
}


#==============================================================================
#'Shannon Wiener Diversity Index
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return Shannon Wiener Diversity Index: -sum(p * ln(p))
#'where:
#'p = the propotion a specific taxon composes of a sample
#'
#'Shannon Wiener diversity is calculated using the vegan package function
#'diversity with the index set to "shannon".
#'@export
shannon <- function(taxa.wide) {
  final.vec <- vegan::diversity(taxa.wide[, 8:ncol(taxa.wide)], "shannon")
  return(final.vec)
}

#==============================================================================
#'Probability of Interspecific Encounter (PIE)
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return Hurlbert's (1971) Probability of Interspecific Encounter (PIE):
#'PIE = (N/(N-1)) * (1 - sum(p^2))
#'where:
#'p = the propotion a specific taxon composes of a sample
#'N = the total abundance of organisms
#'
#'This measurement is equivalent to Simpson's Index but includes a correction
#'factor base on the total abundance of organisms.
#'@export

hurlberts_pie <- function(taxa.wide){
  pie_formula <- function(taxon.df){
    pie_part1 <- rowSums(taxon.df[, 8:ncol(taxa.wide)]) /
      (rowSums(taxon.df[, 8:ncol(taxa.wide)]) - 1)
    pie_part2 <- (taxon.df[, 8:ncol(taxa.wide)] /
                    rowSums(taxon.df[, 8:ncol(taxa.wide)])) ^ 2
    pie_part3 <- 1 - rowSums(pie_part2)
    return(pie_part1 * pie_part3)
  }
  final.vec <- ifelse(rowSums(taxa.wide[, 8:ncol(taxa.wide)]) > 1,
                      pie_formula(taxa.wide), 0)
  return(final.vec)
}

#==============================================================================
#'Abundance
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return A count of the total number of organisms observed.
#'@export

abundance <- function(taxa.wide){
  final.vec <- rowSums(taxa.wide[, 8:ncol(taxa.wide)])
  return(final.vec)
}

#==============================================================================
#'Pielou's Evenness
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@return Pielou's Evenness Index: H'/ ln(S)
#'where:
#'H' = Shannon Weiner Diversity
#'S = taxa richness
#'
#'If species richness is <= 1, this function will return 0.
#'Shannon Wiener diversity is calculated using the vegan package function
#'diversity with the index set to "shannon". Taxa richness is calculated
#'using the vegan package function specnumber.
#'@export

pielou <- function(taxa.wide){
  #Requires vegan package
  richness <- vegan::specnumber(taxa.wide[, 8:ncol(taxa.wide)])
  final.vec <- ifelse(richness > 1,
                      vegan::diversity(taxa.wide[, 8:ncol(taxa.wide)]) /
                        log(richness), 0)
  return(final.vec)
}

#============================================================================
#'EPT Richness Excluding Tolerant Taxa
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@param master = A master taxa list including taxonomic ranks Phylum through
#'the specified taxonomic rank (Family or Genus) and the an
#'associated list of tolerance values. The default is set to the master taxa
#'list included in the BIBI package.  The master taxa list can be viewed with
#'the following script: master.df <- data(master)
#'@param tolerance_value = The name of the column in the master taxon list
#'(specified using the master variable) that contains tolerance values on
#'a scale of 0-10.  Tolerant organisms are classified as organisms with a
#'tolerance value >= 7. The defualt is set to the the BIBI tolerance values,
#'which are tolerance values summarized from multiple sources.
#'@return The number of taxa identified as EPT (Orders: Ephemeroptera,
#' Plecoptera, and Trichoptera), excluding taxa with a tolerance value >= 7.
#'@export

ept_rich_no_tol <- function(long, rank = "FAMILY", master, tolerance_value = "BIBI_TV"){
  wide.df <- wide(long, rank)
  Order <- split(long[, rank], long$ORDER)
  ephem <- unique(Order$EPHEMEROPTERA)
  plecop <- unique(Order$PLECOPTERA)
  trichop <- unique(Order$TRICHOPTERA)
  taxa.list <- c(ephem, plecop, trichop)
  new.df <- wide.df[, names(wide.df) %in% taxa.list]
  new.df <- new.df[, !(names(new.df) %in% "UNIDENTIFIED")]
  
  # Find all of the tolerant taxa
  master$TOLERANCE <- ifelse(master[, tolerance_value] >= 7, "TOLERANT", NA)
  tol <- split(master[, rank], master$TOLERANCE)
  name.list <- as.character(unlist(unique(tol$TOLERANT)))
  # Remove any tolerant taxa from the ept data frame
  no.tol.ept <- data.frame(new.df[, !(names(new.df) %in% name.list)])
  
  # Caculate richness values using vegan
  final.vec <- vegan::specnumber(no.tol.ept[, 8:ncol(no.tol.ept)])
  return(final.vec)
}

#============================================================================
#'EPT Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@return The number of taxa identified as EPT (Orders: Ephemeroptera,
#' Plecoptera, and Trichoptera).
#'@export

rich_ept <- function(long, rank = "FAMILY"){
  ephem <- taxon_richness(long,"EPHEMEROPTERA", "ORDER", rank)
  plecop <- taxon_richness(long,"PLECOPTERA", "ORDER", rank)
  trichop <- taxon_richness(long,"TRICHOPTERA", "ORDER", rank)
  final.vec <- ephem + plecop + trichop
  return(final.vec)
}

#==============================================================================
#'EPT Percent Taxa Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@return The percentage of taxa identified as EPT (Orders: Ephemeroptera,
#'Plecoptera, and Trichoptera).
#'@export

pct_ept_rich <- function(long, rank){
  wide.df <- wide(long, rank)
  Order <- split(long[, rank], long$ORDER)
  ephem <- unique(Order$EPHEMEROPTERA)
  plecop <- unique(Order$PLECOPTERA)
  trichop <- unique(Order$TRICHOPTERA)
  sample.info <- names(wide.df[, 1:7])
  taxa.list <- c(sample.info, ephem, plecop, trichop)
  new.df <- wide.df[, names(wide.df) %in% taxa.list]
  new.df <- new.df[, !(names(new.df) %in% "UNIDENTIFIED")]
  ept.rich <- vegan::specnumber(new.df[, 8:ncol(new.df)])
  total.rich <- vegan::specnumber(wide.df[, 8:ncol(wide.df)])
  final.vec <- (ept.rich / total.rich) * 100
  return(final.vec)
}
#==============================================================================
#'Percentage of EPT Taxa Excluding Tolerant Taxa
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@param master = A master taxa list including taxonomic ranks Phylum through
#'the specified taxonomic rank (Family or Genus) and the an
#'associated list of tolerance values. The default is set to the master taxa
#'list included in the BIBI package.  The master taxa list can be viewed with
#'the following script: master.df <- data(master)
#'@param tolerance_value = The name of the column in the master taxon list
#'(specified using the master variable) that contains tolerance values on
#'a scale of 0-10.  Tolerant organisms are classified as organisms with a
#'tolerance value >= 7.  The defualt is set to the the BIBI tolerance values,
#'which are tolerance values summarized from multiple sources.
#'@return Percent of the assembalge represented by Ephemeroptera, Plecoptera,
#' and Trichoptera (EPT) Familial Richness excluding tolerant EPT taxa
#'@export

pct_ept_rich_no_tol <- function (long, rank, master, tolerance_value = "BIBI_TV") {
  
  wide.df <- wide(long, rank)
  Order <- split(long[, rank], long$ORDER)
  ephem <- unique(Order$EPHEMEROPTERA)
  plecop <- unique(Order$PLECOPTERA)
  trichop <- unique(Order$TRICHOPTERA)
  taxa.list <- c(ephem, plecop, trichop)
  new.df <- wide.df[, names(wide.df) %in% taxa.list]
  new.df <- new.df[, !(names(new.df) %in% "UNIDENTIFIED")]
  
  # Find all of the tolerant taxa
  master$TOLERANCE <- ifelse(master[, tolerance_value] >= 7, "TOLERANT", NA)
  tol <- split(master[, rank], master$TOLERANCE)
  name.list <- as.character(unlist(unique(tol$TOLERANT)))
  # Remove any tolerant taxa from the ept data frame
  no.tol.ept <- data.frame(new.df[, !(names(new.df) %in% name.list)])
  
  # Caculate richness values using vegan
  ept.rich <- vegan::specnumber(no.tol.ept[, 8:ncol(no.tol.ept)])
  total.rich <- vegan::specnumber(wide.df[, 8:ncol(wide.df)])
  final.vec <- (ept.rich / total.rich) * 100
  return(final.vec)
}

#============================================================================
#'COTE Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@return The number of taxa identified as COTE (Orders: Coleoptera, Odonata,
#' Trichoptera, and Ephemeroptera).
#'@export

rich_cote <- function(long, rank = "FAMILY"){
  coleop <- taxon_richness(long,"COLEOPTERA", "ORDER", rank)
  odonat <- taxon_richness(long,"ODONATA", "ORDER", rank)
  trichop <- taxon_richness(long,"TRICHOPTERA", "ORDER", rank)
  ephem <- taxon_richness(long, "EPHEMEROPTERA", "ORDER", rank)
  final.vec <-  coleop + odonat + trichop + ephem
  return(final.vec)
}

#============================================================================
#'POTEC Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. You must
#'sepecify either 'FAMILY' or "GENUS.'
#'@return The number of taxa identified as POTEC (Orders: Plecoptera, Odonata,
#' Trichoptera, Ephemeroptera, and Coleoptera).
#'@export

rich_potec <- function(long, rank = "FAMILY"){
  plecop <- taxon_richness(long,"PLECOPTERA", "ORDER", rank)
  odonat <- taxon_richness(long,"ODONATA", "ORDER", rank)
  trichop <- taxon_richness(long,"TRICHOPTERA", "ORDER", rank)
  ephem <- taxon_richness(long, "EPHEMEROPTERA", "ORDER", rank)
  coleop <- taxon_richness(long,"COLEOPTERA", "ORDER", rank)
  final.vec <-  plecop + odonat + trichop + ephem + coleop
  return(final.vec)
}

#==============================================================================
#'Non-Chironomid and Oligochaet Taxa Richness
#'
#'@param long = Taxonomic counts arrange in a long data format (i.e., each
#'row represents a unique sample and taxon).
#'@param rank = The taxonomic rank used to perform the analysis. This
#'function requires a rank below the Order level taxonomic classification.
#'@return The number of taxa that were not identified as oligochaets
#'(Class: Oligochaeta) and chironomids (Family: Chironomidae).
#'@export


rich_nco <- function(long, rank = "GENUS"){
  
  chiro <- taxon_richness(long, "CHIRONOMIDAE", "FAMILY", rank)
  oligo <- taxon_richness(long, "OLIGOCHAETA", "CLASS", rank)
  chiro_oligo.rich <- chiro + oligo
  #============================================================================
  taxa.wide <- wide(long, rank)
  if (ncol(taxa.wide) < 7) {
    total.rich <- rep(0, nrow(taxa.wide))
  } else {
    if (ncol(taxa.wide) == 7) {
      total.rich <- ifelse(taxa.wide[, 8] > 0, 1, 0)
    } else {
      total.rich <- vegan::specnumber(taxa.wide[, 8:ncol(taxa.wide)])
    }
  }
  #============================================================================
  final.vec <- total.rich - chiro_oligo.rich
  
  return(final.vec)
}

#==============================================================================
#'Effective Number of Taxa
#'
#'@param taxa.wide = Taxonomic counts aggregated at the specific taxonomic
#' classification (e.g., Order, Family, or Genus) in a wide data format.
#'  Use the wide function to prepare the data.
#'@param index = Requires "shannon" (Shannon Wiener Diversity Index) or
#'"invsimpson" (Inverse Simpson's Diversity Index).
#'@return This functions converts the common diversity indices values to
#'effective richness, as outlined by Jost (2006).  Effective richness is a
#'richness value weighted by species evenness; the value represents the taxa
#'richness of a perfectly even assemblage reuquired to produce an equivalent
#'diveristy index value to the observed assemblage. Shannon Wiener Diversity
#'is converted using the exponential function: exp(H')
#'where:
#'H'= Shannon Wiener Diversity
#'
#'Inverse Simpsons Diveristy is converted using:
#'1/(1-(1/D))
#'where:
#'D = vegan Inverse Simpson's Diversity Index  (1/(sum(n/N)^2) )
#'n = the abundance of a particular taxon
#'N = the total abundance of organisms
#'@export
#'

effective_richness <- function(taxa.wide, index){
  
  if(index == "shannon"){
    final.vec <- exp(vegan::diversity(taxa.wide[, 8:ncol(taxa.wide)], "shannon"))
  }
  
  if(index == "invsimpson"){
    final.vec <- 1 / (1 + vegan::diversity(taxa.wide[, 8:ncol(taxa.wide)], "invsimpson"))
  }
  
  return(final.vec)
}

