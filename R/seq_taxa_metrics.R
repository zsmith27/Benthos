#==============================================================================
#The percent of the sample represented by each taxa per taxon level
#==============================================================================
#'The percentage each taxon makes up of a sample
#'
#'@param wide.df = Data in a wide data format.
#'@param taxa.rank = Taxonomic rank (all caps).
#'@param master.df = Master taxa list
#'@return The percent of the sample represented by each taxa per taxon level.
#'@export
#==============================================================================
calc_pct_taxa <- function(wide.df, taxa.rank, master.df){
  if(any(rowSums(wide.df[, 6:ncol(wide.df)]) == 0)) {
    with.counts <- wide.df[rowSums(wide.df[, 6:ncol(wide.df)]) != 0, ]
    with.counts[, 6:ncol(with.counts)] <- (with.counts[, 6:ncol(with.counts)] /
                                     rowSums(with.counts[, 6:ncol(with.counts)])) * 100
    without.counts <- wide.df[rowSums(wide.df[, 6:ncol(wide.df)]) == 0, ]
    
    wide.df <- rbind(with.counts, without.counts)
    
  } else {
    wide.df[, 6:ncol(wide.df)] <- (wide.df[, 6:ncol(wide.df)] /
                                     rowSums(wide.df[, 6:ncol(wide.df)])) * 100
  }
  
  t_list <- unique(toupper(master.df[, taxa.rank]))
  shorten <- wide.df[, colnames(wide.df) %in% t_list]
  cn <- colnames(shorten)
  list_taxa.df <- if(length(cn) == 0){
    list_taxa.df <- data.frame(wide.df[, c("EVENT_ID", "STATION_ID",
                                           "DATE", "SAMPLE_NUMBER",
                                           "AGENCY_CODE")])
    list_taxa.df$NO_MATCH <- 0
    colnames(list_taxa.df) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                "AGENCY_CODE", "NO_MATCH")
    list_taxa.df
  }else{
    list_taxa.df <- data.frame(cbind(wide.df[, 1:5], shorten))
    colnames(list_taxa.df) <- c("EVENT_ID", "STATION_ID", "DATE", "SAMPLE_NUMBER",
                                "AGENCY_CODE", paste("PCT", cn, sep = "_"))
    list_taxa.df
  }
  final.df <- data.frame(list_taxa.df)
  return(data.frame(list_taxa.df[order(list_taxa.df$EVENT_ID),]))
}
#==============================================================================
#'The percentage each taxon makes up of a sample
#'
#'@param long.df = Data in a long data format
#'@param master.df = Master taxa list
#'@return The percent of the sample represented by each taxon per taxonomic level.
#'@export
#==============================================================================

seq_pct_taxa <- function(long.df, master.df){
  print("Sequence %Taxa")
  print("...Phylum Wide")
  phylum <- if("PHYLUM" %in% names(long.df)) wide(long.df, "PHYLUM")
  print("...Subphylum Wide")
  subphylum <- if("SUBPHYLUM" %in% names(long.df)) wide(long.df, "SUBPHYLUM")
  print("...Class Wide")
  class <- if("CLASS" %in% names(long.df)) wide(long.df, "CLASS")
  print("Subclass Wide")
  subclass <- if("SUBCLASS" %in% names(long.df)) wide(long.df, "SUBCLASS")
  print("...Order Wide")
  order <- if("ORDER" %in% names(long.df)) wide(long.df, "ORDER")
  print("...Suborder Wide")
  suborder <- if("SUBORDER" %in% names(long.df)) wide(long.df, "SUBORDER")
  print("...Family Wide")
  family <- if("FAMILY" %in% names(long.df)) wide(long.df, "FAMILY")
  print("...Subfamily Wide")
  subfamily <- if("SUBFAMILY" %in% names(long.df)) wide(long.df, "SUBFAMILY")
  print("...Tribe Wide")
  tribe <- if("TRIBE" %in% names(long.df)) wide(long.df, "TRIBE")
  print("...Genus Wide")
  genus <- if("GENUS" %in% names(long.df)) wide(long.df, "GENUS")
  print("...Species Wide")
  species <- if("SPECIES" %in% names(long.df)) wide(long.df, "SPECIES")
  
  pct_phylum <- if(length(phylum) > 0) calc_pct_taxa(phylum, taxa.rank = "PHYLUM")
  pct_subphylum <- if(length(subphylum) > 0) calc_pct_taxa(subphylum, "SUBPHYLUM")
  pct_class <- if(length(class) > 0) calc_pct_taxa(class, "CLASS")
  pct_subclass <- if(length(subclass) > 0) calc_pct_taxa(subclass, "SUBCLASS")
  pct_order <- if(length(order) > 0) calc_pct_taxa(order, "ORDER")
  pct_suborder <- if(length(suborder) > 0) calc_pct_taxa(suborder, "SUBORDER")
  pct_family <- if(length(family) > 0) calc_pct_taxa(family, "FAMILY")
  pct_subfamily <- if(length(subfamily) > 0) calc_pct_taxa(subfamily, "SUBFAMILY")
  pct_tribe <- if(length(tribe) > 0) calc_pct_taxa(tribe, "TRIBE")
  pct_genus <- if(length(genus) > 0) calc_pct_taxa(genus, "GENUS")
  pct_species <- if(length(species) > 0) calc_pct_taxa(species, "SPECIES")
  
  check_exists <- function(pct_taxa){
    pct_taxa <- if(length(pct_taxa) > 0){
      pct_taxa <- pct_taxa
    } else{
      pct_taxa <- pct_taxa[, 1:5]
    }
  }
  
  checked_pct_phylum <- check_exists(pct_phylum)
  checked_pct_subphylum <- check_exists(pct_subphylum)
  checked_pct_class <- check_exists(pct_class)
  checked_pct_subclass <- check_exists(pct_subclass)
  checked_pct_order <- check_exists(pct_order)
  checked_pct_suborder <- check_exists(pct_suborder)
  checked_pct_family <- check_exists(pct_family)
  checked_pct_subfamily <- check_exists(pct_subfamily)
  checked_pct_tribe <- check_exists(pct_tribe)
  checked_pct_genus <- check_exists(pct_genus)
  checked_pct_species <- check_exists(pct_species)
  
  #comb_all <- cbind(checked_pct_phylum, checked_pct_subphylum,
  #                  checked_pct_class, checked_pct_subclass,
  #                  checked_pct_order, checked_pct_suborder,
  #                  checked_pct_family, checked_pct_subfamily,
  #                  checked_pct_tribe, checked_pct_genus,
  #                  checked_pct_species)
  
  comb_all <- plyr::join_all(list(checked_pct_phylum, checked_pct_subphylum,
                                  checked_pct_class, checked_pct_subclass,
                                  checked_pct_order, checked_pct_suborder,
                                  checked_pct_family, checked_pct_subfamily,
                                  checked_pct_tribe, checked_pct_genus,
                                  checked_pct_species), by =  c("EVENT_ID", "STATION_ID",
                                                                "DATE", "SAMPLE_NUMBER",
                                                                "AGENCY_CODE"),
                             type = "full")
  comb_taxa <- comb_all[, 6:ncol(comb_all)]
  rm.cols <- names(comb_taxa[, colSums(comb_taxa) == 0])
  final.df <- comb_all[, !(names(comb_all) %in% rm.cols)]
  
  return(final.df)
}

#==============================================================================
rich_by_rank <- function(long.df, low.res.rank, high.res.rank) {
  rich.list <- lapply(unique(long[, low.res.rank]), function(x){
    print(paste("...RICH", x, sep = "_"))
    final.list <- taxon_richness(long.df, x, low.res.rank, high.res.rank)
  return(final.list)
  } )
  final.df <- data.frame(do.call(cbind, rich.list))
  names(final.df) <- paste("RICH", unique(long.df[, low.res.rank]), sep = "_")
  return(final.df)
}

seq_taxa_rich <- function(long.df, rank = "GENUS"){
  print("Sequence Taxa Richness")
  
  taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
                 "ORDER", "SUBORDER", "FAMILY", "SUBFAMILY",
                 "TRIBE", "GENUS", "SPECIES")
  last.col <- which(taxa.cols %in% rank) - 1
  
  taxa.final <- taxa.cols[1:last.col]
  #============================================================================
  
  
  phylum.df <- if("PHYLUM" %in% names(long.df)){
    print("...Phylum Richness")
    rich_by_rank(long, "PHYLUM", rank)
  } 
  subphylum.df <- if("SUBPHYLUM" %in% names(long.df)){
    print("...Subphylum Richness")
    rich_by_rank(long, "SUBPHYLUM", rank)
  } 
  class.df <- if("CLASS" %in% names(long.df)){
    print("...Class Richness")
    rich_by_rank(long, "CLASS", rank)
  } 
  subclass.df <- if("SUBCLASS" %in% names(long.df)){
    print("...Subclass Richness")
    rich_by_rank(long, "SUBCLASS", rank)
  } 
  order.df <- if("ORDER" %in% names(long.df)){
    print("...Order Richness")
    rich_by_rank(long, "ORDER", rank)
  } 
  suborder.df <- if("SUBORDER" %in% names(long.df)){
    print("...Suborder Richness")
    rich_by_rank(long, "SUBORDER", rank)
  } 
  family.df <- if("FAMILY" %in% names(long.df)){
    print("...Family Richness")
    rich_by_rank(long, "FAMILY", rank)
  } 
  subfamily.df <- if("SUBFAMILY" %in% names(long.df)){
    print("...Subfamily Richness")
    rich_by_rank(long, "SUBFAMILY", rank)
  } 
  tribe.df <- if("TRIBE" %in% names(long.df)){
    print("...Tribe Richness")
    rich_by_rank(long, "TRIBE", rank)
  } 
  genus.df <- if("GENUS" %in% names(long.df)){
    print("...Genus Richness")
    rich_by_rank(long, "GENUS", rank)
  } 

  taxa.dfs <- taxa.cols[taxa.cols %in% names(long.df)]
  paste0(tolower(taxa.dfs), ".df", collapse = ", ")
  #============================================================================
  # Bind columnes together for the data frames that exist. Data frames that
  # do not exist are first filled with zeros and recieve a name "if "...
  # After the columns are appended, columns containing "if " are removed.
  # MORE ELEGANT WAY TO DO THIS?
  final.df <- cbind(if(exists("phylum.df")) phylum.df else 0,
                    if(exists("subphylum.df")) subphylum.df else 0,
                    if(exists("class.df")) class.df else 0,
                    if(exists("subclass.df")) subclass.df else 0,
                    if(exists("order.df")) order.df else 0,
                    if(exists("suborder.df")) suborder.df else 0,
                    if(exists("family.df")) family.df else 0,
                    if(exists("subfamily.df")) subfamily.df else 0,
                    if(exists("tribe.df")) tribe.df else 0,
                    if(exists("genus.df")) genus.df else 0,
                    if(exists("species.df")) species.df else 0)
  final.df <- final.df[,!grepl("if ", names(final.df))]
  return(final.df)
}
