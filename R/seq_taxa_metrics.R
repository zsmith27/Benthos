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
  if(any(rowSums(wide.df[, 8:ncol(wide.df)]) == 0)) {
    with.counts <- wide.df[rowSums(wide.df[, 8:ncol(wide.df)]) != 0, ]
    with.counts[, 8:ncol(with.counts)] <- (with.counts[, 8:ncol(with.counts)] /
                                     rowSums(with.counts[, 8:ncol(with.counts)])) * 100
    without.counts <- wide.df[rowSums(wide.df[, 8:ncol(wide.df)]) == 0, ]
    
    wide.df <- rbind(with.counts, without.counts)
    
  } else {
    wide.df[, 8:ncol(wide.df)] <- (wide.df[, 8:ncol(wide.df)] /
                                     rowSums(wide.df[, 8:ncol(wide.df)])) * 100
  }
  
  t_list <- unique(toupper(master.df[, taxa.rank]))
  shorten <- wide.df[, colnames(wide.df) %in% t_list]
  cn <- colnames(shorten)
  list_taxa.df <- if(length(cn) == 0){
    list_taxa.df <- data.frame(wide.df[, c("UNIQUE_ID", "STATION_ID",
                                           "AGENCY_CODE", "DATE",
                                           "METHOD", "SAMPLE_NUMBER", 
                                           "CONDITION")])
    list_taxa.df$NO_MATCH <- 0
    colnames(list_taxa.df) <- c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE",
                                "DATE", "METHOD", "SAMPLE_NUMBER",
                                "CONDITION",
                                "NO_MATCH")
    list_taxa.df
  }else{
    list_taxa.df <- data.frame(cbind(wide.df[, 1:7], shorten))
    colnames(list_taxa.df) <- c("UNIQUE_ID", "STATION_ID", "AGENCY_CODE",
                                "DATE", "METHOD", "SAMPLE_NUMBER",
                                "CONDITION",
                                 paste("PCT", cn, sep = "_"))
    list_taxa.df
  }
  final.df <- data.frame(list_taxa.df)
  return(data.frame(list_taxa.df[order(list_taxa.df$UNIQUE_ID),]))
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
  
  pct_phylum <- if(length(phylum) > 8) calc_pct_taxa(phylum, taxa.rank = "PHYLUM", master.df)
  pct_subphylum <- if(length(subphylum) > 8) calc_pct_taxa(subphylum, "SUBPHYLUM", master.df)
  pct_class <- if(length(class) > 8) calc_pct_taxa(class, "CLASS", master.df)
  pct_subclass <- if(length(subclass) > 8) calc_pct_taxa(subclass, "SUBCLASS", master.df)
  pct_order <- if(length(order) > 8) calc_pct_taxa(order, "ORDER", master.df)
  pct_suborder <- if(length(suborder) > 8) calc_pct_taxa(suborder, "SUBORDER", master.df)
  pct_family <- if(length(family) > 8) calc_pct_taxa(family, "FAMILY", master.df)
  pct_subfamily <- if(length(subfamily) > 8) calc_pct_taxa(subfamily, "SUBFAMILY", master.df)
  pct_tribe <- if(length(tribe) > 8) calc_pct_taxa(tribe, "TRIBE", master.df)
  pct_genus <- if(length(genus) > 8) calc_pct_taxa(genus, "GENUS", master.df)
  pct_species <- if(length(species) > 8) calc_pct_taxa(species, "SPECIES", master.df)
  
  check_exists <- function(pct_taxa, long.df){
    if(is.null(pct_taxa)){
      pct_taxa <- unique(long.df[, c("UNIQUE_ID", "STATION_ID",
                                     "AGENCY_CODE","DATE",
                                     "METHOD", "SAMPLE_NUMBER",
                                     "CONDITION")])
    }
    return(pct_taxa)
  }
  
  checked_pct_phylum <- check_exists(pct_phylum, long.df)
  checked_pct_subphylum <- check_exists(pct_subphylum, long.df)
  checked_pct_class <- check_exists(pct_class, long.df)
  checked_pct_subclass <- check_exists(pct_subclass, long.df)
  checked_pct_order <- check_exists(pct_order, long.df)
  checked_pct_suborder <- check_exists(pct_suborder, long.df)
  checked_pct_family <- check_exists(pct_family, long.df)
  checked_pct_subfamily <- check_exists(pct_subfamily, long.df)
  checked_pct_tribe <- check_exists(pct_tribe, long.df)
  checked_pct_genus <- check_exists(pct_genus, long.df)
  checked_pct_species <- check_exists(pct_species, long.df)
  
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
                                  checked_pct_species), by =  c("UNIQUE_ID", "STATION_ID",
                                                                "AGENCY_CODE","DATE",
                                                                "METHOD", "SAMPLE_NUMBER",
                                                                "CONDITION"),
                             type = "full")
  comb_taxa <- comb_all[, 8:ncol(comb_all)]
  rm.cols <- names(comb_taxa[, colSums(comb_taxa) == 0])
  final.df <- comb_all[, !(names(comb_all) %in% rm.cols)]
  
  return(final.df)
}

#==============================================================================
#'Richness of All Taxon
#'
#'@param long.df = Data in a long data format
#'@param low.res.rank = Lower resolution taxonomic rank.
#'@param high.res.rank = Higher resolution taxonomic rank.
#'@param master.df = Master Taxa List.
#'@return Sequences through each taxon at level of the specified low.res.rank
#'and calculates richness based on the observations at the specified
#' high.res.rank.
#'@export
#'
rich_by_rank <- function(long.df, low.res.rank, high.res.rank, master.df) {
  taxa.org <- unique(long.df[, low.res.rank])
  taxa.vec <- taxa.org[taxa.org %in% unique(master.df[, low.res.rank])]
  rich.list <- lapply(taxa.vec, function(x){
    print(paste("...RICH", x, sep = "_"))
    final.list <- taxon_richness(long.df, x, low.res.rank, high.res.rank)
  return(final.list)
  } )
  final.df <- data.frame(do.call(cbind, rich.list))
  names(final.df) <- paste("RICH", taxa.vec, sep = "_")
  return(final.df)
}
#==============================================================================
#'Sequence Through and Calculate the Richness of All Taxon
#'
#'@param long.df = Data in a long data format
#'@param rank = The taxonomic resolution of your dataset.
#'@param master.df = Master Taxa List.
#'@return Sequences through each taxonomic hierarchy and taxon, calculating the
#' richness of tax at the specified taxonomic rank ('''rank''').
#'@export
#'
seq_taxa_rich <- function(long.df, rank = "GENUS", master.df){

    taxa.cols <- c("PHYLUM", "SUBPHYLUM", "CLASS",
                   "SUBCLASS", "ORDER", "SUBORDER",
                   "FAMILY", "SUBFAMILY", "TRIBE",
                   "GENUS", "SPECIES")
    taxa.final <- names(long.df)[which(names(long.df) %in% taxa.cols)]
    taxa.final <- taxa.final[-length(taxa.final)]


  #============================================================================
  
  
  phylum.df <- if("PHYLUM" %in% taxa.final){
    print("...Phylum Richness")
    rich_by_rank(long.df, "PHYLUM", rank, master.df)
  } 
  subphylum.df <- if("SUBPHYLUM" %in% taxa.final){
    print("...Subphylum Richness")
    rich_by_rank(long.df, "SUBPHYLUM", rank, master.df)
  } 
  class.df <- if("CLASS" %in% taxa.final){
    print("...Class Richness")
    rich_by_rank(long.df, "CLASS", rank, master.df)
  } 
  subclass.df <- if("SUBCLASS" %in% taxa.final){
    print("...Subclass Richness")
    rich_by_rank(long.df, "SUBCLASS", rank, master.df)
  } 
  order.df <- if("ORDER" %in% taxa.final){
    print("...Order Richness")
    rich_by_rank(long.df, "ORDER", rank, master.df)
  } 
  suborder.df <- if("SUBORDER" %in% taxa.final){
    print("...Suborder Richness")
    rich_by_rank(long.df, "SUBORDER", rank, master.df)
  } 
  family.df <- if("FAMILY" %in% taxa.final){
    print("...Family Richness")
    rich_by_rank(long.df, "FAMILY", rank, master.df)
  } 
  subfamily.df <- if("SUBFAMILY" %in% taxa.final){
    print("...Subfamily Richness")
    rich_by_rank(long.df, "SUBFAMILY", rank, master.df)
  } 
  tribe.df <- if("TRIBE" %in% taxa.final){
    print("...Tribe Richness")
    rich_by_rank(long.df, "TRIBE", rank, master.df)
  } 
  genus.df <- if("GENUS" %in% taxa.final){
    print("...Genus Richness")
    rich_by_rank(long.df, "GENUS", rank, master.df)
  } 

  #============================================================================
  # Bind columnes together for the data frames that exist. Data frames that
  # do not exist are first filled with zeros and recieve a name "if "...
  # After the columns are appended, columns containing "if " are removed.
  # MORE ELEGANT WAY TO DO THIS?
  num.rows <- nrow(unique(long.df[, 1:7]))
  final.df <- cbind(if(!is.null(phylum.df)) phylum.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(subphylum.df)) subphylum.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(class.df)) class.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(subclass.df)) subclass.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(order.df)) order.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(suborder.df)) suborder.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(family.df)) family.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(subfamily.df)) subfamily.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(tribe.df)) tribe.df else data.frame(rep(NA, num.rows)),
                    if(!is.null(genus.df)) genus.df else data.frame(rep(NA, num.rows)))
  
  final.df <- final.df[, !grepl("rep.NA", names(final.df))]
  #final.df <- final.df[, !grepl("\\.", names(final.df))]
  return(final.df)
}
