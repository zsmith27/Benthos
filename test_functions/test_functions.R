

new_pct_taxon <- function(long.df, unique.id.col, taxon.rank, taxon, count.col) {
  sapply(unique(long.df[, unique.id.col]), function(x){
    sub.df <- long.df[long.df[, unique.id.col] %in% x, ]
    if (!taxon %in% sub.df[, taxon.rank]) {
      0
    } else {
      sum(sub.df[sub.df[, taxon.rank] %in% taxon, count.col]) /
        sum(sub.df[, count.col]) * 100
    }
  })
}
library(dplyr)
unique.id.col <- quo(UNIQUE_ID)
taxon.rank <- quo(GENUS)
taxon.x <- quo(EPEORUS)
count.col <- quo(REPORTING_VALUE)
pct_taxon <- function(long.df, unique.id.col, count.col, taxon.col, taxon) {
  # Prep.
  unique.id.col <- enquo(unique.id.col)
  taxon.col <- enquo(taxon.col)
  count.col <- enquo(count.col)
  #----------------------------------------------------------------------------
  # Calculate the percentage of the specified taxon.
  final.vec <- long.df %>% 
    group_by(!!unique.id.col) %>% 
    summarise(TOTAL = sum(!!count.col),
           INDV = sum(UQ(count.col)[UQ(taxon.col) %in% taxon]),
           PCT = INDV / TOTAL * 100) %>% 
    pull(PCT)
  #----------------------------------------------------------------------------
  return(final.vec)
}

data("onondaga")

long.df <- onondaga




test2 <- long.df %>% 
  distinct(UNIQUE_ID) %>% 
  arrange(UNIQUE_ID) %>% 
  mutate(CHIRO = pct_taxon(long.df, UNIQUE_ID, REPORTING_VALUE, FAMILY, "CHIRONOMIDAE"),
         EPHEM = new_pct_taxon2(long.df, UNIQUE_ID, REPORTING_VALUE, ORDER, c("EPHEMEROPTERA", "PLECOPTERA", "TRICHOPTERA")),
         EPEORUS = new_pct_taxon2(long.df, UNIQUE_ID, REPORTING_VALUE, GENUS, "EPEORUS"),
         HEXAPODA =new_pct_taxon2(long.df, UNIQUE_ID, REPORTING_VALUE, SUBPHYLUM, "HEXAPODA"),
         gold = gold(long.df,
                     unique.id.col = UNIQUE_ID,
                     count.col = REPORTING_VALUE,
                     class.col = CLASS,
                     order.col = ORDER))

test <- onondaga %>% 
  filter(UNIQUE_ID == "CAZ_1_1") %>% 
  mutate(total = sum(REPORTING_VALUE)) %>% 
  group_by(FAMILY) %>% 
  summarise(PCT = sum(REPORTING_VALUE) / unique(total) * 100)
  

