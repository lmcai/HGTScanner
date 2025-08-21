library(taxizedb)
library(dplyr)
library(stringr)
library(pbapply)  # for progress bar
library(rentrez)
library(tidyr)

# Load local NCBI database
db_load_ncbi()

# Read species list
species <- readLines("species.txt")

get_family_order <- function(sp_name) {
  # Extract genus
  genus <- word(sp_name, 1)
  
  # Query genus in local NCBI database
  res <- tryCatch(classification(genus, db="ncbi")[[1]], error=function(e) NULL)
  
  if (!is.null(res) && "rank" %in% colnames(res)) {
    family <- res$name[res$rank == "family"]
    order  <- res$name[res$rank == "order"]
    
    if (length(family) == 0) family <- NA
    if (length(order) == 0)  order  <- NA
    
    return(data.frame(Species=sp_name,
                      Genus=genus,
                      Family=family,
                      Order=order,
                      stringsAsFactors = FALSE))
  } else {
    return(data.frame(Species=sp_name,
                      Genus=genus,
                      Family=NA,
                      Order=NA,
                      stringsAsFactors = FALSE))
  }
}

# Use pbapply::pblapply to show progress
results <- pbapply::pblapply(species, get_family_order)

# Combine results
taxonomy_df <- do.call(rbind, results)

# Save to CSV
write.csv(taxonomy_df, "species_genus_family_order.csv", row.names=FALSE)

#FOr species without family and order assignment in NCBI:
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("wcornwell/taxonlookup")
library(taxonlookup)
res=lookup_table(species$V1, by_species=TRUE)



#################
#Mitochondria
install.packages("taxizedb")
library(purrr)
library(taxizedb)
db_download_ncbi()
src <- src_ncbi()
get_taxid_from_acc <- function(acc) {
  tryCatch({
    s <- rentrez::entrez_summary(db = "nuccore", id = acc)
    return(s$taxid)
  }, error = function(e) {
    message("No record for accession: ", acc)
    return(NA)   # return NA instead of error
  })
}

acc_list=read.csv('ncbi_id.list')

taxids <- pbapply::pblapply(a, get_taxid_from_acc)
taxids <- as.character(unlist(taxids))
names(taxids) <- a
tax_class <- classification(taxids, db="ncbi", src=src)
# remove NULL or non-data.frame entries
tax_class_clean <- keep(tax_class, ~is.data.frame(.x))
tax_df <- bind_rows(tax_class_clean, .id="taxid") %>%
  left_join(
    data.frame(taxid=as.character(taxids),
               accession=names(taxids),
               stringsAsFactors = FALSE),
    by="taxid"
  )
  
tax_wide <- tax_df %>%
  select(accession, taxid, rank, name) %>%
  distinct() %>%
  pivot_wider(names_from = rank, values_from = name)

tax_flat <- tax_wide
tax_flat[] <- lapply(tax_flat, function(x) {
  if (is.list(x)) sapply(x, toString) else x
})

write.csv(tax_flat, "test.csv", row.names = FALSE)