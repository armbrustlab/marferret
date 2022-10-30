
### 1. Aggregation of taxonomic and functional annotations

### this will go under Stephen's best_pfam.py script


# merged pfam data saved out here:
PFAM_DAT="${MARFERRET_DIR}/data/pfam/MarFERReT.v1.Pfam_dat.csv"
# with the following fields:
"ref_id,aa_id,pfam,pfam_eval,pfam_score"

output file (pretty large file with the best pfam annotations for all entries):
PFAM_DAT="${MARFERRET_DIR}/data/pfam/MarFERReT.v1.Pfam_dat.csv"
with the following fields:
"ref_id,aa_id,pfam,pfam_eval,pfam_score"

PFAM_DAT="${MARFERRET_DIR}/data/pfam/MarFERReT.v1.Pfam_dat.csv"

### 2. Calculation of CTG inventory in transcript bins

#### Core Transcribed Genes - MarFERReT v1.0 ####
library(dplyr)
library(reshape2)
library(ggplot2)

# Define local MarFERReT dir:
mft_dir = "/Users/rgroussman/Dropbox/Armbrust/EukRefDB/"

# Load entry metadata:
mft_dat = read.csv(paste0(mft_dir,"data/MarFERReT.entry_metadata.v1.csv"), stringsAsFactors = FALSE)
# Load the taxonomic relationships file (taxa_csv_):
TAXA_CSV = read.csv(paste0(mft_dir,"data/ncbi/MarFERReT.v1.taxa.csv"), stringsAsFactors = FALSE)
# Load in Pfam annotation data::
pfam_dat = read.csv(paste0(mft_dir,"data/pfam/MarFERReT.v1.Pfam_dat.csv"), stringsAsFactors = FALSE, header=TRUE)

# Group and summarize pfam_dat by ref_id and pfam to get number of contigs (n_contigs) and reduce dataframe size.
pfam_dat2 = pfam_dat %>% group_by(ref_id,pfam) %>% summarize(n=n())

# Join in metadata to join the tax_ID and name fields by ref_id:
pfam_dat2 = left_join(pfam_dat2, select(mft_dat, ref_id, marferret_name, tax_id, data_type), by="ref_id")

# Filter this dataframe to include only transcriptome entries,
# (TSA or SAT), removing genomic entries:
	# Keep only rows where data_type == "TSA" or "SAT"
pfam_dat2 = pfam_dat2 %>% filter(data_type %in% c("TSA","SAT"))

# Reduce the df further by restricting to entries with species nodes
# (annotated to species-level or below)
pfam_dat2 = left_join(pfam_dat2, select(TAXA_CSV, tax_id, species), by="tax_id")
# Collapse subspecies and multiple entries within species taxids into species:
pfam_dat2 = pfam_dat2 %>% group_by(species, pfam) %>% summarize(n=sum(n), marferret_name=first(marferret_name)) %>% rename(tax_id = species)

### Assign the high-level lineage
# These lineages are spread across multiple phylogenetic ranks
# and contain nested taxa; the order of assignment matters.
pfam_dat2 = left_join(pfam_dat2, select(TAXA_CSV, tax_id, class, phylum, superkingdom, superkingdom_, superkingdom__, superkingdom___), by="tax_id")
pfam_dat2 = pfam_dat2 %>% mutate(lineage = case_when(class == '2864' ~ 'Dinophyceae',
  phylum == '2836' ~ 'Bacillariophyta',  phylum == '5878' ~ 'Ciliophora', phylum == '2830' ~ 'Haptophyta',
  phylum == '3041' ~ 'Chlorophyta', superkingdom__ == '543769' ~ 'Rhizaria',
	superkingdom___ == '2696291' ~ 'Ochrophyta', superkingdom_ == '554915' ~ 'Amoebozoa', superkingdom_ == '33154' ~ 'Opisthokonta',
  superkingdom == '2759' ~ 'Eukaryota', TRUE ~ 'Other'))

# Get the unique subset of tax_id and pfam_id from here::
uniq_species_and_pfams <- pfam_dat2 %>% select(tax_id, pfam_id) %>% unique()
dim(uniq_species_and_pfams) # [1] 1169450

# Prepare to create the presence/absence matrix;
# Positive tax_id pfam_id combinations have '1'
uniq_species_and_pfams$presence = 1
# cast the presence/absence matrix using pfam_id and tax_id;
# this fills out the presence/absence matrix with '0' values for absent tax_id/pfam_id pairs
pfam_presence_by_species.tsa = dcast(uniq_species_and_pfams, pfam_id ~ tax_id, value.var = 'presence', fill = 0 )
# Save out the file (optional)
write.csv(pfam_presence_by_species.tsa, paste0(mft_dir,"data/pfam/MarFERReT_v1.pfam_by_species.tsa.csv"), row.names = FALSE, quote = FALSE)

# Get a vector of pfam_id
v_pfam <- pfam_presence_by_species.tsa$pfam_id
# Transpose the dataframe, removing the 'pfam_id' for a
species_by_pfam_presence <- t(pfam_presence_by_species.tsa[,-1])
colnames(species_by_pfam_presence) <- v_pfam # set column names
# convert to dataframe
species_by_pfam_presence <- as.data.frame(species_by_pfam_presence)
# get rowsums (number of pfam identified in a species bin)
species_by_pfam_presence$rowsum = rowSums(species_by_pfam_presence)
# set 'tax_id' column from rownames
species_by_pfam_presence$tax_id = rownames(species_by_pfam_presence)

#### HISTOGRAM OF N PFAM PER SPECIES TAX_ID ####
# Plot a histogram to show the number of Pfams in species bins
ggplot(species_by_pfam_presence, aes(rowsum)) +
  geom_histogram() + theme_bw() + geom_vline(xintercept=1000,col="red",linetype="dashed") +
  ylab("Count") + xlab("Number of Pfam functions found in species")

#### SELECT SPECIES WITH >1000 PFAMS ####

# Filter the species*pfam dataframe to species with at least 1000 Pfams.
# Get a list of the species tax_ids with >= 1000 Pfam functions:
species_1kpfams = species_by_pfam_presence %>% filter(rowsum >= 1000 & tax_id != "NA") %>% rownames() # get a vector of the tax_ids with >=900 knums

# Return to the previous dataframe of pfam_id*tax_id presence:
# remove NA columns:
pfam_presence_by_species.tsa = pfam_presence_by_species.tsa %>% select(-"NA")
rownames(pfam_presence_by_species.tsa) = pfam_presence_by_species.tsa$pfam_id

# then keep only the columns (species) with at least 1000 Pfams:
pfam_presence_by_species.tsa = pfam_presence_by_species.tsa %>% select(one_of(species_1kpfams))
dim(pfam_presence_by_species.tsa) # [1] 11247   368

# Assign lineage to the entra metadata file:
mft_dat = mft_dat %>% left_join(ref_dat_taxid, by="tax_id")

mft_dat = left_join(mft_dat, select(TAXA_CSV, tax_id, species, class, phylum, superkingdom, superkingdom_, superkingdom__, superkingdom___), by="tax_id")

mft_dat = mft_dat %>% mutate(lineage = case_when(class == '2864' ~ 'Dinophyceae',
  phylum == '2836' ~ 'Bacillariophyta',  phylum == '5878' ~ 'Ciliophora', phylum == '2830' ~ 'Haptophyta',
  phylum == '3041' ~ 'Chlorophyta', superkingdom__ == '543769' ~ 'Rhizaria',
	superkingdom___ == '2696291' ~ 'Ochrophyta', superkingdom_ == '554915' ~ 'Amoebozoa', superkingdom_ == '33154' ~ 'Opisthokonta',
  superkingdom == '2759' ~ 'Eukaryota', TRUE ~ 'Other'))


#### ITERATE THROUGH LINEAGES TO GENERATE CTG CATALOG ####

# Now we use these lineages to iterate through pfam_presence_by_species.tsa and generate a combined df:

# initialize the  df:
ctg_dat = data.frame(matrix(ncol = 4, nrow = 0), stringsAsFactors=FALSE) %>% mutate(across(everything(), as.character))
colnames(ctg_dat)[1:4] = c("lineage","n_species","pfam_id","frequency")

# for each lineage:
for (this_lineage in lineages) {
  # get vector of species tax_id in lineage:
  lineage_species = mft_dat %>% filter(lineage == this_lineage) %>% ungroup() %>% pull(species)
  # subset pfam_presence_by_species.tsa select species within the lineage:
  pfam_sub = pfam_presence_by_species.tsa %>% select(one_of(as.character(lineage_species)))
  # Get the rowsum (number of species tax_id each Pfam is observed in)
  pfam_sub$rowsum = rowSums(pfam_sub)
  #	Divide the number of Pfams observed by the number of species to get a frequency of observation:
  pfam_sub$frequency = pfam_sub$rowsum / (ncol(pfam_sub)-1) # number of columns minus one (rowsum) = number of taxids

  # Generate a figure: Pfam frequency by species
  xlabel = paste("Proportion of lineage ", this_lineage, sep="")
  pfam_histo <- ggplot(pfam_sub, aes(frequency)) +
    geom_histogram(binwidth = 0.025) + ggtitle("") + xlab(xlabel) + ylab("Frequency") +
    geom_vline(xintercept=0.95, color="red") + theme_bw()
  pfam_histo
  out_image = paste0("Pfam_frequency.", this_lineage, ".png")
  ggsave(out_image, width = 5, height = 4, units = "in")

  #	keep the pfam_id that are found in at least 95% of the species of this group;
	# these are the core transcribed genes of this lineage:
  core_pfams_subset <- pfam_sub[pfam_sub$frequency >= 0.95,] %>% rownames()
  # prepare output:
  pfam_sub = pfam_sub[pfam_sub$frequency >= 0.95,]
  pfam_sub$pfam_id = rownames(pfam_sub)
  pfam_sub = pfam_sub %>% relocate(c("pfam_id", "rowsum", "frequency"))
  out_csv = paste0(mft_dir,"data/ctg/MarFERReT.pfam_presence.", this_lineage, ".ctg_catalog.v1.csv")
  write.csv(pfam_sub, out_csv, row.names = FALSE, quote = FALSE)

  # Join to the larger table, ctg_dat:
  pfam_sub$lineage = this_lineage
  ctg_dat_sub = pfam_sub %>% select(lineage, n_species = rowsum, pfam_id, frequency) %>% mutate(across(everything(), as.character))
  ctg_dat = bind_rows(ctg_dat, ctg_dat_sub)
}

#### all-Eukaryota CTGs ####

# Now generate the CTGs for the set of eukaryota as a whole:
this_lineage = "Eukaryota"

# Do the same as above but for all the eukaryote transcriptomes:
pfam_sub = pfam_presence_by_species.tsa
# Calculate the core pfams
pfam_sub$rowsum = rowSums(pfam_sub)
#	Divide the number of Pfams observed by the number of species to get a frequency of observation:
pfam_sub$frequency = pfam_sub$rowsum / (ncol(pfam_sub)-1) # number of columns minus one (rowsum) = number of taxids

# Generate a figure: Pfam frequency by species
xlabel = paste("Proportion of lineage ", this_lineage, sep="")
pfam_histo <- ggplot(pfam_sub, aes(frequency)) +
  geom_histogram(binwidth = 0.025) + ggtitle("") + xlab(xlabel) + ylab("Frequency") +
  geom_vline(xintercept=0.95, color="red") + theme_bw()
pfam_histo
out_image = paste0("/Users/rgroussman/Dropbox/Armbrust/EukRefDB/figs/ctg/supp/Pfam_frequency.", this_lineage, ".png")
ggsave(out_image, width = 5, height = 4, units = "in")

#	keep the pfam_id that are found in at least 95% of the species of this group;
# these are the core transcribed genes of this lineage:
core_pfams_subset <- pfam_sub[pfam_sub$frequency >= 0.95,] %>% rownames()
length(core_pfams_subset)

# prepare output:
pfam_sub = pfam_sub[pfam_sub$frequency >= 0.95,]
pfam_sub$pfam_id = rownames(pfam_sub)
pfam_sub = pfam_sub %>% relocate(c("pfam_id", "rowsum", "frequency"))
out_csv = paste0(mft_dir,"data/ctg/MarFERReT.pfam_presence.", this_lineage, ".ctg_catalog.v1.csv")
write.csv(pfam_sub, out_csv, row.names = FALSE, quote = FALSE)

# Join to the larger table, ctg_dat:
pfam_sub$lineage = this_lineage
ctg_dat_sub = pfam_sub %>% select(lineage, n_species = rowsum, pfam_id, frequency) %>% mutate(across(everything(), as.character))
ctg_dat = bind_rows(ctg_dat, ctg_dat_sub)

#### OUTPUT CTG CATALOG ####

# write the list of core Pfams and their frequencies:
write.csv(ctg_dat, paste0(mft_dir,"data/ctg/MarFERReT.core_genes.v1.csv"), row.names = FALSE, quote = FALSE)
