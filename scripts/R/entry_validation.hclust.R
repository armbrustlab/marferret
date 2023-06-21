# AUTHOR: Ryan Groussman

# Purpose: Tree-based visualization of protein-coding content
# in candidate entries for the MarFERReT reference sequence library. 
# This code inputs hmmsearch-based Pfam protein family annotations 
# along with metadata and taxonomy information to conduct hierarchical clustering
# of entries through a binary matrix of Pfam protein-coding content. 

#### Dependencies and Data ####

# load package dependencies:
library("ggplot2")
library("tidyverse")
library("reshape2")
library("ggdendro")
library("tidytree")
library("ggtree")
library("ape")

## load data files:
# define the local MarFERReT dir:
mft_dir = "/home/user/data/marferret/" # replace with local installation directory

# load entry curation metadata:
mft_dat = read.csv(paste0(mft_dir,"data/MarFERReT.v1.entry_curation.csv"))

# load the NCBI taxonomy relationship file:
TAXA_CSV = read.csv(paste0(mft_dir,"data/MarFERReT.v1.taxa.csv"), stringsAsFactors = FALSE)

# This block assigns a 'lineage' value and requires some fine-tuning as the
# lineage values are chosen from different taxonomic rank columns
# assign lineage to taxa.csv:
TAXA_CSV = TAXA_CSV %>% mutate(lineage = case_when(class == '2864' ~ 'Dinophyceae',
      class == '3027' ~ 'Cryptophyceae',
      phylum == '2836' ~ 'Bacillariophyta',
      phylum == '5878' ~ 'Ciliophora',
      phylum == '2830' ~ 'Haptophyta',
      phylum == '3041' ~ 'Chlorophyta',
      phylum == '2763' ~ 'Rhodophyta',
      superkingdom___ == '2696291' ~ 'Ochrophyta',
      superkingdom__ == '543769' ~ 'Rhizaria',
      superkingdom__ == '33630' ~ 'Alveolata',
      superkingdom__ == '33634' ~ 'Stramenopiles',
      superkingdom_ == '2611352' ~ 'Discoba',
      superkingdom_ == '554915' ~ 'Amoebozoa',
      superkingdom_ == '33154' ~ 'Opisthokonta',
      TRUE ~ 'Other'))

# Load the entry Pfam annotations for all entries, unclustered
aa_seqs_pfam_df = read_csv("data/MarFERReT.v1.aa_seqs.Pfam34.annotations.csv.gz")
dim(aa_seqs_pfam_df) # [1] 12554711        4
aa_seqs_pfam_df$pfam_id %>% unique() %>% length() # 12549

# separate the handle to get the entry_id, then bring in the tax_id from metadata.csv
aa_seqs_pfam_df = aa_seqs_pfam_df %>% 
  separate(taxid, into = "ref_id", sep = "_", remove = TRUE, extra="drop")

## generate a presence/absence table from pfam_df
# get the unique subset of entriesand pfam_ids from here:
ref_pfam_df <- aa_seqs_pfam_df %>% select(ref_id, pfam_id) %>% unique()
dim(ref_pfam_df) # [1] 2502744       2

# set a 'presence' column for detected pairs of entries and pfams:
ref_pfam_df$presence = 1
# cast the presence/absence matrix using pfam_id and entry (ref_id)
ref_pfam_df_wide = dcast(ref_pfam_df, ref_id ~ pfam_id, value.var = 'presence', fill = 0 )
# At this point expect 899 rows and 12550 columns:
dim(ref_pfam_df_wide) # [1]   899 12550

# The format should look like this:
head(ref_pfam_df_wide)[1:5]
  # ref_id PF00001.23 PF00002.26 PF00003.24 PF00004.31
  # 1      1          1          0          1          1
  # 2     10          1          0          1          1
  # 3    100          0          0          1          1
  # 4    101          0          0          0          1
  # 5    102          0          0          0          1
  # 6    103          0          1          0          1

# Save out the file at this point if desired:
write.csv(ref_pfam_df_wide, paste0(mft_dir,"DataDescriptor/data/MarFERReT_v1.aa_seqs.ref_pfam_df.csv"), row.names = FALSE, quote = FALSE)


# Create ref_counts dataframe
ref_counts = ref_pfam_df_wide %>%
  rowwise() %>%
  mutate(counts = sum(c_across(starts_with("PF")))) %>%
  select(ref_id, counts)

dim(ref_counts) # [1] 899   2

# Histogram plots for showing the distribution of total Pfam counts in entries
ref_hist = ggplot(ref_counts, aes(counts)) +
  geom_histogram() + theme_bw() + #geom_vline(xintercept=1000,col="red",linetype="dashed") +
  ylab("Number of entries in bin") + xlab("N Pfams found in entry")

# Save plot to PNG file:
file_out = (paste0(mft_dir,"DataDescriptor/figs/taxid_pfam/mft.sci_dat.candidate_pfam_counts.v1.png"))
ggsave(file_out, width = 4, height = 3, units = "in")

# Create pfamid_counts dataframe
pfamid_counts <- ref_pfam_df_wide %>%
  select(starts_with("PF")) %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "pfam_id", values_to = "counts")

dim(pfamid_counts) # [1] 12549     2

# Histogram plots for showing the distribution of references possessing a given Pfam
pfam_hist = ggplot(pfamid_counts, aes(counts)) +
  geom_histogram() + theme_bw() + #geom_vline(xintercept=1000,col="red",linetype="dashed") +
  ylab("Number of Pfams in bin") + xlab("N of refs with a Pfam")

# Save plot to PNG file:
file_out = (paste0(mft_dir,"DataDescriptor/figs/taxid_pfam/mft.sci_dat.pfam_candidate_counts.v1.png"))
ggsave(file_out, width = 4, height = 3, units = "in")

# create the initial matrix, leaving out the ref_id column:
refid_pfam_matrix <- as.matrix(ref_pfam_df_wide %>% select(-ref_id))
# set the ref_id as the colname
rownames(refid_pfam_matrix) = ref_pfam_df_wide$ref_id

# construct the binary distance matrix from the initial matrix:
d = dist(x = refid_pfam_matrix, method = "binary")

# if desired, write out the dist object as a matrix:
write.csv(as.matrix(d), paste0(mft_dir,"DataDescriptor/data/pfam/899_candidate.dist.csv"), row.names = FALSE, quote = FALSE)


# conduct complete hierarchical clustering on the binary distance matrix
hc = hclust(d)
# Call:
#   hclust(d = d)
# 
# Cluster method   : complete 
# Distance         : binary 
# Number of objects: 899 

# convert the hclust object to a tree object
ref_mft_tree = as.phylo(hc)
# Phylogenetic tree with 899 tips and 898 internal nodes.


# Convert the tree object to a tidytree object using the as_tibble() function from the tidytree package
ref_mft_tree_df <- as_tibble(ref_mft_tree)


# add in taxid from metadata:
ref_mft_tree_df = ref_mft_tree_df %>% left_join(
  select(mft_dat,candidate_id,tax_id)%>% mutate(label = as.character(candidate_id)),by="label")


# add in taxid, name and lineage information by 'tax_id':
ref_mft_tree_df = ref_mft_tree_df %>% left_join(
  select(TAXA_CSV,tax_id,rank,tax_name,lineage),by="tax_id")

# convert mft_tree_df back to tree data:
ref_mft_tree2 = ref_mft_tree_df %>% as.treedata() 
# Phylogenetic tree with 899 tips and 898 internal nodes.

# Declare the colors for each lineage in the figure
fig_cols = c("Other"= "#999999",
             "Amoebozoa"= "#88419d",
             "Alveolata" = "#80b1d3",
             "Opisthokonta"="#e41a1c",
             "Discoba" = "#9ebcda", 
             "Cryptophyceae" = "#8dd3c7", "Rhodophyta" = "#fb8072",
             "Stramenopiles" = "#fee391",
             "Chlorophyta"= "#4daf4a",
             "Haptophyta"= "#f781bf",
             "Ochrophyta"= "#fe9929",
             "Rhizaria" = "#b15928",
             "Dinophyceae"="#045a8d",
             "Bacillariophyta"= "#cc4c02",
             "Ciliophora" = "#984ea3")

# Create a smaller version of the tree plot with labels hidden:
ref_mft_pfam_tree <- ggtree(ref_mft_tree2, layout="rectangular", aes(color=lineage)) + 
  scale_color_manual("Lineage", values=fig_cols) + 
  guides(color=guide_legend(title="Lineage")) + 
  theme(legend.position = 'left', legend.background = element_rect()) 
# save the unlabeled tree out as a PNG
file_out = (paste0(mft_dir,"DataDescriptor/figs/taxid_pfam/draft/mft.sci_dat.pfam_counts.entry_tree1.png"))
ggsave(file_out, width = 4, height = 5, units = "in")

#### ADD FLAGS FOR 899 CANDIDATE ENTRIES ####

# follow these steps to add the flagging data to the plot from the entry_curation.csv

# For the labeled tree, alter the df to create legible labels for 
# full tree visualization, setting tax_name + tax_id as the label
tree_df = ref_mft_tree_df %>% mutate(label = case_when(!is.na(tax_name) ~ paste0(label, '_', tax_name), TRUE ~ label))

# load in the flag data from the curation csv
cur_dat = read_csv("data/MarFERReT.v1.entry_curation.csv")

# set FLAGGED = True if flagged, False otherwise:
cur_dat = cur_dat %>% mutate(FLAGGED = case_when(
  is.na(FLAG) ~ NA,
  !is.na(FLAG) ~ TRUE
))

# Investigate number of flags
cur_dat %>% filter(FLAGGED==TRUE) %>% dim() # [1] 115  18
cur_dat %>% filter(FLAGGED==FALSE) %>% dim() # [1] 0  18
cur_dat %>% filter(is.na(FLAGGED)) %>% dim() # [1] 787  18

## do for each of the three flag types:
# outer circle: tree identification
cur_dat = cur_dat %>% mutate(HCLUST_FLAGGED = case_when(
  is.na(HCLUST_FLAG) ~ NA, !is.na(HCLUST_FLAG) ~ TRUE))
cur_dat %>% filter(HCLUST_FLAGGED==TRUE) %>% dim() # [1] 54
cur_dat %>% filter(is.na(HCLUST_FLAGGED)) %>% dim() # [1] 848

# inner circle: vlierberghe
cur_dat = cur_dat %>% mutate(MMETSP_FLAGGED = case_when(
  is.na(FLAG_VanVlierberghe) ~ NA, !is.na(FLAG_VanVlierberghe) ~ TRUE))
cur_dat %>% filter(MMETSP_FLAGGED==TRUE) %>% dim() # [1] 82
cur_dat %>% filter(is.na(MMETSP_FLAGGED)) %>% dim() # [1] 820

# Lasek-Nesselquist (ciliates only)
cur_dat = cur_dat %>% mutate(CILIATE_FLAGGED = case_when(
  is.na(`FLAG_Lasek-Nesselquist`) ~ NA, !is.na(`FLAG_Lasek-Nesselquist`) ~ TRUE))
cur_dat %>% filter(CILIATE_FLAGGED==TRUE) %>% dim() # [1] 18
cur_dat %>% filter(is.na(CILIATE_FLAGGED)) %>% dim() # [1] 884

# we will use this to filter whether or not an entry gets a labeled

# join flag data to tree table
tree_df = tree_df %>% left_join(
  select(cur_dat,ref_id,candidate_id,FLAGGED,HCLUST_FLAGGED,MMETSP_FLAGGED,CILIATE_FLAGGED),by="candidate_id")

# convert table back to ggtree object
flag_tree = tree_df %>% as.treedata()

# Create a version of the plot with all points labeled:
flag_tree_plt <- ggtree(flag_tree, layout="circular", aes(color=lineage)) + 
  # scale lineage color
  scale_color_manual("Lineage", values=fig_cols) +
  # Add HCLUST_FLAGGED layer
  geom_point(aes(shape = HCLUST_FLAGGED), na.rm = TRUE, fill = '#008B8B', stroke=0.01, size = 1.5, position = position_nudge(x = .01)) +
  # Add MMETSP_FLAGGED layer
  geom_point(aes(shape = MMETSP_FLAGGED), na.rm = TRUE, fill = '#FF6347', stroke=0.01, size = 1.5, position = position_nudge(x = .025)) +
  # Add CILIATE_FLAGGED layer
  geom_point(aes(shape = CILIATE_FLAGGED), na.rm = TRUE, fill = "#FFD700", stroke=0.01, size = 1.5, position = position_nudge(x = .04)) +
  scale_shape_manual(values=21) +
  
  guides(color=guide_legend(title="Lineage")) + 
  theme(legend.position = 'right', legend.background = element_rect()) 

flag_tree_plt

# save out PNG or PDF:
file_out = (paste0(mft_dir,"DataDescriptor/figs/taxid_pfam/mft.sci_dat.pfam_counts.entry_tree_flagged.v2.png"))
ggsave(file_out, width = 6, height = 6, units = "in")
file_out = (paste0(mft_dir,"DataDescriptor/figs/taxid_pfam/mft.sci_dat.pfam_counts.entry_tree_flagged.v2.pdf"))
ggsave(file_out, width = 6, height = 6, units = "in")

# create a variant of the figure for the legend block:

# Create a small data frame to hold the legend information
df_leg <- data.frame(
  x = numeric(0), 
  y = numeric(0), 
  ValidationSource = factor(c(), levels = c("hclust", "VanVlierberghe", "Lasek-Nesselquist"))
)

# Create the ggplot object
legend_plot = ggplot(df_leg, aes(x=x, y=y, color=ValidationSource)) +
  geom_blank() +
  scale_color_manual(name = "Validation source",
                     values = c("hclust" = "#008B8B", "VanVlierberghe" = "#FF6347", "Lasek-Nesselquist" = "#FFD700")) +
  theme(legend.position = "bottom")

legend_plot
