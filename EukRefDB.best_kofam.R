#############
## 01/06/22 #
#############

# implemented in 
# "Function_DB.EukRefDB_2021.log.sh" #### Best-scoring Pfam ####

library(dplyr)

setwd("/mnt/nfs/projects/ryan/EukRefDB_2021/annotations/Pfam_34.0/combined/")
#setwd("/Users/rgroussman/data/EukRefDB/temp_test/best_kofam")

# set list from 1 to 891:
ref_ids = 1:891
#ref_ids = 1:10
# remove 665 (dead reference)
ref_ids = ref_ids[ref_ids!= 665]

# iterate through ref_ids:
for (ref_id in ref_ids) {
  print(ref_id) 
  # open ref_id(n).hmm_out.csv
  filename = paste(c("ref_id",toString(ref_id),".hmm_out.csv"),collapse="")
  outfile = paste(c("ref_id",toString(ref_id),".best_pfam.hmm_out.csv"),collapse="")
  pfam_dat = read.csv(file=filename, sep=',', header=TRUE)
  # group by aa_id and slice by maximum knum_score (it is the pfam score in this instance)
  best_pfam_dat <- pfam_dat %>% group_by(aa_id) %>% slice(which.max(knum_score))
  # save out results
  write.csv(best_pfam_dat, outfile, row.names = FALSE, quote = FALSE)
  
  # we have the best pfam per FRAME, now let's group further and get the best 
  # pfam per CONTIG. we'll slice off the aa_id to get nt_id:
  
  best_pfam_dat$nt_id = substr(best_pfam_dat$aa_id,1,nchar(as.character(best_pfam_dat$aa_id))-2)
  best_pfam_dat2 <- best_pfam_dat %>% group_by(nt_id) %>% slice(which.max(knum_score))
  
  # write out results:
  outfile_contig = paste(c("ref_id",toString(ref_id),".best_pfam_nt.hmm_out.csv"),collapse="")
  write.csv(best_pfam_dat2, outfile_contig, row.names = FALSE, quote = FALSE)
  
} 




#

