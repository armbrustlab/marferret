
# Collecting sequence data from source

A total of 902 candidate entry sequences were retrieved for validation and potential inclusion in the final MarFERReT build from four major contributors and a handful of smaller sequencing projects. We do not host the primary raw data in the final MarFERReT data products; the original links, raw filenames, and other source reference information for final data are available in the [MarFERReT metadata file](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/data/MarFERReT.v1.metadata.csv) and for all candidate entries in the [entry and taxID curation table](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/data/MarFERReT.v1.entry_curation.csv).

Below are instructions for downloading primary sequence data for candidate entries via command line (bash) or manual (graphic user interface) approaches. For replicating the MarFERReT build process from scratch; source sequences should be downloaded to the `data/source_seqs/` subdirectories unter the main MarFERReT directory locally. 

## Marine Microbial Eukaryote Sequence Project

The Marine Microbial Eukaryote Sequence Project (MMETSP, Keeling et al., 2014) is the largest contributor to MarFERReT. We collected the MMETSP sequences as peptide translations for 678 transcriptomes from Version 2 of the MMETSP re-assemblies ([Johnson et al., 2018](https://doi.org/10.5281/zenodo.740440)).

Download peptide translations of transcriptome assemblies from Zenodo using wget:

```
wget https://zenodo.org/record/3247846/files/mmetsp_dib_trinity2.2.0_pep_zenodo.tar.gz
```

Unpack the tarball:

```
tar -xvf mmetsp_dib_trinity2.2.0_pep_zenodo.tar.gz
```

## Zenodo repositories

A total of 15 candidate entries were retrieved from Zenodo repositories for three separate sequencing projects. 

Ten diatom transcriptomes from Guajardo et al., 2021 were downloaded from the associated [Zenodo repository](https://doi.org/10.5281/zenodo.4591037). Download nucleotide transcriptome assemblies from Zenodo using wget:

```
wget https://zenodo.org/record/4591037/files/ThTSP-01_Minidiscusspinulatus-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-02_RCC4665-Minidiscusvariabilis-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-03_RCC4660-Minidiscuscomicus-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-04_RCC4219-Thalassiosirasp-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-05_RCC4593-Thalassiosiraminima-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-06_RCC4590-Minidiscussp-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-07_RCC4582-Minidiscussp-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-08_RCC4606-Thalassiosirasp-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-09_RCC4583-Thalassiosiraminima-Trinity.fasta
wget https://zenodo.org/record/4591037/files/ThTSP-10_RCC4584-Minidiscussp-Trinity.fasta
```
Two haptophyte and two diatom transcriptomes from environmental isolates (Lambert et al., 2022) were retrieved from their [Zenodo repository](https://zenodo.org/record/7336407):

```
wget https://zenodo.org/record/7336407/files/ALTEMP.filtered.Trinity.fasta.gz
wget https://zenodo.org/record/7336407/files/AT125.filtered.Trinity.fasta.gz
wget https://zenodo.org/record/7336407/files/KBHA01.filtered.Trinity.fasta.gz
wget https://zenodo.org/record/7336407/files/PnB2.filtered.Trinity.fasta.gz
```
One haptophyte transcriptome was retrieved from this [Zenodo repository](https://zenodo.org/record/1476491):

```
wget https://zenodo.org/record/1476491/files/phaeocystisglobosa_euk_seqs.fasta
```

## Seeleuthner et al. SAGs (2018)

Predicted proteins from eight single-cell ampified genomes of Ochrophyte were downloaded from [genoscope](https://www.genoscope.cns.fr/tara/).

Download predicted proteins using wget:
```
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/m3a.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/m3f.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/m4a1.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/m4a2.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/m4c.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/m4e.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/chryh1.faa.gz
wget https://www.genoscope.cns.fr/tara/localdata/data/SAG-v1/chryh2.faa.gz
```
Unzip the compressed files

```
gunzip *.faa.gz
```

## NCBI Genbank downloads

Eighteen transcriptomes and 23 predicted protein gene models from genomes were gathered from [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/). FASTA files were downloaded manually from Genbank; the original links, accession numbers and raw sequence filenames used are listed in the [MarFERReT metadata file](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/data/MarFERReT.v1.metadata.csv).

## Roscoff Culture Collection

A total of 41 nucleotide transcriptomes were initially gathered from [METBdb](http://metdb.sb-roscoff.fr/metdb/).
Download nucleotide transcriptomes using wget:

```
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00017-ankylochrysis-sp-malift191.5pg1-paired/Trinity/RCC-METDB_00017-ankylochrysis-sp-malift191.5pg1-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00030-aureococcus-anophagefferens-ccmp1784-paired/Trinity/RCC-METDB_00030-aureococcus-anophagefferens-ccmp1784-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00044-bigelowiella-natans-bl_33-3-paired/Trinity/RCC-METDB_00044-bigelowiella-natans-bl_33-3-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00045-bigelowiella-sp-bl_34-paired/Trinity/RCC-METDB_00045-bigelowiella-sp-bl_34-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00069-chloroparvula-sp_b2-indianocean_1-1-paired/Trinity/RCC-METDB_00069-chloroparvula-sp_b2-indianocean_1-1-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00070-chloroparvula-sp_b3-biosope_46c4s-paired/Trinity/RCC-METDB_00070-chloroparvula-sp_b3-biosope_46c4s-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00074-chloropicon-maureeniae-lrf_5-paired/Trinity/RCC-METDB_00074-chloropicon-maureeniae-lrf_5-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00082-chloropicon-sp-ccmp2113-paired/Trinity/RCC-METDB_00082-chloropicon-sp-ccmp2113-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00086-chroodactylon-ramnosum-sag103.79-paired/Trinity/RCC-METDB_00086-chroodactylon-ramnosum-sag103.79-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00146-emiliania-huxleyi-ac675syn.os2-paired/Trinity/RCC-METDB_00146-emiliania-huxleyi-ac675syn.os2-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00150-emiliania-huxleyi-boum77-paired/Trinity/RCC-METDB_00150-emiliania-huxleyi-boum77-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00151-emiliania-huxleyi-chc2-paired/Trinity/RCC-METDB_00151-emiliania-huxleyi-chc2-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00180-gephyrocapsa-mullerae-chc184-paired/Trinity/RCC-METDB_00180-gephyrocapsa-mullerae-chc184-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00185-glossomastix-chrysoplasta-ccmp1537-paired/Trinity/RCC-METDB_00185-glossomastix-chrysoplasta-ccmp1537-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00191-guinardia-flaccida-ra120910p1e1-paired/Trinity/RCC-METDB_00191-guinardia-flaccida-ra120910p1e1-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00198-haslea-sp-aka_5Ap52.f2-paired/Trinity/RCC-METDB_00198-haslea-sp-aka_5ap52.f2-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00246-mesopedinella-arctica-prosope_2-paired/Trinity/RCC-METDB_00246-mesopedinella-arctica-prosope_2-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00250-micromonas-commoda-biosope_46b4s-paired/Trinity/RCC-METDB_00250-micromonas-commoda-biosope_46b4s-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00262-minidiscus-sp-ra080513-10-paired/Trinity/RCC-METDB_00262-minidiscus-sp-ra080513-10-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00269-navicula-sp-ra120910p2a1-paired/Trinity/RCC-METDB_00269-navicula-sp-ra120910p2a1-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00278-ochromonas-triangulata-chr25-paired/Trinity/RCC-METDB_00278-ochromonas-triangulata-chr25-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00299-partenskyella-glossopodia-prosope_99-paired/Trinity/RCC-METDB_00299-partenskyella-glossopodia-prosope_99-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00309-pelagomonas-calceolata-ccmp1214-paired/Trinity/RCC-METDB_00309-pelagomonas-calceolata-ccmp1214-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00310-pelagomonas-calceolata-boum81-4-paired/Trinity/RCC-METDB_00310-pelagomonas-calceolata-boum81-4-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00311-pelagomonas-calceolata-indianocean_17-1-paired/Trinity/RCC-METDB_00311-pelagomonas-calceolata-indianocean_17-1-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00312-pelagomonas-calceolata-biosope_111c1e-paired/Trinity/RCC-METDB_00312-pelagomonas-calceolata-biosope_111c1e-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00316-pelagophyceae_xxx-sp-akap71.e10-paired/Trinity/RCC-METDB_00316-pelagophyceae_xxx-sp-akap71.e10-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00325-phaeocystis-antarctica-ccmp1374-paired/Trinity/RCC-METDB_00325-phaeocystis-antarctica-ccmp1374-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00326-phaeocystis-antarctica-ccmp1871-paired/Trinity/RCC-METDB_00326-phaeocystis-antarctica-ccmp1871-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00327-phaeocystis-cordata-ccmp2495-paired/Trinity/RCC-METDB_00327-phaeocystis-cordata-ccmp2495-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00329-phaeocystis-globosa-rikz-graves-paired/Trinity/RCC-METDB_00329-phaeocystis-globosa-rikz-graves-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00330-phaeocystis-rex-he001206-d2-b1-paired/Trinity/RCC-METDB_00330-phaeocystis-rex-he001206-d2-b1-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00334-phaeocystis-sp-nies1396-paired/Trinity/RCC-METDB_00334-phaeocystis-sp-nies1396-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00335-phaeocystis-sp-biosope202-paired/Trinity/RCC-METDB_00335-phaeocystis-sp-biosope202-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00346-pleurosigma-sp-ra120910p1b3-paired/Trinity/RCC-METDB_00346-pleurosigma-sp-ra120910p1b3-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00408-shionodiscus-bioculatus-malinaft99.3pg3bis-paired/Trinity/RCC-METDB_00408-shionodiscus-bioculatus-malinaft99.3pg3bis-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00442-synedra-sp-malinas93p44.d11-paired/Trinity/RCC-METDB_00442-synedra-sp-malinas93p44.d11-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00475-triparma-pacifica-oli94sch-paired/Trinity/RCC-METDB_00475-triparma-pacifica-oli94sch-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00476-triparma-pacifica-min129-20maa-paired/Trinity/RCC-METDB_00476-triparma-pacifica-min129-20maa-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00485-phaeocystis-sp-biosope226fl2-3-paired/Trinity/RCC-METDB_00485-phaeocystis-sp-biosope226fl2-3-paired.fasta
wget http://metdb.sb-roscoff.fr/metdb/download/download/rcc-metdb_00336-phaeocystis-sp-biosope219fl2-3-paired/Trinity/RCC-METDB_00336-phaeocystis-sp-biosope219fl2-3-paired.fasta
```

## JGI PHYCOCOSM

A total of 117 predicted-protein gene models were gathered from [JGI Phycocosm](https://phycocosm.jgi.doe.gov/phycocosm/home).

Protein gene models were downloaded manually from Phycocosm as compressed FASTA files (aa.fasta.gz). 

## datadryad.org

Two single-cell amplified transcriptomes (SATs) were manually downloaded from a repository at datadryad.org [here](https://datadryad.org/stash/dataset/doi:10.5061/dryad.pg4f4qrk0) in a single data file `doi_10.5061_dryad.pg4f4qrk0__v9.zip` containing two FASTA files; `Abedinium_dasypus_DSEL18_54_rnaspades.fasta` and `Abedinium_folium_DICHO17_03_rnaspades.fasta`. 



