2023-11-20  Ryan Groussman, PhD <rgrous83@uw.edu>

## MarFERReT protein sequence reference library changelog

### v1.1.1

**Improvements to MarFERReT assembly scripts**
- Incorporated inclusion/exclusion flag into build procedure to integrate QC process
- Added check to see if translated protein sequence already exists
- Changed scripts to pull database version from input metadata file
- Streamlined the internal file naming process
  - Important output files are now generated with the same version ID as the input metadata.csv
- The /data/source_seqs/ directory for input source sequences is now pre-generated

**Added test data set and demo tutorial**
- Added a small data set to test the MarFERRet build process under /versions/v1_demo/
  - Seven-entry metadata.csv file (MarFERReT.v1_demo.metadata.csv)
  - Compressed tarball with small set of FASTA sequences (demo_mini_entries.tar.gz)
  - Demo tutorial script (marferret.v1_demo.tutorial.sh)



