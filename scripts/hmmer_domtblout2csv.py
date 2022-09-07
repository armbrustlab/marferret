#!/usr/bin/env python


# hmmer_tab2csv.py


"""

Input: HMMER3 'tabular' output files
#                                                                             --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                      accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#              ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
S23C1_TRINITY_DN1574590_c3_g2_i1_1 -          1-cysPrx_C           PF10417.8      9e-15   58.5   0.0   1.5e-14   57.8   0.0   1.4   1   0   0   1   1   1   1 len167
S06C1_TRINITY_DN2124915_c0_g3_i2_1 -          1-cysPrx_C           PF10417.8    4.7e-13   53.0   0.1   9.4e-13   52.0   0.1   1.5   1   0   0   1   1   1   1 len226
S20C1_TRINITY_DN1612005_c1_g3_i4_1 -          1-cysPrx_C           PF10417.8    1.1e-12   51.8   0.0     2e-12   51.0   0.0   1.5   1   0   0   1   1   1   1 len167

Input: HMMER3 'domtblout' file:
#                                                                                           --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name                       accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#               ------------------- ---------- ----- -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
S08C1_TRINITY_DN1312418_c0_g1_i2_1  -            460 K00001               -            306   1.3e-70  249.5   0.2   1   1   7.6e-74   1.6e-70  249.1   0.2     2   305    51   418    46   419 0.90 len413
S17C1_TRINITY_DN2335443_c3_g1_i4_2  -            495 K00001               -            306   1.5e-70  249.2   0.2   1   1   8.8e-74   1.9e-70  248.9   0.2     5   306    69   436    53   436 0.89 len437


Output: A CSV file
	# Extract the wanted data from each line:
	for line in pfam_data:
		line_elts = line.split()
		contig_id = line_elts[0]
		pfam_name = line_elts[2]
		pfam_id = line_elts[3]
		pfam_eval = line_elts[4]
		pfam_score = line_elts[5]

	outline = contig_id + "," + pfam_name + "," + pfam_id + "," + pfam_score + "," + pfam_eval + "\n"


"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--tab", help="Input hmm_out.tab", type=str)
parser.add_argument("-o", "--out_file", help="Specify the name of the outfile", type=str)

args = parser.parse_args()

# open files:
pfam_data = open(args.tab, 'r')
outfile = open(args.out_file,'w')


def parse_line (line):

	line_elts = line.split()
	contig_id = line_elts[0]
	knum = line_elts[3]
	knum_eval = line_elts[6]
	knum_score = line_elts[7]

	outline = contig_id + "," + knum + "," + knum_eval + "," + knum_score + "\n"
	outfile.write(outline)

# Extract the wanted data from each line:
# header: "contig_id, pfam_name, pfam_id, pfam_score, pfam_eval"

# write out header to outfile:
header="aa_id,knum,knum_eval,knum_score\n"
outfile.write(header)

for line in pfam_data:
	if line.startswith("#") == False:
		parse_line(line)
