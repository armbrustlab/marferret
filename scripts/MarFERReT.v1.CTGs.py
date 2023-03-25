import pandas as pd
from pathlib import Path
from plotnine import ggplot, aes, geom_histogram, geom_vline, ggtitle, xlab, ylab, theme_bw

def calculate_lineage(df):
    conditions = [
        df['class'] == '2864',
        df['phylum'] == '2836',
        df['phylum'] == '5878',
        df['phylum'] == '2830',
        df['phylum'] == '3041',
        df['superkingdom__'] == '543769',
        df['superkingdom__'] == '33630',
        df['superkingdom___'] == '2696291',
        df['superkingdom_'] == '554915',
        df['superkingdom_'] == '33154',
        df['superkingdom'] == '2759'
    ]

    choices = [
        'Dinophyceae',
        'Bacillariophyta',
        'Ciliophora',
        'Haptophyta',
        'Chlorophyta',
        'Rhizaria',
        'Alveolata',
        'Ochrophyta',
        'Amoebozoa',
        'Opisthokonta',
        'Eukaryota'
    ]

    df['lineage'] = pd.select(conditions, choices, default='Other')
    return df

def process_lineage(lineage, mft_dat, pfam_presence_by_species_tsa):
    print(lineage)

    lineage_species = mft_dat[mft_dat['lineage'] == lineage]['species'].unique()
    pfam_sub = pfam_presence_by_species_tsa[lineage_species]

    pfam_sub['rowsum'] = pfam_sub.sum(axis=1)
    pfam_sub['frequency'] = pfam_sub['rowsum'] / (pfam_sub.shape[1] - 1)

    core_pfams_subset = pfam_sub[pfam_sub['frequency'] >= 0.95].index

    pfam_sub['pfam_id'] = pfam_sub.index
    pfam_sub['lineage'] = lineage
    ctg_dat_sub = pfam_sub[['lineage', 'rowsum', 'pfam_id', 'frequency']].rename(columns={'rowsum': 'n_species'})

    out_csv = mft_dir / f"data/ctg/MarFERReT.pfam_presence.{lineage}.ctg_catalog.v1.csv"
    pfam_sub.to_csv(out_csv, index=False)

    return ctg_dat_sub

mft_dir = Path("/Users/rgroussman/Dropbox/Armbrust/EukRefDB/SciData_submission/")
data_dir = mft_dir / "data"
pfam_dir = data_dir / "pfam"

mft_dat = pd.read_csv(data_dir / "MarFERReT.v1.metadata.csv")
TAXA_CSV = pd.read_csv(data_dir / "MarFERReT.v1.taxa.csv")
pfam_dat = pd.read_csv(pfam_dir / "MarFERReT.v1.entry_pfam_sums.csv")

unique_pfams = pfam_dat['pfam_name'].nunique()

pfam_dat2 = (
    pfam_dat
    .merge(mft_dat[['ref_id', 'marferret_name', 'tax_id', 'data_type']], on='ref_id', how='left')
    .query('data_type in ["TSA", "SAT"]')
    .merge(TAXA_CSV[['tax_id', 'species']], on='tax_id', how='left')
    .drop(['ref_id', 'marferret_name', 'tax_id', 'data_type'], axis=1)
    .pivot_table(index='pfam_name', columns='species', values='entry_count', fill_value=0)
)

mft_dat = calculate_lineage(mft_dat)

lineages = [
    'Dinophyceae',
    'Bacillariophyta',
    'Ciliophora',
    'Haptophyta',
    'Chlorophyta',
    'Rhizaria',
    'Alveolata',
    'Ochrophyta',
    'Amoebozoa',
    'Opisthokonta',
    'Eukaryota',
    'Other'
]

ctg_dat = pd.concat([process_lineage(lineage, mft_dat, pfam_dat2) for lineage in lineages])
ctg_dat.to_csv(mft_dir / "data/ctg/MarFERReT.pfam_presence.ctg_catalog.v1.csv", index=False)

plot_data = ctg_dat.query('lineage != "Eukaryota" and lineage != "Other"')

(ggplot(plot_data, aes(x='n_species'))
 + geom_histogram(binwidth=1)
 + geom_vline(aes(xintercept=unique_pfams * 0.95), linetype='dashed', color='red')
 + ggtitle("Core PFAMs by lineage")
 + xlab("Number of species with PFAM")
 + ylab("Frequency")
 + theme_bw()
).save(mft_dir / "figures/ctg/MarFERReT.pfam_presence.ctg_catalog.v1.pdf", dpi=300)
