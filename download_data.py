from pathlib import Path
from utils import fetch_file

# ==============================================================================================================
# Settings
# ==============================================================================================================
main_path = Path.cwd()
refs_dir = main_path / 'data/refs'
dbs_dir = main_path / 'data/dbs'

ensembl_version = 116

# ==============================================================================================================
# Download files
# ==============================================================================================================
refs_urls = {
    f'Homo_sapiens.GRCh38.{ensembl_version}.chr.gtf.gz': f'https://ftp.ensembl.org/pub/release-{ensembl_version}/gtf/homo_sapiens/Homo_sapiens.GRCh38.{ensembl_version}.chr.gtf.gz',
}

dbs_urls = {
    'variant_summary.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz',
    f'Homo_sapiens.GRCh38.{ensembl_version}.refseq.tsv.gz': f'https://ftp.ensembl.org/pub/release-{ensembl_version}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{ensembl_version}.refseq.tsv.gz',
    f'Homo_sapiens.GRCh38.{ensembl_version}.uniprot.tsv.gz': f'https://ftp.ensembl.org/pub/release-{ensembl_version}/tsv/homo_sapiens/Homo_sapiens.GRCh38.{ensembl_version}.uniprot.tsv.gz',
    'MANE.GRCh38.v1.5.summary.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz',
    'HUMAN_9606_idmapping.dat.gz': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz',
    'AlphaMissense_aa_substitutions.tsv.gz': 'https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz'
}

for url in refs_urls:
    fetch_file(refs_urls[url], refs_dir)

for url in dbs_urls:
    fetch_file(dbs_urls[url], dbs_dir)

# ==============================================================================================================
# Prepare data
# ==============================================================================================================

# --------------------------------------------------------------------------------------------------------------
# Prepare genes_x_phenotypes_matrix
# --------------------------------------------------------------------------------------------------------------