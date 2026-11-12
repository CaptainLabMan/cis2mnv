import re
import gzip
from Bio.SeqUtils import seq3
import pandas as pd


def extract_protein_position(value):
    '''N109S -> 109
    p.Asn109Ser -> 109
    p.Gly1046_Arg1048del -> 1046
    p.Val110delinsArgLeu -> 110'''

    match = re.search(r'\d+', value)
    if match is None:
        raise ValueError(f'Protein position not found in {value!r}')
    return int(match.group())


def to_hgvsp(value):
    '''N109S -> p.Asn109Ser'''
    return f'p.{seq3(value[0])}{value[1:-1]}{seq3(value[-1])}'


def prepare_alphamissense(input_filepath, output_filepath):
    transitions = {'pathogenic': 'P', 'ambiguous': 'A', 'benign': 'B'}

    with gzip.open(input_filepath, 'rt') as input_file, gzip.open(output_filepath, 'wt') as output_file:
        for row in input_file:
            row = row.rstrip('\r\n')
            if row and not row.startswith('#'):
                splitted_row = row.split('\t')
                if row.startswith('uniprot_id'):
                    splitted_row.insert(1, 'protein_position')
                    output_file.write('\t'.join(splitted_row) + '\n')
                else:
                    uniprot_id, protein_variant, am_pathogenicity, am_class = splitted_row
                    protein_position = extract_protein_position(protein_variant)
                    prepared_values = [uniprot_id, protein_position, to_hgvsp(protein_variant), round(float(am_pathogenicity), 2), transitions[am_class]]
                    output_file.write('\t'.join(str(value) for value in prepared_values) + '\n')

    return output_filepath


def prepare_clinvar(input_filpath, output_filepath, assembly='GRCh38'):
    '''assemblies available: GRCh37 / GRCh38'''

    clinvar_df = pd.read_csv(input_filpath, sep='\t', usecols=['Name', 'ClinicalSignificance', 'Assembly'], dtype={'ClinicalSignificance': 'category', 'Assembly': 'category'})

    clinvar_df_filtered = (clinvar_df[clinvar_df['Assembly'] == assembly]
                       .replace('-', pd.NA)
                       .dropna(subset=['Name', 'ClinicalSignificance'], how='any'))

    clinvar_df_filtered[['HGVSc', 'HGVSp']] = clinvar_df_filtered['Name'].str.split(' ', n=1, expand=True).reindex(columns=[0, 1])
    clinvar_df_filtered[['transcript', 'HGVSc']] = clinvar_df_filtered['HGVSc'].str.split(':', n=1, expand=True).reindex(columns=[0, 1])
    clinvar_df_filtered[['transcript', 'gene_name']] = clinvar_df_filtered['transcript'].str.split('(', n=1, expand=True).reindex(columns=[0, 1])
    clinvar_df_filtered[['transcript', 'transcript_version']] = clinvar_df_filtered['transcript'].str.split('.', n=1, expand=True).reindex(columns=[0, 1])
    clinvar_df_filtered['HGVSp'] = clinvar_df_filtered['HGVSp'].str.strip('()')
    clinvar_df_filtered['gene_name'] = clinvar_df_filtered['gene_name'].str.strip(')')

    clinvar_df_short = (clinvar_df_filtered[['gene_name', 'transcript', 'transcript_version', 'HGVSc', 'HGVSp', 'ClinicalSignificance']]
                    .dropna(how='any', axis=0)
                    .query('HGVSp.str.startswith("p.") and ~HGVSp.str.contains("fs|Ter|=|[?]") and ~HGVSp.str.endswith("del")', engine='python'))
    
    clinvar_df_short = (clinvar_df_short.groupby(['transcript', 'HGVSp'], as_index=False)['ClinicalSignificance']
                             .agg(lambda values: '; '.join(dict.fromkeys(
                                item.strip()
                                for value in values.astype(str)
                                for item in value.split(';')
                                if item.strip()
                            ))))

    clinvar_df_short['protein_position'] = (clinvar_df_short['HGVSp']
                                            .apply(extract_protein_position)
                                            .astype('Int64'))

    transitions = {
    'Uncertain significance': 'U',
    'Conflicting classifications of pathogenicity': 'C',
    'Likely benign': 'LB',
    'Likely pathogenic': 'LP',
    'Benign': 'B',
    'Pathogenic': 'P',
    'not provided': 'NP',
    'other': 'O',
    'drug response': 'DR',
    'no classification for the single variant': 'NC',
    'risk factor': 'RF',
    'association': 'As',
    'Affects': 'Af',
    }

    def replace_clin_sig(string, replace_dict):
        for key, value in replace_dict.items():
            string = string.replace(key, value)
        return string

    clinvar_df_short['ClinicalSignificance'] = clinvar_df_short['ClinicalSignificance'].apply(replace_clin_sig, args=(transitions, ))

    clinvar_df_short = clinvar_df_short[['transcript', 'protein_position', 'HGVSp', 'ClinicalSignificance']]

    clinvar_df_short.to_csv(output_filepath, sep='\t', compression='gzip', index=False)

    return output_filepath

    
