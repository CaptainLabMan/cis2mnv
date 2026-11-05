import subprocess
import gzip
from contextlib import nullcontext
import pysam
import pandas as pd
import pyranges1 as pr


def filter_gtf_by_features(input_gtf_filepath: str,
                           features: list | set,
                           add_chr: bool = False,
                           remove_chr: bool = False,
                           output_gtf_filepath: str | None = None) -> str | pd.DataFrame:
    
    if add_chr and remove_chr:
        raise ValueError("add_chr and remove_chr cannot both be True")

    def chr_transition(splitted_row):
        if add_chr:
            if not splitted_row[0].startswith('chr'):
                splitted_row[0] = 'chr' + splitted_row[0]
            if splitted_row[0] == 'chrMT':
                splitted_row[0] = 'chrM'
        elif remove_chr:
            splitted_row[0] = splitted_row[0].removeprefix('chr')
            if splitted_row[0] == 'M':
                splitted_row[0] = 'MT'
        return splitted_row

    features = set(features)

    with gzip.open(input_gtf_filepath, 'rt', encoding="utf-8") as input_gtf, (pysam.BGZFile(output_gtf_filepath, 'w') if output_gtf_filepath else nullcontext()) as output_gtf:
        selected_rows = []
        for row in input_gtf:
            if row.startswith('#'):
                if output_gtf_filepath:
                    output_gtf.write(row.encode('utf-8'))
            else:
                splitted_row = row.rstrip('\r\n').split('\t')
                feature = splitted_row[2]
                if feature in features:
                    if add_chr or remove_chr:
                        splitted_row = chr_transition(splitted_row)
                    if output_gtf_filepath:
                        output_gtf.write(('\t'.join(splitted_row) + '\n').encode("utf-8"))
                    else:
                        selected_rows.append(splitted_row)

    if output_gtf_filepath:
        return output_gtf_filepath
    return pd.DataFrame(
        selected_rows,
        columns=[
            'Chromosome',
            'Source',
            'Feature',
            'Start',
            'End',
            'Score',
            'Strand',
            'Frame',
            'Attributes',
        ],
    ).astype({'Start': 'int64', 'End': 'int64'})


def keep_only_main_gtf_transcripts(gtf_input: str | pr.PyRanges) -> pr.PyRanges:
    # ==============================================================================================================
    # Оставляем только 1 транскрипт на каждый ген (gene_id).
    # Исключением являются гены для которых есть MANE_Plus_Clinical, у них может быть два транскрипта.
    # Возвращает pr содержащий строки только с транскриптами.
    # ==============================================================================================================

    # --------------------------------------------------------------------------------------------------------------
    # Подгружаем pr содержащий только транскрипты.
    # (Можно загрузить весь и фильтрануть: Feature == "transcript")
    # --------------------------------------------------------------------------------------------------------------
    if isinstance(gtf_input, str):
        gtf_pr_transcripts = pr.read_gtf(gtf_input, duplicate_attr=True)
    elif isinstance(gtf_input, pr.PyRanges):
        gtf_pr_transcripts = gtf_input
    else:
        raise TypeError('gtf_input must be a file path string or a PyRanges object')
    
    gtf_pr_transcripts = gtf_pr_transcripts[gtf_pr_transcripts['Feature'] == 'transcript'].copy()

    # --------------------------------------------------------------------------------------------------------------
    # Превращение тэгов в список и разбиение на разные строки (взрыв/explode).
    # --------------------------------------------------------------------------------------------------------------
    gtf_pr_transcripts['tag'] = gtf_pr_transcripts['tag'].str.split(',')
    gtf_pr_transcripts_exploded = gtf_pr_transcripts.explode('tag')

    # --------------------------------------------------------------------------------------------------------------
    # Задаем и присваиваем приоритет каждому транскрипту по тэгу и биотипу транскриптов для дальнейшей сортировки с их участием
    # --------------------------------------------------------------------------------------------------------------
    tag_priority = {
        'MANE_Select': 0,
        'Ensembl_canonical': 1,
        'RefSeq Select': 1,
        'gencode_primary': 2,
        'gencode_basic': 3,
        'CCDS': 4,
    }

    transcript_biotype_priority = {
        'protein_coding': 0,
    }

    gtf_pr_transcripts_exploded['tag_priority'] = gtf_pr_transcripts_exploded['tag'].map(tag_priority).astype('Int64')
    gtf_pr_transcripts_exploded['transcript_biotype_priority'] = gtf_pr_transcripts_exploded['transcript_biotype'].map(transcript_biotype_priority).astype('Int64')

    # --------------------------------------------------------------------------------------------------------------
    # Сортируем
    # --------------------------------------------------------------------------------------------------------------
    gtf_pr_transcripts_exploded = gtf_pr_transcripts_exploded.sort_values(['gene_id', 'tag_priority', 'transcript_support_level', 'transcript_biotype_priority'], na_position='last')

    # --------------------------------------------------------------------------------------------------------------
    # Оставляем только первый транскрипт для каждого гена (gene_id)
    # --------------------------------------------------------------------------------------------------------------
    gtf_pr_transcripts_singe_transcript = gtf_pr_transcripts_exploded.drop_duplicates(subset=['gene_id'], keep='first')

    # --------------------------------------------------------------------------------------------------------------
    # Собираем MANE_Plus_Clinical транскрипты
    # --------------------------------------------------------------------------------------------------------------
    gtf_pr_transcripts_mane_plus_clinical = gtf_pr_transcripts_exploded.query('tag == "MANE_Plus_Clinical"')

    # --------------------------------------------------------------------------------------------------------------
    # Объединяем датафреймы с основными транскриптами и транскриптами MANE_Plus_Clinical.
    # Убираем дубликаты на тот случай если MANE_Plus_Clinical уже является основным.
    # --------------------------------------------------------------------------------------------------------------
    gtf_pr_main_transcripts = pd.concat([gtf_pr_transcripts_singe_transcript, gtf_pr_transcripts_mane_plus_clinical], ignore_index=True)
    gtf_pr_main_transcripts = gtf_pr_main_transcripts.drop_duplicates(subset=['gene_id', 'transcript_id'])
    gtf_pr_main_transcripts = gtf_pr_main_transcripts.drop(['tag_priority', 'transcript_biotype_priority'], axis=1)

    return gtf_pr_main_transcripts


def filter_gtf_by_transcripts(input_gtf_filepath, output_gtf_filepath, transcripts_set) -> str:
    '''Keeps every gene (including exon, cds etc.) only by 1 transcript (+ MANE Clinical Plus if available)'''
    transcripts_set = set(transcripts_set)
    with gzip.open(input_gtf_filepath, 'rt', encoding="utf-8") as input_gtf, pysam.BGZFile(output_gtf_filepath, 'w') as output_gtf:
        for row in input_gtf:
            if row.startswith('#'):
                output_gtf.write(row.encode("utf-8"))
            else:
                splitted_row = row.split('\t')
                feature = splitted_row[2]
                if feature == 'gene':
                    output_gtf.write(row.encode('utf-8'))
                else:
                    attributes = splitted_row[8].strip().rstrip(';').split('; ')
                    attribute_dict = {}
                    for attribute in attributes:
                        key, value = attribute.split(' ', 1)
                        attribute_dict[key] = value.strip('"')
                    if attribute_dict.get('transcript_id') in transcripts_set:
                        output_gtf.write(row.encode('utf-8'))
    return output_gtf_filepath


def sort_gtf(gtf_input_filepath, output_gtf_filepath):
    command = f'''gunzip -c "{gtf_input_filepath}" \
        | grep -v '^#' \
        | LC_ALL=C sort -t $'\t' -k1,1 -k4,4n -k5,5n \
        | bgzip -c > "{output_gtf_filepath}"'''
    subprocess.run(command, shell=True, check=True)


def gz_to_bgz(input_filepath, output_filepath):
    command = f'gunzip -c "{input_filepath}" | bgzip > "{output_filepath}"'
    subprocess.run(command, shell=True, check=True)


def index_gtf(gtf_input_filepath):
    command = f'tabix -p gff {gtf_input_filepath}'
    subprocess.run(command, shell=True, check=True)


def calculate_cds_pos(gtf_table: pr.PyRanges) -> pr.PyRanges:
    '''Calculates CDS (+ stop_codon) positions by chromosome positions.
    Applicable only for CDS and stop_codon features. GTF should be filetered before usage.'''
    strand_map = {'+': 1, '-': -1}

    result = (
        gtf_table.copy()
        .assign(
            length=lambda df: df['End'] - df['Start'],
            weighed_start=lambda df: (
                df['Start'] * df['Strand'].map(strand_map).astype('int64')
            ),
        )
        .sort_values(['transcript_id', 'weighed_start', 'End'])
        .assign(
            cds_end=lambda df: df.groupby('transcript_id')['length'].cumsum(),
            cds_start=lambda df: df['cds_end'] - df['length'] + 1,
        )
    )

    return pr.PyRanges(result)