import subprocess
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from urllib.parse import urlsplit
from urllib.request import urlretrieve
import pandas as pd
import pyranges1 as pr
import pysam


def fetch_file(link: str, output_dir: str, filename: str = None) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    filename = filename or PurePosixPath(urlsplit(link).path).name
    if not filename:
        raise ValueError("Filename doesn't specified")

    output_path = output_dir / filename

    if output_path.is_file():
        return output_path

    with TemporaryDirectory(dir=output_dir) as temp_dir:
        temp_path = Path(temp_dir) / 'download'
        urlretrieve(link, temp_path)
        temp_path.replace(output_path)

    return output_path


def index_bam(bam_filepath):
    bam_filepath = Path(bam_filepath)

    index_paths = [
        Path(str(bam_filepath) + ".bai"),
        bam_filepath.with_suffix(".bai"),
        Path(str(bam_filepath) + ".csi"),
    ]

    if not any(path.is_file() for path in index_paths):
        subprocess.run(
            ["samtools", "index", str(bam_filepath)],
            check=True,
        )


def calculate_aa_pos(chrom: str, pos: int, cds_table_df: pd.DataFrame):
    pos_zero_based = int(pos) - 1

    # Находим интервалы, содержащие позицию.
    matched = cds_table_df.loc[
        cds_table_df['Chromosome'].eq(chrom)
        & cds_table_df['Start'].le(pos_zero_based)
        & cds_table_df['End'].gt(pos_zero_based)
    ]

    # Рассчитываем смещение с учётом цепи.
    # .where(условие, другое_значение) означает: если условие истинно — оставить исходное значение, иначе взять другое.
    offset = (pos_zero_based - matched['Start']).where(
        matched['Strand'].eq('+'),
        matched['End'] - 1 - pos_zero_based,
    )

    # Получаем позиции CDS и белка.
    matched = matched.assign(
        variant_cds_pos=matched['cds_start'] + offset,
        variant_protein_pos=lambda df: (df['variant_cds_pos'] - 1) // 3 + 1,
    )

    results = list(
        matched[
            ['transcript_id', 'variant_protein_pos']
        ].itertuples(index=False, name=None)
    )

    return results


def create_mnv_candidate_regions(
    vcf_filepath,
    bed_filepath=None,
    max_region_length=3,
    merge_overlapped_regions=0,
    only_pass_variants=1,
    cds_table_filepath=None
):
    cds_table_df = (pd.read_parquet(cds_table_filepath) if cds_table_filepath is not None else None)
    if cds_table_df is not None:
        cds_table_df = cds_table_df.astype({
            'Chromosome': str,
            'Start': 'int64',
            'End': 'int64',
            'cds_start': 'int64',
        })
    
    variants = []
    with pysam.VariantFile(vcf_filepath) as vcf:
        if bed_filepath is None:
            variant_records = vcf
        else:
            intervals = pd.read_csv(
                bed_filepath,
                sep="\t",
                header=None,
                usecols=range(3),
                names=["Chromosome", "Start", "End"],
                dtype={"Chromosome": str, "Start": int, "End": int},
            )

            variant_records = (
                variant
                for chrom, start, end in intervals.itertuples(
                    index=False, name=None
                )
                for variant in vcf.fetch(chrom, start, end)
            )

        for variant in variant_records:
            if only_pass_variants and "PASS" not in variant.filter:
                continue

            if variant.stop - variant.start > max_region_length:
                continue

            variants.append((variant.chrom, variant.start, variant.stop))

    df = (
        pd.DataFrame(variants, columns=["Chromosome", "Start", "End"])
        .drop_duplicates()
        .sort_values(["Chromosome", "Start"])
    )

    regions = []

    for chrom, group in df.groupby("Chromosome", sort=False):
        records = list(
            group[["Start", "End"]].itertuples(index=False, name=None)
        )
        last_end = -1

        for i, (start, end) in enumerate(records):
            count = 1
            limit = start + max_region_length

            for j in range(i + 1, len(records)):
                next_start, next_end = records[j]

                if next_start == start:
                    continue

                if next_start >= limit:
                    break

                if next_end > limit:
                    continue

                # Проверяем совпадение кодона, если передана таблица CDS.
                if cds_table_df is not None:
                    variant_codons = set(
                        calculate_aa_pos(chrom, start + 1, cds_table_df)
                    )
                    next_variant_codons = set(
                        calculate_aa_pos(chrom, next_start + 1, cds_table_df)
                    )

                    if not (variant_codons & next_variant_codons):
                        continue

                end = max(end, next_end)
                count += 1

            # Keep overlapping windows, but skip contained ones.
            if count >= 2 and end > last_end:
                regions.append((chrom, start, end))
                last_end = end

    regions_pr = pr.PyRanges(pd.DataFrame(
        regions,
        columns=["Chromosome", "Start", "End"],
    ))

    if merge_overlapped_regions:
        regions_pr = (regions_pr
                      .merge_overlaps()
                      .sort_ranges())

    return regions_pr