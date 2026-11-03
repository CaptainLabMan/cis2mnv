# cis2mnv

mamba env create -f environment.yaml  

for Linux:  
    mamba env create -f environment.vep.yaml  
for MacOS:  
    CONDA_SUBDIR=osx-64 mamba env create -f environment.vep.yaml  

mamba env create -f environment.snpeff.yaml   

mamba run -n cis2mnv-vep.venv vep  

mamba run -n cis2mnv-snpeff.venv snpeff  
mamba run -n cis2mnv-snpeff.venv snpeff download -v hg38  

pyranges1  
pandas  
freebayes
samtools
tabix
bcftools
ensembl-vep
genozip

https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files

https://console.cloud.google.com/storage/browser/dm_alphamissense;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false
https://huggingface.co/spaces/ntranoslab/esm_variants/tree/main

https://pcingola.github.io/SnpEff/snpeff/human_genomes/  