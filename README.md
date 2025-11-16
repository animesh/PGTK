# PGTK - ProteoGenomics of Myeloma Cell Lines

## Overview

PGTK (ProteoGenomics TK) is a comprehensive multi-omics study of seven in-house Myeloma patient-derived cell lines (TK) with a more primary cell-like phenotype. This repository contains information and resources related to the integrated proteogenomics analysis of these cell lines.

## Publication

The detailed methodology and findings of this study are published in:

**Article:** [Cancers - MDPI](https://www.mdpi.com/2072-6694/16/23/3963)

Please cite this publication when using data or results from this project.

## Data Availability

### Transcriptome Data

The RNA-seq transcriptome data for the seven Myeloma cell lines is publicly available from NCBI BioProject:

- **NCBI BioProject ID:** [PRJNA1176350](https://www.ncbi.nlm.nih.gov/bioproject/1176350)
- **Access:** [https://www.ncbi.nlm.nih.gov/bioproject/1176350](https://www.ncbi.nlm.nih.gov/bioproject/1176350)

#### trancriptome https://www.ncbi.nlm.nih.gov/bioproject/1176350 download using https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit 

```
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
tar xvzf sratoolkit.3.2.1-ubuntu64.tar.gz
SRATOOL_DIR="$PWD$/sratoolkit.3.2.1-ubuntu64"; mkdir -p ./fasterq_tmp; curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1176350&usehistory=y" | grep -oP '(?<=<QueryKey>)\d+|(?<=<WebEnv>)[^<]+' | { read qk; read we; curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$qk&WebEnv=$we&rettype=runinfo&retmode=text"; } | tail -n +2 | cut -d',' -f1 | xargs -P4 -I{} "$SRATOOL_DIR/bin/fasterq-dump" {} --split-files -e 4 -p -O . --temp ./fasterq_tmp
```


### Proteomics Data

The mass spectrometry proteomics data is available through the PRIDE Archive:

- **PRIDE Project 1:** [PXD033531](https://www.ebi.ac.uk/pride/archive/projects/PXD033531)
  - Access: [https://www.ebi.ac.uk/pride/archive/projects/PXD033531](https://www.ebi.ac.uk/pride/archive/projects/PXD033531)

- **PRIDE Project 2:** [PXD033510](https://www.ebi.ac.uk/pride/archive/projects/PXD033510/)
  - Access: [https://www.ebi.ac.uk/pride/archive/projects/PXD033510/](https://www.ebi.ac.uk/pride/archive/projects/PXD033510/)

  #### proteome https://www.ebi.ac.uk/pride/archive/projects/PXD033531 and https://www.ebi.ac.uk/pride/archive/projects/PXD033510/ dowload recursively from the ftp

```
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033531/"
wget -r "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/11/PXD033510/"
```


## Project Description

This study focuses on seven patient-derived Myeloma cell lines (TK) that exhibit a more primary cell-like phenotype, making them valuable models for studying Multiple Myeloma biology. The integrated proteogenomics approach combines:

- **Transcriptomics:** Gene expression profiling through RNA-seq
- **Proteomics:** Protein abundance measurements through mass spectrometry

This multi-omics integration enables comprehensive characterization of the molecular landscape of these Myeloma cell lines.

## Citation

If you use this data or findings in your research, please cite:

```
[Citation details from the MDPI Cancers publication]
DOI: Available at https://www.mdpi.com/2072-6694/16/23/3963
```

## License

See [LICENSE](LICENSE) file for details.

## Contact

me 🤓