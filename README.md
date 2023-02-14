# apoe-ad-age-atlas
This repository contains the code needed to reproduce all results and plots in Millet, Ledo et al. Code snippets are separated into files mirroring the organization of figures in the manuscript.

**Fig_1_and_Fig_S1.R**: code needed to generate Seurat structures for the 10wk/20wk/96wk atlas as well as microglial subclustering, plus code for all indicated figures (summary plots, DGE analysis, and preparation of code for CellRank). For convenience, the .rds files for those structures and their original CellRanger outputs are included as well. (Raw sequencing data is uploaded to GEO, accession number XXXXXXXXX.)

> *Accompanying code: CellRank.ipynb*: This provides the Python notebook used to run CellRank on the data produced by the above R script.

**Fig_2_and_Fig_S2A.R**: code needed to run SCENIC, CellPhoneDB, and Compass analysis.

> *Accompanying code: SCENIC.ipynb*: This provides the Python notebook used to run SCENIC on the data produced by the above R script.

**Fig_S2B-D.R**: code needed to analyze scUTRquant output. scUTRquant was run exactly per [suggested documentation](https://github.com/Mayrlab/scUTRquant) using the Snakemake pipeline, with bx_whitelist: "extdata/bxs/3M-february-2018.txt", tech: "10xv3", strand: "--fr-stranded", min_umis: 1000.


**Fig_3_and_Fig_S3.R**: code needed to generate Signac structures for the 60wk multiome sample as well as microglial subclustering, plus code for all indicated figures (summary plots, footprinting, LDA modeling, SCENIC+, CellPhoneDB, and Compass). For convenience, the .rds files for those structures and their original CellRanger outputs are included as well. (Raw sequencing data is uploaded to GEO, accession number XXXXXXXX.)

> *Accompanying code: SCENICplus.ipynb & MIRA.ipynb*: These provide the Python notebooks used to run SCENIC+ and MIRA on the data produced by the above R script.

**Fig_4.R**: code needed to anchor integrate the listed previously published human and mouse single-cell datasets with the data generated in this study and to produce the relevant plots in the manuscript. 
