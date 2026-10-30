# The egg-egg and cell-cell vaccination arms response in PBMC scRNA-seq

Gene browser: [index.html](index.html)

## Software

- Cell Ranger 10.1.0, Apptainer container cellranger_10.1.0.sif
- /cluster/VAST/wanx-general-lab/Nina/Mandy_sc/Working/software/cellranger_10.1.0.sif
- /cluster/VAST/wanx-general-lab/Nina/Mandy_sc/Data/GR0389

## Run Cell Ranger

- Genome reference: ASM1176430v1.1 (NCBI 102), cellranger mkref, MusPutFur_ASM1176430v1_1_ncbi102, 09_mkref_ferret_ncbi.sbatch
- V(D)J reference: Walsh et al. 2025 + IMGT, /cluster/VAST/wanx-general-lab/Nina/Mandy_sc/Working/reference/MusPutFur_vdj_imgt/fasta/regions.fa
- Run [cellranger multi](https://www.10xgenomics.com/support/software/cell-ranger/8.0/analysis/running-pipelines/cr-5p-multi) on both pools, 10_multi_ocm_ncbi.sbatch
- Output per sample, Cell_ranger.zip
  - Count matrix: sample_filtered_feature_bc_matrix
  - B cell receptors: vdj_b/filtered_contig_annotations.csv, airr_rearrangement.tsv
  - T cell receptors: vdj_t/filtered_contig_annotations.csv

### Reference genome

- v1_ensembl: old, Ensembl 116. No IG constant genes, no MT genes
- v2_ncbi: current, NCBI 102. All 6 IGH constant genes, 13 MT genes

### Gene names

- Replace LOC gene symbols with human gene symbols, ferret_LOC_to_symbol_map

## Downstream analysis, Downstream_analysis.zip

- Files: [Analysis](Analysis)
- All the functions are stored in _setup.R
- Step 0: step_00_preprocessing.Rmd, step_00_preprocessing.html
- Step 1: step_10_annotation.Rmd, find the annotation manually first, step_10_annotation.html, step_15_subclustering.html
- Step 2: step_20_cellproportion.Rmd, step_20_cellproportion.html
- Step 3: step_30_bcr_tcr_detection.Rmd, step_30_bcr_tcr_detection.html
- Step 4: step_40_bcr_plasmablast.Rmd, step_40_bcr_plasmablast.html
- Step 5: step_50_T_proportion.Rmd, step_50_T_proportion.html
