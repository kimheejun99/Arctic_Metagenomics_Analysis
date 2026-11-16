# Arctic Metagenome Pipeline — Stage 2: Comparative Analysis

This repository covers Stage 2 of a two-stage Arctic metagenome analysis.
Stage 1 (raw data processing): https://github.com/kimheejun99/Arctic_Metagenomics_Pipeline

## Repository Structure
- 01_taxonomy/       : Taxonomic profiling (Kraken2)
- - 02_diversity/      : Alpha & Beta diversity analysis
  - - 03_ARG/            : ARG quantification (ARGs-OAP)
    - - 04_VF/             : Virulence factor analysis (VFDB)
      - - 05_methanogenesis/ : Methanogenesis KEGG module analysis
        - - 06_host_tracing/   : Read-level host tracing
          - - envs/              : Conda environment files
           
            - ## Data
            - - NCBI BioProject: PRJNA268009
              - - Samples: 49 peat soil metagenomes (AL: 36, PERM: 13)
               
                - ## Citation
                - Kim, H. (2026). Comparative Metagenomic Analysis of Active Layer and
                - Permafrost Peat Soils in Utqiagvik, Alaska.
                - Kangwon National University, MS Thesis.
