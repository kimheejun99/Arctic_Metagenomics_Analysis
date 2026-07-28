# Arctic Metagenome Analysis (Stage 2)

This is the second-stage analysis pipeline that aggregates, compares, and visualizes the per-sample outputs (Kraken2 report, VFDB annotation, ARGs-OAP results, etc.) produced by the primary circum-Arctic metagenome pipeline ([arctic-metagenome-pipeline](https://github.com/kimheejun99/Arctic_Metagenomics_Pipeline)) **across multiple samples**.

## Pipeline overview

```
sample_coords.txt (+ sample_list.txt) + metadata.txt
   │
   ├─ Step0_make_googlemap.py    Coordinate-based Google Maps visualization (only when a coordinates file exists)
   ├─ Step1_organize_samples.sh  Organizes flat outputs into folders by type
   ├─ Step2_vf_arg_plots.py      Overall VF/ARG averages + group comparisons (bar graphs, heatmaps)
   ├─ Step3_host_tracing.py      Traces the host phylum of reads with detected VF/ARG (joined with Kraken2)
   └─ Step4_merge_kreports.py    Merges kreports -> ASV/taxonomy/metadata tables
```

## Repository structure

```
.
├── environment.yml           # conda environment (ArcticAnalysis)
├── .gitignore
├── run_pipeline.sh           # Orchestrates the full pipeline
├── Step0_make_googlemap.py
├── Step1_organize_samples.sh
├── Step2_vf_arg_plots.py
├── Step3_host_tracing.py
├── Step4_merge_kreports.py
└── sample_coords.txt.example
```

## Requirements

Unlike the primary pipeline, no bioinformatics-specific binaries such as DIAMOND or Kraken2 are required. It runs purely on bash and Python (pandas, matplotlib).

```bash
conda env create -f environment.yml
conda activate ArcticAnalysis
```

| Tool | Purpose |
|---|---|
| Python 3.10 | Steps 0, 2, 3, 4 |
| pandas | tsv/kreport parsing, group aggregation |
| matplotlib | Step2 bar graphs/heatmaps, Step3 stacked bar |
| bash | run_pipeline.sh, Step1 |

## Input files

| File | Required | Description |
|---|---|---|
| `sample_coords.txt` | Required (tab-separated, header required: `Run latitude longitude`) | Used for Step0's Google Map. The Run column also serves as the sample list to process. |
| `sample_list.txt` | Optional | Auto-generated from the Run column of `sample_coords.txt` if absent. Only prepare this manually if you want to run a specific subset. |
| `metadata.txt` | Optional (tab-separated: `sample-id group`) | If absent, all samples are treated as a single "All" group. |

See `sample_coords.txt.example` for a template (rename it to `sample_coords.txt` and fill in your actual data when using it).

### Per-sample source files (outputs of the primary pipeline, placed flat in the current folder)

Before running Step1, files must be scattered in the current directory as follows:

```
./ (current location)
├── sample_coords.txt
├── (sample_list.txt)
├── metadata.txt
│
├── SRR0000001_1.fastq.gz
├── SRR0000001_2.fastq.gz
├── SRR0000001.kreport
├── SRR0000001.kraken.out
├── SRR0000001_annotated_parsed.tsv   <- requires a Category column
├── SRR0000001_args_output/
│   ├── normalized_cell.type.txt
│   ├── normalized_cell.subtype.txt
│   ├── normalized_cell.gene.txt
│   └── blastout.filtered.txt
│
└── (other samples follow the same pattern)
```

After running Step1, files are organized into folders by type (`./kraken/`, `./vfdb/`, `./arg/`); fastq.gz files are left in place and not moved. Steps 2-4 reference these organized folders.

## External database

Steps 3 and 4 require the NCBI taxdump (`nodes.dmp`, `names.dmp`) to convert taxid to phylum/lineage.

- This is likely already included in your existing Kraken2 DB (e.g., PlusPFP-16) as `taxonomy/nodes.dmp` and `taxonomy/names.dmp`.
- - If not, download `taxdump.tar.gz` from https://ftp.ncbi.nlm.nih.gov/pub/taxdump/ and extract it.
 
  - ## Usage
 
  - ```bash
    conda activate ArcticAnalysis
    bash run_pipeline.sh sample_coords.txt sample_list.txt metadata.txt
    # or without sample_list.txt:
    bash run_pipeline.sh sample_coords.txt "" metadata.txt
    ```

    Values that can be overridden via environment variables:

    ```bash
    GOOGLE_MAPS_API_KEY=your_key \
    TAXDUMP=/mnt/program/db/plusPFP_16/taxonomy \
    JOBS=8 \
    bash run_pipeline.sh sample_coords.txt sample_list.txt metadata.txt
    ```

    | Environment variable | Default | Description |
    |---|---|---|
    | `GOOGLE_MAPS_API_KEY` | `YOUR_GOOGLE_MAPS_API_KEY` (placeholder) | A valid key is required for the Step0 map to actually load |
    | `TAXDUMP` | `/mnt/program/db/plusPFP_16` | Parent path to nodes.dmp/names.dmp for Step3/4 |
    | `JOBS` | `6` | Number of parallel processes for Step3 |

    ## Outputs

    | Step | Output |
    |---|---|
    | Step0 | `sample_map.html` |
    | Step1 | `./kraken/`, `./vfdb/`, `./arg/` |
    | Step2 | `./plots/VF_total_category_barplot.png`, `VF_group_heatmap.png`, `ARG_total_type_barplot_*.png`, `ARG_group_heatmap_*.png`, etc. |
    | Step3 | `./host_tracing/VF_host_tracing_by_group.png`, `ARG_host_tracing_by_group.png`, raw tsv |
    | Step4 | `./merged_tables/1_ASV_table.txt`, `2_taxonomy_table.txt`, `3_metadata.txt` |

    ## Notes

    - `Step0_make_googlemap.py` uses the Google Maps JavaScript API, so you need to obtain a Maps JavaScript API key from the Google Cloud Console to actually view the map. Running it without a key will still create the HTML file, but the map will not load.
    - - All input files use a `.txt` extension for consistency. In practice they are tab-separated text, and the code parses them directly using `sep="\t"` (or automatic delimiter detection, `sep=None`, for Step0) rather than relying on the extension, so the extension itself has no effect on functionality.
      - 
