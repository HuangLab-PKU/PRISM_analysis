# PRISM Analysis

This repository contains the analysis code and data for the PRISM (Probe-based RNA Imaging for Spatial Transcriptomics) paper. PRISM is a probe-based spatial transcriptomics technology that enables high-resolution, multi-gene simultaneous detection.

## Project Overview

PRISM technology achieves simultaneous detection of multiple gene expression patterns in tissue sections through the design of specific probes combined with fluorescence in situ hybridization (FISH) and imaging technology. This analysis covers various application scenarios from 2D to 3D, including mouse brain, mouse embryo, and human hepatocellular carcinoma samples.

## Folder Structure and Corresponding Paper Figures

### Methodology/ - Methodological Validation (Fig. S14)
Contains comparative analysis of PRISM technology with other spatial transcriptomics technologies:
- **gene_pattern_comparison/**: Gene expression pattern comparison with MERFISH, MERScope, BGI Stereo-seq and other technologies
- Validates the accuracy and reliability of PRISM technology

### MouseBrain2D/ - Mouse Brain 2D Analysis (Fig. 1, S16)
2D spatial transcriptomics analysis of mouse brain sections:
- **PRISM_cell_typing_and_analysis.ipynb**: Cell type identification and spatial analysis
- **N_cluster_scanpy.py**: Clustering number optimization analysis
- Demonstrates the application of PRISM in 2D tissue sections

### MouseEmbryo/ - Mouse Embryo Analysis (Fig. 2, 3)
Developmental analysis of mouse embryos at E12.5-E13.5 stages:
- **Cell_typing_Embryo.ipynb**: Embryonic cell type identification and spatial distribution analysis
- Studies gene expression patterns during embryonic development

### HCC2D/ - Human Hepatocellular Carcinoma 2D Analysis (Fig. 4, S30, S31, S33, S34)
2D spatial transcriptomics analysis of human hepatocellular carcinoma samples:
- **cell_typing.ipynb**: Cancer cell and immune cell type identification
- **interaction_graph.ipynb**: Cell-cell interaction network analysis
- **nhood_enrichment.ipynb**: Neighborhood enrichment analysis
- **GASTON_demo.ipynb**: GASTON algorithm demonstration
- **STAGATE_region_define.ipynb**: Spatial region definition
- Reveals cellular composition and spatial organization of tumor microenvironment

### HCC3Dpseudo/ - Human Hepatocellular Carcinoma Pseudo-3D Analysis (Fig. 5, S35, S36, S37, S38)
3D spatial analysis reconstructed from 2D sections:
- **cell_typing.ipynb**: 3D cell type identification
- **STAGATE_analysis.ipynb**: STAGATE spatial analysis
- **CAF_sep.ipynb**: Cancer-associated fibroblast separation analysis
- **interaction_graph.ipynb**: 3D cell interaction analysis
- **roi_variation.ipynb**: Region of interest variation analysis
- **cell_projection.ipynb**: Cell projection analysis
- **nhood_enrichment.ipynb**: Neighborhood enrichment analysis
- Constructs 3D spatial structure model of tumors

### MouseBrain3D/ - Mouse Brain 3D Analysis (Fig. 6, S41, S42, S43)
3D spatial transcriptomics analysis of different brain regions in mouse brain:
- **CTX_cell_typing_and_analysis.ipynb**: Cortical region analysis
- **HP_cell_typing_and_analysis.ipynb**: Hippocampal region analysis
- **HY_cell_typing_and_analysis.ipynb**: Hypothalamic region analysis
- **TH_cell_typing_and_analysis.ipynb**: Thalamic region analysis
- **subcellular_analysis.ipynb**: Subcellular localization analysis
- **rna_distribution_cal.ipynb**: RNA distribution calculation
- Reveals 3D gene expression patterns in different brain regions

## Technical Features

- **High Resolution**: Single-cell level resolution
- **Multi-gene Detection**: Simultaneous detection of 30 genes
- **Spatial Information Preservation**: Complete preservation of cell spatial position information
- **3D Reconstruction**: Support for 3D spatial structure reconstruction from 2D sections
- **Algorithm Integration**: Integration of multiple spatial analysis algorithms (STAGATE, GASTON, DECIPHER, etc.)

## Environment Requirements

Main dependencies:
- scanpy
- pandas
- numpy
- matplotlib
- seaborn
- scikit-image
- scipy

For detailed dependencies, please refer to `requirements.txt`

## Usage Instructions

1. Ensure all necessary dependencies are installed
2. Navigate to the corresponding folder based on specific analysis needs
3. Follow the instructions in the notebook to run the analysis pipeline
4. Pay attention to modifying data paths and parameter configurations

## Citation

If you use this analysis code, please cite the relevant PRISM paper.

## Contact

For questions or suggestions, please contact the project maintainers.
