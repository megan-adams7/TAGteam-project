# TAGteam-project

This repository contains code, data, and documentation for a research project analyzing the conservation and distribution of TAGteam motifs in the upstream regions of genes across 300+ *Drosophila* species.

## Project Overview

This project explores the abundance, distribution, and evolutionary conservation of TAGteam motifs—short DNA sequences involved in early gene regulation—across over 300 Drosophila species. TAGteam motifs are known to play a key role in the maternal-to-zygotic transition, and their presence upstream of genes can influence developmental gene expression.

Using custom computational pipelines written in Python and Bash, we process genomic annotation files to extract upstream regions of genes, scan for canonical and variant TAGteam motifs, and analyze their frequency and positional patterns. The project aims to identify which genes are consistently associated with high TAGteam density, detect lineage-specific motif losses or gains, and assess broader phylogenetic trends.

Results are visualized using R to highlight gene- and species-specific differences, contributing insights into the evolution of regulatory networks and developmental timing across the Drosophila genus.

## Tools and Technologies
- **Languages**: Python, Bash, R
- **Libraries**: `bedops`, `bedtools`, `seqkit`, `ggplot2`, `gggenomes`
- **Data**: GFF/GTF gene annotations, FASTA genome files for >300 *Drosophila* species
- **Environment**: UNIX-based command line, RStudio, VS Code

## Features

- motif scanning
- multi-species analysis
- unbiased motif search
- data-visualizations

## Repository Structure

<pre> ``` TAGteam-project/ │ ├── gtf_files # gtf files for Drosophila species ├── zipped_fasta # fasta files for Drosophila species ├── gene_name # Working directory │ ├── main.sh # Main script │ ├── vis.Rmd # Visualization pipeline │ ├── table.csv # empty csv table │ ├── gene_name_protein.fasta # sample gene protein sequence fasta file └── README.md ``` </pre>

## Sample results
![TAGteam Motif Distribution](/sna_plot.png)