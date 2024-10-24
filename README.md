README for Microbial Community Invasion Simulations and Sequencing Data

Overview

This repository contains MATLAB code and R scripts used in the study "Collective dynamical regimes predict invasion success and impacts in microbial communities." The code enables the simulation of microbial community dynamics, including pH-mediated interactions, as well as the processing of sequencing data related to community compositions. Below is a description of the key scripts, how to run them, and the necessary working environment.

File Descriptions

1. LV_Invade.m
This script implements a generalized Lotka-Volterra (gLV) model to simulate microbial community dynamics with species interactions. It includes the function LV_compute_invasion, which calculates invasion outcomes for the invader species, allowing you to evaluate its effects on the resident community.

Dependencies: There are no external function dependencies, as LV_compute_invasion is fully defined within this script.
Input: Species interaction matrices and species pool size.
Output: Time series data of species abundances and invasion success.

2. pH_Invade.m
This script extends the gLV model by introducing pH-mediated interactions, where pH impacts the growth rates and interaction strengths between species. It simulates how varying pH conditions influence invasion dynamics.

Dependencies: As with LV_Invade.m, all necessary functions are self-contained.
Input: Similar to LV_Invade.m, with additional parameters related to pH effects.
Output: Time series data reflecting changes in community composition due to pH-mediated interactions.

3. Hu-Merging-SeqWorkflow.R
This R script processes raw sequencing data and generates amplicon sequence variants (ASVs) using the DADA2 pipeline. It is used to analyze microbial community compositions before and after invasions.

Input: FASTQ sequencing files.
Output: ASV tables and taxonomic classifications.
Working Environment

MATLAB Setup:
MATLAB Version: The scripts were developed and tested on MATLAB 2024.
Required Toolboxes: Ensure that the Statistics and Machine Learning Toolbox and Optimization Toolbox are installed, as they may be required for data fitting and simulations.
Function Dependencies: All necessary functions, including LV_compute_invasion, are defined within the provided scripts. No external functions are required.
R Setup:
R Version: The R script requires R version 3.6 or later.
R Packages: The following R packages must be installed:
DADA2: for processing sequencing data.
ggplot2: for visualizing results.
phyloseq: for downstream analysis of ASVs.
Additional Resources: Please refer to the DADA2 documentation for more details on setting up your working environment.
Running the Code

MATLAB:
Load the relevant interaction matrices and species pool sizes into MATLAB.
Run LV_Invade.m for simulations using the generalized Lotka-Volterra model. The script includes the LV_compute_invasion function, so no additional files are needed.
Optionally, run pH_Invade.m for simulations incorporating pH-mediated effects on species interactions.
R:
Load the sequencing data (FASTQ files) and sample sheet into your R environment.
Run Hu-Merging-SeqWorkflow.R to generate ASVs and assign taxonomic identities using the SILVA database (version 132).
Use phyloseq for further community composition analysis.
Dryad Data Repository

The raw sequencing data used in this study is available on Dryad at the following link:
DOI: 10.5061/dryad.8gtht76xz

This repository contains all data needed to reproduce the analyses presented in the manuscript, including the processed ASV tables and community metadata.
