# cilia-ad-ppi
# Ciliary Genes and Alzheimer's Disease: Protein-Protein Interaction (PPI) Analysis
This repository contains the code and analysis for the paper "Analysis of Ciliary Genes Reveals Novel Protein-Protein InteractionsTherapeutic Targets in Alzheimer's Disease". 
The study explores the role of ciliary genes in Alzheimer's Disease (AD) progression and identifies potential therapeutic targets by predicting novel protein-protein interactions (PPIs) and performing gene-drug interaction analysis.

## Table of Contents
- [Overview](#overview)
- [Data Sources](#data-sources)
- [Code Description](#code-description)
  

  ## Overview
  The goal of this study is to investigate the contribution of ciliary genes to Alzheimer's disease, with a specific focus on:
    - Identifying differentially expressed ciliary genes in AD using public datasets.
    - Predicting novel protein-protein interactions (PPIs) within the ciliary interactome.
    - Performing pathway and gene-drug interaction analysis to identify potential therapeutic targets, including drug repurposing opportunities.
  ## Data Sources
  The following datasets are used in this study:
    - CiliaCarta Compendium – A comprehensive list of ciliary genes.
    - GWAS (Genome-Wide Association Studies) – Alzheimer's disease-related genes
    - Human Protein Reference Database (HPRD) & BioGRID – Protein-protein interactions (PPIs).
    - Drug Gene Interaction Database (DGIdb) – Gene-drug interactions.
    - Gene Expression Omnibus (GEO) Series:
          - GSE48350
          - GSE29378
          - GSE28146
          -  GSE53987 (Schizophrenia dataset for comparison).
  ## Code Description
  1. PPI Prediction Script
  The main script (ppi_prediction.py) includings the following steps
  - Load ciliary gene data: Data is cleaned processed from from the source files
  - Interaction analysis: Known and novel interactors of ciliary genes are identified using PPI datasets
  - Novel interaction prediction: A A Random Forest model predicts novel PPIs based on known datasets
 


