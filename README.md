# cilia-gene
# Analysis of Ciliary Genes and Their Role in Alzheimer's Disease
This repository contains the code and analysis for the paper "Analysis of Ciliary Genes Reveals Novel Protein-Protein Interactions Therapeutic Targets in Alzheimer's Disease". 
The study explores the role of ciliary genes in Alzheimer's Disease (AD) progression and identifies potential therapeutic targets by predicting novel protein-protein interactions (PPIs) and performing gene-drug interaction analysis.

## Table of Contents
- [Overview](#overview)
- [Data Sources](#data-sources)
- [Code Description](#code-description)
- [LICENSE](#license)
  
  

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
          - GSE53987 (Schizophrenia dataset for comparison).
  
  In addition to the primary data, I have also generated new datasets through various data handling and processing steps, such as combining,
  filtering, and transforming the original data. These processed datasets are uploaded separately and reflect the data used in the final analysis.
  ## Code Description
  1. PPI Prediction Script
  2. Protein and Drug Interaction 
  3. Comparison of Interaction Scores between Cilia and Non-Cilia Genes with Target Drugs
  4. Gene Expression Analysis Based on Age 
  5. Gene Expression Analysis Based on Gender

  Note: Data handling, cleaning, normalization, and transformation code is excluded from this file
  ## LICENSE
  MIT License

  Copyright (c) [2024] [Thi Ngo]

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

