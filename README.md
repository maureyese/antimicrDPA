# Antimicrobial Activity of DPA in exposure to skin microbiome from DPAHelix project

This repository corresponds to one of the five computational biology analyses of the DPAHelix project in the GOGEC Competition of 2026.

> **Objective**: Analyze imbalance of skin microbiome under DPA-based sunscreen mixture.

This analysis is based on a network pharmacology approach. To predict the impact of DPA as a main active sunscreen component on the skin microbiome, we employed a network pharmacology-based framework. We first identified core genera of the skin microbiome and curated a list of their essential genes using the Database of Essential Genes (DEG). The tridimensional structures of these protein targets were retrieved from the AlphaFold Database. A blind molecular docking screening was performed with DPA as the ligand to identify potential blinding interaction. The resulting high-confidence binding events were used to construct a ligand-protein interaction network. Finally, we performed functional enrichment analysis on the targeted proteins to elucidate the biological processes and pathways most susceptible to modulation by DPA. This approach provides a hypothesis for its potential antimicrobial and dysbiosis-inducing effects.

We worked using VS Code, along with the WSL extension (Ubuntu 22.04) as our development environmnet (Only-Windows users).

WSL Installation: https://learn.microsoft.com/en-us/windows/wsl/install
VS Code Installation: https://code.visualstudio.com/docs/setup/windows#_install-vs-code-on-windows
WSL extension on VS Code: https://code.visualstudio.com/docs/remote/wsl

## Key resources

1. DEG database: http://origin.tubic.org/deg/public/index.php
_(Finish setting resources)_

## Activity 1. Identify key skin microbiome on different skin surface zones.

_(under development...)_

## Activity 2. Identify representatitve essential genes from main genera

- Identify potential extracellular transporters or receptors: https://bioinformatics.ysu.edu/tools/subcell.html

_(under development...)_

## Activity 3. Perform blind docking screening using DPA

- Set a threeshold to validate potential binding affinities.

_(under development...)_

## Activity 4. Generate ligand-protein network

_(under development...)_

## Activity 5. Perform functional enrichment on critical targets