# sars2nterm
Code Repository for SARS-CoV-2 N-terminomics data from Meyer et al. 2020. This repository contains all scripts and data files needed to process and reproduce the figures from the [Meyer et al. bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2020.09.16.297945v1). 

This repository contains the latest version of the code, if you need to access specific versions associated with the different preprint releases, please use the relevent release.

# Abstract:
Many events in viral infection are regulated through post-translational modification. In the case of positive-sense RNA viruses such as SARS-CoV-2 which encodes two viral proteases, proteolytic cleavage of both viral and cellular proteins is crucial for productive viral infection. Through this work we identified multiple viral proteins cleaved during the infection timecourse, including crucial antigens such as Spike and N, as well as cellular proteins which are cleaved at positions matching the consensus sequences for cleavage by the viral proteases. We further showed that siRNA inhibition and inhibitor treatment targeting these cellular proteins inhibited viral replication and the production of infectious virus, suggesting avenues for treatments for COVID-19 disease.

# In detail: 
We studied proteolysis in the context of SARS-CoV-2 infection of ACE2-A549 and Vero E6 cells, infected at a multiplicity of infection of 1, with samples harvested at 0,6,12 and 24h post-infection. Neo-N-termini generated by protease activity and lysine residues were labelled at the protein level with TMTpro reagents prior to tryptic digestion. For unenriched samples, the trypsin-digested material was then fractionated offline, and analysed by LC-MS/MS for the study of total protein abundance over the infection timecourse. For enriched samples, unblocked neo-N-termini generated by tryptic digested were labelled with undecanal and depleted, yielding a sample enriched for blocked neo-N-termini (TMTpro, N-terminal acetylation, and pyroGlu). The enriched samples were similarly fractionated offline and analysed by LC-MS/MS. 

These scripts take and process MaxQuant output files, as well as analysing and plotting qRT-PCR Ct values, and cell viability and inhibitor data. Only the minimal proteomic data files are provided here that are required for these scripts. The full MaxQuant output files are available through the PRIDE repository (PXD021145, PXD021152, PXD021153, PXD021154, PXD021402 & PXD021403 - pre-dataset release, reviewer usenames and passwords are available to access this data in the preprint).

# To use these scripts:
- All scripts access a global variable 'path', which should be the folder you download this archive to, you will need to update path to the directory where you have saved this data.

# Scripts 
Note: all scripts are self-contained. 

- Fig_TimecoursePlots_20200923.m Analyses and produces the plots for Figure 1.
- figure2_ViralNtermini_20200923.m Analyses and produces the plots and tables for Figure 2, and Tables S1, S2
- EnrichedCellularNterm_20200923.m - Analyses and produces the plots and tables for Figure 4, and Tables S4, S5
- siRNA_Fig_20200923.m - Analyses and produces the plots for Figure 5
- Inhibitors_20200923.m - Analyses/Fits inhibitor data and produces the plots for Figure 6, S7 and S8
- ntermQC_20200923.m - This analyses the proteomics data and produces the plots for Figure S1.
- SFig_SupViralNtermini_20200923.m - This analyses the proteomics data and produces the plots for Figure S2.
- siRNA_cytotox_20200923.m - This plots cytotoxicity data for the siRNA experiment for Fig S5.
- siRNA_KO_efficiency_20200923.m - This analyses qRT-PCR data to calculate knockown efficiency for the siRNA experiment, shown in Fig S6.
- STable_AllViralPeptides20200923.m - This collates all the viral peptides identified accross the various experiments, generating table S3
- ntermQC_ProteaseInhibitor.m - Analyses and produces the plots for Figure S3 (Protease inhibitor HUNTER experiment).
-ProteaseInhibitorHunter.m - Analyses and produces the plots for Figure 3 and S3 (Protease inhibitor HUNTER experiment).

# Other files:
- Minimal data files required for the analysis are in /data/ to permit all the scripts to generate the data as seen in the preprint - for the full data see the PRIDE repositorys
- The experimental design files necessary to reorder the TMT labels to account for label randomisation are in /data/ and are 'SARS2a549tmtlabelling_20200507.csv' and 'Verotmtlabelling_20200511.csv' respectively.
- Protein annotations for A549 and Vero indicating signal and transit peptides for annotating N-termini results in /data/. human_acc_to_gn_signaltransit.csv and vero_acc_to_gn_signaltransit.csv respectively.
- 6x6p.cxv chimeraX command file is included to reproduce the structure visualisation in Figure 2. See /other.
- Custom .fasta file for SARS2 proteins, including those identified by Riboseq (Finkel et al. 2020). This is 'SARS2_Custom_20200518.fasta' and in /data/ 

Code requires Matlab 2019b and the Statistics and Machine Learning toolbox. Some of the functions e.g. tiledlayout were only introduced in R2019b, the scripts will work in earlier versions but will require some adjustment (e.g. tiledlayout -> subplot).

Code: Ed Emmott (Github/Twitter: @edemmott, Email: e.emmott@liverpool.ac.uk), University of Liverpool, 2020.
