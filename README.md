# Diversity-Stability Relationships in Metacommunities

This repository contains code to replicate the analysis of 

Wisnoski, N.I., R. Andrade, M.C.N. Castorani, C.P. Catano, A. Compagnoni, T. Lamy, N.K. Lany, L. Marazzi, S. Record, A.C. Smith, C.M. Swan, J.D. Tonkin, N.M. Voelker, P.L. Zarnetske, E.R. Sokol. Diversity-stability relationships become decoupled across spatial scales: a synthesis of organism and ecosystem types. Ecology. 

The project is part of a broader synthesis working group funded by NSF Long-Term Ecological Research (LTER) Network and supported by NCEAS. The main respository for the broader working group can be found here: https://github.com/sokole/ltermetacommunities.

To replicate our analysis run the code in the following order: 

1. `0_fetch_data.R` = This is a preliminary step needed to download the datasets from publicly available data repositories, mostly the Environmental Data Initiative (EDI) and Ecological Archives (hosted on Figshare). Running this code will download and save formatted datasets into the subdirectory `data` in this repository. 

2. `1_analyze_variability.R` = This script reads in the saved datasets and computes a number of diversity and variability metrics on the data, then saves the output into the `results` folder. 

3. `2_dsr_analysis.R` = This script computes relationships between diversity and variability at different spatial scales. It writes a summary table of results to the `results` folder. 

4. `3_sem.R` = This script performs a structural equation model relating diversity and variability at different spatial scales of the analysis. 


