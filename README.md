Title: Deeper Caribbean reef fish communities exhibit increasing taxonomic and functional distinctiveness in dominance structure over a nine-year period

Authors: James S. Boon*, Sally A. Keith, Dan A. Exton, Erika Gress, Dominic A. Andradi-Brown and Richard Field*

*School of Geography, University of Nottingham, Nottingham, UK

Corresponding author: James S. Boon

Accepted for Publication in Coral Reefs

This repository contains the data and code required to undertake the analysis presented in the paper

* The folder "species_abundance_and_trait_data" contains two datasets. The first is called "abundance_sites_alpha", which contains the primary abundance data used in the analysis. Each column is a site/depth/period combination, rows are species, and inputs are their individual abundances. The other file is "species_traits", which contains the species traits we used in the functional analysis.

*The folder "inext_outputs" contains 3 folders. Each contain .csv files with the outputs from the inext analysis which i used to make the plots seen in the paper. I found it tricky to extract the data i needed from the inext outputs neatly in r, so i just uploaded it into a csv file and added the seperate site, depth and time period columns i needed to plot the results. The "inext_outputs_alpha" folder contains two files with the outputs from the alpha diversity analysis, one for the taxonomic ("alpha_diversity_output_sites_depthscombined") and another for the functional ("functional_alpha_data_qo_revised") data. The "inext_outputs_beta" folder contains two files with the outputs from the beta diversity analysis, one for the taxonomic ("taxonomic_beta_groupeddepth_byyear") and another for the functional ("functional_beta_groupeddepth_byyear") data. Finally, the "inext_outputs_heatmaps" folder contains several files with the outputs from the beta diversity analysis for each individual site and is the data we used to make the heatmaps. 
