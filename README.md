# "Seed plant families with diverse mycorrhizal states have higher diversification rates"

Codes and data for Mujica et al. Link to the preprint will be available soon.

## Instructions to reproduce the manuscript

1. Set the working directory to the root of the repository
2. Run the following commands:
> library("rmarkdown")
> rmarkdown::render("./manuscript/mujica_etal.Rmd")
> rmarkdown::render("./manuscript/supp_mat.Rmd")
3. PDFs will be available in the ./'manuscript' folder

Should you want to re-run the analysis, run the following replacing 'filename' by the corresponding script:
> source("./R/filename.R")
but keep in mind that some of the analyses take a few seconds/minutes to run, so be patient.

Some of the figures needed to be adapted using the vector graphics editor [inkscape](https://inkscape.org/). The SVG files are available within the './output' folder, along with the raw pdf versions of the figures straight from R.

## Abstract

Most of plant species have mycorrhizas, which can be classified in four types: Arbuscular (AM), Ecto (EM), Orchid (OM) and Ericoid Mycorrhiza (ER). Since the AM ancestral state, some plant lineages have switched partner (EM, OM and ER) or lost the association (NM). Evolutionary transitions to a novel mycorrhizal state (MS) might allow plant lineages to access new resources, enhancing diversification rates. However, some clades are not restricted to one MS, and this variability might promote diversification. Here, we address the relationship between MS diversity and seed plant diversification. Using the Fungal-root database, which compiles plant species and their MS, we assigned a single MS to each plant family, calculated the MS heterogeneity and estimated their diversification rates using the method-of-moments. Our results showed higher diversification rates in families with mixed MS, and a positive relationship between MS heterogeneity and diversification rates, which suggests that MS lability promotes plant diversification.
