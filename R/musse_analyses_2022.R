library("diversitree")
library("tidyverse")
library("foreach")
library("doMC")
library("BAMMtools")
library("V.PhyloMaker")
library("taxonlookup")
library("caper")
library("ape")
library("phangorn")
library("phytools")
library("geiger")
library("picante")

## GBIF.org (02 June 2022) GBIF Occurrence Download  https://doi.org/10.15468/dl.beduwc

here::i_am("R/musse_analyses_2022.R")

## Importing trees
fulltree.mol <- GBOTB.extended
fulltree.full <- read.tree(here::here("data/smith_brown_2018_trees/ALLOTB.tre"))

## Importing mycorrhizal data
species.data <- read.csv(here::here("data/full_database_Soudzilovskaia_2019.csv"), na.strings = "NA")
genus.data <- read.csv(here::here("data/genus_database_Soudzilovskaia.csv"), na.strings = "NA")

## Pruning trees to keep only species with data
tree.mol.species <- drop.tip(fulltree.mol, tip = fulltree.mol$tip.label[!(fulltree.mol$tip.label %in% species.data$species)])
tree.full.species <- drop.tip(fulltree.full, tip = fulltree.full$tip.label[!(fulltree.full$tip.label %in% species.data$species)])
tree.full.species <- multi2di(tree.full.species)
tree.full.species <- force.ultrametric(tree.full.species, method = "extend")

## ## Showing the distribution of species with data in the phylogenies
## png(filename = here::here("output/sb_moltree_with_data_species_level.png"), width = 15, height = 15, units = "in", res = 300)
## #pdf(file = here::here("output/sb_moltree_with_data_species_level.pdf"), width = 15, height = 15)
## plot(fulltree.mol, type = "fan", show.tip.label = FALSE, edge.col = "grey")
## tiplabels(tip = match(tree.mol.species$tip.label, fulltree.mol$tip.label), pch = ".", col = "red", bg = "red", cex = 1, offset = 5)
## dev.off()

## png(file = here::here("output/sb_fulltree_with_data_species_level.png"), width = 15, height = 15, units = "in", res = 300)
## #pdf(file = here::here("output/sb_fulltree_with_data_species_level.pdf"), width = 15, height = 15)
## plot(fulltree.full, type = "fan", show.tip.label = FALSE, edge.col = "grey")
## tiplabels(tip = match(tree.full.species$tip.label, fulltree.full$tip.label), pch = ".", col = "red", bg = "red", cex = 1, offset = 5)
## dev.off()

## Reclassifying mycorrhizal types according to our previous criteria
species.data$mycorrhiza.type[which(species.data$curator_remark_1_comment == "probably AM")] <- "AM"
species.data$curator_remark_1_comment[which(species.data$curator_remark_1_comment == "probably AM")] <- NA
species.data$mycorrhiza.type[which(species.data$curator_remark_1_comment == "probably ErM")] <- "ER"
species.data$curator_remark_1_comment[which(species.data$curator_remark_1_comment == "probably ErM")] <- NA
species.data$mycorrhiza.type[which(species.data$curator_remark_1_comment == "probably NM")] <- "NM"
species.data$curator_remark_1_comment[which(species.data$curator_remark_1_comment == "probably NM")] <- NA
species.data$mycorrhiza.type[which(species.data$curator_remark_1_comment == "probably OM")] <- "OM"
species.data$curator_remark_1_comment[which(species.data$curator_remark_1_comment == "probably OM")] <- NA

reclass.myco <- data.frame(
    original = unique(species.data$mycorrhiza.type),
    new = c("OM", "AM", "EM", "EM", "NM", "AM", "AM", "MIX", NA, NA, "NM", "ER", "NV", "ER", "MIX", "MIX", "OTHER")
)

species.data$new.myco <- reclass.myco$new[match(species.data$mycorrhiza.type, reclass.myco$original)]

## Removing species without phylogenetic data from the datasets
species.data.clean.mol <- species.data[species.data$species %in% tree.mol.species$tip.label,]
species.data.clean.full <- species.data[species.data$species %in% tree.full.species$tip.label,]

## Aggregating types per species
agg.data.mol <- aggregate(species.data.clean.mol$new.myco, by = list(species.data.clean.mol$species), FUN = unique)
agg.data.full <- aggregate(species.data.clean.full$new.myco, by = list(species.data.clean.full$species), FUN = unique)

## Recoding species with multiple types to be MIX
agg.data.mol[sapply(agg.data.mol[,2], function(x){length(x) > 1}), 2] <- "MIX"
agg.data.full[sapply(agg.data.full[,2], function(x){length(x) > 1}), 2] <- "MIX"

agg.data.mol[,2] <- unlist(agg.data.mol[,2])
agg.data.full[,2] <- unlist(agg.data.full[,2])

names(agg.data.mol) <- c("species", "trait")
names(agg.data.full) <- c("species", "trait")

agg.data.mol <- agg.data.mol[!is.na(agg.data.mol$trait),]
agg.data.full <- agg.data.full[!is.na(agg.data.full$trait),]

agg.data.mol <- agg.data.mol[agg.data.mol$trait != "OTHER",]
agg.data.full <- agg.data.full[agg.data.full$trait != "OTHER",]

## Removing species with OM and ER
agg.data.mol.noomer <- agg.data.mol[!(agg.data.mol$trait %in% c("OM", "ER")),]
agg.data.full.noomer <- agg.data.full[!(agg.data.full$trait %in% c("OM", "ER")),]

agg.data.mol$trait.num <- match(agg.data.mol$trait, sort(unique(agg.data.mol$trait)))
agg.data.full$trait.num <- match(agg.data.full$trait, sort(unique(agg.data.full$trait)))

agg.data.mol.noomer$trait.num <- match(agg.data.mol.noomer$trait, sort(unique(agg.data.mol.noomer$trait)))
agg.data.full.noomer$trait.num <- match(agg.data.full.noomer$trait, sort(unique(agg.data.full.noomer$trait)))

## write.tree(tree.mol.species, file = here::here("output/musse/species_data/moltree.tre"))
## write.tree(tree.full.species, file = here::here("output/musse/species_data/fulltree.tre"))

## write.table(agg.data.mol[, -2], file = here::here("output/musse/species_data/moldata.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
## write.table(agg.data.full[, -2], file = here::here("output/musse/species_data/fulldata.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

## write.table(agg.data.mol.noomer[, -2], file = here::here("output/musse/species_data/moldata_noomer.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
## write.table(agg.data.full.noomer[, -2], file = here::here("output/musse/species_data/fulldata_noomer.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

### ANALYSES SHOULD BE ALSO DONE WITH THE EXTERNAL SCRIPT FOR THE BAYESIAN IMPLEMENTATION


#######################################
#### Preparing dataset using genus data
#######################################

genus.mol <- data.frame(species = fulltree.mol$tip.label,
                        genus = sapply(strsplit(fulltree.mol$tip.label, split = "_", fixed = TRUE), "[[", 1))
genus.full <- data.frame(species = fulltree.full$tip.label,
                        genus = sapply(strsplit(fulltree.full$tip.label, split = "_", fixed = TRUE), "[[", 1))

genus.mol$trait <- genus.data$type[match(genus.mol$genus, genus.data$genus)]
genus.full$trait <- genus.data$type[match(genus.full$genus, genus.data$genus)]

genus.mol$trait[which(!is.na(match(genus.mol$species, agg.data.mol$species)))] <- na.omit(agg.data.mol$trait[match(genus.mol$species, agg.data.mol$species)])
genus.full$trait[which(!is.na(match(genus.full$species, agg.data.full$species)))] <- na.omit(agg.data.full$trait[match(genus.full$species, agg.data.full$species)])

reclass.myco.genus <- data.frame(
    original = unique(genus.full$trait),
    new = c("AM", NA, "EM", "MIX", "UNC", "NM", "MIX", "MIX", "ER", "EM", "MIX", "NM", "MIX", "OM", "ER", "MIX", "OTHER")
    )

genus.mol$new.trait <- reclass.myco.genus$new[match(genus.mol$trait, reclass.myco.genus$original)]
genus.full$new.trait <- reclass.myco.genus$new[match(genus.full$trait, reclass.myco.genus$original)]

genus.mol <- genus.mol[!is.na(genus.mol$trait),]
genus.full <- genus.full[!is.na(genus.full$trait),]

genus.mol.noomer <- genus.mol[!(genus.mol$new.trait %in% c("OM", "ER", "OTHER", "UNC")), ]
genus.full.noomer <- genus.full[!(genus.full$new.trait %in% c("OM", "ER", "OTHER", "UNC")), ]

genus.mol$trait.num <- match(genus.mol$new.trait, sort(unique(genus.mol$new.trait)))
genus.full$trait.num <- match(genus.full$new.trait, sort(unique(genus.full$new.trait)))

genus.mol.noomer$trait.num <- match(genus.mol.noomer$new.trait, sort(unique(genus.mol.noomer$new.trait)))
genus.full.noomer$trait.num <- match(genus.full.noomer$new.trait, sort(unique(genus.full.noomer$new.trait)))

## write.table(genus.mol[, c(1, 5)], file = here::here("output/musse/genus_data/moldata.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
## write.table(genus.full[, c(1, 5)], file = here::here("output/musse/genus_data/fulldata.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

## write.table(genus.mol.noomer[, c(1, 5)], file = here::here("output/musse/genus_data/moldata_noomer.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
## write.table(genus.full.noomer[, c(1, 5)], file = here::here("output/musse/genus_data/fulldata_noomer.txt"), quote = FALSE, row.names = FALSE, sep = "\t")


## MuSSE - Standard implementation

### Species level data

tree.mol.species.musse <- drop.tip(tree.mol.species, tree.mol.species$tip.label[!(tree.mol.species$tip.label %in% agg.data.mol[,1])])
musse.mol.sp <- make.musse(tree.mol.species.musse, setNames(agg.data.mol[,3], agg.data.mol[,1]), k = length(unique(agg.data.mol[, 3])), sampling.f = Ntip(tree.mol.species)/Ntip(fulltree.full))

tree.full.species.musse <- drop.tip(tree.full.species, tree.full.species$tip.label[!(tree.full.species$tip.label %in% agg.data.full[,1])])
musse.full.sp <- make.musse(tree.full.species.musse, setNames(agg.data.full[,3], agg.data.full[,1]), k = 6, sampling.f = Ntip(tree.full.species)/Ntip(fulltree.full))

starting.pars.mol.species.musse <- starting.point.musse(tree.mol.species.musse, k = length(unique(agg.data.mol[,3])))
starting.pars.full.species.musse <- starting.point.musse(tree.full.species.musse, k = length(unique(agg.data.full[,3])))


tree.mol.species.musse.noomer <- drop.tip(tree.mol.species, tree.mol.species$tip.label[!(tree.mol.species$tip.label %in% agg.data.mol.noomer[,1])])
musse.mol.sp.noomer <- make.musse(tree.mol.species.musse.noomer, setNames(agg.data.mol.noomer[,3], agg.data.mol.noomer[,1]), k = length(unique(agg.data.mol.noomer[,3])), sampling.f = Ntip(tree.mol.species.musse.noomer)/Ntip(fulltree.full))

tree.full.species.musse.noomer <- drop.tip(tree.full.species, tree.full.species$tip.label[!(tree.full.species$tip.label %in% agg.data.full.noomer[,1])])
musse.full.sp.noomer <- make.musse(tree.full.species.musse.noomer, setNames(agg.data.full.noomer[,3], agg.data.full.noomer[,1]), k = length(unique(agg.data.full.noomer[,3])), sampling.f = Ntip(tree.full.species.musse.noomer)/Ntip(fulltree.full))

starting.pars.mol.species.musse.noomer <- starting.point.musse(tree.mol.species.musse.noomer, k = length(unique(agg.data.mol.noomer[,3])))
starting.pars.full.species.musse.noomer <- starting.point.musse(tree.full.species.musse.noomer, k = length(unique(agg.data.full.noomer[,3])))

musse.mol.species.noomer <- find.mle(musse.mol.sp.noomer, starting.pars.mol.species.musse.noomer, method = "subplex")
saveRDS(musse.mol.species.noomer, file = here::here("output/musse_mol_species_noomer_ML.RDS"))

musse.full.species.noomer <- find.mle(musse.full.sp.noomer, starting.pars.full.species.musse.noomer, method = "subplex")
saveRDS(musse.full.species.noomer, file = here::here("output/musse_full_species_noomer_ML.RDS"))

## musse.mol.species <- find.mle(musse.mol.sp, starting.pars.mol.species.musse, method = "subplex")
## saveRDS(musse.mol.species, file = here::here("output/musse_mol_species_ML.RDS"))

## musse.full.species <- find.mle(musse.full.sp, starting.pars.full.species.musse, method = "subplex")
## saveRDS(musse.full.species, file = here::here("output/musse_full_species_ML.RDS"))

### Genus level data

fulltree.genus <- multi2di(fulltree.mol)

tree.mol.genus.musse <- drop.tip(fulltree.genus, fulltree.genus$tip.label[!(fulltree.genus$tip.label %in% genus.mol$species)])
tree.mol.genus.musse.noomer <- drop.tip(fulltree.genus, fulltree.genus$tip.label[!(fulltree.genus$tip.label %in% genus.mol.noomer$species)])

musse.mol.genus.noomer <- make.musse(tree.mol.genus.musse.noomer, setNames(genus.mol.noomer$trait.num, genus.mol.noomer$species), k = length(unique(genus.mol.noomer$trait.num)), sampling.f = Ntip(tree.mol.genus.musse.noomer)/Ntip(fulltree.full))
starting.pars.mol.genus.musse.noomer <- starting.point.musse(tree.mol.genus.musse.noomer, k = length(unique(genus.mol.noomer$trait.num)))
musse.mol.genus.noomer <- find.mle(musse.mol.genus.noomer, starting.pars.mol.genus.musse.noomer, method = "subplex")
saveRDS(musse.mol.genus.noomer, file = here::here("output/musse_mol_genus_noomer_ML.RDS"))

musse.mol.genus <- make.musse(tree.mol.genus.musse, setNames(genus.mol$trait.num, genus.mol$species), k = length(unique(genus.mol$trait.num)), sampling.f = Ntip(tree.mol.genus.musse)/Ntip(fulltree.full))
starting.pars.mol.genus.musse <- starting.point.musse(tree.mol.genus.musse, k = length(unique(genus.mol$trait.num)))
musse.mol.genus <- find.mle(musse.mol.genus, starting.pars.mol.genus.musse, method = "subplex")
saveRDS(musse.mol.genus, file = here::here("output/musse_mol_genus_ML.RDS"))
