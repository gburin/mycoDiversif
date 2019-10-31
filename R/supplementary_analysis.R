###
library("caper")
library("tidyverse")
library("cowplot")
library("RColorBrewer")
library("geiger")
library("phytools")
library("ggrepel")

#######################
### Analysis per genera including all families
#######################

family.data.gen <- read.csv("./data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

age.data <- read.csv("./data/data_all_families.csv", sep = ";")

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.raw.perc == 1),]

## Checking which families have less than 5% or 8spp in the dataset
which(family.data.gen$perc.sp.unk > 0.95)
which((family.data.gen$rich - family.data.gen$UNK) <= 8)

family.data.gen$stem.age <- age.data$stem.age[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e05 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.5)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae"), family.data.gen$family),]

## Checking representativeness of each family

sum(family.data.gen$NA.perc > 0.95)
sum((family.data.gen$rich - family.data.gen$UNK) < 8)

## min(((family.data.gen$rich - family.data.gen$UNK)/family.data.gen$rich))

## We will not remove any family due to poor sampling since all families have either at least 58.3% of total richness or a minimum of 8 species sampled

family.data.gen <- family.data.gen
family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
data.pgls <- comparative.data(tree.pruned, family.data.gen, names.col = "family")

phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r05 <- pgls(r.e05 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")
phylosig.age <- pgls(stem.age ~ 1, data = data.pgls, lambda = "ML")
phylosig.rich <- pgls(rich ~ 1, data = data.pgls, lambda = "ML")

## Fitting PGLS
mod.r0 <- caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
mod.r05 <- caper::pgls(r.e05 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r05)$param[2], NULL))
mod.r09 <- caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r09)$param[2], NULL))

## Fitting standard linear models
lm.r0 <- lm(r.e0 ~ shannon, data = family.data.gen)
lm.r05 <- lm(r.e05 ~ shannon, data = family.data.gen)
lm.r09 <- lm(r.e09 ~ shannon, data = family.data.gen)

## phylANOVA
data.aov <- family.data.gen[-which(is.na(match(family.data.gen$family, fulltree$tip.label))), ]
data.aov <- data.aov[-which(data.aov$type.50 == "ER"), ]
data.aov <- data.aov[-which(is.na(data.aov$r.e0)),]

phyaov.r0 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
phyaov.r05 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e05, data.aov$family))
phyaov.r09 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

## Age vs rich
pgls.age.sh <- caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.age)$param[2], NULL))
pgls.rich.sh <- caper::pgls(rich ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.rich)$param[2], NULL))

lm.age.sh <- lm(stem.age ~ shannon, data = family.data.gen)
lm.rich.sh <- lm(rich ~ shannon, data = family.data.gen)


## Appending results to original data frame for plotting
family.data.gen$lm.r0[which(!is.na(family.data.gen$r.e0))] <- fitted(lm.r0)
family.data.gen$lm.r05[which(!is.na(family.data.gen$r.e0))] <- fitted(lm.r05)
family.data.gen$lm.r09[which(!is.na(family.data.gen$r.e0))] <- fitted(lm.r09)

family.data.gen$pgls.r0 <- fitted(mod.r0)[,1][match(family.data.gen$family, names(fitted(mod.r0)[,1]))]
family.data.gen$pgls.r05 <- fitted(mod.r05)[,1][match(family.data.gen$family, names(fitted(mod.r05)[,1]))]
family.data.gen$pgls.r09 <- fitted(mod.r09)[,1][match(family.data.gen$family, names(fitted(mod.r09)[,1]))]

save.image(file = "./output/suppdata_genus_full.RData")


######################
##Analysis per species
######################

family.data <- read.csv("./data/family_data_full.csv", stringsAsFactors = FALSE)

## Filtering for families with more than 5% or 8 species

family.data.clean <- family.data[family.data$samp.sp >= 8 | family.data$perc.rich >= 0.05,]

## Identifying the families used in the study for later plotting
family.data$study <- "no"
family.data$study[family.data$samp.sp >= 8 | family.data$perc.rich >= 0.05] <- "yes"

fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Removing Orchidaceae and Ericaceae from analyses
family.data.clean <- family.data.clean[-match(c("Orchidaceae", "Ericaceae"), family.data.clean$family),]

## Dropping families without reliable mycorrhizal data
tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.clean$family))])

## Calculating diversification rates using the method-of-moments (Magallón & Sanderson, 2001)
family.data.clean$r.e0.stem <- bd.ms(time = family.data.clean$stem.age, n = family.data.clean$global.rich, crown = FALSE, epsilon = 0)
family.data.clean$r.e05.stem <- bd.ms(time = family.data.clean$stem.age, n = family.data.clean$global.rich, crown = FALSE, epsilon = 0.5)
family.data.clean$r.e09.stem <- bd.ms(time = family.data.clean$stem.age, n = family.data.clean$global.rich, crown = FALSE, epsilon = 0.9)

family.data.clean$shannon <- vegan::diversity(family.data.clean[, c(4, 5, 6, 8, 9)], index = "shannon")

## Preparing data for PGLS
data.pgls <- comparative.data(tree.pruned, family.data.clean[, -16], names.col = "family")

## Testing for correlation between Stem age and Species Richness (both considering or not the phylogenetic structure)

stem.rich.phylo <- pgls(global.rich ~ stem.age, data = data.pgls, lambda = "ML")
stem.rich.std <- lm(global.rich ~ stem.age, data = family.data.clean)


stem.age.vasc <- max(branching.times(fulltree.vasc)) - findMRCA(fulltree.vasc, c(fulltree.vasc$tip.label[621:622]), type = "height")

r.vasc.stem.r0 <- bd.ms(time = stem.age.vasc, n = sum(family.data$global.rich), crown = FALSE, epsilon = 0)
r.vasc.stem.r05 <- bd.ms(time = stem.age.vasc, n = sum(family.data$global.rich), crown = FALSE, epsilon = 0.5)
r.vasc.stem.r0.9 <- bd.ms(time = stem.age.vasc, n = sum(family.data$global.rich), crown = FALSE, epsilon = 0.9)


### Calculating the the 95% confidence intervals of species richness based on global r, epsilon and time for vascular plants.

## Stem age
### Generating table for plotting

limits.vasc.stem <- data.frame(
    time = seq(1, 300, by = 1),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0, epsilon = 0)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r05, epsilon = 0.5)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0.9, epsilon = 0.9))
    )
names(limits.vasc.stem) <- c("age", "lb.0", "ub.0", "lb.05", "ub.05", "lb.09", "ub.09")



## Calculating Pagel's lambda for each epsilon value
phylosig.r0 <- pgls(r.e0.stem ~ 1, data = data.pgls, lambda = "ML")
phylosig.r05 <- pgls(r.e05.stem ~ 1, data = data.pgls, lambda = "ML")
phylosig.r09 <- pgls(r.e09.stem ~ 1, data = data.pgls, lambda = "ML")

## Fitting PGLS
mod.r0 <- caper::pgls(r.e0.stem ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
mod.r05 <- caper::pgls(r.e05.stem ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r05)$param[2], NULL))
mod.r09 <- caper::pgls(r.e09.stem ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r09)$param[2], NULL))

## Fitting standard linear models
lm.r0 <- lm(r.e0.stem ~ shannon, data = family.data.clean)
lm.r05 <- lm(r.e05.stem ~ shannon, data = family.data.clean)
lm.r09 <- lm(r.e09.stem ~ shannon, data = family.data.clean)

## phylANOVA
data.aov <- family.data.clean[-which(is.na(match(family.data.clean$family, tree.pruned$tip.label))), ]
phyaov.r0 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0.stem, data.aov$family))
phyaov.r05 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e05.stem, data.aov$family))
phyaov.r09 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09.stem, data.aov$family))

## Parametric ANOVA
aov.r0 <- aov(r.e0.stem ~ type.60, data = family.data.clean)
aov.r05 <- aov(r.e05.stem ~ type.60, data = family.data.clean)
aov.r09 <- aov(r.e09.stem ~ type.60, data = family.data.clean)

## Appending results to original data frame for plotting
family.data.clean$lm.r0 <- fitted(lm.r0)
family.data.clean$lm.r05 <- fitted(lm.r05)
family.data.clean$lm.r09 <- fitted(lm.r09)

family.data.clean$pgls.r0 <- fitted(mod.r0)[,1][match(family.data.clean$family, names(fitted(mod.r0)[,1]))]
family.data.clean$pgls.r05 <- fitted(mod.r05)[,1][match(family.data.clean$family, names(fitted(mod.r05)[,1]))]
family.data.clean$pgls.r09 <- fitted(mod.r09)[,1][match(family.data.clean$family, names(fitted(mod.r09)[,1]))]


## Age and Richness as a function of Shannon index 

pgls.stem.age.sh <- caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = "ML")
pgls.rich.sh <- caper::pgls(global.rich ~ shannon, data = data.pgls, lambda = "ML")

lm.stem.age.sh <- lm(stem.age ~ shannon, data = family.data.clean)
lm.rich.sh <- lm(global.rich ~ shannon, data = family.data.clean)

save.image("./output/suppdata_species_noremarks.RData")






######################
##Analysis per species using records with remarks 'probably XX'
######################

family.data.rem <- read.csv("./data/family_data_full_incl_remarks.csv", stringsAsFactors = FALSE)

## Filtering for families with more than 5% or 8 species

family.data.clean.rem <- family.data.rem[family.data.rem$samp.sp >= 8 | family.data.rem$perc.rich >= 0.05,]

## Identifying the families used in the study for later plotting
family.data.rem$study <- "no"
family.data.rem$study[family.data.rem$samp.sp >= 8 | family.data.rem$perc.rich >= 0.05] <- "yes"

fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Removing Orchidaceae and Ericaceae from analyses
family.data.clean.rem <- family.data.clean.rem[-match(c("Orchidaceae", "Ericaceae"), family.data.clean.rem$family),]

## Dropping families without reliable mycorrhizal data
tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.clean.rem$family))])

## Calculating diversification rates using the method-of-moments (Magallón & Sanderson, 2001)
family.data.clean.rem$r.e0.stem <- bd.ms(time = family.data.clean.rem$stem.age, n = family.data.clean.rem$global.rich, crown = FALSE, epsilon = 0)
family.data.clean.rem$r.e05.stem <- bd.ms(time = family.data.clean.rem$stem.age, n = family.data.clean.rem$global.rich, crown = FALSE, epsilon = 0.5)
family.data.clean.rem$r.e09.stem <- bd.ms(time = family.data.clean.rem$stem.age, n = family.data.clean.rem$global.rich, crown = FALSE, epsilon = 0.9)

family.data.clean.rem$shannon <- vegan::diversity(family.data.clean.rem[, c(4, 5, 6, 8, 9)], index = "shannon")

## Preparing data for PGLS
data.pgls.rem <- comparative.data(tree.pruned, family.data.clean.rem, names.col = "family")

## Testing for correlation between Stem age and Species Richness (both considering or not the phylogenetic structure)

stem.rich.phylo.rem <- pgls(global.rich ~ stem.age, data = data.pgls.rem, lambda = "ML")

stem.rich.std.rem <- lm(global.rich ~ stem.age, data = family.data.clean.rem)


r.vasc.stem.r0.rem <- bd.ms(time = stem.age.vasc, n = sum(family.data$global.rich), crown = FALSE, epsilon = 0)
r.vasc.stem.r05.rem <- bd.ms(time = stem.age.vasc, n = sum(family.data$global.rich), crown = FALSE, epsilon = 0.5)
r.vasc.stem.r0.9.rem <- bd.ms(time = stem.age.vasc, n = sum(family.data$global.rich), crown = FALSE, epsilon = 0.9)


### Calculating the the 95% confidence intervals of species richness based on global r, epsilon and time for vascular plants.

## Stem age
### Generating table for plotting

limits.vasc.stem.rem <- data.frame(
    time = seq(1, 300, by = 1),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0.rem, epsilon = 0)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r05.rem, epsilon = 0.5)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0.9.rem, epsilon = 0.9))
    )
names(limits.vasc.stem.rem) <- c("age", "lb.0", "ub.0", "lb.05", "ub.05", "lb.09", "ub.09")



## Calculating Pagel's lambda for each epsilon value
phylosig.r0.rem <- pgls(r.e0.stem ~ 1, data = data.pgls.rem, lambda = "ML")
phylosig.r05.rem <- pgls(r.e05.stem ~ 1, data = data.pgls.rem, lambda = "ML")
phylosig.r09.rem <- pgls(r.e09.stem ~ 1, data = data.pgls.rem, lambda = "ML")

## Fitting PGLS
mod.r0.rem <- caper::pgls(r.e0.stem ~ shannon, data = data.pgls.rem, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
mod.r05.rem <- caper::pgls(r.e05.stem ~ shannon, data = data.pgls.rem, lambda = setNames(summary(phylosig.r05)$param[2], NULL))
mod.r09.rem <- caper::pgls(r.e09.stem ~ shannon, data = data.pgls.rem, lambda = setNames(summary(phylosig.r09)$param[2], NULL))

## Fitting standard linear models
lm.r0.rem <- lm(r.e0.stem ~ shannon, data = family.data.clean.rem)
lm.r05.rem <- lm(r.e05.stem ~ shannon, data = family.data.clean.rem)
lm.r09.rem <- lm(r.e09.stem ~ shannon, data = family.data.clean.rem)

## phylANOVA
data.aov.rem <- family.data.clean.rem[-which(is.na(match(family.data.clean.rem$family, tree.pruned$tip.label))), ]
phyaov.r0.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0.stem, data.aov$family))
phyaov.r05.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e05.stem, data.aov$family))
phyaov.r09.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09.stem, data.aov$family))

## Parametric ANOVA
aov.r0.rem <- aov(r.e0.stem ~ type.60, data = family.data.clean.rem)
aov.r05.rem <- aov(r.e05.stem ~ type.60, data = family.data.clean.rem)
aov.r09.rem <- aov(r.e09.stem ~ type.60, data = family.data.clean.rem)

## Appending results to original data frame for plotting
family.data.clean.rem$lm.r0[which(!is.na(family.data.clean.rem$r.e0.stem))] <- fitted(lm.r0.rem)
family.data.clean.rem$lm.r05[which(!is.na(family.data.clean.rem$r.e05.stem))] <- fitted(lm.r05.rem)
family.data.clean.rem$lm.r09[which(!is.na(family.data.clean.rem$r.e09.stem))] <- fitted(lm.r09.rem)

family.data.clean.rem$pgls.r0 <- fitted(mod.r0)[,1][match(family.data.clean.rem$family, names(fitted(mod.r0.rem)[,1]))]
family.data.clean.rem$pgls.r05 <- fitted(mod.r05)[,1][match(family.data.clean.rem$family, names(fitted(mod.r05.rem)[,1]))]
family.data.clean.rem$pgls.r09 <- fitted(mod.r09)[,1][match(family.data.clean.rem$family, names(fitted(mod.r09.rem)[,1]))]


## Age and Richness as a function of Shannon index 

pgls.stem.age.sh.rem <- caper::pgls(stem.age ~ shannon, data = data.pgls.rem, lambda = "ML")
pgls.rich.sh.rem <- caper::pgls(global.rich ~ shannon, data = data.pgls.rem, lambda = "ML")

lm.stem.age.sh.rem <- lm(stem.age ~ shannon, data = family.data.clean.rem)
lm.rich.sh.rem <- lm(global.rich ~ shannon, data = family.data.clean.rem)

save.image("./output/suppdata_species_withremarks.RData")
