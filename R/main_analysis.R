###
library("caper")
library("tidyverse")
library("cowplot")
library("RColorBrewer")
library("geiger")
library("phytools")
library("ggrepel")

#######################
### Analysis per genera
#######################

family.data.gen <- read.csv("./data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"
family.data.gen$perc.sp.unk <- family.data.gen$UNK/family.data.gen$rich

age.data <- read.csv("./data/data_all_families.csv", sep = ";")

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

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

min(((family.data.gen$rich - family.data.gen$UNK)/family.data.gen$rich))

## We will not remove any family due to poor sampling since all families have either at least 58.3% of total richness or a minimum of 8 species sampled

family.data.gen.valid <- family.data.gen
family.data.gen.valid$shannon.valid <- vegan::diversity(family.data.gen.valid[, c(8, 9, 10, 12, 13)])

tree.pruned.valid <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen.valid$family))])
data.pgls.valid <- comparative.data(tree.pruned.valid, family.data.gen.valid, names.col = "family")

phylosig.valid.r0 <- pgls(r.e0 ~ 1, data = data.pgls.valid, lambda = "ML")
phylosig.valid.r05 <- pgls(r.e05 ~ 1, data = data.pgls.valid, lambda = "ML")
phylosig.valid.r09 <- pgls(r.e09 ~ 1, data = data.pgls.valid, lambda = "ML")

## Fitting PGLS
mod.valid.r0 <- caper::pgls(r.e0 ~ shannon.valid, data = data.pgls.valid, lambda = setNames(summary(phylosig.valid.r0)$param[2], NULL))
mod.valid.r05 <- caper::pgls(r.e05 ~ shannon.valid, data = data.pgls.valid, lambda = setNames(summary(phylosig.valid.r05)$param[2], NULL))
mod.valid.r09 <- caper::pgls(r.e09 ~ shannon.valid, data = data.pgls.valid, lambda = setNames(summary(phylosig.valid.r09)$param[2], NULL))

## Fitting standard linear models
lm.valid.r0 <- lm(r.e0 ~ shannon.valid, data = family.data.gen.valid)
lm.valid.r05 <- lm(r.e05 ~ shannon.valid, data = family.data.gen.valid)
lm.valid.r09 <- lm(r.e09 ~ shannon.valid, data = family.data.gen.valid)

## phylANOVA
data.aov.valid <- family.data.gen.valid[-which(is.na(match(family.data.gen.valid$family, fulltree$tip.label))), ]
data.aov.valid <- data.aov.valid[-which(data.aov.valid$type.50 == "ER"), ]
data.aov.valid <- data.aov.valid[-which(is.na(data.aov.valid$r.e0)),]

phyaov.valid.r0 <- phylANOVA(drop.tip(tree.pruned.valid, tip = tree.pruned.valid$tip.label[which(is.na(match(tree.pruned.valid$tip.label, data.aov.valid$family)))]), x = setNames(data.aov.valid$type.60.valid, data.aov.valid$family), y = setNames(data.aov.valid$r.e0, data.aov.valid$family))
phyaov.valid.r05 <- phylANOVA(drop.tip(tree.pruned.valid, tip = tree.pruned.valid$tip.label[which(is.na(match(tree.pruned.valid$tip.label, data.aov.valid$family)))]), x = setNames(data.aov.valid$type.60.valid, data.aov.valid$family), y = setNames(data.aov.valid$r.e05, data.aov.valid$family))
phyaov.valid.r09 <- phylANOVA(drop.tip(tree.pruned.valid, tip = tree.pruned.valid$tip.label[which(is.na(match(tree.pruned.valid$tip.label, data.aov.valid$family)))]), x = setNames(data.aov.valid$type.60.valid, data.aov.valid$family), y = setNames(data.aov.valid$r.e09, data.aov.valid$family))


## Age vs rich
pgls.valid.age.sh <- caper::pgls(stem.age ~ shannon.valid, data = data.pgls.valid, lambda = "ML")
pgls.valid.rich.sh <- caper::pgls(rich ~ shannon.valid, data = data.pgls.valid, lambda = "ML")

lm.valid.age.sh <- lm(stem.age ~ shannon.valid, data = family.data.gen.valid)
lm.valid.rich.sh <- lm(rich ~ shannon.valid, data = family.data.gen.valid)


## Appending results to original data frame for plotting
family.data.gen.valid$lm.valid.r0[which(!is.na(family.data.gen.valid$r.e0))] <- fitted(lm.valid.r0)
family.data.gen.valid$lm.valid.r05[which(!is.na(family.data.gen.valid$r.e0))] <- fitted(lm.valid.r05)
family.data.gen.valid$lm.valid.r09[which(!is.na(family.data.gen.valid$r.e0))] <- fitted(lm.valid.r09)

family.data.gen.valid$pgls.r0 <- fitted(mod.valid.r0)[,1][match(family.data.gen.valid$family, names(fitted(mod.valid.r0)[,1]))]
family.data.gen.valid$pgls.r05 <- fitted(mod.valid.r05)[,1][match(family.data.gen.valid$family, names(fitted(mod.valid.r05)[,1]))]
family.data.gen.valid$pgls.r09 <- fitted(mod.valid.r09)[,1][match(family.data.gen.valid$family, names(fitted(mod.valid.r09)[,1]))]


