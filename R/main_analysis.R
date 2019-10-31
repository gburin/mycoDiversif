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

age.data <- read.csv("./data/data_all_families.csv", sep = ";")

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.raw.perc == 1),]

family.data.gen$stem.age <- age.data$stem.age[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e05 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.5)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

## write.table(family.data.gen, file = "./data/family_data_genus_final.csv", sep = ",", quote = FALSE, row.names = FALSE)

family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae"), family.data.gen$family),]

## min(((family.data.gen$rich - family.data.gen$UNK)/family.data.gen$rich))

## We will not remove any family due to poor sampling since all families have either at least 58.3% of total richness or a minimum of 8 species sampled

tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
data.pgls <- comparative.data(tree.pruned, family.data.gen[family.data.gen$MIX.raw.perc != 1,], names.col = "family")

phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r05 <- pgls(r.e05 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")
phylosig.rich <- pgls(rich ~ 1, data = data.pgls, lambda = "ML")
phylosig.age <- pgls(stem.age ~ 1, data = data.pgls, lambda = "ML")

## Fitting PGLS excluding families with != 100% MIX
mod.r0 <- caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
mod.r05 <- caper::pgls(r.e05 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r05)$param[2], NULL))
mod.r09 <- caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r09)$param[2], NULL))

## Fitting standard linear models excluding families with != 100% MIX
lm.r0 <- lm(r.e0 ~ shannon, data = family.data.gen[family.data.gen$MIX.raw.perc != 1,])
lm.r05 <- lm(r.e05 ~ shannon, data = family.data.gen[family.data.gen$MIX.raw.perc != 1,])
lm.r09 <- lm(r.e09 ~ shannon, data = family.data.gen[family.data.gen$MIX.raw.perc != 1,])

## phylANOVA excluding families with != 100% MIX
data.aov <- family.data.gen[-which(is.na(match(family.data.gen$family, fulltree$tip.label))), ]
data.aov <- data.aov[-which(is.na(data.aov$r.e0)),]
data.aov <- data.aov[which(data.aov$MIX.raw.perc != 1), ]
data.aov <- data.aov[-which(data.aov$type.60 == "ER"), ]

phyaov.r0 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
phyaov.r05 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e05, data.aov$family))
phyaov.r09 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

## ANOVA excluding families with != 100% MIX
aov.r0 <- aov(r.e0 ~ type.60, data = data.aov)
aov.r05 <- aov(r.e05 ~ type.60, data = data.aov)
aov.r09 <- aov(r.e09 ~ type.60, data = data.aov)

TukeyHSD(aov.r0)
TukeyHSD(aov.r05)
TukeyHSD(aov.r09)

## Age vs rich
pgls.age.sh <- caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.age)$param[2], NULL))
pgls.rich.sh <- caper::pgls(rich ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.rich)$param[2], NULL))

lm.age.sh <- lm(stem.age ~ shannon, data = family.data.gen)
lm.rich.sh <- lm(rich ~ shannon, data = family.data.gen)

## Appending results to original data frame for plotting
family.data.gen <- family.data.gen[family.data.gen$MIX.raw.perc != 1, ]
family.data.gen$lm.r0[which(!is.na(family.data.gen$r.e0))] <- fitted(lm.r0)
family.data.gen$lm.r05[which(!is.na(family.data.gen$r.e0))] <- fitted(lm.r05)
family.data.gen$lm.r09[which(!is.na(family.data.gen$r.e0))] <- fitted(lm.r09)

family.data.gen$pgls.r0 <- fitted(mod.r0)[,1][match(family.data.gen$family, names(fitted(mod.r0)[,1]))]
family.data.gen$pgls.r05 <- fitted(mod.r05)[,1][match(family.data.gen$family, names(fitted(mod.r05)[,1]))]
family.data.gen$pgls.r09 <- fitted(mod.r09)[,1][match(family.data.gen$family, names(fitted(mod.r09)[,1]))]


save.image(file = "./output/main_results.RData")
