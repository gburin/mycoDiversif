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
load("../output/results_20perc.RData")

age.data <- read.csv("../data/data_all_families.csv", sep = ";")
family.data.gen <- read.csv("../data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

## Importing tree
fulltree <- read.tree("../data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.raw.perc == 1),]

family.data.gen$stem.age <- age.data$stem.age[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e05 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.5)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae", "Diapensiaceae"), family.data.gen$family),]

tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
data.pgls <- comparative.data(tree.pruned, family.data.gen, names.col = "family")





######################
##Analysis per species
######################
rm(list = ls())
family.data <- read.csv("../data/family_data_full.csv", stringsAsFactors = FALSE)

## Filtering for families with more than 5% or 8 species

family.data.clean <- family.data[family.data$samp.sp >= 8 | family.data$perc.rich >= 0.05,]

## Identifying the families used in the study for later plotting
family.data$study <- "no"
family.data$study[family.data$samp.sp >= 8 | family.data$perc.rich >= 0.05] <- "yes"

fulltree <- read.tree("../data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

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
data.aov <- data.aov[data.aov$type.50 != "ER", ]
data.aov <- data.aov[data.aov$type.60 != "ER", ]
data.aov <- data.aov[data.aov$type.80 != "ER", ]
data.aov <- data.aov[data.aov$type.100 != "ER", ]


### Threshold 50%
phyaov.r0.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$r.e0.stem, data.aov$family))
phyaov.r09.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$r.e09.stem, data.aov$family))

### Threshold 60%
phyaov.r0.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0.stem, data.aov$family))
phyaov.r09.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09.stem, data.aov$family))

### Threshold 80%
phyaov.r0.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$r.e0.stem, data.aov$family))
phyaov.r09.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$r.e09.stem, data.aov$family))

### Threshold 100%
phyaov.r0.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$r.e0.stem, data.aov$family))
phyaov.r09.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov$family))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$r.e09.stem, data.aov$family))

## Parametric ANOVA

### Threshold 50%
aov.r0.50 <- aov(r.e0.stem ~ type.50, data = family.data.clean)
aov.r09.50 <- aov(r.e09.stem ~ type.50, data = family.data.clean)

### Threshold 60%
aov.r0.60 <- aov(r.e0.stem ~ type.60, data = family.data.clean)
aov.r09.60 <- aov(r.e09.stem ~ type.60, data = family.data.clean)

### Threshold 80%
aov.r0.80 <- aov(r.e0.stem ~ type.80, data = family.data.clean)
aov.r09.80 <- aov(r.e09.stem ~ type.80, data = family.data.clean)

### Threshold 100%
aov.r0.100 <- aov(r.e0.stem ~ type.100, data = family.data.clean)
aov.r09.100 <- aov(r.e09.stem ~ type.100, data = family.data.clean)

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

save.image("../output/suppdata_species_noremarks.RData")






######################
##Analysis per species using records with remarks 'probably XX'
######################

rm(list = ls())
family.data.rem <- read.csv("../data/family_data_full_incl_remarks.csv", stringsAsFactors = FALSE)

## Filtering for families with more than 5% or 8 species

family.data.clean.rem <- family.data.rem[family.data.rem$samp.sp >= 8 | family.data.rem$perc.rich >= 0.05,]

## Identifying the families used in the study for later plotting
family.data.rem$study <- "no"
family.data.rem$study[family.data.rem$samp.sp >= 8 | family.data.rem$perc.rich >= 0.05] <- "yes"

fulltree <- read.tree("../data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

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


## Calculating Pagel's lambda for each epsilon value
phylosig.r0.rem <- pgls(r.e0.stem ~ 1, data = data.pgls.rem, lambda = "ML")
phylosig.r05.rem <- pgls(r.e05.stem ~ 1, data = data.pgls.rem, lambda = "ML")
phylosig.r09.rem <- pgls(r.e09.stem ~ 1, data = data.pgls.rem, lambda = "ML")

## Fitting PGLS
mod.r0.rem <- caper::pgls(r.e0.stem ~ shannon, data = data.pgls.rem, lambda = setNames(summary(phylosig.r0.rem)$param[2], NULL))
mod.r05.rem <- caper::pgls(r.e05.stem ~ shannon, data = data.pgls.rem, lambda = setNames(summary(phylosig.r05.rem)$param[2], NULL))
mod.r09.rem <- caper::pgls(r.e09.stem ~ shannon, data = data.pgls.rem, lambda = setNames(summary(phylosig.r09.rem)$param[2], NULL))

## Fitting standard linear models
lm.r0.rem <- lm(r.e0.stem ~ shannon, data = family.data.clean.rem)
lm.r05.rem <- lm(r.e05.stem ~ shannon, data = family.data.clean.rem)
lm.r09.rem <- lm(r.e09.stem ~ shannon, data = family.data.clean.rem)



## phylANOVA
data.aov.rem <- family.data.clean.rem[-which(is.na(match(family.data.clean.rem$family, tree.pruned$tip.label))), ]
data.aov.rem <- data.aov.rem[data.aov.rem$type.50 != "ER", ]
data.aov.rem <- data.aov.rem[data.aov.rem$type.60 != "ER", ]
data.aov.rem <- data.aov.rem[data.aov.rem$type.80 != "ER", ]
data.aov.rem <- data.aov.rem[data.aov.rem$type.100 != "ER", ]

### Threshold 50%
phyaov.r0.50.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.50, data.aov.rem$family), y = setNames(data.aov.rem$r.e0.stem, data.aov.rem$family))
phyaov.r09.50.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.50, data.aov.rem$family), y = setNames(data.aov.rem$r.e09.stem, data.aov.rem$family))

### Threshold 60%
phyaov.r0.60.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.60, data.aov.rem$family), y = setNames(data.aov.rem$r.e0.stem, data.aov.rem$family))
phyaov.r09.60.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.60, data.aov.rem$family), y = setNames(data.aov.rem$r.e09.stem, data.aov.rem$family))

### Threshold 80%
phyaov.r0.80.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.80, data.aov.rem$family), y = setNames(data.aov.rem$r.e0.stem, data.aov.rem$family))
phyaov.r09.80.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.80, data.aov.rem$family), y = setNames(data.aov.rem$r.e09.stem, data.aov.rem$family))

### Threshold 100%
phyaov.r0.100.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.100, data.aov.rem$family), y = setNames(data.aov.rem$r.e0.stem, data.aov.rem$family))
phyaov.r09.100.rem <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[is.na(match(tree.pruned$tip.label, data.aov.rem$family))]), x = setNames(data.aov.rem$type.100, data.aov.rem$family), y = setNames(data.aov.rem$r.e09.stem, data.aov.rem$family))

## Parametric ANOVA

### Threshold 50%
aov.r0.50.rem <- aov(r.e0.stem ~ type.50, data = family.data.clean.rem)
aov.r09.50.rem <- aov(r.e09.stem ~ type.50, data = family.data.clean.rem)

### Threshold 60%
aov.r0.60.rem <- aov(r.e0.stem ~ type.60, data = family.data.clean.rem)
aov.r09.60.rem <- aov(r.e09.stem ~ type.60, data = family.data.clean.rem)

### Threshold 80%
aov.r0.80.rem <- aov(r.e0.stem ~ type.80, data = family.data.clean.rem)
aov.r09.80.rem <- aov(r.e09.stem ~ type.80, data = family.data.clean.rem)

### Threshold 100%
aov.r0.100.rem <- aov(r.e0.stem ~ type.100, data = family.data.clean.rem)
aov.r09.100.rem <- aov(r.e09.stem ~ type.100, data = family.data.clean.rem)


## Appending results to original data frame for plotting
family.data.clean.rem$lm.r0[which(!is.na(family.data.clean.rem$r.e0.stem))] <- fitted(lm.r0.rem)
family.data.clean.rem$lm.r05[which(!is.na(family.data.clean.rem$r.e05.stem))] <- fitted(lm.r05.rem)
family.data.clean.rem$lm.r09[which(!is.na(family.data.clean.rem$r.e09.stem))] <- fitted(lm.r09.rem)

family.data.clean.rem$pgls.r0 <- fitted(mod.r0.rem)[,1][match(family.data.clean.rem$family, names(fitted(mod.r0.rem)[,1]))]
family.data.clean.rem$pgls.r05 <- fitted(mod.r05.rem)[,1][match(family.data.clean.rem$family, names(fitted(mod.r05.rem)[,1]))]
family.data.clean.rem$pgls.r09 <- fitted(mod.r09.rem)[,1][match(family.data.clean.rem$family, names(fitted(mod.r09.rem)[,1]))]


## Age and Richness as a function of Shannon index 

pgls.stem.age.sh.rem <- caper::pgls(stem.age ~ shannon, data = data.pgls.rem, lambda = "ML")
pgls.rich.sh.rem <- caper::pgls(global.rich ~ shannon, data = data.pgls.rem, lambda = "ML")

lm.stem.age.sh.rem <- lm(stem.age ~ shannon, data = family.data.clean.rem)
lm.rich.sh.rem <- lm(global.rich ~ shannon, data = family.data.clean.rem)

save.image("../output/suppdata_species_withremarks.RData")


## Replicating analysis with Harris & Davies' phylogeny

rm(list = ls())

tree <- read.tree("../data/harris_davies_phylo.txt")
tree$tip.label <- str_to_title(tree$tip.label)

data.match <- setNames(read.csv("../data/harris_davies_data.csv", header = FALSE, stringsAsFactors = FALSE), c("node", "family"))
data.match$family <- str_to_title(data.match$family)

rates <- read.csv("../data/harris_davies_rates.csv")

rates$clade <- data.match$family[1:nrow(rates)]
rates <- rates[order(rates$clade),]

data.muj <- read.csv("../data/family_data_genus_final.csv")

data.muj <- cbind(data.muj, rates[match(data.muj$family, rates$clade), c("r.e0", "r.e09")])
names(data.muj)[36:37] <- c("new.phylo.r.e0", "new.phylo.r.e09")

data.muj <- data.muj[-match(c("Orchidaceae", "Ericaceae"), data.muj$family), ]

tree.pruned <- drop.tip(tree, tip = tree$tip.label[is.na(match(tree$tip.label, data.muj$family))])

data.pgls <- comparative.data(tree.pruned, data.muj[data.muj$MIX.raw.perc != 1, ], names.col = "family")

phylosig.new.phylo.r0 <- pgls(new.phylo.r.e0 ~ 1, data.pgls, lambda = "ML")
phylosig.new.phylo.r09 <- pgls(new.phylo.r.e09 ~ 1, data.pgls, lambda = "ML")

mod.np.r0 <- caper::pgls(new.phylo.r.e0 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.new.phylo.r0)$param[2], NULL))
mod.np.r09 <- caper::pgls(new.phylo.r.e09 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.new.phylo.r0)$param[2], NULL))

lm.np.r0 <- lm(new.phylo.r.e0 ~ shannon, data = data.muj[data.muj$MIX.raw.perc != 1,])
lm.np.r09 <- lm(new.phylo.r.e09 ~ shannon, data = data.muj[data.muj$MIX.raw.perc != 1,])

data.aov <- data.muj[-which(is.na(match(data.muj$family, tree$tip.label))), ]
data.aov <- data.aov[which(data.aov$MIX.raw.perc != 1), ]
data.aov <- data.aov[-which(data.aov$type.60 == "ER"), ]
data.aov$type.60 <- as.character(data.aov$type.60)

phyaov.np.r0.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$new.phylo.r.e0, data.aov$family))

phyaov.np.r09.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$new.phylo.r.e09, data.aov$family))

aov.np.r0.60 <- aov(new.phylo.r.e0 ~ type.60, data = data.aov)
aov.np.r09.60 <- aov(new.phylo.r.e09 ~ type.60, data = data.aov)

# TukeyHSD(aov.np.r0.60)
# TukeyHSD(aov.np.r09.60)

data.muj.plot <- data.muj[data.muj$MIX.raw.perc != 1, ]
data.muj.plot$np.lm.r0[which(!is.na(data.muj.plot$new.phylo.r.e0))] <- fitted(lm.np.r0)
data.muj.plot$np.lm.r09[which(!is.na(data.muj.plot$new.phylo.r.e09))] <- fitted(lm.np.r09)

data.muj.plot$np.plgs.r0 <- fitted(mod.np.r0)[,1][match(data.muj.plot$family, names(fitted(mod.np.r0)[,1]))]
data.muj.plot$np.plgs.r09 <- fitted(mod.np.r09)[,1][match(data.muj.plot$family, names(fitted(mod.np.r09)[,1]))]

fulltree <- read.tree("../data/fam_tree_family_full.tre")

cophylo.myco <- cophylo(tree, fulltree)

save.image(file = "../output/harris_davies_results.RData")
