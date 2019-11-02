library("tidyverse")
library("ape")
library("caper")
library("geiger")
library("phytools")
library("foreach")
library("doMC")
library("plyr")
load("./output/sampled_datasets_20perc.RData")

## Importing data
family.data.gen <- read.csv("./data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"
family.data.gen <- family.data.gen[family.data.gen$MIX.raw.perc != 1, ]

age.data <- read.csv("./data/data_all_families.csv", sep = ";")

family.data.gen$stem.age <- age.data$stem.age[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, sampled.datasets[[1]]$family))])

data.pgls <- comparative.data(tree.pruned, sampled.datasets[[1]], names.col = "family")

phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")

foo <- function(x){
    data.pgls <- comparative.data(tree.pruned, sampled.datasets[[x]], names.col = "family")
    ## PGLS
    mod.r0 <- caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
    mod.r09 <- caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
    ## LM
    lm.r0 <- lm(r.e0 ~ shannon, data = data.pgls$data)
    lm.r09 <- lm(r.e09 ~ shannon, data = data.pgls$data)
    ## phyANOVA
    data.aov <- data.pgls$data
    data.aov <- data.aov[data.aov$type.50 != "ER", ]
    data.aov <- data.aov[data.aov$type.60 != "ER", ]
    data.aov <- data.aov[data.aov$type.80 != "ER", ]
    data.aov <- data.aov[data.aov$type.100 != "ER", ]
    phyaov.r0.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.50, rownames(data.aov)), y = setNames(data.aov$r.e0, rownames(data.aov)))
    phyaov.r09.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.50, rownames(data.aov)), y = setNames(data.aov$r.e09, rownames(data.aov)))
    phyaov.r0.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.60, rownames(data.aov)), y = setNames(data.aov$r.e0, rownames(data.aov)))
    phyaov.r09.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.60, rownames(data.aov)), y = setNames(data.aov$r.e09, rownames(data.aov)))
    phyaov.r0.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.80, rownames(data.aov)), y = setNames(data.aov$r.e0, rownames(data.aov)))
    phyaov.r09.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.80, rownames(data.aov)), y = setNames(data.aov$r.e09, rownames(data.aov)))
    phyaov.r0.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.100, rownames(data.aov)), y = setNames(data.aov$r.e0, rownames(data.aov)))
    phyaov.r09.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, rownames(data.aov))))]), x = setNames(data.aov$type.100, rownames(data.aov)), y = setNames(data.aov$r.e09, rownames(data.aov)))
    ## Standard ANOVA
    aov.r0.50 <- aov(r.e0 ~ type.50, data = data.aov)
    aov.r09.50 <- aov(r.e09 ~ type.50, data = data.aov)
    aov.r0.60 <- aov(r.e0 ~ type.60, data = data.aov)
    aov.r09.60 <- aov(r.e09 ~ type.60, data = data.aov)
    aov.r0.80 <- aov(r.e0 ~ type.80, data = data.aov)
    aov.r09.80 <- aov(r.e09 ~ type.80, data = data.aov)
    aov.r0.100 <- aov(r.e0 ~ type.100, data = data.aov)
    aov.r09.100 <- aov(r.e09 ~ type.100, data = data.aov)
    res.reg <- data.frame(
        epsilon = c(0, 0.9, 0, 0.9),
        model = c("PGLS", "PGLS", "LM", "LM"),
        pvalue = c(summary(mod.r0)$coefficients[2, 4],
                   summary(mod.r09)$coefficients[2, 4],
                   summary(lm.r0)$coefficients[2, 4],
                   summary(lm.r09)$coefficients[2, 4]),
        R2 = c(summary(mod.r0)$r.squared,
               summary(mod.r09)$r.squared,
               summary(lm.r0)$r.squared,
               summary(lm.r09)$r.squared)
        )
    res.aov <- data.frame(
        threshold = c(50, 60, 80, 100, 50, 60, 80, 100),
        model = rep(c("phy", "std"), each = 4),
        pvalue.r0 = c(phyaov.r0.50$Pf,
                      phyaov.r0.60$Pf,
                      phyaov.r0.80$Pf,
                      phyaov.r0.100$Pf,
                      summary(aov.r0.50)[[1]][1, 5],
                      summary(aov.r0.60)[[1]][1, 5],
                      summary(aov.r0.80)[[1]][1, 5],
                      summary(aov.r0.100)[[1]][1, 5]
                      ),
        pvalue.r09 = c(phyaov.r09.50$Pf,
                      phyaov.r09.60$Pf,
                      phyaov.r09.80$Pf,
                      phyaov.r09.100$Pf,
                      summary(aov.r09.50)[[1]][1, 5],
                      summary(aov.r09.60)[[1]][1, 5],
                      summary(aov.r09.80)[[1]][1, 5],
                      summary(aov.r09.100)[[1]][1, 5]
                      )
    )
    return(list(res.reg, res.aov))
}

results.sampled <- lapply(1:50, foo)

pvalue.reg <- ldply(results.sampled, function(x){x[[1]]})
pvalue.aov <- ldply(results.sampled, function(x){x[[2]]})

save(pvalue.reg, pvalue.aov, file = "./output/results_20perc.RData")
