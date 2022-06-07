###
library("caper")
library("tidyverse")
library("cowplot")
library("RColorBrewer")
library("geiger")
library("phytools")
#library("ggrepel")
library("foreach")
library("doMC")
library("reshape2")
library("plyr")

here::i_am("R/main_analysis_zanne.R")

#######################
### Analysis per genera
#######################

age.data <- read.csv(here::here("data/data_all_families.csv"), sep = ";")
fulltree <- read.tree(here::here("data/fam_tree_family_full.tre"))
fulltree$node.label <- NULL

nrep = length(list.files(here::here("output/simulated_datasets/")))

main.analysis <- function(x, age, fulltree){
    print(paste0("Replica ", x, " of ", length(list.files(here::here("output/simulated_datasets/")))))
    family.data.gen <- read.csv(here::here(paste0("output/simulated_datasets/random_data_", sprintf("%05d", x), ".csv")), stringsAsFactors = FALSE)
    family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
    family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

    ## Removing families with unknown mycorrhizal type
    family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

    family.data.gen$stem.age <- age$stem.age[match(family.data.gen$family, age$familia)]
    family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
    family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

    family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

    family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae", "Diapensiaceae"), family.data.gen$family),]

    tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
    data.pgls <- comparative.data(tree.pruned, family.data.gen, names.col = "family")

    ## phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.rich <- pgls(rich ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.age <- pgls(stem.age ~ 1, data = data.pgls, lambda = "ML")

    ## Fitting PGLS excluding families with != 100% MIX
    mod.r0 <- caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = "ML")
    mod.r09 <- caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = "ML")

    ## Fitting standard linear models excluding families with != 100% MIX
    lm.r0 <- lm(r.e0 ~ shannon, data = family.data.gen)
    lm.r09 <- lm(r.e09 ~ shannon, data = family.data.gen)

    ## phylANOVA excluding families with != 100% MIX
    data.aov <- family.data.gen[-which(is.na(match(family.data.gen$family, fulltree$tip.label))), ]
    data.aov <- data.aov[-which(is.na(data.aov$r.e0)),]
    
    ## Thresholds
### 50%
    phyaov.r0.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
    phyaov.r09.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

### 60%
    phyaov.r0.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
    phyaov.r09.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

### 80%
    phyaov.r0.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
    phyaov.r09.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

### 100%
    phyaov.r0.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
    phyaov.r09.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

    ## ANOVA excluding families with != 100% MIX
    ## Thresholds
### 50%
    aov.r0.50 <- aov(r.e0 ~ type.50, data = data.aov)
    aov.r09.50 <- aov(r.e09 ~ type.50, data = data.aov)

### 60%
    aov.r0.60 <- aov(r.e0 ~ type.60, data = data.aov)
    aov.r09.60 <- aov(r.e09 ~ type.60, data = data.aov)

### 80%
    aov.r0.80 <- aov(r.e0 ~ type.80, data = data.aov)
    aov.r09.80 <- aov(r.e09 ~ type.80, data = data.aov)

### 100%
    aov.r0.100 <- aov(r.e0 ~ type.100, data = data.aov)
    aov.r09.100 <- aov(r.e09 ~ type.100, data = data.aov)

    ## Age vs rich
    pgls.age.sh <- caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = "ML")
    pgls.rich.sh <- caper::pgls(rich ~ shannon, data = data.pgls, lambda = "ML")

    lm.age.sh <- lm(stem.age ~ shannon, data = family.data.gen)
    lm.rich.sh <- lm(rich ~ shannon, data = family.data.gen)

    results <- data.frame(pgls.int.r0 = coef(mod.r0)[1],
                          pgls.slope.r0 = coef(mod.r0)[2],
                          pgls.r2.r0 = summary(mod.r0)$r.squared,
                          pgls.pvalue.r0 = summary(mod.r0)$coefficients[2, 4],
                          phyaov.pvalue.r0.50 = phyaov.r0.50$Pf,
                          phyaov.pvalue.r0.60 = phyaov.r0.60$Pf,
                          phyaov.pvalue.r0.80 = phyaov.r0.80$Pf,
                          phyaov.pvalue.r0.100 = phyaov.r0.100$Pf,
                          lm.int.r0 = coef(lm.r0)[1],
                          lm.slope.r0 = coef(lm.r0)[2],
                          lm.r2.r0 = summary(lm.r0)$r.squared,
                          lm.pvalue.r0 = summary(lm.r0)$coefficients[2, 4],
                          aov.pvalue.r0.50 = summary(aov.r0.50)[[1]][1, 5],
                          aov.pvalue.r0.60 = summary(aov.r0.60)[[1]][1, 5],
                          aov.pvalue.r0.80 = summary(aov.r0.80)[[1]][1, 5],
                          aov.pvalue.r0.100 = summary(aov.r0.100)[[1]][1, 5],
                          pgls.int.r09 = coef(mod.r09)[1],
                          pgls.slope.r09 = coef(mod.r09)[2],
                          pgls.r2.r09 = summary(mod.r09)$r.squared,
                          pgls.pvalue.r09 = summary(mod.r09)$coefficients[2, 4],
                          phyaov.pvalue.r09.50 = phyaov.r09.50$Pf,
                          phyaov.pvalue.r09.60 = phyaov.r09.60$Pf,
                          phyaov.pvalue.r09.80 = phyaov.r09.80$Pf,
                          phyaov.pvalue.r09.100 = phyaov.r09.100$Pf,
                          lm.int.r09 = coef(lm.r09)[1],
                          lm.slope.r09 = coef(lm.r09)[2],
                          lm.r2.r09 = summary(lm.r09)$r.squared,
                          lm.pvalue.r09 = summary(lm.r09)$coefficients[2, 4],
                          aov.pvalue.r09.50 = summary(aov.r09.50)[[1]][1, 5],
                          aov.pvalue.r09.60 = summary(aov.r09.60)[[1]][1, 5],
                          aov.pvalue.r09.80 = summary(aov.r09.80)[[1]][1, 5],
                          aov.pvalue.r09.100 = summary(aov.r09.100)[[1]][1, 5],
                          pgls.int.age.sh = coef(pgls.age.sh)[1],
                          pgls.slope.age.sh = coef(pgls.age.sh)[2],
                          pgls.r2.age.sh = summary(pgls.age.sh)$r.squared,
                          pgls.pvalue.age.sh = summary(pgls.age.sh)$coefficients[2, 4],
                          pgls.int.rich.sh = coef(pgls.rich.sh)[1],
                          pgls.slope.rich.sh = coef(pgls.rich.sh)[2],
                          pgls.r2.rich.sh = summary(pgls.rich.sh)$r.squared,
                          pgls.pvalue.rich.sh = summary(pgls.rich.sh)$coefficients[2, 4],
                          lm.int.age.sh = coef(lm.age.sh)[1],
                          lm.slope.age.sh = coef(lm.age.sh)[2],
                          lm.r2.age.sh = summary(lm.age.sh)$r.squared,
                          lm.pvalue.age.sh = summary(lm.age.sh)$coefficients[2, 4],
                          lm.int.rich.sh = coef(lm.rich.sh)[1],
                          lm.slope.rich.sh = coef(lm.rich.sh)[2],
                          lm.r2.rich.sh = summary(lm.rich.sh)$r.squared,
                          lm.pvalue.rich.sh = summary(lm.rich.sh)$coefficients[2, 4], row.names = NULL)
    
    ## return(results)
    return(list(results,
                phyaov.r0.50 = phyaov.r0.50$Pt,
                phyaov.r09.50 = phyaov.r09.50$Pt,
                phyaov.r0.60 = phyaov.r0.60$Pt,
                phyaov.r09.60 = phyaov.r09.60$Pt,
                phyaov.r0.80 = phyaov.r0.80$Pt,
                phyaov.r09.80 = phyaov.r09.80$Pt,
                phyaov.r0.100 = phyaov.r0.100$Pt,
                phyaov.r09.100 = phyaov.r09.100$Pt,
                aov.r0.50 = TukeyHSD(aov.r0.50)$type.50,
                aov.r09.50 = TukeyHSD(aov.r09.50)$type.50,
                aov.r0.60 = TukeyHSD(aov.r0.60)$type.60,
                aov.r09.60 = TukeyHSD(aov.r09.60)$type.60,
                aov.r0.80 = TukeyHSD(aov.r0.80)$type.80,
                aov.r09.80 = TukeyHSD(aov.r09.80)$type.80,
                aov.r0.100 = TukeyHSD(aov.r0.100)$type.100,
                aov.r09.100 = TukeyHSD(aov.r09.100)$type.100
                )
           )
                
    
}

registerDoMC(50)

## main.results <- ldply(1:nrep, main.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE)

main.results <- llply(1:nrep, main.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE)

## write.table(main.results, file = "./output/main_results_zanne_2021.csv", sep = ",", quote = FALSE, row.names = FALSE)

saveRDS(main.results, file = here::here("output/main_results_zanne_2021.RDS"))
