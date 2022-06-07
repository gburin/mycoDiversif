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

#######################
### Analysis per genera
#######################

age.data <- read.csv("../data/data_all_families.csv", sep = ";")
fulltree <- read.tree("../data/fam_tree_family_full.tre")
fulltree$node.label <- NULL

here::i_am("R/side_analysis_zanne_raboskybenson_epsilon.R")

nrep <- sort(sample(1:10000, 1000))

side.analysis <- function(x, age, fulltree, perc){
    print(paste0("Replica ", x, " of ", length(list.files("../output/simulated_datasets/"))))
    family.data.gen <- read.csv(paste0("../output/simulated_datasets/random_data_", sprintf("%05d", x), ".csv"), stringsAsFactors = FALSE)
    family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
    family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

    ## Removing families with unknown mycorrhizal type
    family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

    family.data.gen$stem.age <- age$stem.age[match(family.data.gen$family, age$familia)]
    family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
    family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)
    family.sub <- sort(sample(1:nrow(family.data.gen), ceiling((perc * nrow(family.data.gen))/100)))
    family.data.gen$r.e0[family.sub] <- sapply(family.sub, function(x){bd.ms(time = family.data.gen$stem.age[x], n = family.data.gen$rich[x], crown = FALSE, epsilon = runif(1, 0, 0.9999))})
    family.data.gen$r.e09[family.sub] <- sapply(family.sub, function(x){bd.ms(time = family.data.gen$stem.age[x], n = family.data.gen$rich[x], crown = FALSE, epsilon = runif(1, 0, 0.9999))})

    family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

    family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae", "Diapensiaceae"), family.data.gen$family),]
    family.data.gen <- family.data.gen[!is.na(family.data.gen$r.e0),]

    tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
    ## Sampling families which will be set to having r = 0
    data.pgls <- comparative.data(tree.pruned, family.data.gen, names.col = "family")

    ## phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.rich <- pgls(rich ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.age <- pgls(stem.age ~ 1, data = data.pgls, lambda = "ML")

    ## Fitting PGLS excluding families with != 100% MIX
    mod.r0 <- tryCatch(tryCatch(caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = "ML"), error = caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = "ML", bounds = list(lambda = c(1e-6, 0.8)))), error = function(x){NA})
    mod.r09 <- tryCatch(tryCatch(caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = "ML"), error = caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = "ML", bounds = list(lambda = c(1e-6, 0.8)))), error = function(x){NA})

    ## Fitting standard linear models excluding families with != 100% MIX
    lm.r0 <- lm(r.e0 ~ shannon, data = family.data.gen)
    lm.r09 <- lm(r.e09 ~ shannon, data = family.data.gen)

    ## Age vs rich
    pgls.age.sh <- tryCatch(tryCatch(caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = "ML"), error = caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = "ML", bounds = list(lambda = c(1e-6, 0.8)))), error = function(x){NA})
    pgls.rich.sh <- tryCatch(tryCatch(caper::pgls(rich ~ shannon, data = data.pgls, lambda = "ML"), error = caper::pgls(rich ~ shannon, data = data.pgls, lambda = "ML", bounds = list(lambda = c(1e-6, 0.8)))), error = function(x){NA})

    lm.age.sh <- lm(stem.age ~ shannon, data = family.data.gen)
    lm.rich.sh <- lm(rich ~ shannon, data = family.data.gen)

    ## phylANOVA excluding families with != 100% MIX
    data.aov <- family.data.gen[which(!is.na(match(family.data.gen$family, fulltree$tip.label))), ]
    #data.aov <- data.aov[-which(is.na(data.aov$r.e0)),]
    
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

    results <- data.frame(pgls.int.r0 = tryCatch(coef(mod.r0)[1], error = function(x){NA}),
                          pgls.slope.r0 = tryCatch(coef(mod.r0)[2], error = function(x){NA}),
                          pgls.r2.r0 = tryCatch(summary(mod.r0)$r.squared, error = function(x){NA}),
                          pgls.pvalue.r0 = tryCatch(summary(mod.r0)$coefficients[2, 4], error = function(x){NA}),
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
                          pgls.int.r09 = tryCatch(coef(mod.r09)[1], error = function(x){NA}),
                          pgls.slope.r09 = tryCatch(coef(mod.r09)[2], error = function(x){NA}),
                          pgls.r2.r09 = tryCatch(summary(mod.r09)$r.squared, error = function(x){NA}),
                          pgls.pvalue.r09 = tryCatch(summary(mod.r09)$coefficients[2, 4], error = function(x){NA}),
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
                          pgls.int.age.sh = tryCatch(coef(pgls.age.sh)[1], error = function(x){NA}),
                          pgls.slope.age.sh = tryCatch(coef(pgls.age.sh)[2], error = function(x){NA}),
                          pgls.r2.age.sh = tryCatch(summary(pgls.age.sh)$r.squared, error = function(x){NA}),
                          pgls.pvalue.age.sh = tryCatch(summary(pgls.age.sh)$coefficients[2, 4], error = function(x){NA}),
                          pgls.int.rich.sh = tryCatch(coef(pgls.rich.sh)[1], error = function(x){NA}),
                          pgls.slope.rich.sh = tryCatch(coef(pgls.rich.sh)[2], error = function(x){NA}),
                          pgls.r2.rich.sh = tryCatch(summary(pgls.rich.sh)$r.squared, error = function(x){NA}),
                          pgls.pvalue.rich.sh = tryCatch(summary(pgls.rich.sh)$coefficients[2, 4], error = function(x){NA}),
                          lm.int.age.sh = coef(lm.age.sh)[1],
                          lm.slope.age.sh = coef(lm.age.sh)[2],
                          lm.r2.age.sh = summary(lm.age.sh)$r.squared,
                          lm.pvalue.age.sh = summary(lm.age.sh)$coefficients[2, 4],
                          lm.int.rich.sh = coef(lm.rich.sh)[1],
                          lm.slope.rich.sh = coef(lm.rich.sh)[2],
                          lm.r2.rich.sh = summary(lm.rich.sh)$r.squared,
                          lm.pvalue.rich.sh = summary(lm.rich.sh)$coefficients[2, 4], row.names = NULL)
    
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

registerDoMC(56)

side.results.5 <- llply(nrep, side.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE, perc = 5)

#write.table(side.results.5$results, file = "../output/fit_data_random_datasets_side_rabosky_benson_epsilon_5perc.csv", sep = ",", quote = FALSE, row.names = FALSE)

save(side.results.5, file = "../output/side_analysis_rabosky_benson_epsilon_5perc.RData")

side.results.20 <- ldply(nrep, side.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE, perc = 20)

#write.table(side.results$results, file = "../output/fit_data_random_datasets_side_rabosky_benson_epsilon_20perc.csv", sep = ",", quote = FALSE, row.names = FALSE)

save(side.results.20, file = "../output/side_analysis_rabosky_benson_epsilon_20perc.RData")

side.results.50 <- llply(nrep, side.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE, perc = 50)

#write.table(side.results$results, file = "../output/fit_data_random_datasets_side_rabosky_benson.csv", sep = ",", quote = FALSE, row.names = FALSE)

save(side.results.50, file = "../output/side_analysis_rabosky_benson_epsilon_50perc.RData")
