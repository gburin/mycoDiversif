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

here::i_am("R/main_analysis_igeatanen.R")

######################
### Analysis per genus
######################

bamm.data <- read.csv(here::here("output/bamm_rates.csv"), sep = ",")
fulltree <- read.tree(here::here("data/fam_tree_family_full.tre"))
fulltree$node.label <- NULL

bamm.data$netdiv <- bamm.data$lambda - bamm.data$mu
age.data <- as.data.frame(aggregate(bamm.data[, c("lambda", "mu", "netdiv")], by = list(bamm.data$family), FUN = mean))
names(age.data)[1] <- "family"

nrep = length(list.files(here::here("output/simulated_datasets/")))

main.analysis <- function(x, age, fulltree){
    print(paste0("Replica ", x))
    family.data.gen <- read.csv(here::here(paste0("output/20perc_data/", list.files("./output/20perc_data/")[x])), stringsAsFactors = FALSE)
    ## family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
    ## family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

    #family.data.gen$stem.age <- age$stem.age[match(family.data.gen$family, age$familia)]
    #family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
    #family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)
    family.data.gen$lambda <- age$lambda[match(family.data.gen$family, age$family)]
    family.data.gen$mu <- age$mu[match(family.data.gen$family, age$family)]
    family.data.gen$netdiv <- age$netdiv[match(family.data.gen$family, age$family)]
    
    family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

    tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
    data.pgls <- comparative.data(tree.pruned, family.data.gen, names.col = "family")

    phylosig.lambda <- pgls(lambda ~ 1, data = data.pgls, lambda = "ML")
    phylosig.mu <- pgls(mu ~ 1, data = data.pgls, lambda = "ML")
    phylosig.netdiv <- pgls(netdiv ~ 1, data = data.pgls, lambda = "ML")
    phylosig.rich <- pgls(rich ~ 1, data = data.pgls, lambda = "ML")
    ## phylosig.age <- pgls(stem.age ~ 1, data = data.pgls, lambda = "ML")

    ## Fitting PGLS excluding families with != 100% MIX
    mod.lambda <- caper::pgls(lambda ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.lambda)$param[2], NULL))
    mod.mu <- caper::pgls(mu ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.mu)$param[2], NULL))
    mod.netdiv <- caper::pgls(netdiv ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.netdiv)$param[2], NULL))

    ## Fitting standard linear models excluding families with != 100% MIX
    lm.lambda <- lm(lambda ~ shannon, data = family.data.gen)
    lm.mu <- lm(mu ~ shannon, data = family.data.gen)
    lm.netdiv <- lm(netdiv ~ shannon, data = family.data.gen)

    ## phylANOVA excluding families with != 100% MIX
    data.aov <- family.data.gen[-which(is.na(match(family.data.gen$family, fulltree$tip.label))), ]
    data.aov <- data.aov[-which(is.na(data.aov$lambda)),]
    
    ## Thresholds
## ### 50%
##     phyaov.lambda.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$lambda, data.aov$family))
##     phyaov.mu.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$mu, data.aov$family))
##     phyaov.netdiv.50 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.50, data.aov$family), y = setNames(data.aov$netdiv, data.aov$family))

## ### 60%
##     phyaov.lambda.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$lambda, data.aov$family))
##     phyaov.mu.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$mu, data.aov$family))
##     phyaov.netdiv.60 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$netdiv, data.aov$family))

### 80%
    phyaov.lambda.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$lambda, data.aov$family))
    phyaov.mu.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$mu, data.aov$family))
    phyaov.netdiv.80 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.80, data.aov$family), y = setNames(data.aov$netdiv, data.aov$family))

### 90%
    phyaov.lambda.90 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.90, data.aov$family), y = setNames(data.aov$lambda, data.aov$family))
    phyaov.mu.90 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.90, data.aov$family), y = setNames(data.aov$mu, data.aov$family))
    phyaov.netdiv.90 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.90, data.aov$family), y = setNames(data.aov$netdiv, data.aov$family))

### 100%
    phyaov.lambda.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$lambda, data.aov$family))
    phyaov.mu.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$mu, data.aov$family))
    phyaov.netdiv.100 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.100, data.aov$family), y = setNames(data.aov$netdiv, data.aov$family))

    ## ANOVA excluding families with != 100% MIX
    ## Thresholds
## ### 50%
##     aov.lambda.50 <- aov(lambda ~ type.50, data = data.aov)
##     aov.mu.50 <- aov(mu ~ type.50, data = data.aov)
##     aov.netdiv.50 <- aov(netdiv ~ type.50, data = data.aov)

## ### 60%
##     aov.lambda.60 <- aov(lambda ~ type.60, data = data.aov)
##     aov.mu.60 <- aov(mu ~ type.60, data = data.aov)
##     aov.netdiv.60 <- aov(netdiv ~ type.60, data = data.aov)

### 80%
    aov.lambda.80 <- aov(lambda ~ type.80, data = data.aov)
    aov.mu.80 <- aov(mu ~ type.80, data = data.aov)
    aov.netdiv.80 <- aov(netdiv ~ type.80, data = data.aov)

### 90%
    aov.lambda.90 <- aov(lambda ~ type.90, data = data.aov)
    aov.mu.90 <- aov(mu ~ type.90, data = data.aov)
    aov.netdiv.90 <- aov(netdiv ~ type.90, data = data.aov)

### 100%
    aov.lambda.100 <- aov(lambda ~ type.100, data = data.aov)
    aov.mu.100 <- aov(mu ~ type.100, data = data.aov)
    aov.netdiv.100 <- aov(netdiv ~ type.100, data = data.aov)

    ## Age vs rich
    ## pgls.age.sh <- caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.age)$param[2], NULL))
    pgls.rich.sh <- caper::pgls(rich ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.rich)$param[2], NULL))

    ## lm.age.sh <- lm(stem.age ~ shannon, data = family.data.gen)
    lm.rich.sh <- lm(rich ~ shannon, data = family.data.gen)

    results <- data.frame(pgls.int.lambda = coef(mod.lambda)[1],
                          pgls.slope.lambda = coef(mod.lambda)[2],
                          pgls.r2.lambda = summary(mod.lambda)$r.squared,
                          pgls.pvalue.lambda = summary(mod.lambda)$coefficients[2, 4],
                          #phyaov.pvalue.lambda.50 = phyaov.lambda.50$Pf,
                          #phyaov.pvalue.lambda.60 = phyaov.lambda.60$Pf,
                          phyaov.pvalue.lambda.80 = phyaov.lambda.80$Pf,
                          phyaov.pvalue.lambda.90 = phyaov.lambda.90$Pf,
                          phyaov.pvalue.lambda.100 = phyaov.lambda.100$Pf,
                          lm.int.lambda = coef(lm.lambda)[1],
                          lm.slope.lambda = coef(lm.lambda)[2],
                          lm.r2.lambda = summary(lm.lambda)$r.squared,
                          lm.pvalue.lambda = summary(lm.lambda)$coefficients[2, 4],
                          aov.pvalue.lambda.50 = summary(aov.lambda.50)[[1]][1, 5],
                          #aov.pvalue.lambda.60 = summary(aov.lambda.60)[[1]][1, 5],
                          aov.pvalue.lambda.80 = summary(aov.lambda.80)[[1]][1, 5],
                          aov.pvalue.lambda.90 = summary(aov.lambda.90)[[1]][1, 5],
                          aov.pvalue.lambda.100 = summary(aov.lambda.100)[[1]][1, 5],
                          pgls.int.mu = coef(mod.mu)[1],
                          pgls.slope.mu = coef(mod.mu)[2],
                          pgls.r2.mu = summary(mod.mu)$r.squared,
                          pgls.pvalue.mu = summary(mod.mu)$coefficients[2, 4],
                          #phyaov.pvalue.mu.50 = phyaov.mu.50$Pf,
                          #phyaov.pvalue.mu.60 = phyaov.mu.60$Pf,
                          phyaov.pvalue.mu.80 = phyaov.mu.80$Pf,
                          phyaov.pvalue.mu.90 = phyaov.mu.90$Pf,
                          phyaov.pvalue.mu.100 = phyaov.mu.100$Pf,
                          lm.int.mu = coef(lm.mu)[1],
                          lm.slope.mu = coef(lm.mu)[2],
                          lm.r2.mu = summary(lm.mu)$r.squared,
                          lm.pvalue.mu = summary(lm.mu)$coefficients[2, 4],
                          #aov.pvalue.mu.50 = summary(aov.mu.50)[[1]][1, 5],
                          #aov.pvalue.mu.60 = summary(aov.mu.60)[[1]][1, 5],
                          aov.pvalue.mu.80 = summary(aov.mu.80)[[1]][1, 5],
                          aov.pvalue.mu.90 = summary(aov.mu.90)[[1]][1, 5],
                          aov.pvalue.mu.100 = summary(aov.mu.100)[[1]][1, 5],
                          pgls.int.netdiv = coef(mod.netdiv)[1],
                          pgls.slope.netdiv = coef(mod.netdiv)[2],
                          pgls.r2.netdiv = summary(mod.netdiv)$r.squared,
                          pgls.pvalue.netdiv = summary(mod.netdiv)$coefficients[2, 4],
                          #phyaov.pvalue.netdiv.50 = phyaov.netdiv.50$Pf,
                          #phyaov.pvalue.netdiv.60 = phyaov.netdiv.60$Pf,
                          phyaov.pvalue.netdiv.80 = phyaov.netdiv.80$Pf,
                          phyaov.pvalue.netdiv.90 = phyaov.netdiv.90$Pf,
                          phyaov.pvalue.netdiv.100 = phyaov.netdiv.100$Pf,
                          lm.int.netdiv = coef(lm.netdiv)[1],
                          lm.slope.netdiv = coef(lm.netdiv)[2],
                          lm.r2.netdiv = summary(lm.netdiv)$r.squared,
                          lm.pvalue.netdiv = summary(lm.netdiv)$coefficients[2, 4],
                          #aov.pvalue.netdiv.50 = summary(aov.netdiv.50)[[1]][1, 5],
                          #aov.pvalue.netdiv.60 = summary(aov.netdiv.60)[[1]][1, 5],
                          aov.pvalue.netdiv.80 = summary(aov.netdiv.80)[[1]][1, 5],
                          aov.pvalue.netdiv.90 = summary(aov.netdiv.90)[[1]][1, 5],
                          aov.pvalue.netdiv.100 = summary(aov.netdiv.100)[[1]][1, 5],
                          #pgls.int.age.sh = coef(pgls.age.sh)[1],
                          #pgls.slope.age.sh = coef(pgls.age.sh)[2],
                          #pgls.r2.age.sh = summary(pgls.age.sh)$r.squared,
                          #pgls.pvalue.age.sh = summary(pgls.age.sh)$coefficients[2, 4],
                          pgls.int.rich.sh = coef(pgls.rich.sh)[1],
                          pgls.slope.rich.sh = coef(pgls.rich.sh)[2],
                          pgls.r2.rich.sh = summary(pgls.rich.sh)$r.squared,
                          pgls.pvalue.rich.sh = summary(pgls.rich.sh)$coefficients[2, 4],
                          #lm.int.age.sh = coef(lm.age.sh)[1],
                          #lm.slope.age.sh = coef(lm.age.sh)[2],
                          #lm.r2.age.sh = summary(lm.age.sh)$r.squared,
                          #lm.pvalue.age.sh = summary(lm.age.sh)$coefficients[2, 4],
                          lm.int.rich.sh = coef(lm.rich.sh)[1],
                          lm.slope.rich.sh = coef(lm.rich.sh)[2],
                          lm.r2.rich.sh = summary(lm.rich.sh)$r.squared,
                          lm.pvalue.rich.sh = summary(lm.rich.sh)$coefficients[2, 4], row.names = NULL)
    
    return(results)       
}    


registerDoMC(50)

## main.results <- ldply(1:nrep, main.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE)

main.results <- llply(1:nrep, main.analysis, age = age.data, fulltree = fulltree, .parallel = TRUE)

## write.table(main.results, file = "./output/main_results_zanne_2021.csv", sep = ",", quote = FALSE, row.names = FALSE)

saveRDS(main.results, file = here::here("output/main_results_igeatanen_2025.RDS"))
