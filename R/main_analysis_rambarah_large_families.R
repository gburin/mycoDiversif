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
RB.tree.RC.complete <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_RC_complete_MCCv_2.tre")
RB.tree.UC.complete <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_UC_complete_MCCv_2.tre")
RB.tree.CC.complete <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_CC_complete_MCCv_2.tre")

RB.tree.RC.conservative <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_RC_conservative_MCCv_2.tre")
RB.tree.UC.conservative <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_UC_conservative_MCCv_2.tre")
RB.tree.CC.conservative <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_CC_conservative_MCCv_2.tre")

## Extracting families from tip information

RB.tip.list <- RB.tree.RC.complete$tip.label
RB.family.list <- unique(setNames(sapply(RB.tip.list, function(x){strsplit(x, split = "_")[[1]][2]}), NULL))

## Replacing species name by the family, and dropping repeated tips
RB.tree.RC.complete.pruned <- drop.tip(RB.tree.RC.complete, tip = grep(RB.family.list[1], RB.tree.RC.complete$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.RC.complete.pruned <- drop.tip(RB.tree.RC.complete.pruned, tip = grep(RB.family.list[i], RB.tree.RC.complete.pruned$tip.label)[-1])
}
RB.tree.RC.complete.pruned$tip.label <- setNames(sapply(RB.tree.RC.complete.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)

RB.tree.UC.complete.pruned <- drop.tip(RB.tree.UC.complete, tip = grep(RB.family.list[1], RB.tree.UC.complete$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.UC.complete.pruned <- drop.tip(RB.tree.UC.complete.pruned, tip = grep(RB.family.list[i], RB.tree.UC.complete.pruned$tip.label)[-1])
}
RB.tree.UC.complete.pruned$tip.label <- setNames(sapply(RB.tree.UC.complete.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)

RB.tree.CC.complete.pruned <- drop.tip(RB.tree.CC.complete, tip = grep(RB.family.list[1], RB.tree.CC.complete$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.CC.complete.pruned <- drop.tip(RB.tree.CC.complete.pruned, tip = grep(RB.family.list[i], RB.tree.CC.complete.pruned$tip.label)[-1])
}
RB.tree.CC.complete.pruned$tip.label <- setNames(sapply(RB.tree.CC.complete.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)



RB.tree.RC.conservative.pruned <- drop.tip(RB.tree.RC.conservative, tip = grep(RB.family.list[1], RB.tree.RC.conservative$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.RC.conservative.pruned <- drop.tip(RB.tree.RC.conservative.pruned, tip = grep(RB.family.list[i], RB.tree.RC.conservative.pruned$tip.label)[-1])
}
RB.tree.RC.conservative.pruned$tip.label <- setNames(sapply(RB.tree.RC.conservative.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)

RB.tree.UC.conservative.pruned <- drop.tip(RB.tree.UC.conservative, tip = grep(RB.family.list[1], RB.tree.UC.conservative$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.UC.conservative.pruned <- drop.tip(RB.tree.UC.conservative.pruned, tip = grep(RB.family.list[i], RB.tree.UC.conservative.pruned$tip.label)[-1])
}
RB.tree.UC.conservative.pruned$tip.label <- setNames(sapply(RB.tree.UC.conservative.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)

RB.tree.CC.conservative.pruned <- drop.tip(RB.tree.CC.conservative, tip = grep(RB.family.list[1], RB.tree.CC.conservative$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.CC.conservative.pruned <- drop.tip(RB.tree.CC.conservative.pruned, tip = grep(RB.family.list[i], RB.tree.CC.conservative.pruned$tip.label)[-1])
}
RB.tree.CC.conservative.pruned$tip.label <- setNames(sapply(RB.tree.CC.conservative.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)


rb.RC.complete.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_RC_complete.csv", as.is = TRUE)
rb.UC.complete.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_UC_complete.csv", as.is = TRUE)
rb.CC.complete.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_CC_complete.csv", as.is = TRUE)

rb.RC.conservative.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_RC_conservative.csv", as.is = TRUE)
rb.UC.conservative.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_UC_conservative.csv", as.is = TRUE)
rb.CC.conservative.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_CC_conservative.csv", as.is = TRUE)


## Missing families from Ramírez-Barahona et al. 2020
age.data$familia[is.na(match(age.data$familia, rb.RC.complete.ages$Family))]

## Adding ages from Ramírez-Barahona et al. 2020
age.data$rb.RC.complete.stem <- rb.RC.complete.ages$Stem_BEAST[match(age.data$familia, rb.RC.complete.ages$Family)]
age.data$rb.UC.complete.stem <- rb.RC.complete.ages$Stem_BEAST[match(age.data$familia, rb.UC.complete.ages$Family)]
age.data$rb.CC.complete.stem <- rb.RC.complete.ages$Stem_BEAST[match(age.data$familia, rb.CC.complete.ages$Family)]

age.data$rb.RC.conservative.stem <- rb.RC.conservative.ages$Stem_BEAST[match(age.data$familia, rb.RC.conservative.ages$Family)]
age.data$rb.UC.conservative.stem <- rb.RC.conservative.ages$Stem_BEAST[match(age.data$familia, rb.UC.conservative.ages$Family)]
age.data$rb.CC.conservative.stem <- rb.RC.conservative.ages$Stem_BEAST[match(age.data$familia, rb.CC.conservative.ages$Family)]


nrep = length(list.files("../output/simulated_datasets/"))


### Calculating phylogenetic signal to be used in PGLS/phylANOVA only once to speed up the script

family.data.gen <- read.csv("../output/simulated_datasets/random_data_00001.csv", stringsAsFactors = FALSE)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"
## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

calibs <- c("RC.complete", "UC.complete", "CC.complete", "RC.conservative", "UC.conservative", "CC.conservative")

for(i in 1:length(calibs)){
    family.data.gen[, calibs[i]] <- age.data[match(family.data.gen$family, age.data$familia), paste0("rb.", calibs[i], ".stem")]
    family.data.gen[, paste0(calibs[i], ".r.e0")] <- bd.ms(time = family.data.gen[, calibs[i]], n = family.data.gen$rich, crown = FALSE, epsilon = 0)
    family.data.gen[, paste0(calibs[i], ".r.e09")] <- bd.ms(time = family.data.gen[, calibs[i]], n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)
}

family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae", "Diapensiaceae"), family.data.gen$family),]
family.data.gen <- family.data.gen[family.data.gen$rich >= 20,]

for(i in 1:length(calibs)){
    assign(paste0("tree.pruned.", calibs[i]),
           drop.tip(get(paste0("RB.tree.", calibs[i], ".pruned")),
                    tip = get(paste0("RB.tree.", calibs[i], ".pruned"))$tip.label[is.na(match(get(paste0("RB.tree.", calibs[i], ".pruned"))$tip.label, family.data.gen$family))]))
    assign(paste0("data.pgls.", calibs[i]), comparative.data(get(paste0("tree.pruned.", calibs[i])), family.data.gen, names.col = "family"))
    assign(paste0("phylosig.r0.", calibs[i]), phylosig(get(paste0("data.pgls.", calibs[i]))$phy, setNames(get(paste0("data.pgls.", calibs[i]))$data[, paste0(calibs[i], ".r.e0")], rownames(get(paste0("data.pgls.", calibs[i]))$data)), method = "lambda", test = TRUE))
    assign(paste0("phylosig.r09.", calibs[i]), phylosig(get(paste0("data.pgls.", calibs[i]))$phy, setNames(get(paste0("data.pgls.", calibs[i]))$data[, paste0(calibs[i], ".r.e09")], rownames(get(paste0("data.pgls.", calibs[i]))$data)), method = "lambda", test = TRUE))
    assign(paste0("phylosig.rich.", calibs[i]), phylosig(get(paste0("data.pgls.", calibs[i]))$phy, setNames(get(paste0("data.pgls.", calibs[i]))$data[, "rich"], rownames(get(paste0("data.pgls.", calibs[i]))$data)), method = "lambda", test = TRUE))
    assign(paste0("phylosig.age.", calibs[i]), phylosig(get(paste0("data.pgls.", calibs[i]))$phy, setNames(get(paste0("data.pgls.", calibs[i]))$data[, calibs[i]], rownames(get(paste0("data.pgls.", calibs[i]))$data)), method = "lambda", test = TRUE))
}

main.analysis <- function(x, age, fulltree, calib){
    print(paste0("Replica ", x, " of ", length(list.files("../output/simulated_datasets/"))))
    family.data.gen <- read.csv(paste0("../output/simulated_datasets/random_data_", sprintf("%05d", x), ".csv"), stringsAsFactors = FALSE)
    family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
    family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"
    ## Removing families with unknown mycorrhizal type
    family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

    family.data.gen$stem.age <- age.data[match(family.data.gen$family, age$familia), paste0("rb.", calib, ".stem")]
    family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
    family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

    family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

    family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae", "Diapensiaceae"), family.data.gen$family),]

    tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
    data.pgls <- comparative.data(tree.pruned, family.data.gen, names.col = "family")

    ## Fitting PGLS excluding families with != 100% MIX
    mod.r0 <- caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = get(paste0("phylosig.r0.", calib))$lambda)
    mod.r09 <- caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = get(paste0("phylosig.r09.", calib))$lambda)

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
    pgls.age.sh <- caper::pgls(stem.age ~ shannon, data = data.pgls, lambda = get(paste0("phylosig.rich.", calib))$lambda)
    pgls.rich.sh <- caper::pgls(rich ~ shannon, data = data.pgls, lambda = get(paste0("phylosig.age.", calib))$lambda)

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

results.RC.complete.large <- llply(1:nrep, main.analysis, age = age.data, fulltree = RB.tree.RC.complete.pruned, calib = "RC.complete", .parallel = TRUE)
save(results.RC.complete.large, file = "../output/RB_results_RC_complete_large.RData")
rm(results.RC.complete)
results.UC.complete.large <- llply(1:nrep, main.analysis, age = age.data, fulltree = RB.tree.UC.complete.pruned, calib = "UC.complete", .parallel = TRUE)
save(results.UC.complete.large, file = "../output/RB_results_UC_complete_large.RData")
rm(results.UC.complete)
results.CC.complete.large <- llply(1:nrep, main.analysis, age = age.data, fulltree = RB.tree.CC.complete.pruned, calib = "CC.complete", .parallel = TRUE)
save(results.CC.complete.large, file = "../output/RB_results_CC_complete_large.RData")
rm(results.CC.complete)

results.RC.conservative.large <- llply(1:nrep, main.analysis, age = age.data, fulltree = RB.tree.RC.conservative.pruned, calib = "RC.conservative", .parallel = TRUE)
save(results.RC.conservative, file = "../output/RB_results_RC_conservative_large.RData")
rm(results.RC.conservative)
results.UC.conservative.large <- llply(1:nrep, main.analysis, age = age.data, fulltree = RB.tree.UC.conservative.pruned, calib = "UC.conservative", .parallel = TRUE)
save(results.UC.conservative, file = "../output/RB_results_UC_conservative_large.RData")
rm(results.UC.conservative)
results.CC.conservative.large <- llply(1:nrep, main.analysis, age = age.data, fulltree = RB.tree.CC.conservative.pruned, calib = "CC.conservative", .parallel = TRUE)
save(results.CC.conservative, file = "../output/RB_results_CC_conservative_large.RData")
rm(results.CC.conservative)

## write.table(main.results, file = "./output/fit_data_random_datasets.csv", sep = ",", quote = FALSE, row.names = FALSE)

## save(main.results, file = "./output/anova_posthoc_tables.RData")
