library("tidyverse")
library("ape")
library("caper")
library("geiger")
library("phytools")
library("foreach")
library("doMC")
library("plyr")

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

registerDoMC(50)

## Generating dummy data frame to replace with sampled data - Original type was assigned to each cell proportional to the number of species per type in the original dataset

myco.sampling <- function(x){
    expanded.data <- data.frame(
        family = rep(family.data.gen$family, times = apply(round(family.data.gen[, 3:7]), 1, sum)),
        orig.type = unlist(apply(family.data.gen[, 3:7], 1, function(x){rep(gsub("rich.sh.", "", names(x)), times = round(x))})),
        stringsAsFactors = FALSE
    )

    ## Sampling cells that will be changed
    sp.to.change <- sort(sample(1:nrow(expanded.data), ceiling(0.2 * nrow(expanded.data))))

    expanded.data$new.type <- NA
    expanded.data$new.type[-sp.to.change] <- as.character(expanded.data$orig.type[-sp.to.change])

    ## Substituting sampled lines for new type (excluding original type)
    for(i in 1:length(sp.to.change)){
        print(paste(i, "of", length(sp.to.change), sep = " "))
        expanded.data$new.type[sp.to.change[i]] <- sample(unique(expanded.data$orig.type)[-which(unique(expanded.data$orig.type) == expanded.data$orig.type[sp.to.change[i]])], 1)
    }

    data.sampled <- data.frame(
        family = unique(expanded.data$family),
        type = aggregate(expanded.data$new.type, by = list(expanded.data$family), FUN = table)$x
    )

    names(data.sampled) <- c("family", "AM", "EM", "ER", "NM", "OM")

    data.sampled <- cbind(data.sampled, data.sampled[, 2:6]/apply(data.sampled[,2:6], 1, sum))
    names(data.sampled)[7:11] <- c("AM.perc", "EM.perc", "ER.perc", "NM.perc", "OM.perc")

    mico.classificator <- function(x, thresh = 0.5){
        if(any(sapply(x, is.nan))){
            return("UNK")
        } else if(any(x == 1)){
            return(gsub(".valid", "", gsub(".perc", "", names(x)[which(x == 1)])))
        } else if(sum(x > thresh) != 0){
            return(gsub(".valid", "", gsub(".perc", "", names(x)[which.max(x)])))
        } else {
            return("MIX")
        }
    }

    data.sampled$type.50 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.5)
    data.sampled$type.60 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.6)
    data.sampled$type.80 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.8)
    data.sampled$type.100 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.999999)

    data.sampled$shannon <- vegan::diversity(data.sampled[, 2:6])
    data.sampled$r.e0 <- family.data.gen$r.e0[match(data.sampled$family, family.data.gen$family)]
    data.sampled$r.e09 <- family.data.gen$r.e09[match(data.sampled$family, family.data.gen$family)]
    return(data.sampled)
}

sampled.datasets <- llply(as.list(1:50), .fun = myco.sampling, .parallel = TRUE)

save(sampled.datasets, file = "./output/sampled_datasets_20perc.RData")

tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, data.sampled$family))])
data.pgls <- comparative.data(tree.pruned, data.sampled, names.col = "family")

phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")

## PGLS
mod.r0 <- caper::pgls(r.e0 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))
mod.r09 <- caper::pgls(r.e09 ~ shannon, data = data.pgls, lambda = setNames(summary(phylosig.r0)$param[2], NULL))

## LM
lm.r0 <- lm(r.e0 ~ shannon, data = data.sampled)
lm.r09 <- lm(r.e09 ~ shannon, data = data.sampled)

## phyANOVA
data.aov <- data.sampled[-which(is.na(match(data.sampled$family, fulltree$tip.label))), ]
data.aov <- data.aov[-which(is.na(data.aov$r.e0)),]
data.aov <- data.aov[-which(data.aov$type.60 == "ER"), ]

phyaov.r0 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e0, data.aov$family))
phyaov.r09 <- phylANOVA(drop.tip(tree.pruned, tip = tree.pruned$tip.label[which(is.na(match(tree.pruned$tip.label, data.aov$family)))]), x = setNames(data.aov$type.60, data.aov$family), y = setNames(data.aov$r.e09, data.aov$family))

