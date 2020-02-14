library("tidyverse")
library("ape")
library("caper")
library("geiger")
library("phytools")
library("foreach")
library("doMC")
library("plyr")

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

mico.amostrator <- function(x){
    types <- c("AM", "EM", "NM", "OM", "ER")
    results <- setNames(sample(types[-match(x, types)], 1), NULL)
    return(results)
}


## Importing data
age.data <- read.csv("./data/data_all_families.csv", sep = ";")
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL

nrep <-  length(list.files("./output/simulated_datasets/"))

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
registerDoMC(2)

## Generating dummy data frame to replace with sampled data - Original type was assigned to each cell proportional to the number of species per type in the original dataset. We did it for 1000 different, randomly sampled datasets from the main analysis.

myco.sampling <- function(x, age, fulltree, reps){
    print(paste0("Replica ", x))
    family.data.gen <- read.csv(paste0("./output/simulated_datasets/random_data_", x, ".csv"), stringsAsFactors = FALSE)
    family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
    family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"
    ## Removing families with unknown mycorrhizal type
    family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]
    family.data.gen$stem.age <- age$stem.age[match(family.data.gen$family, age$familia)]
    family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
    family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)
    family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])
    family.data.gen <- family.data.gen[-match(c("Orchidaceae", "Ericaceae", "Diapensiaceae"), family.data.gen$family),]

    for(j in 1:reps){
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
        new.types <- setNames(sapply(expanded.data$orig.type[sp.to.change], FUN = mico.amostrator), NULL)
        expanded.data$new.type[sp.to.change] <- new.types
        data.sampled <- data.frame(
            family = unique(expanded.data$family),
            type = aggregate(expanded.data$new.type, by = list(expanded.data$family), FUN = table)$x
        )
        names(data.sampled) <- c("family", "AM", "EM", "ER", "NM", "OM")
        data.sampled <- cbind(data.sampled, data.sampled[, 2:6]/apply(data.sampled[,2:6], 1, sum))
        names(data.sampled)[7:11] <- c("AM.perc", "EM.perc", "ER.perc", "NM.perc", "OM.perc")
        data.sampled$type.50 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.5)
        data.sampled$type.60 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.6)
        data.sampled$type.80 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.8)
        data.sampled$type.100 <- apply(data.sampled[, 7:11], 1, mico.classificator, thresh = 0.999999)

        data.sampled$shannon <- vegan::diversity(data.sampled[, 2:6])
        data.sampled$rich <- family.data.gen$rich[match(data.sampled$family, family.data.gen$family)]
        data.sampled$r.e0 <- family.data.gen$r.e0[match(data.sampled$family, family.data.gen$family)]
        data.sampled$r.e09 <- family.data.gen$r.e09[match(data.sampled$family, family.data.gen$family)]
        write.table(data.sampled, file = paste0("./output/20perc_data/20perc_data_tree", x, "_replica", j, ".csv"), quote = FALSE, sep = ",", row.names = FALSE)
        }
    }

llply(as.list(1:100), .fun = myco.sampling, age = age.data, fulltree = fulltree, reps = 10, .parallel = TRUE)
