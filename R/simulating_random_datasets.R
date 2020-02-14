###
library("caper")
library("tidyverse")
library("cowplot")
library("RColorBrewer")
library("geiger")
library("phytools")
library("ggrepel")
library("foreach")
library("doMC")
library("reshape2")
library("plyr")

mico.classificator <- function(x, thresh = 0.5){
    if(any(sapply(x, is.nan))){
        return("UNK")
    } else if(any(x == 1)){
        return(gsub(".perc", "", names(x)[which(x == 1)]))
    } else if(sum(x > thresh) != 0){
        return(gsub(".perc", "", names(x)[which.max(x)]))
    } else {
        return("MIX")
    }
}

weigh.type <- function(x){data <- x
    types.raw <- rep(data$mico.new, times = data$richness)
    types.agg <- table(types.raw)
    return(types.agg)
}

##############################
### Generating random datasets
##############################

nrep <- 10000

genus.soudz <- read.csv("./data/richness_mico_genera_soudz_mujica.csv", stringsAsFactors = FALSE)[, 2:5]

genus.soudz <- genus.soudz[-which(genus.soudz$genus == "Thysanotus"),]

conversion <- data.frame(raw = unique(genus.soudz$mico),
           new = c("MIX", "UNK", "AM", "NM", "MIX", "UNK", "EM", "OM", "NM", "UNK", "UNK", "ER", "UNK", "UNK")
           )

genus.soudz$mico.new <- conversion$new[match(genus.soudz$mico, conversion$raw)]

foo <- function(replica, data){
    print(replica)
    genus.soudz.sh <- data.frame(
        genus = data$genus,
        family = data$family,
        AM = rep(0, nrow(data)),
        EM = rep(0, nrow(data)),
        NM = rep(0, nrow(data)),
        OM = rep(0, nrow(data)),
        ER = rep(0, nrow(data)),
        UNK = rep(0, nrow(data))
    )
    mix.indices <- which(data$mico.new == "MIX")
    genus.soudz.sh$AM[which(data$mico.new == "AM")] <- data$richness[which(data$mico.new == "AM")]
    genus.soudz.sh$NM[which(data$mico.new == "NM")] <- data$richness[which(data$mico.new == "NM")]
    genus.soudz.sh$EM[which(data$mico.new == "EM")] <- data$richness[which(data$mico.new == "EM")]
    genus.soudz.sh$OM[which(data$mico.new == "OM")] <- data$richness[which(data$mico.new == "OM")]
    genus.soudz.sh$ER[which(data$mico.new == "ER")] <- data$richness[which(data$mico.new == "ER")]
    genus.soudz.sh$UNK[which(data$mico.new == "UNK")] <- data$richness[which(data$mico.new == "UNK")]
    for(i in 1:length(mix.indices)){
        if(data$mico[mix.indices[i]] == "NM-AM"){
            prop.samp <- sample(0:data$rich[mix.indices[i]], 1)
            prop.samp <- c(prop.samp, data$rich[mix.indices[i]] - prop.samp)
            genus.soudz.sh$NM[mix.indices[i]] <- prop.samp[1]
            genus.soudz.sh$AM[mix.indices[i]] <- prop.samp[2]
        } else if(data$mico[mix.indices[i]] == "EcM-AM"){
            prop.samp <- sample(0:data$rich[mix.indices[i]], 1)
            prop.samp <- c(prop.samp, data$rich[mix.indices[i]] - prop.samp)
            genus.soudz.sh$EM[mix.indices[i]] <- prop.samp[1]
            genus.soudz.sh$AM[mix.indices[i]] <- prop.samp[2]
        }
        }

    rich.per.type.sh <- dcast(aggregate(melt(genus.soudz.sh)$value, by = list(melt(genus.soudz.sh)$family, melt(genus.soudz.sh)$variable), FUN = sum), Group.1 ~ Group.2)

    family.data <- data.frame(
        family = unique(data$family),
        rich = aggregate(data$richness, by = list(data$family), FUN = sum)$x,
        rich.per.type.sh[, 2:7]
    )

    family.data <- cbind(family.data, family.data[, 3:8]/apply(family.data[, 3:8], 1, sum))
    names(family.data)[9:14] <- paste0(names(family.data)[3:8], ".perc")

    family.data$type.50 <- apply(family.data[, 9:14], 1, mico.classificator, thresh = 0.5)
    family.data$type.60 <- apply(family.data[, 9:14], 1, mico.classificator, thresh = 0.6)
    family.data$type.80 <- apply(family.data[, 9:14], 1, mico.classificator, thresh = 0.8)
    family.data$type.100 <- apply(family.data[, 9:14], 1, mico.classificator, thresh = 0.99999999)

    write.table(family.data, file = paste0("./output/simulated_datasets/random_data_", replica, ".csv"), quote = FALSE, row.names = FALSE, sep = ",")
}

registerDoMC(3)

sim.datasets <- llply(1:nrep, foo, data = genus.soudz, .parallel = TRUE)
