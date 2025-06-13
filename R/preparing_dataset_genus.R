library("reshape2")
here::i_am("R/preparing_dataset_genus.R")

genus.soudz <- read.csv(here::here("data/richness_mico_genera_soudz_mujica.csv"), stringsAsFactors = FALSE)[, 2:5]

genus.soudz <- genus.soudz[-which(genus.soudz$genus == "Thysanotus"),]

conversion <- data.frame(raw = unique(genus.soudz$mico),
           new = c("MIX", "UNK", "AM", "NM", "MIX", "UNK", "EM", "OM", "NM", "MIX", "MIX", "ER", "MIX", "MIX")
           )

genus.soudz$mico.new <- conversion$new[match(genus.soudz$mico, conversion$raw)]

genus.soudz.sh <- data.frame(
    genus = genus.soudz$genus,
    family = genus.soudz$family,
    AM = rep(0, nrow(genus.soudz)),
    EM = rep(0, nrow(genus.soudz)),
    NM = rep(0, nrow(genus.soudz)),
    OM = rep(0, nrow(genus.soudz)),
    ER = rep(0, nrow(genus.soudz))
)

genus.soudz.sh$AM[which(genus.soudz$mico.new == "AM")] <- genus.soudz$richness[which(genus.soudz$mico.new == "AM")]
genus.soudz.sh$NM[which(genus.soudz$mico.new == "NM")] <- genus.soudz$richness[which(genus.soudz$mico.new == "NM")]
genus.soudz.sh$EM[which(genus.soudz$mico.new == "EM")] <- genus.soudz$richness[which(genus.soudz$mico.new == "EM")]
genus.soudz.sh$OM[which(genus.soudz$mico.new == "OM")] <- genus.soudz$richness[which(genus.soudz$mico.new == "OM")]
genus.soudz.sh$ER[which(genus.soudz$mico.new == "ER")] <- genus.soudz$richness[which(genus.soudz$mico.new == "ER")]


genus.soudz.sh$AM[which(genus.soudz$mico == "NM-AM")] <- genus.soudz.sh$AM[which(genus.soudz$mico == "NM-AM")] + (0.5 * genus.soudz$richness[which(genus.soudz$mico == "NM-AM")])
genus.soudz.sh$AM[which(genus.soudz$mico == "EcM-AM")] <- genus.soudz.sh$AM[which(genus.soudz$mico == "EcM-AM")] + (0.5 * genus.soudz$richness[which(genus.soudz$mico == "EcM-AM")])

genus.soudz.sh$EM[which(genus.soudz$mico == "EcM-AM")] <- genus.soudz.sh$EM[which(genus.soudz$mico == "EcM-AM")] + (0.5 * genus.soudz$richness[which(genus.soudz$mico == "EcM-AM")])

genus.soudz.sh$NM[which(genus.soudz$mico == "NM-AM")] <- genus.soudz.sh$NM[which(genus.soudz$mico == "NM-AM")] + (0.5 * genus.soudz$richness[which(genus.soudz$mico == "NM-AM")])

rich.per.type.sh <- dcast(aggregate(melt(genus.soudz.sh)$value, by = list(melt(genus.soudz.sh)$family, melt(genus.soudz.sh)$variable), FUN = sum), Group.1 ~ Group.2)

family.data <- data.frame(
    family = unique(genus.soudz$family),
    rich = aggregate(genus.soudz$richness, by = list(genus.soudz$family), FUN = sum)$x,
    rich.per.type.sh[, 2:6]
    )

family.data <- cbind(family.data, family.data[, 3:7]/apply(family.data[, 3:7], 1, sum))
names(family.data)[8:12] <- paste0(names(family.data)[3:7], ".perc")

weigh.type <- function(x){
    data <- x
    types.raw <- rep(data$mico.new, times = data$richness)
    types.agg <- table(types.raw)
    return(types.agg)
}

family.data[, 13:16] <- NA
for(i in 1:length(unique(family.data$family))){
    tmp <- unlist(weigh.type(genus.soudz[genus.soudz$family == unique(genus.soudz$family)[i],]))
    tmp <- setNames(tmp[match(c("AM", "MIX", "NM", "UNK"), names(tmp))], c("AM", "MIX", "NM", "UNK"))
    family.data[i, 13:16] <- tmp
}

names(family.data)[13:16] <- paste0(names(weigh.type(genus.soudz[genus.soudz$family == unique(genus.soudz$family)[1],])), ".raw")

family.data <- cbind(family.data, setNames(family.data[, 13:16] / apply(family.data[, 13:16], 1, FUN = sum), paste0(names(family.data)[13:16], ".perc")))




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



## family.data$type.50 <- apply(family.data[, 8:12], 1, mico.classificator, thresh = 0.5)
## family.data$type.60 <- apply(family.data[, 8:12], 1, mico.classificator, thresh = 0.6)
family.data$type.80 <- apply(family.data[, 8:12], 1, mico.classificator, thresh = 0.8)
family.data$type.90 <- apply(family.data[, 8:12], 1, mico.classificator, thresh = 0.9)
family.data$type.100 <- apply(family.data[, 8:12], 1, mico.classificator, thresh = 0.99999999)

write.table(family.data, here::here("data/family_data_genus_classif.csv"), sep = ",", quote = FALSE, row.names = FALSE)
