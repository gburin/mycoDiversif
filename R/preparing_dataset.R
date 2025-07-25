## library("taxize")

## ################
## raw.data <- read.csv("./data/full_database_Soudzilovskaia_2019.csv", stringsAsFactors = FALSE)
## ## raw.data$species <- gsub(".*\ $", "", raw.data$species, perl = TRUE)
## ## raw.data$species <- gsub("  ", " ", raw.data$species)
## ## raw.data$species <- gsub(" ", "_", raw.data$species)
## raw.data$genus <- sapply(raw.data$species, function(x){strsplit(x, split = "_")[[1]][1]})
## ## raw.data$genus[grep("\n", raw.data$genus)] <- c("Achatocarpus", "Guapira")
## ## raw.data$species <- gsub("\n", "_", raw.data$species)

## full_genera <- read.csv("./data/Spermatophyta_Genera_edit.csv", stringsAsFactors = FALSE)

## ## Checking how many species with remarks by any author
## sum(apply(raw.data[, c(7, 9)], 1, function(x){any(!is.na(x))}))

## ## Cleaning raw data
## data.clean <- raw.data[which(apply(raw.data[, c(7, 9)], 1, function(x){any(!is.na(x))}) == FALSE), c("genus", "species", "mycorrhiza.type", "number.of.plants.studied")]
## data.clean$family <- full_genera$family[match(data.clean$genus, full_genera$genus)]
## data.clean$family[grep(".*aceae$", x = data.clean$genus)] <- data.clean$genus[grep(".*aceae$", x = data.clean$genus)]

## missing.genera <- sort(unique(data.clean$genus[is.na(data.clean$family)]))

## missing.family <- vector(mode = "character", length = length(missing.genera))

## for(i in 1:length(missing.genera)){
##     print(paste(i, "of", length(missing.genera)))
##     missing.family[i] <- ipni_search(genus = missing.genera[i])$family[1]
## }

## ## Fixing missing data by hand
## missing.genera[is.na(missing.family)]
## missing.family[which(missing.genera == "Acrosticum")] <- "Pteridaceae"
## missing.family[which(missing.genera == "Conocephalum")] <- "Conocephalaceae"
## missing.family[which(missing.genera == "Cratoxilon")] <- "Hypericaceae"
## missing.family[which(missing.genera == "Dendroceros")] <- "Dendrocerotaceae"
## missing.family[which(missing.genera == "Deplauchea")] <- "Bignoniaceae"
## missing.family[which(missing.genera == "Folioceros")] <- "Anthocerotaceae"
## missing.family[which(missing.genera == "Lachenalii")] <- "Pteridaceae"
## #missing.family[which(missing.genera == "Leiosporoceros")] <- "Leiosporocerotaceae"
## missing.family[which(missing.genera == "Livistonia")] <- "Arecaceae"
## missing.family[which(missing.genera == "Lunularia")] <- "Lunulariaceae"
## missing.family[which(missing.genera == "Megaceros")] <- "Dendrocerotaceae"
## missing.family[which(missing.genera == "Myrovernix")] <- "Asteraceae"
## missing.family[which(missing.genera == "Nornanbya")] <- "Arecaceae"
## missing.family[which(missing.genera == "Nothoceros")] <- "Dendrocerotaceae"
## missing.family[which(missing.genera == "Notothylas")] <- "Notothyladaceae"
## missing.family[which(missing.genera == "Pellia")] <- "Pelliaceae"
## missing.family[which(missing.genera == "Pentaphylla")] <- "Pteridaceae"
## missing.family[which(missing.genera == "Phaeoceros")] <- "Notothyladaceae"
## missing.family[which(missing.genera == "Phaeomegaceros")] <- "Dendrocerotaceae"
## missing.family[which(missing.genera == "Qulllaja")] <- "Quillajaceae"


## data.clean$family[is.na(data.clean$family)] <- missing.family[match(data.clean$genus[is.na(data.clean$family)], missing.genera)]
## data.clean$number.of.plants.studied <- as.integer(data.clean$number.of.plants.studied)

## write.table(data.clean, file = "./data/soudz_w_family.csv", sep = "\t", quote = FALSE, row.names = FALSE)




### Includings species for which the remarks point towards a probable type

library("taxize")
here::i_am("R/preparing_dataset.R")

################
raw.data <- read.csv(here::here("data/full_database_Soudzilovskaia_2019.csv"), stringsAsFactors = FALSE)
## raw.data$species <- gsub(".*\ $", "", raw.data$species, perl = TRUE)
## raw.data$species <- gsub("  ", " ", raw.data$species)
## raw.data$species <- gsub(" ", "_", raw.data$species)
raw.data$genus <- sapply(raw.data$species, function(x){strsplit(x, split = "_")[[1]][1]})
## raw.data$genus[grep("\n", raw.data$genus)] <- c("Achatocarpus", "Guapira")
## raw.data$species <- gsub("\n", "_", raw.data$species)

full_genera <- read.csv(here::here("data/Spermatophyta_Genera.csv"), stringsAsFactors = FALSE)

## Checking how many species with remarks by any author
sum(apply(raw.data[, c(7, 9)], 1, function(x){any(!is.na(x))}))

raw.data$mycorrhiza.type[which(raw.data$curator_remark_1_comment == "probably AM")] <- "AM"
raw.data$curator_remark_1_comment[which(raw.data$curator_remark_1_comment == "probably AM")] <- NA
raw.data$mycorrhiza.type[which(raw.data$curator_remark_1_comment == "probably ErM")] <- "ER"
raw.data$curator_remark_1_comment[which(raw.data$curator_remark_1_comment == "probably ErM")] <- NA
raw.data$mycorrhiza.type[which(raw.data$curator_remark_1_comment == "probably NM")] <- "NM"
raw.data$curator_remark_1_comment[which(raw.data$curator_remark_1_comment == "probably NM")] <- NA
raw.data$mycorrhiza.type[which(raw.data$curator_remark_1_comment == "probably OM")] <- "OM"
raw.data$curator_remark_1_comment[which(raw.data$curator_remark_1_comment == "probably OM")] <- NA

## Cleaning raw data
data.clean <- raw.data[which(apply(raw.data[, c(7, 9)], 1, function(x){any(!is.na(x))}) == FALSE), c("genus", "species", "mycorrhiza.type", "number.of.plants.studied")]
data.clean$family <- full_genera$family[match(data.clean$genus, full_genera$genus)]
data.clean$family[grep(".*aceae$", x = data.clean$genus)] <- data.clean$genus[grep(".*aceae$", x = data.clean$genus)]

missing.genera <- sort(unique(data.clean$genus[is.na(data.clean$family)]))

missing.family <- vector(mode = "character", length = length(missing.genera))

for(i in 1:length(missing.genera)){
    print(paste(i, "of", length(missing.genera)))
    tmp <- tryCatch(ipni_search(genus = missing.genera[i]), error = function(x){NA})
    if(dim(tmp)[1] != 0){
        missing.family[i] <- tmp$family[1]
        }
}

## Fixing missing data by hand
missing.genera[missing.family == ""]
missing.family[which(missing.genera == "Acrosticum")] <- "Pteridaceae"
missing.family[which(missing.genera == "Conocephalum")] <- "Conocephalaceae"
missing.family[which(missing.genera == "Cratoxilon")] <- "Hypericaceae"
missing.family[which(missing.genera == "Dendroceros")] <- "Dendrocerotaceae"
missing.family[which(missing.genera == "Deplauchea")] <- "Bignoniaceae"
missing.family[which(missing.genera == "Folioceros")] <- "Anthocerotaceae"
missing.family[which(missing.genera == "Lachenalii")] <- "Pteridaceae"
#missing.family[which(missing.genera == "Leiosporoceros")] <- "Leiosporocerotaceae"
missing.family[which(missing.genera == "Livistonia")] <- "Arecaceae"
#missing.family[which(missing.genera == "Lunularia")] <- "Lunulariaceae"
missing.family[which(missing.genera == "Megaceros")] <- "Dendrocerotaceae"
#missing.family[which(missing.genera == "Myrovernix")] <- "Asteraceae"
missing.family[which(missing.genera == "Nornanbya")] <- "Arecaceae"
missing.family[which(missing.genera == "Nothoceros")] <- "Dendrocerotaceae"
missing.family[which(missing.genera == "Notothylas")] <- "Notothyladaceae"
missing.family[which(missing.genera == "Pellia")] <- "Pelliaceae"
missing.family[which(missing.genera == "Petalostemum")] <- "Fabaceae"
missing.family[which(missing.genera == "Pentaphylla")] <- "Pteridaceae"
missing.family[which(missing.genera == "Phaeoceros")] <- "Notothyladaceae"
missing.family[which(missing.genera == "Phaeomegaceros")] <- "Dendrocerotaceae"
missing.family[which(missing.genera == "Qulllaja")] <- "Quillajaceae"


data.clean$family[is.na(data.clean$family)] <- missing.family[match(data.clean$genus[is.na(data.clean$family)], missing.genera)]
## data.clean$number.of.plants.studied <- as.integer(data.clean$number.of.plants.studied)

write.table(data.clean, file = here::here("data/soudz_w_family_incl_remarks.csv"), sep = "\t", quote = FALSE, row.names = FALSE)





##### Analysis

fulldata.fam <- read.csv(here::here("data/soudz_w_family_incl_remarks.csv"), sep = "\t", stringsAsFactors = FALSE, row.names = NULL)

fulldata <- fulldata.fam[-which(sapply(strsplit(fulldata.fam$species, split = "_", fixed = TRUE), length) == 1), ]


## First approach: using classification straight from Soudzilovskaia et al. 2019

raw.types <- unique(fulldata$mycorrhiza.type)
new.types <- c("OM", "AM", "EM", "EM", "AM", "AM", "MIX", "NM", "NM", "OTHER", "ER", "NV", "NM", "MIX", "OTHER")

fulldata$myc.type.mujica <- new.types[match(fulldata$mycorrhiza.type, raw.types)]

## Removing non-vascular plants and the 6 species classified as "other"
fulldata <- fulldata[fulldata$myc.type.mujica != "NV", ]
fulldata <- fulldata[fulldata$myc.type.mujica != "OTHER", ]

families.tax.probs <- read.csv(here::here("data/families_with_taxonomic_issues.csv"), sep = ",")

fulldata$family[match(families.tax.probs$species, fulldata$species)] <- as.character(families.tax.probs$new.class)
fulldata$family[fulldata$family == "Salviniaceae"] <- "FERN"

## Removing Ferns and Mosses
fulldata <- fulldata[fulldata$family != "FERN",]
fulldata <- fulldata[fulldata$family != "moss",]


#### FIX TABULATION BELOW


## Tabulating data per species
fulldata.tab <- setNames(as.data.frame(matrix(NA, nrow = length(unique(fulldata$species)), ncol = 7)), c("species", "AM", "EM", "ER", "MIX", "NM", "OM"))

for(i in 1:length(unique(fulldata$species))){
    print(i)
    tmp <- table(fulldata$myc.type.mujica[fulldata$species == unique(fulldata$species)[i]])
    fulldata.tab[i,] <- c(unique(fulldata$species)[i], tmp[match(c("AM", "EM", "ER", "MIX", "NM", "OM"), names(tmp))])
}

fulldata.tab <- plyr::llply(1:length(unique(fulldata$species)), function(x){table(fulldata$myc.type.mujica[which(fulldata$species == unique(fulldata$species)[x])])})
names(fulldata.tab) <- c("species", "AM", "EM", "ER", "MIX", "NM", "OM")
fulldata.tab$genus <- fulldata$genus[match(fulldata.tab$species, fulldata$species)]
fulldata.tab$family <- fulldata$family[match(fulldata.tab$species, fulldata$species)]

## Reorganizing columns, just because :D
fulldata.tab <- fulldata.tab[, c(1, 8:9, 2:7)]
fulldata.tab[, c(4:9)] <- apply(fulldata.tab[, c(4:9)], 2, as.integer)
fulldata.tab[is.na(fulldata.tab)] <- 0

## If there are two types, but one of them is NM, then we assume that the type is the other, "more defined"

foo.tiago <- function(x){
    perc <- x/sum(x)
    if(any(perc == 1)){
        return(names(perc)[which(perc == 1)])
    } else if(perc["NM"] > 0 & perc["NM"] < 1){
        n.myc <- sum(perc > 0)
        if(n.myc == 2){
            temp <- perc[-which(names(perc) == "NM")]
            return(names(temp)[which(temp != 0)])
        } else {
            return("MIX")
        }
    } else {
        return("MIX")
    }
}


fulldata.tab$type <- apply(fulldata.tab[, 4:9], 1, foo.tiago)

## Summarizing data

genus.per.family <- aggregate(fulldata.tab$genus, by = list(fulldata.tab$family), FUN = function(x){length(unique(x))})
sp.per.family <- aggregate(fulldata.tab$species, by = list(fulldata.tab$family), FUN = function(x){length(unique(x))})

family.data <- setNames(sp.per.family, c("family", "samp.sp"))
family.data$samp.ge <- genus.per.family$x
family.type <- as.data.frame(matrix(NA, nrow = dim(family.data)[1], ncol = length(unique(names(unlist(aggregate(fulldata.tab$type, by = list(fulldata.tab$family), FUN = table)[,2])))) + 1))
names(family.type) <- c("family", sort(unique(names(unlist(aggregate(fulldata.tab$type, by = list(fulldata.tab$family), FUN = table)[,2])))))

for(i in 1:length(family.data$family)){
    print(i)
    tmp <- table(fulldata.tab$type[fulldata.tab$family == family.data$family[i]])
    family.type[i,] <- c(family.data$family[i], tmp[match(names(family.type)[2:7], names(tmp))])
    }
family.type[, c(2:7)] <- apply(family.type[, c(2:7)], 2, as.integer)
family.type[is.na(family.type)] <- 0

family.data <- cbind(family.data,
                     family.type[, 2:7]
                     )


## Checking if species richness values are correct
sum(family.data$samp.sp == apply(family.data[, 4:9], 1, sum)) == dim(family.data)[1]

foo.family <- function(x, thresh){
    perc <- x/sum(x)
    if(any(perc == 1)){
        return(names(perc)[which(perc == 1)])
    } else if(any(perc > thresh)){
        return(names(perc)[which(perc > thresh)])
    } else {
        return("MIX")
    }
}

## family.data$type.50 <- apply(family.data[, 4:9], 1, foo.family, thresh = 0.5)
## family.data$type.60 <- apply(family.data[, 4:9], 1, foo.family, thresh = 0.6)
family.data$type.80 <- apply(family.data[, 4:9], 1, foo.family, thresh = 0.8)
family.data$type.90 <- apply(family.data[, 4:9], 1, foo.family, thresh = 0.9)
family.data$type.100 <- apply(family.data[, 4:9], 1, foo.family, thresh = 1)


## Importing global richness per family
data.global <- read.csv("./data/data_all_families.csv", sep = ";")

family.data$global.rich <- data.global$nro_especies[match(family.data$family, data.global$familia)]

family.data$global.rich[is.na(family.data$global.rich)] <- 124

family.data$perc.rich <- family.data$samp.sp/family.data$global.rich

family.data <- cbind(family.data,
                     data.global[match(family.data$family, data.global$familia), c("crown.age", "stem.age")]
)

family.data <- family.data[!is.na(family.data$crown.age),]

write.table(family.data, file = here::here("data/family_data_full.csv"))
