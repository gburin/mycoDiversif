library("tidyverse")
library("ape")

family.data.gen <- read.csv("./data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"
family.data.gen$perc.sp.unk <- family.data.gen$UNK/family.data.gen$rich

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]


expanded.data <- data.frame(
    family = rep(family.data.gen$family, times = family.data.gen$rich),
    orig.type = rep(as.character(family.data.gen$type.60.valid), times = family.data.gen$rich),
    stringsAsFactors = FALSE
)

sp.to.change <- sort(sample(1:nrow(expanded.data), ceiling(0.2 * nrow(expanded.data))))

expanded.data$new.type <- NA
expanded.data$new.type[-sp.to.change] <- as.character(expanded.data$orig.type[-sp.to.change])

for(i in 1:length(sp.to.change)){
    expanded.data$new.type[sp.to.change[i]] <- sample(unique(expanded.data$orig.type)[-which(unique(expanded.data$orig.type) == expanded.data$orig.type[sp.to.change[i]])], 1)
}

