library("ape")
library("tidyverse")
library("RColorBrewer")
library("geiger")
library("phytools")
library("cowplot")
library("ggrepel")
library("caper")
library("patchwork")
library("plyr")
library("reshape2")
load("../output/anova_posthoc_tables.RData")

## Importing results from random datasets
fullresults <- read.csv("../output/fit_data_random_datasets.csv")

matrix.mutator <- function(x){
    temp <- matrix(NA, 4, 4, dimnames = list(c("AM", "EM", "MIX", "NM"), c("AM", "EM", "MIX", "NM")))
    diag(temp) <- 1
    temp[lower.tri(temp)] <- round(x[,4], 3)
    temp[upper.tri(temp)] <- t(temp)[upper.tri(temp)]
    return(temp)
}

formatted.results <- lapply(main.results, function(y){
    res <- vector(mode = "list", length = 16)
    res[1:8] <- y[1:8]
    res[9:16] <- lapply(y[9:16], matrix.mutator)
    names(res) <- names(y)
    return(res)
})
        
for(i in 1:16){
    temp <- ldply(formatted.results, function(x){x[[i]][upper.tri(x[[i]], diag = TRUE)] <- NA; return(x[[i]][, -4])})
    assign(names(main.results[[1]])[i], melt(data.frame(
                                            temp,
                                            type.2 = rep(c("AM", "EM", "MIX", "NM"), 10000))))
}

r0.50 <- na.omit(data.frame(rbind(phyaov.r0.50, aov.r0.50), class = rep(c("phy", "lm"), each = 120000)))
r0.60 <- na.omit(data.frame(rbind(phyaov.r0.60, aov.r0.60), class = rep(c("phy", "lm"), each = 120000)))
r0.80 <- na.omit(data.frame(rbind(phyaov.r0.80, aov.r0.80), class = rep(c("phy", "lm"), each = 120000)))
r0.100 <- na.omit(data.frame(rbind(phyaov.r0.100, aov.r0.100), class = rep(c("phy", "lm"), each = 120000)))

r09.50 <- na.omit(data.frame(rbind(phyaov.r09.50, aov.r09.50), class = rep(c("phy", "lm"), each = 120000)))
r09.60 <- na.omit(data.frame(rbind(phyaov.r09.60, aov.r09.60), class = rep(c("phy", "lm"), each = 120000)))
r09.80 <- na.omit(data.frame(rbind(phyaov.r09.80, aov.r09.80), class = rep(c("phy", "lm"), each = 120000)))
r09.100 <- na.omit(data.frame(rbind(phyaov.r09.100, aov.r09.100), class = rep(c("phy", "lm"), each = 120000)))


posthoc.r0.50 <-
    ggplot() +
    geom_histogram(data = subset(r0.50, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r0.50, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
    #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable, drop = TRUE) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))


posthoc.r0.60 <-
    ggplot() +
    geom_histogram(data = subset(r0.60, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r0.60, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
    #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))


posthoc.r0.80 <-
    ggplot() +
    geom_histogram(data = subset(r0.80, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r0.80, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
    #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

posthoc.r0.100 <-
    ggplot() +
    geom_histogram(data = subset(r0.100, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r0.100, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
    #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))




posthoc.r09.50 <-
    ggplot() +
    geom_histogram(data = subset(r09.50, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r09.50, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
                                        #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))


posthoc.r09.60 <-
    ggplot() +
    geom_histogram(data = subset(r09.60, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r09.60, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
                                        #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

posthoc.r09.80 <-
    ggplot() +
    geom_histogram(data = subset(r09.80, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r09.80, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
                                        #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

posthoc.r09.100 <-
    ggplot() +
    geom_histogram(data = subset(r09.100, class == "lm"), aes(x = value, y = ..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(data = subset(r09.100, class == "phy"), aes(x = value, y = -..density.., fill = class), colour = NA, alpha = 0.5, bins = 50) +
    geom_vline(xintercept = 0.05, linetype = "dashed") +
    #geom_hline(yintercept = 0, size = 0.5) +
                                        #coord_flip() +
    scale_colour_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    scale_fill_brewer(palette = "Set1", direction = -1, labels = c("Linear", "Phylogenetic")) +
    labs(x = "p-value", y = "Density", colour = "Model Type", fill = "Model Type", title = "C") +
    facet_grid(type.2 ~ variable) +
    theme_cowplot() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))




## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Importing table with rates for plotting
family.data.gen <- read.csv("./output/simulated_datasets/random_data_09986.csv", stringsAsFactors = FALSE)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

age.data <- read.csv("./data/data_all_families.csv", sep = ";")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

family.data.gen$stem.age <- age.data$stem.age[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e05 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.5)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)
family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

## Calculating expected limits from all vascular plants
stem.age.vasc <- max(branching.times(fulltree.vasc)) - findMRCA(fulltree.vasc, c(fulltree.vasc$tip.label[621:622]), type = "height")

r.vasc.stem.r0 <- bd.ms(time = stem.age.vasc, n = sum(family.data.gen$rich), crown = FALSE, epsilon = 0)
r.vasc.stem.r05 <- bd.ms(time = stem.age.vasc, n = sum(family.data.gen$rich), crown = FALSE, epsilon = 0.5)
r.vasc.stem.r0.9 <- bd.ms(time = stem.age.vasc, n = sum(family.data.gen$rich), crown = FALSE, epsilon = 0.9)





box.r0.50 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.50, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r0.50 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.50), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.50), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r0.50 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.50), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.50), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))



box.r0.60 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.60, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r0.60 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.60), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.60), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r0.60 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.60), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.60), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))



box.r0.80 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.80, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r0.80 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.80), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.80), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r0.80 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.80), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.80), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))



box.r0.100 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.100, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r0.100 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.100), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.100), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r0.100 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r0.100), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r0.100), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))







box.r09.50 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.50, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r09.50 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r09.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.50), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.50), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r09.50 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r09.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r09.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.50), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.50), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))



box.r09.60 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.60, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r09.60 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r09.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r09.60 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r09.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r09.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))



box.r09.80 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.80, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r09.80 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r09.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.80), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.80), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r09.80 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r09.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r09.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.80), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.80), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))



box.r09.100 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.100, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")


pvalue.phyanova.r09.100 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r09.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.100), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.100), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value", title = "B") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))

pvalue.anova.r09.100 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r09.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    geom_histogram(aes(x = aov.pvalue.r09.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.100), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.100), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 10))


save.image("../output/anova_boxplots.RData")


(box.r0.50) + (((pvalue.phyanova.r0.50 + pvalue.anova.r0.50) / posthoc.r0.50) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))

(box.r0.60) + (((pvalue.phyanova.r0.60 + pvalue.anova.r0.60) / posthoc.r0.60) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))

(box.r0.80) + (((pvalue.phyanova.r0.60 + pvalue.anova.r0.60) / posthoc.r0.80) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))

(box.r0.100) + (((pvalue.phyanova.r0.60 + pvalue.anova.r0.60) / posthoc.r0.100) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))




(box.r09.50) + (((pvalue.phyanova.r09.50 + pvalue.anova.r09.50) / posthoc.r09.50) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))

(box.r09.60) + (((pvalue.phyanova.r09.60 + pvalue.anova.r09.60) / posthoc.r09.60) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))

## ggsave(filename = "./output/figs/anova_results_main.pdf", width = 8, height = 4)

(box.r09.80) + (((pvalue.phyanova.r09.60 + pvalue.anova.r09.60) / posthoc.r09.80) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))

(box.r09.100) + (((pvalue.phyanova.r09.60 + pvalue.anova.r09.60) / posthoc.r09.100) + plot_layout(heights = c(3, 7))) +
    plot_layout(widths = c(6, 4))
