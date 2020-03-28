library("patchwork")
library("ggplot2")
library("cowplot")
library("RColorBrewer")
library("geiger")
library("phytools")

fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

## Main analyses - epsilon = 0
fullresults <- read.csv("../output/fit_data_random_datasets.csv")

## Importing table with rates for plotting
family.data.gen <- read.csv("../output/simulated_datasets/random_data_09986.csv", stringsAsFactors = FALSE)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

age.data <- read.csv("../data/data_all_families.csv", sep = ";")

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
r.vasc.stem.r0.9 <- bd.ms(time = stem.age.vasc, n = sum(family.data.gen$rich), crown = FALSE, epsilon = 0.9)

limits.vasc.stem <- data.frame(
    time = seq(1, 300, by = 1),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0, epsilon = 0)),
    #t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r05, epsilon = 0.5)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0.9, epsilon = 0.9))
    )
names(limits.vasc.stem) <- c("age", "lb.0", "ub.0", "lb.09", "ub.09")

save.image("../output/background_rates_plots.RData")

## Simple plot - Raw data
scatter.mtdi.r0 <-
    ggplot(fullresults) +
    geom_abline(mapping = aes(intercept = pgls.int.r0, slope = pgls.slope.r0, colour = brewer.pal(3, "Set1")[1]), show.legend = TRUE, alpha = 1) +
    geom_abline(mapping = aes(intercept = lm.int.r0, slope = lm.slope.r0, colour = brewer.pal(3, "Set1")[2]), show.legend = TRUE, alpha = 1) +
    geom_point(data = na.omit(family.data.gen), aes(x = shannon, y = r.e0), size = 2, alpha = 0.5) +
    labs(x = "Mycorrhizal State Diversity Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme_cowplot() +
    theme(legend.position = "top")

r2.pgls.r0 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.r2.r0, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.r2.r0), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2))+
    theme_cowplot()

pvalue.pgls.r0 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = log10(pgls.pvalue.r0), y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = log10(median(fullresults$pgls.pvalue.r0)), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    scale_x_continuous(breaks = c(-13, -12, -11, -10, -9), labels = c(bquote('10' ^-13), bquote('10' ^-12), bquote('10' ^-11), bquote('10' ^-10), bquote('10' ^-9))) +
    theme_cowplot()

r2.lm.r0 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.r2.r0, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.r2.r0), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2)) +
    theme_cowplot()

pvalue.lm.r0 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = log10(lm.pvalue.r0), y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = log10(median(fullresults$lm.pvalue.r0)), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    scale_x_continuous(breaks = c(-16, -14, -12, -10), labels = c(bquote('10' ^-16), bquote('10' ^-14), bquote('10' ^-12), bquote('10' ^-10))) +
    theme_cowplot()

scatter.mtdi.r0 | ((r2.lm.r0 + pvalue.lm.r0) / (r2.pgls.r0 + pvalue.pgls.r0))

save.image("../output/scatter_r0.RData")


### Boxplots

rm(list = ls())

fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

## Main analyses - epsilon = 0
fullresults <- read.csv("../output/fit_data_random_datasets.csv")

## Importing table with rates for plotting
family.data.gen <- read.csv("../output/simulated_datasets/random_data_09986.csv", stringsAsFactors = FALSE)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

age.data <- read.csv("../data/data_all_families.csv", sep = ";")

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
r.vasc.stem.r0.9 <- bd.ms(time = stem.age.vasc, n = sum(family.data.gen$rich), crown = FALSE, epsilon = 0.9)

## Epsilon = 0 for 50% threshold
box.r0.50 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.50, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

pvalue.anovas.r0.50 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r0.pvalue.50 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r0.50 + (pvalue.anovas.r0.50 / anovas.r0.pvalue.50)

## Epsilon = 0.9 for 50% threshold
box.r09.50 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.50, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.50, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.50, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))

pvalue.anovas.r09.50 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r09.pvalue.50 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.50, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r09.50, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r09.50 + (pvalue.anovas.r09.50 / anovas.r09.pvalue.50)




## Epsilon = 0 for 60% threshold
box.r0.60 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.60, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

pvalue.anovas.r0.60 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r0.pvalue.60 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r0.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r0.60 + (pvalue.anovas.r0.60 / anovas.r0.pvalue.60)


## Epsilon = 0 for 80% threshold
box.r0.80 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.80, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

pvalue.anovas.r0.80 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r0.pvalue.80 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r0.80 + (pvalue.anovas.r0.80 / anovas.r0.pvalue.80)

## Epsilon = 0.9 for 80% threshold
box.r09.80 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.80, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.80, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.80, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))

pvalue.anovas.r09.80 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r09.pvalue.80 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.80, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r09.80, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r09.80 + (pvalue.anovas.r09.80 / anovas.r09.pvalue.80)


## Epsilon = 0 for 100% threshold
box.r0.100 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.100, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

pvalue.anovas.r0.100 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r0.pvalue.100 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r0.100 + (pvalue.anovas.r0.100 / anovas.r0.pvalue.100)

## Epsilon = 0.9 for 100% threshold
box.r09.100 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.100, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
    geom_jitter(aes(x = factor(type.100, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.100, levels = c("AM", "EM", "NM", "MIX"))), size = 3, alpha = 0.5) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))

pvalue.anovas.r09.100 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    #geom_histogram(aes(x = aov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

anovas.r09.pvalue.100 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r0.100, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r09.100, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

box.r09.100 + (pvalue.anovas.r09.100 / anovas.r09.pvalue.100)

save.image("../output/boxplots_supp.RData")



### Analysis including species with remarks "Probably XX"

rm(list = ls())
load("../output/suppdata_species_withremarks.RData")

## Simple plot - Raw data
raw.r0.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = r.e0.stem), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r0.rem)[1], sl = coef(lm.r0.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r0.rem)[1], sl = coef(mod.r0.rem)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    labs(x = "Mycorrhizae Type Diversity Index", y = "Diversification Rate", col = "Model Type") +
    theme_cowplot() +
    theme(legend.position = "top") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1)

raw.r09.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = r.e09.stem), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r09.rem)[1], sl = coef(lm.r09.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r09.rem)[1], sl = coef(mod.r09.rem)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    labs(x = "Mycorrhizae Type Diversity Index", y = "Diversification Rate", col = "Model Type") +
    theme_cowplot() +
    theme(legend.position = "top") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1)


## Individual plots
stem.age.sh.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = stem.age), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.stem.age.sh.rem)[1], sl = coef(lm.stem.age.sh.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(pgls.stem.age.sh.rem)[1], sl = coef(pgls.stem.age.sh.rem)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, ,colour = col), size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Family Age", col = element_blank()) +
    theme_cowplot() +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme(legend.position = "top")

rich.sh.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = global.rich), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.rich.sh.rem)[1], sl = coef(lm.rich.sh.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(pgls.rich.sh.rem)[1], sl = coef(pgls.rich.sh.rem)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, ,colour = col), size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Species Richness", col = element_blank()) +
    theme_cowplot() +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme(legend.position = "top")


plot_grid(raw.r0.rem,
          stem.age.sh.rem,
          raw.r09.rem,
          rich.sh.rem,
          ncol = 2,
          align = 'hv',
          axis = 'tlbr',
          labels = letters[1:4]
          )

### Boxplots

family.data.clean.rem$type.50 <- as.character(family.data.clean.rem$type.50)
family.data.clean.rem$type.60 <- as.character(family.data.clean.rem$type.60)
family.data.clean.rem$type.80 <- as.character(family.data.clean.rem$type.80)
family.data.clean.rem$type.100 <- as.character(family.data.clean.rem$type.100)

box.r0.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e0.stem, colour = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e0.stem), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")) +
    labs(x = "Mycorrhizal Type", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

box.r09.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e09.stem, colour = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e09.stem), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")) +
    labs(x = "Mycorrhizal Type", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))


plot_grid(box.r0.rem,
          box.r09.rem,
          nrow = 2,
          align = 'v',
          labels = c("a", "b")
          )

save.image("../output/scatter_boxplots_remarks.RData")



## ggsave(filename = "../output/supp_figs/boxplots_netdiv_myctype_sp_rem.pdf", width = 11, height = 7, units = "in")


### 20% analysis

rm(list = ls())
load("../output/results_20perc.RData")

## r_epsilon = 0
hist.pgls.r2.r0 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = pgls.r2.r0), fill = "red", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(pgls.r2.r0)), linetype = "dashed", size = 2, colour = "red") +
    theme_cowplot()

hist.pgls.pvalue.r0 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = log10(pgls.pvalue.r0)), fill = "red", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(log10(pgls.pvalue.r0 + 1e-16))), linetype = "dashed", size = 2, colour = "red") +
    scale_x_continuous(breaks = c(-16, -14, -12, -10, -8, -6), labels = c(bquote('10' ^-16), bquote('10' ^-14), bquote('10' ^-12), bquote('10' ^-10), bquote('10' ^-8), bquote('10' ^-6))) +
    labs(x = "p-value", y = element_blank()) +
    theme_cowplot()

hist.lm.r2.r0 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = lm.r2.r0), fill = "blue", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(lm.r2.r0)), linetype = "dashed", size = 2, colour = "blue") +
    theme_cowplot()

hist.lm.pvalue.r0 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = log10(lm.pvalue.r0)), fill = "blue", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(log10(lm.pvalue.r0))), linetype = "dashed", size = 2, colour = "blue") +
    scale_x_continuous(breaks = c(-12, -10, -8, -6, -4), labels = c(bquote('10' ^-12), bquote('10' ^-10), bquote('10' ^-8), bquote('10' ^-6), bquote('10' ^-4))) +
    labs(x = "p-value", y = element_blank()) +
    theme_cowplot()



hist.pgls.slope.r0 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = pgls.slope.r0), fill = "red", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(pgls.slope.r0)), linetype = "dashed", size = 2, colour = "red") +
    labs(x = bquote("R", ^2), y = element_blank()) +
    theme_cowplot()

hist.lm.slope.r0 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = lm.slope.r0), fill = "blue", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(lm.slope.r0)), linetype = "dashed", size = 2, colour = "blue") +
    labs(x = "slope", y = element_blank()) +
    theme_cowplot()


(hist.lm.pvalue.r0 + hist.lm.r2.r0 + hist.lm.slope.r0) / (hist.pgls.pvalue.r0 + hist.pgls.r2.r0 + hist.pgls.slope.r0)



## r_epslion = 0.9
hist.pgls.r2.r09 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = pgls.r2.r09), fill = "red", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(pgls.r2.r09)), linetype = "dashed", size = 2, colour = "red") +
    labs(x = bquote("R", ^2), y = element_blank()) +
    theme_cowplot()

hist.pgls.pvalue.r09 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = log10(pgls.pvalue.r09)), fill = "red", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(log10(pgls.pvalue.r09 + 1e-16))), linetype = "dashed", size = 2, colour = "red") +
    scale_x_continuous(breaks = c(-12.5, -10, -7.5, -5), labels = c(bquote('10' ^-12.5), bquote('10' ^-10), bquote('10' ^-7.5), bquote('10' ^-5))) +
    labs(x = "p-value", y = element_blank()) +
    theme_cowplot()

hist.lm.r2.r09 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = lm.r2.r09), fill = "blue", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(lm.r2.r09)), linetype = "dashed", size = 2, colour = "blue") +
    labs(x = bquote("R", ^2), y = element_blank()) +
    theme_cowplot()

hist.lm.pvalue.r09 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = log10(lm.pvalue.r09)), fill = "blue", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(log10(lm.pvalue.r09))), linetype = "dashed", size = 2, colour = "blue") +
    scale_x_continuous(breaks = c(-12, -10, -8, -6, -4), labels = c(bquote('10' ^-12), bquote('10' ^-10), bquote('10' ^-8), bquote('10' ^-6), bquote('10' ^-4))) +
    labs(x = "p-value", y = element_blank()) +
    theme_cowplot()



hist.pgls.slope.r09 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = pgls.slope.r09), fill = "red", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(pgls.slope.r09)), linetype = "dashed", size = 2, colour = "red") +
    labs(x = "slope", y = element_blank()) +
    theme_cowplot()

hist.lm.slope.r09 <-
    ggplot(data = random.results) +
    geom_histogram(aes(x = lm.slope.r09), fill = "blue", colour = NA, alpha = 0.3, bins = 100) +
    geom_vline(aes(xintercept = mean(lm.slope.r09)), linetype = "dashed", size = 2, colour = "blue") +
    labs(x = "slope", y = element_blank()) +
    theme_cowplot()


(hist.lm.pvalue.r09 + hist.lm.r2.r09 + hist.lm.slope.r09) / (hist.pgls.pvalue.r09 + hist.pgls.r2.r09 + hist.pgls.slope.r09)

save.image("../output/hist_20perc.RData")
