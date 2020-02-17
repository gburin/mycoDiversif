library("ape")
library("tidyverse")
library("RColorBrewer")
library("geiger")
library("phytools")
library("cowplot")
library("ggrepel")
library("caper")
library("patchwork")

## Importing tree
fulltree <- read.tree("./data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("./data/Vascular_Plants_rooted.dated.tre")

## Importing results from random datasets
fullresults <- read.csv("./output/fit_data_random_datasets.csv")

## Importing table with rates for plotting
family.data.gen <- read.csv("./data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

age.data <- read.csv("./data/data_all_families.csv", sep = ";")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.raw.perc == 1),]

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

limits.vasc.stem <- data.frame(
    time = seq(1, 300, by = 1),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0, epsilon = 0)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r05, epsilon = 0.5)),
    t(mapply(stem.limits, seq(1, 300, by = 1), r = r.vasc.stem.r0.9, epsilon = 0.9))
    )
names(limits.vasc.stem) <- c("age", "lb.0", "ub.0", "lb.05", "ub.05", "lb.09", "ub.09")

## Stem age
ggplot(data = limits.vasc.stem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = family.data.gen, aes(x = stem.age, y = rich, colour = type.60), size = 2.5) +
    geom_text_repel(data = family.data.gen, aes(x = stem.age, y = rich, label = family, colour = type.60)) +
    scale_y_log10() +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal State")

ggsave(filename = "./output/figs/magsand_stem_labeled.pdf")

# Without labels
ggplot(data = limits.vasc.stem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = family.data.gen[family.data.gen$type.60 != "ER",], aes(x = stem.age, y = rich, colour = type.60), size = 2.5) +
    scale_y_log10() +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal State")

ggsave(filename = "./output/figs/magsand_stem_nolabel.pdf")


## Scatter - MTDI vs. Net Diversification

## Using 50/50 MIX division just to represent points

scatter.mtdi.r09 <-
    ggplot(fullresults) +
    geom_abline(mapping = aes(intercept = pgls.int.r09, slope = pgls.slope.r09), colour = brewer.pal(3, "Set1")[1], show.legend = TRUE, alpha = 0.01) +
    geom_abline(mapping = aes(intercept = lm.int.r09, slope = lm.slope.r09), colour = brewer.pal(3, "Set1")[2], show.legend = TRUE, alpha = 0.01) +
    geom_point(data = na.omit(family.data.gen), aes(x = shannon, y = r.e09), size = 2, alpha = 0.5) +
    labs(x = "Mycorrhizal State Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme_cowplot() +
    theme(legend.position = "top")

r2.pgls.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.r2.r09, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.r2.r09), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2))+
    theme_cowplot()

pvalue.pgls.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.pvalue.r09, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.pvalue.r09), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

r2.lm.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.r2.r09, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.r2.r09), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2)) +
    theme_cowplot()

pvalue.lm.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.pvalue.r09, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.pvalue.r09), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

scatter.mtdi.r09 | ((r2.lm.r09 + pvalue.lm.r09) / (r2.pgls.r09 + pvalue.pgls.r09))

## Age vs Shannon
scatter.age.sh <-
    ggplot(fullresults) +
    geom_abline(mapping = aes(intercept = pgls.int.age.sh, slope = pgls.slope.age.sh), colour = brewer.pal(3, "Set1")[1], show.legend = TRUE, alpha = 0.01) +
    geom_abline(mapping = aes(intercept = lm.int.age.sh, slope = lm.slope.age.sh), colour = brewer.pal(3, "Set1")[2], show.legend = TRUE, alpha = 0.01) +
    geom_point(data = na.omit(family.data.gen), aes(x = shannon, y = stem.age), size = 2, alpha = 0.5) +
    labs(x = "Mycorrhizal State Shannon Index", y = "Stem Age", col = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme_cowplot() +
    theme(legend.position = "top")

r2.pgls.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.r2.age.sh, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.r2.age.sh), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2))+
    theme_cowplot()

pvalue.pgls.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.pvalue.age.sh, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.pvalue.age.sh), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

r2.lm.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.r2.age.sh, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.r2.age.sh), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2)) +
    theme_cowplot()

pvalue.lm.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.pvalue.age.sh, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.pvalue.age.sh), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

scatter.age.sh | ((r2.pgls.age.sh + pvalue.pgls.age.sh) / (r2.lm.age.sh + pvalue.lm.age.sh))


## Richness vs Shannon
scatter.rich.sh <-
    ggplot(fullresults) +
    geom_abline(mapping = aes(intercept = pgls.int.rich.sh, slope = pgls.slope.rich.sh), colour = brewer.pal(3, "Set1")[1], show.legend = TRUE, alpha = 0.01) +
    geom_abline(mapping = aes(intercept = lm.int.rich.sh, slope = lm.slope.rich.sh), colour = brewer.pal(3, "Set1")[2], show.legend = TRUE, alpha = 0.01) +
    geom_point(data = na.omit(family.data.gen), aes(x = shannon, y = rich), size = 2, alpha = 0.5) +
    labs(x = "Mycorrhizal State Shannon Index", y = "Stem Rich", col = "Model Type") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    coord_trans(y = "log10") +
    scale_y_continuous(breaks = c(10, 50, 100, 500, 1000, 10000, 20000, 30000)) +
    theme_cowplot() +
    theme(legend.position = "top")

r2.pgls.rich.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.r2.rich.sh, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.r2.rich.sh), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2))+
    theme_cowplot()

pvalue.pgls.rich.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.pvalue.rich.sh, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.pvalue.rich.sh), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

r2.lm.rich.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.r2.rich.sh, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.r2.rich.sh), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = bquote('R' ^2)) +
    theme_cowplot()

pvalue.lm.rich.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.pvalue.rich.sh, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.pvalue.rich.sh), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = element_blank(), x = "p-value") +
    theme_cowplot()

scatter.rich.sh | ((r2.pgls.rich.sh + pvalue.pgls.rich.sh) / (r2.lm.rich.sh + pvalue.lm.rich.sh))

ggsave(filename = "./output/figs/scatterplots_lm_pgls_stem.pdf", width = 5.5, height = 13)



### Boxplots
## We aggregate mean rates per type for each random dataset

family.data.gen$type.50 <- as.character(family.data.gen$type.50)
family.data.gen$type.60 <- as.character(family.data.gen$type.60)
family.data.gen$type.80 <- as.character(family.data.gen$type.80)
family.data.gen$type.100 <- as.character(family.data.gen$type.100)
family.data.gen <- family.data.gen[family.data.gen$type.60 != "ER",]

box.r0 <-
    ggplot(data = family.data.gen[-which(is.na(family.data.gen$r.e0)),]) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e0, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e0), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

box.r09 <-
    ggplot(data = family.data.gen[-which(is.na(family.data.gen$r.e09)),]) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme_cowplot() +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))

ggsave(filename = "./output/figs/boxplots_netdiv_myctype_r09.pdf", width = 11, height = 5.5, units = "in")


plot_grid(box.r0,
          box.r09,
          nrow = 2,
          align = 'v',
          labels = c("a", "b")
          )

ggsave(filename = "./output/figs/boxplots_netdiv_myctype.pdf", width = 11, height = 7, units = "in")

