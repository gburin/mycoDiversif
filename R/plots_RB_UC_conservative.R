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

age.data <- read.csv("../data/data_all_families.csv", sep = ";")
RB.tree.UC.conservative <- read.nexus("../data/ramirez_barahona_data/ramirez_barahona_UC_conservative_MCCv_2.tre")

## Extracting families from tip information

RB.tip.list <- RB.tree.UC.conservative$tip.label
RB.family.list <- unique(setNames(sapply(RB.tip.list, function(x){strsplit(x, split = "_")[[1]][2]}), NULL))

## Replacing species name by the family, and dropping repeated tips
RB.tree.UC.conservative.pruned <- drop.tip(RB.tree.UC.conservative, tip = grep(RB.family.list[1], RB.tree.UC.conservative$tip.label)[-1])
for(i in 2:length(RB.family.list)){
    RB.tree.UC.conservative.pruned <- drop.tip(RB.tree.UC.conservative.pruned, tip = grep(RB.family.list[i], RB.tree.UC.conservative.pruned$tip.label)[-1])
}
RB.tree.UC.conservative.pruned$tip.label <- setNames(sapply(RB.tree.UC.conservative.pruned$tip.label, function(x){strsplit(x, split = "_")[[1]][2]}), NULL)

rb.UC.conservative.ages <- read.csv("../data/ramirez_barahona_data/ramirez_barahona_Ages_UC_conservative.csv", as.is = TRUE)

## Missing families from Ramírez-Barahona et al. 2020
age.data$familia[is.na(match(age.data$familia, rb.UC.conservative.ages$Family))]

## Adding ages from Ramírez-Barahona et al. 2020
age.data$rb.UC.conservative.stem <- rb.UC.conservative.ages$Stem_BEAST[match(age.data$familia, rb.UC.conservative.ages$Family)]

## Importing results from random datasets
load("../output/RB_results_UC_conservative.RData")
nrep <- 10000
fullresults <- ldply(1:nrep, function(x){results.UC.conservative[[x]][[1]]})

## Importing table with rates for plotting
family.data.gen <- read.csv("../output/simulated_datasets/random_data_09986.csv", stringsAsFactors = FALSE)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.perc == 1),]

family.data.gen$stem.age <- age.data$rb.UC.conservative.stem[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e05 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.5)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)
family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

## Calculating expected limits from all vascular plants

fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

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
    scale_y_log10(breaks = c(10, 1000, 10000)) +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal State")

ggsave(filename = "../output/figs/magsand_stem_labeled.pdf")

# Without labels
ggplot(data = limits.vasc.stem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = family.data.gen, aes(x = stem.age, y = rich, colour = type.60), size = 2.5) +
    scale_y_log10(breaks = c(10, 1000, 10000)) +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal State")

ggsave(filename = "../output/figs/magsand_stem_nolabel.pdf")


## Scatter - MTDI vs. Net Diversification

## Using replica 9986 (with R2 closer to the median) just to represent points

scatter.mtdi.r09 <-
    ggplot(fullresults) +
    geom_abline(mapping = aes(intercept = pgls.int.r09, slope = pgls.slope.r09, colour = brewer.pal(3, "Set1")[1]), show.legend = TRUE, alpha = 1) +
    geom_abline(mapping = aes(intercept = lm.int.r09, slope = lm.slope.r09, colour = brewer.pal(3, "Set1")[2]), show.legend = TRUE, alpha = 1) +
    geom_point(data = na.omit(family.data.gen), aes(x = shannon, y = r.e09), size = 2, alpha = 0.5) +
    labs(x = "Mycorrhizal State Diversity Index", y = "Diversification Rate", col = "Model Type", title = "A") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme_cowplot() +
    theme(legend.position = "top")

r2.lm.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.r2.r09), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.r2.r09), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = bquote('R' ^2), title = "B") +
    theme_cowplot()

pvalue.lm.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = log10(lm.pvalue.r09)), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = log10(median(fullresults$lm.pvalue.r09)), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = "p-value (log10)") +
    #scale_x_continuous(breaks = c(-14, -12, -10), labels = c(bquote('10' ^-14), bquote('10' ^-12), bquote('10' ^-10))) +
    theme_cowplot()

r2.pgls.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.r2.r09), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.r2.r09), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = bquote('R' ^2), title = "C")+
    theme_cowplot()

pvalue.pgls.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = log10(pgls.pvalue.r09)), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = log10(median(fullresults$pgls.pvalue.r09)), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = "p-value (log10)") +
    #scale_x_continuous(breaks = c(-12, -10, -8), labels = c(bquote('10' ^-12), bquote('10' ^-10), bquote('10' ^-8))) +
    theme_cowplot()


scatter.mtdi.r09 | ((r2.lm.r09 + pvalue.lm.r09) / (r2.pgls.r09 + pvalue.pgls.r09))

ggsave(filename = "../output/figs/scatterplots_lm_pgls_stem_r09_UC_conservative.pdf", height = 5.5, width = 13)

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
    geom_histogram(aes(x = pgls.r2.age.sh), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.r2.age.sh), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(x = bquote('R' ^2), y = "Frequency") +
    theme_cowplot()

pvalue.pgls.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = pgls.pvalue.age.sh), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$pgls.pvalue.age.sh), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = "p-value") +
    theme_cowplot()

r2.lm.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.r2.age.sh), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.r2.age.sh), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = bquote('R' ^2)) +
    theme_cowplot()

pvalue.lm.age.sh <-
    ggplot(fullresults) +
    geom_histogram(aes(x = lm.pvalue.age.sh), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    geom_vline(xintercept = median(fullresults$lm.pvalue.age.sh), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = "p-value") +
    theme_cowplot()

scatter.age.sh | ((r2.pgls.age.sh + pvalue.pgls.age.sh) / (r2.lm.age.sh + pvalue.lm.age.sh))

ggsave(filename = "../output/figs/scatterplots_lm_pgls_stem_age_UC_conservative.pdf", height = 5.5, width = 13)


## Richness vs Shannon
scatter.rich.sh <-
    ggplot(fullresults) +
    geom_abline(mapping = aes(intercept = pgls.int.rich.sh, slope = pgls.slope.rich.sh), colour = brewer.pal(3, "Set1")[1], show.legend = TRUE, alpha = 0.01) +
    geom_abline(mapping = aes(intercept = lm.int.rich.sh, slope = lm.slope.rich.sh), colour = brewer.pal(3, "Set1")[2], show.legend = TRUE, alpha = 0.01) +
    geom_point(data = na.omit(family.data.gen), aes(x = shannon, y = rich), size = 2, alpha = 0.5) +
    labs(x = "Mycorrhizal State Shannon Index", y = "Species Richness per Family", col = "Model Type") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    #coord_trans(y = "log10") +
    #scale_y_continuous(breaks = c(10, 50, 100, 500, 1000, 10000, 20000, 30000)) +
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

ggsave(filename = "../output/figs/scatterplots_lm_pgls_stem_rich_UC_conservative.pdf", height = 5.5, width = 13)


### Boxplots
## We aggregate mean rates per type for each random dataset

age.data$r.e0 <- geiger::bd.ms(time = age.data$rb.UC.conservative.stem, n = age.data$nro_especies, crown = FALSE, epsilon = 0)
age.data$r.e09 <- geiger::bd.ms(time = age.data$rb.UC.conservative.stem, n = age.data$nro_especies, crown = FALSE, epsilon = 0.9)

random.data <- lapply(paste0("../output/simulated_datasets/", list.files("../output/simulated_datasets/")), read.csv, stringsAsFactors = FALSE)

random.types <- data.frame(
    family = age.data$familia,
    r.e0 = age.data$r.e0,
    r.e09 = age.data$r.e09
)

for(i in 1:length(random.data)){
    random.types <- cbind(random.types,
                          random.data[[i]]$type.60
                          )
}

## Removing families that are 100% UNK
random.types <- random.types[-which(random.types[,4] == "UNK"),]

anova.data.r0 <- plyr::ldply(random.types[4:10003], function(x){aggregate(random.types$r.e0, by = list(x), FUN = mean)})
anova.data.r09 <- plyr::ldply(random.types[4:10003], function(x){aggregate(random.types$r.e09, by = list(x), FUN = mean)})


## box.r09 <-
##     ggplot(data = anova.data.r09[is.na(match(anova.data.r09$Group.1, c("OM", "ER"))),]) +
##     geom_boxplot(aes(x = factor(Group.1, levels = c("AM", "EM", "NM", "MIX")), y = x, colour = factor(Group.1, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, outlier.alpha = 1, outlier.shape = 3, size = 1) +
##     geom_point(aes(x = factor(Group.1, levels = c("AM", "EM", "NM", "MIX")), y = x, colour = factor(Group.1, levels = c("AM", "EM", "NM", "MIX"))), alpha = 0.05, position = position_jitterdodge(jitter.width = 2.2, dodge.width = 1)) +
##     stat_summary(aes(x = factor(Group.1, levels = c("AM", "EM", "NM", "MIX")), y = x, colour = factor(Group.1, levels = c("AM", "EM", "NM", "MIX"))), fun.y = mean, geom = "point", colour = "black", size = 2) +
##     theme_cowplot() +
##     theme(legend.position = "none") +
##     scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
##     labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))

box.r09 <-
    ggplot(data = family.data.gen[is.na(match(family.data.gen$type.60, c("OM", "ER"))),]) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX"))), fill = NA, size = 1.5, outlier.shape = 3, alpha = 0.5) +
    geom_jitter(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), y = r.e09, colour = factor(type.60, levels = c("AM", "EM", "NM", "MIX")), size = shannon), alpha = 0.5) +
    labs(x = "Mycorrhizal State", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"), title = "A") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[-4]) +
    theme_cowplot() +
    theme(legend.position = "none")

ggsave(filename = "../output/figs/boxplots_anova_myctype_r09_UC_conservative.pdf", box.r09, width = 11, height = 5.5, units = "in")

pvalue.phyanova.r09 <-
    ggplot(fullresults) +
    geom_histogram(aes(x = phyaov.pvalue.r09.60), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 50) +
    #geom_histogram(aes(x = aov.pvalue.r09.60, y = ..density..), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 100) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = "p-value", title = "B") +
    xlim(-0.005, 0.2) + 
    theme_cowplot()

pvalue.anova.r09 <-
    ggplot(fullresults) +
    #geom_histogram(aes(x = phyaov.pvalue.r09.60, y = ..density..), fill = brewer.pal(3, "Set1")[1], colour = NA, alpha = 0.5, bins = 100) +
    geom_histogram(aes(x = aov.pvalue.r09.60), fill = brewer.pal(3, "Set1")[2], colour = NA, alpha = 0.5, bins = 50) +
    #geom_vline(xintercept = median(fullresults$phyaov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[1], linetype = "dashed", size = 1.5) +
    #geom_vline(xintercept = median(fullresults$aov.pvalue.r09.60), colour = brewer.pal(3, "Set1")[2], linetype = "dashed", size = 1.5) +
    geom_vline(xintercept = 0.05, colour = "black", linetype = "dashed", size = 1.5) +
    labs(y = "Frequency", x = "p-value", title = "C") +
    #xlim(-0.005, 0.2) + 
    theme_cowplot()

(box.r09) / (pvalue.phyanova.r09 + pvalue.anova.r09)

ggsave(filename = "../output/figs/box_hist_anova_myctype_r09_UC_conservative.pdf", width = 11, height = 9, units = "in")


anovas.r09.pvalue <-
    ggplot(fullresults) +
    geom_jitter(aes(x = factor("phyANOVA"), y = phyaov.pvalue.r09.60), alpha = 0.1, colour = brewer.pal(3, "Set1")[1], width = 0.2) +
    geom_jitter(aes(x = factor("ANOVA"), y = aov.pvalue.r09.60), alpha = 0.1, colour = brewer.pal(3, "Set1")[2], width = 0.2) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    labs(x = "Model", y = "p-value") +
    theme_cowplot()

ggsave(filename = "../output/figs/jitter_pvalue_anovas_r09_UC_conservative.pdf", anovas.r09.pvalue, width = 11, height = 5.5, units = "in")


