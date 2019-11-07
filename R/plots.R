library("ape")
library("tidyverse")
library("RColorBrewer")
library("geiger")
library("phytools")
library("cowplot")
load("./output/main_results.RData")

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
    geom_text_repel(data = family.data.gen, aes(x = stem.age, y = rich, label = family, colour = type.60.valid)) +
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
mtdi.r0 <-
    ggplot(data = na.omit(family.data.gen)) +
    geom_point(aes(x = shannon, y = r.e0), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r0)[1], sl = coef(lm.r0)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, linetype = "dashed", show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r0)[1], sl = coef(mod.r0)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    labs(x = "Mycorrhizal Type Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme_cowplot() +
    theme(legend.position = "top")

mtdi.r09 <-
    ggplot(data = na.omit(family.data.gen)) +
    geom_point(aes(x = shannon, y = r.e0), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r09)[1], sl = coef(lm.r09)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, linetype = "dashed", show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r09)[1], sl = coef(mod.r09)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    labs(x = "Mycorrhizal Type Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"), direction = -1) +
    theme_cowplot() +
    theme(legend.position = "top")


stem.age.sh <-
    ggplot(data = na.omit(family.data.gen)) +
    geom_point(aes(x = shannon, y = stem.age), size = 2) +
    geom_abline(data = data.frame(int = coef(pgls.age.sh)[1], sl = coef(pgls.age.sh)[2], col = "a"), mapping = aes(intercept = int, slope = sl, ,colour = col),  size = 1.5) +
    geom_abline(data = data.frame(int = coef(lm.age.sh)[1], sl = coef(lm.age.sh)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, colour = col), linetype = "dashed", size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Family Age", colour = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("PGLS", "Linear Model"))+
    theme_cowplot() +
    theme(legend.position = "top") 

rich.sh <-
    ggplot(data = na.omit(family.data.gen)) +
    geom_point(aes(x = shannon, y = rich), size = 2) +
    geom_abline(data = data.frame(int = coef(pgls.rich.sh)[1], sl = coef(pgls.rich.sh)[2], col = "a"), mapping = aes(intercept = int, slope = sl, ,colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(lm.rich.sh)[1], sl = coef(lm.rich.sh)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, colour = col), linetype = "dashed", size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Species Richness", colour = "Model Type") +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("PGLS", "Linear Model")) +
    theme_cowplot() +
    theme(legend.position = "top")


plot_grid(mtdi.r0,
          stem.age.sh,
          mtdi.r09,
          rich.sh,
          ncol = 2,
          align = 'hv',
          axis = 'tlbr',
          labels = letters[1:4]
          )

ggsave(filename = "./output/figs/scatterplots_lm_pgls_stem.pdf", width = 11, height = 7)



### Boxplots

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


plot_grid(box.r0,
          box.r09,
          nrow = 2,
          align = 'v',
          labels = c("a", "b")
          )

ggsave(filename = "./output/figs/boxplots_netdiv_myctype.pdf", width = 11, height = 7, units = "in")

