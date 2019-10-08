source("./R/supplementary_analysis.R")

## Stem age
ggplot(data = limits.vasc.stem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = subset(family.data, study == "yes"), aes(x = stem.age, y = global.rich, size = study, colour = type.60), size = 2.5) +
    geom_point(data = subset(family.data, study == "no"), aes(x = stem.age, y = global.rich, size = study), size = 1, colour = "black", alpha = 0.15) +
    geom_text_repel(data = subset(family.data, study == "yes"), aes(x = stem.age, y = global.rich, label = family, colour = type.60)) +
    scale_y_log10() +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    #scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Dark2"))(nrow(fulldata))) +
    theme(legend.position = "bottom") +
    #scale_alpha_manual(values = c(0.3, 1)) +
    #scale_size_manual(values = c(1, 2.5)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal Type")

ggsave(filename = "./output/supp_figs/magsand_stem_labeled_sp.pdf")

# Without labels
ggplot(data = limits.vasc.stem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = subset(family.data, study == "yes"), aes(x = stem.age, y = global.rich, size = study, colour = type.60), size = 2.5) +
    geom_point(data = subset(family.data, study == "no"), aes(x = stem.age, y = global.rich, size = study), size = 1, colour = "black", alpha = 0.15) +
    #geom_text_repel(data = subset(family.data, study == "yes"), aes(x = stem.age, y = global.rich, label = family, colour = type.60)) +
    scale_y_log10() +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    #scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Dark2"))(nrow(fulldata))) +
    theme(legend.position = "bottom") +
    #scale_alpha_manual(values = c(0.3, 1)) +
    #scale_size_manual(values = c(1, 2.5)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal Type")

ggsave(filename = "./output/supp_figs/magsand_stem_nolabel_sp.pdf")


## Simple plot - Raw data
raw.r0 <-
    ggplot(data = na.omit(family.data.clean)) +
    geom_point(aes(x = shannon, y = r.e0.stem), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r0)[1], sl = coef(lm.r0)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r0)[1], sl = coef(mod.r0)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, linetype = "dashed", show.legend = TRUE) +
    labs(x = "Mycorrhizae Type Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    theme(legend.position = "top") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"))

raw.r09 <-
    ggplot(data = na.omit(family.data.clean)) +
    geom_point(aes(x = shannon, y = r.e09.stem), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r09)[1], sl = coef(lm.r09)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r09)[1], sl = coef(mod.r09)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, linetype = "dashed", show.legend = TRUE) +
    labs(x = "Mycorrhizae Type Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    theme(legend.position = "top") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"))


## Individual plots
stem.age.sh <-
    ggplot(data = na.omit(family.data.clean)) +
    geom_point(aes(x = shannon, y = stem.age), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.stem.age.sh)[1], sl = coef(lm.stem.age.sh)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(pgls.stem.age.sh)[1], sl = coef(pgls.stem.age.sh)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, ,colour = col), linetype = "dashed", size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Family Age", col = element_blank()) +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS")) +
    theme(legend.position = "top")

rich.sh <-
    ggplot(data = na.omit(family.data.clean)) +
    geom_point(aes(x = shannon, y = global.rich), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.rich.sh)[1], sl = coef(lm.rich.sh)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(pgls.rich.sh)[1], sl = coef(pgls.rich.sh)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, ,colour = col), linetype = "dashed", size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Species Richness", col = element_blank()) +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS")) +
    theme(legend.position = "top")


plot_grid(raw.r0,
          stem.age.sh,
          raw.r09,
          rich.sh,
          ncol = 2,
          align = 'hv',
          axis = 'tlbr',
          labels = letters[1:4]
          )

ggsave(filename = "./output/supp_figs/scatterplots_lm_pgls_stem_sp.pdf")



### Boxplots

family.data.clean$type.50 <- as.character(family.data.clean$type.50)
family.data.clean$type.60 <- as.character(family.data.clean$type.60)
family.data.clean$type.80 <- as.character(family.data.clean$type.80)
family.data.clean$type.100 <- as.character(family.data.clean$type.100)

box.r0 <-
    ggplot(data = na.omit(family.data.clean)) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e0.stem, colour = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e0.stem), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal Type", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

box.r09 <-
    ggplot(data = na.omit(family.data.clean)) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e09.stem, colour = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e09.stem), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")[c(1:3, 5)]) +
    labs(x = "Mycorrhizal Type", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))


plot_grid(box.r0,
          box.r09,
          nrow = 2,
          align = 'v',
          labels = c("a", "b")
          )

ggsave(filename = "./output/supp_figs/boxplots_netdiv_myctype_sp.pdf", width = 11, height = 7, units = "in")




### Analysis including species with remarks "Probably XX"

## Stem age
ggplot(data = limits.vasc.stem.rem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = subset(family.data.rem, study == "yes"), aes(x = stem.age, y = global.rich, size = study, colour = type.60), size = 2.5) +
    geom_point(data = subset(family.data.rem, study == "no"), aes(x = stem.age, y = global.rich, size = study), size = 1, colour = "black", alpha = 0.15) +
    geom_text_repel(data = subset(family.data.rem, study == "yes"), aes(x = stem.age, y = global.rich, label = family, colour = type.60)) +
    scale_y_log10() +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    #scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Dark2"))(nrow(fulldata))) +
    theme(legend.position = "bottom") +
    #scale_alpha_manual(values = c(0.3, 1)) +
    #scale_size_manual(values = c(1, 2.5)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal Type")

ggsave(filename = "./output/supp_figs/magsand_stem_labeled_sp_rem.pdf")

# Without labels
ggplot(data = limits.vasc.stem.rem) +
    geom_line(aes(x = age, y = lb.0)) +
    geom_line(aes(x = age, y = ub.0)) +
    geom_line(aes(x = age, y = lb.09), linetype = "dashed") +
    geom_line(aes(x = age, y = ub.09), linetype = "dashed") +
    geom_point(data = subset(family.data.rem, study == "yes"), aes(x = stem.age, y = global.rich, size = study, colour = type.60), size = 2.5) +
    geom_point(data = subset(family.data.rem, study == "no"), aes(x = stem.age, y = global.rich, size = study), size = 1, colour = "black", alpha = 0.15) +
    #geom_text_repel(data = subset(family.data, study == "yes"), aes(x = stem.age, y = global.rich, label = family, colour = type.60)) +
    scale_y_log10() +
    labs(x = "Age of Clade (MY)", y = "Number of Species") +
    #scale_colour_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Dark2"))(nrow(fulldata))) +
    theme(legend.position = "bottom") +
    #scale_alpha_manual(values = c(0.3, 1)) +
    #scale_size_manual(values = c(1, 2.5)) +
    scale_colour_brewer(palette = "Dark2") +
    labs(colour = "Mycorrhizal Type")

ggsave(filename = "./output/supp_figs/magsand_stem_nolabel_sp_rem.pdf")


## Simple plot - Raw data
raw.r0.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = r.e0.stem), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r0.rem)[1], sl = coef(lm.r0.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r0.rem)[1], sl = coef(mod.r0.rem)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, linetype = "dashed", show.legend = TRUE) +
    labs(x = "Mycorrhizae Type Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    theme(legend.position = "top") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"))

raw.r09.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = r.e09.stem), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.r09.rem)[1], sl = coef(lm.r09.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, show.legend = TRUE) +
    geom_abline(data = data.frame(int = coef(mod.r09.rem)[1], sl = coef(mod.r09.rem)[2], col = "b"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5, linetype = "dashed", show.legend = TRUE) +
    labs(x = "Mycorrhizae Type Shannon Index", y = "Diversification Rate", col = "Model Type") +
    #xlim(0, 1) +
    theme(legend.position = "top") +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS"))


## Individual plots
stem.age.sh.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = stem.age), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.stem.age.sh.rem)[1], sl = coef(lm.stem.age.sh.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(pgls.stem.age.sh.rem)[1], sl = coef(pgls.stem.age.sh.rem)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, ,colour = col), linetype = "dashed", size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Family Age", col = element_blank()) +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS")) +
    theme(legend.position = "top")

rich.sh.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = shannon, y = global.rich), size = 2) +
    geom_abline(data = data.frame(int = coef(lm.rich.sh.rem)[1], sl = coef(lm.rich.sh.rem)[2], col = "a"), mapping = aes(intercept = int, slope = sl, colour = col), size = 1.5) +
    geom_abline(data = data.frame(int = coef(pgls.rich.sh.rem)[1], sl = coef(pgls.rich.sh.rem)[2], col = "blue"), mapping = aes(intercept = int, slope = sl, ,colour = col), linetype = "dashed", size = 1.5) +
    labs(x = "Mycorrhizal Type Diversity Index", y = "Species Richness", col = element_blank()) +
    #xlim(0, 1) +
    scale_colour_brewer(palette = "Set1", labels = c("Linear Model", "PGLS")) +
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

ggsave(filename = "./output/supp_figs/scatterplots_lm_pgls_stem_sp_rem.pdf")



### Boxplots

family.data.clean.rem$type.50 <- as.character(family.data.clean.rem$type.50)
family.data.clean.rem$type.60 <- as.character(family.data.clean.rem$type.60)
family.data.clean.rem$type.80 <- as.character(family.data.clean.rem$type.80)
family.data.clean.rem$type.100 <- as.character(family.data.clean.rem$type.100)

box.r0.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e0.stem, colour = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e0.stem), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")) +
    labs(x = "Mycorrhizal Type", y = expression("Diversification Rate ("~epsilon~" = 0 )"))

box.r09.rem <-
    ggplot(data = na.omit(family.data.clean.rem)) +
    geom_point(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e09.stem, colour = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), size = shannon), alpha = 0.75, position = position_jitterdodge(jitter.width = 0.75, dodge.width = 1)) +
    geom_boxplot(aes(x = factor(type.60, levels = c("AM", "EM", "NM", "ER", "MIX")), y = r.e09.stem), fill = NA, colour = "darkgrey", outlier.alpha = 1) +
    theme(legend.position = "none") +
    scale_colour_manual(values = brewer.pal(5, "Set1")) +
    labs(x = "Mycorrhizal Type", y = expression("Diversification Rate ("~epsilon~" = 0.9 )"))


plot_grid(box.r0.rem,
          box.r09.rem,
          nrow = 2,
          align = 'v',
          labels = c("a", "b")
          )

ggsave(filename = "./output/supp_figs/boxplots_netdiv_myctype_sp_rem.pdf", width = 11, height = 7, units = "in")
