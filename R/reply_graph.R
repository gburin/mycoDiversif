library("caper")
library("tidyverse")
library("cowplot")
library("RColorBrewer")
library("geiger")
library("phytools")
library("ggrepel")
library("plyr")
library("doMC")
library("patchwork")

#######################
### Analysis per genera
#######################

family.data.gen <- read.csv("../data/family_data_genus_classif.csv", stringsAsFactors = FALSE, row.names = NULL)
family.data.gen$family[family.data.gen$family == "Leguminosae"] <- "Fabaceae"
family.data.gen$family[family.data.gen$family == "Compositae"] <- "Asteraceae"

age.data <- read.csv("../data/data_all_families.csv", sep = ";")

## Importing tree
fulltree <- read.tree("../data/fam_tree_family_full.tre")
fulltree$node.label <- NULL
fulltree.vasc <- read.tree("../data/Vascular_Plants_rooted.dated.tre")

## Removing families with unknown mycorrhizal type
family.data.gen <- family.data.gen[-which(family.data.gen$UNK.raw.perc == 1),]

family.data.gen$stem.age <- age.data$stem.age[match(family.data.gen$family, age.data$familia)]
family.data.gen$r.e0 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0)
family.data.gen$r.e05 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.5)
family.data.gen$r.e09 <- bd.ms(time = family.data.gen$stem.age, n = family.data.gen$rich, crown = FALSE, epsilon = 0.9)

family.data.gen$shannon <- vegan::diversity(family.data.gen[, 3:7])

tree.pruned <- drop.tip(fulltree, tip = fulltree$tip.label[is.na(match(fulltree$tip.label, family.data.gen$family))])
data.pgls <- comparative.data(tree.pruned, family.data.gen[family.data.gen$MIX.raw.perc != 1,], names.col = "family")

phylosig.r0 <- pgls(r.e0 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r05 <- pgls(r.e05 ~ 1, data = data.pgls, lambda = "ML")
phylosig.r09 <- pgls(r.e09 ~ 1, data = data.pgls, lambda = "ML")
phylosig.rich <- pgls(rich ~ 1, data = data.pgls, lambda = "ML")
phylosig.age <- pgls(stem.age ~ 1, data = data.pgls, lambda = "ML")

nsim <- 10000

## sim.random <- llply(1:nsim,
##                     function(x){data.frame(
##                                     family = family.data.gen$family,
##                                     rich = family.data.gen$rich,
##                                     age = family.data.gen$stem.age,
##                                     AM = 0,
##                                     EM = 0,
##                                     ER = 0,
##                                     NM = 0,
##                                     OM = 0,
##                                     r.e0 = family.data.gen$r.e0,
##                                     r.e09 = family.data.gen$r.e09
##                                 )
##                     }
##                     )

## for(i in 1:10000){
##     print(i)
##     for(j in 1:nrow(sim.random[[i]])){
##         sampled.types <- table(sample(4:8, sim.random[[i]]$rich[j], replace = TRUE))
##         sim.random[[i]][j, as.integer(names(sampled.types))] <- sampled.types
##     }
##     sim.random[[i]]$shannon <- vegan::diversity(sim.random[[i]][, 4:8])
## }

## #registerDoMC(3)

## corr.random <- ldply(1:length(sim.random), function(x){
##     print(x)
##     comp.data <- comparative.data(tree.pruned, sim.random[[x]], names.col = "family")
##     pgls.rich <- caper::pgls(shannon ~ rich, data = comp.data, lambda = setNames(summary(phylosig.rich)$param[2], NULL))
##     pgls.r0 <- caper::pgls(shannon ~ r.e0, data = comp.data, lambda = setNames(summary(phylosig.rich)$param[2], NULL))
##     pgls.r09 <- caper::pgls(shannon ~ r.e09, data = comp.data, lambda = setNames(summary(phylosig.rich)$param[2], NULL))
##     lm.rich <- lm(shannon ~ rich, data = sim.random[[x]])
##     lm.r0 <- lm(shannon ~ r.e0, data = sim.random[[x]])
##     lm.r09 <- lm(shannon ~ r.e09, data = sim.random[[x]])
##     results <- data.frame(pgls.int.rich = coef(pgls.rich)[1],
##                           pgls.slope.rich = coef(pgls.rich)[2],
##                           pgls.r2.rich = summary(pgls.rich)$r.squared,
##                           pgls.pvalue.rich = summary(pgls.rich)$coefficients[2, 4],
##                           pgls.int.r0 = coef(pgls.r0)[1],
##                           pgls.slope.r0 = coef(pgls.r0)[2],
##                           pgls.r2.r0 = summary(pgls.r0)$r.squared,
##                           pgls.pvalue.r0 = summary(pgls.r0)$coefficients[2, 4],
##                           pgls.int.r09 = coef(pgls.r09)[1],
##                           pgls.slope.r09 = coef(pgls.r09)[2],
##                           pgls.r2.r09 = summary(pgls.r09)$r.squared,
##                           pgls.pvalue.r09 = summary(pgls.r09)$coefficients[2, 4],
##                           lm.int.rich = coef(lm.rich)[1],
##                           lm.slope.rich = coef(lm.rich)[2],
##                           lm.r2.rich = summary(lm.rich)$r.squared,
##                           lm.pvalue.rich = summary(lm.rich)$coefficients[2, 4],
##                           lm.int.r0 = coef(lm.r0)[1],
##                           lm.slope.r0 = coef(lm.r0)[2],
##                           lm.r2.r0 = summary(lm.r0)$r.squared,
##                           lm.pvalue.r0 = summary(lm.r0)$coefficients[2, 4],
##                           lm.int.r09 = coef(lm.r09)[1],
##                           lm.slope.r09 = coef(lm.r09)[2],
##                           lm.r2.r09 = summary(lm.r09)$r.squared,
##                           lm.pvalue.r09 = summary(lm.r09)$coefficients[2, 4])
## }, .progress = "text"#, .parallel = TRUE
## )

## save(sim.random, corr.random, file = "../output/simulations_reply_2020.RData")

load("../output/simulations_reply_2020.RData")

scatter.rich <- ggplot(sim.random[[1]]) +
    geom_point(aes(x = rich, y = shannon), alpha = 0.3) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness") +
    theme_cowplot()

scatter.r0 <- ggplot(sim.random[[1]]) +
    geom_point(aes(x = r.e0, y = shannon), alpha = 0.3) +
    labs(y = "Mycorrhizal State Diversity Index", x = expression("Diversification Rate ("~epsilon~" = 0 )")) +
    theme_cowplot()

scatter.r09 <- ggplot(sim.random[[1]]) +
    geom_point(aes(x = r.e09, y = shannon), alpha = 0.3) +
    labs(y = "Mycorrhizal State Diversity Index", x = expression("Diversification Rate ("~epsilon~" = 0.9 )")) +
    theme_cowplot()

scatter.rich + scatter.r0 + scatter.r09

ggsave(filename = "../output/scatter_rich_r0_r09_reply.pdf", width = 13, height = 6)

main.res <- read.csv("../output/fit_data_random_datasets.csv", as.is = TRUE)

rich.r2 <- data.frame(
    r2 = c(main.res$lm.r2.rich.sh, main.res$pgls.r2.rich.sh, corr.random$lm.r2.rich, corr.random$pgls.r2.rich),
    pvalue = c(main.res$lm.pvalue.rich.sh, main.res$pgls.pvalue.rich.sh, corr.random$lm.pvalue.rich, corr.random$pgls.pvalue.rich),
    type = rep(rep(c("lm", "PGLS"), each = 10000), 2),
    analysis = rep(c("main", "random"), each = 20000)
)
    
rich.r2.plot <-
    ggplot(rich.r2) +
    geom_histogram(aes(x = r2, fill = analysis, colour = analysis), alpha = 0.5) +
    labs(x = "R2", y = "Frequency", colour = "Analysis", fill = "Analysis") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(legend.position = "none") +
    facet_wrap(~ type)

rich.pvalue.plot <-
    ggplot(rich.r2) +
    geom_histogram(aes(x = pvalue, fill = analysis, colour = analysis), alpha = 0.5) +
    labs(x = "p-value", y = "Frequency", colour = "Analysis", fill = "Analysis") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    facet_wrap(~ type)

(rich.r2.plot / rich.pvalue.plot)

ggsave("../output/random_rich_shannon_r2_reply.pdf")

r0.r2 <- data.frame(
    r2 = c(main.res$lm.r2.r0, main.res$pgls.r2.r0, corr.random$lm.r2.r0, corr.random$pgls.r2.r0),
    pvalue = c(main.res$lm.pvalue.r0, main.res$pgls.pvalue.r0, corr.random$lm.pvalue.r0, corr.random$pgls.pvalue.r0),
    type = rep(rep(c("lm", "PGLS"), each = 10000), 2),
    analysis = rep(c("main", "random"), each = 20000)
)

r0.r2.plot <-
    ggplot(r0.r2) +
    geom_histogram(aes(x = r2, fill = analysis, colour = analysis), alpha = 0.5) +
    labs(x = "R2", y = "Frequency", colour = "Analysis", fill = "Analysis") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(legend.position = "none") +
    facet_wrap(~ type)

r0.pvalue.plot <-
    ggplot(r0.r2) +
    geom_histogram(aes(x = pvalue, fill = analysis, colour = analysis), alpha = 0.5) +
    labs(x = "p-value", y = "Frequency", colour = "Analysis", fill = "Analysis") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    facet_wrap(~ type)

(r0.r2.plot / r0.pvalue.plot)

ggsave("../output/random_r0_shannon_r2_reply.pdf")


r09.r2 <- data.frame(
    r2 = c(main.res$lm.r2.r09, main.res$pgls.r2.r09, corr.random$lm.r2.r09, corr.random$pgls.r2.r09),
    pvalue = c(main.res$lm.pvalue.r09, main.res$pgls.pvalue.r09, corr.random$lm.pvalue.r09, corr.random$pgls.pvalue.r09),
    type = rep(rep(c("lm", "PGLS"), each = 10000), 2),
    analysis = rep(c("main", "random"), each = 20000)
)

r09.r2.plot <-
    ggplot(r09.r2) +
    geom_histogram(aes(x = r2, fill = analysis, colour = analysis), alpha = 0.5) +
    labs(x = "R2", y = "Frequency", colour = "Analysis", fill = "Analysis") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(legend.position = "none") +
    facet_wrap(~ type)

r09.pvalue.plot <-
    ggplot(r09.r2) +
    geom_histogram(aes(x = pvalue, fill = analysis, colour = analysis), alpha = 0.5) +
    labs(x = "p-value", y = "Frequency", colour = "Analysis", fill = "Analysis") +
    scale_fill_brewer(palette = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(legend.position = "bottom") +
    facet_wrap(~ type)

(r09.r2.plot / r09.pvalue.plot)

ggsave("../output/random_r09_shannon_r2_reply.pdf")





### Increasing the number of categories

sim.random.50 <- data.frame(
    family = family.data.gen$family,
    rich = family.data.gen$rich,
    age = family.data.gen$stem.age,
    replicate(50, rep(0, nrow(family.data.gen))),
    r.e0 = family.data.gen$r.e0,
    r.e09 = family.data.gen$r.e09
)
                    

for(i in 1:nrow(sim.random.50)){
    print(i)
                                        #for(j in 1:nrow(sim.random[[i]])){
    sampled.types <- table(sample(4:53, sim.random.50$rich[i], replace = TRUE))
    sim.random.50[i, as.integer(names(sampled.types))] <- sampled.types
                                        #}
    sim.random.50$shannon <- vegan::diversity(sim.random.50[, 4:53])
}


sim.random.500 <- data.frame(
    family = family.data.gen$family,
    rich = family.data.gen$rich,
    age = family.data.gen$stem.age,
    replicate(500, rep(0, nrow(family.data.gen))),
    r.e0 = family.data.gen$r.e0,
    r.e09 = family.data.gen$r.e09
)
                    

for(i in 1:nrow(sim.random.500)){
    print(i)
                                        #for(j in 1:nrow(sim.random[[i]])){
    sampled.types <- table(sample(4:503, sim.random.500$rich[i], replace = TRUE))
    sim.random.500[i, as.integer(names(sampled.types))] <- sampled.types
                                        #}
    sim.random.500$shannon <- vegan::diversity(sim.random.500[, 4:503])
}


sim.random.5000 <- data.frame(
    family = family.data.gen$family,
    rich = family.data.gen$rich,
    age = family.data.gen$stem.age,
    replicate(5000, rep(0, nrow(family.data.gen))),
    r.e0 = family.data.gen$r.e0,
    r.e09 = family.data.gen$r.e09
)
                    

for(i in 1:nrow(sim.random.5000)){
    print(i)
                                        #for(j in 1:nrow(sim.random[[i]])){
    sampled.types <- table(sample(4:5003, sim.random.5000$rich[i], replace = TRUE))
    sim.random.5000[i, as.integer(names(sampled.types))] <- sampled.types
                                        #}
    sim.random.5000$shannon <- vegan::diversity(sim.random.5000[, 4:5003])
}


sim.random.5.plot <- ggplot(sim.random[[1]]) +
    geom_point(aes(x = rich, y = shannon), alpha = 0.3) +
    #geom_smooth(aes(x = rich, y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness", main = "5 Categories (Original)") +
    theme_cowplot()

sim.random.50.plot <- ggplot(sim.random.50) +
    geom_point(aes(x = rich, y = shannon), alpha = 0.3) +
    #geom_smooth(aes(x = rich, y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness", main = "50 Categories") +
    theme_cowplot()

sim.random.500.plot <- ggplot(sim.random.500) +
    geom_point(aes(x = rich, y = shannon), alpha = 0.3) +
    #geom_smooth(aes(x = rich, y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness", main = "500 Categories") +
    theme_cowplot()

sim.random.5000.plot <- ggplot(sim.random.5000) +
    geom_point(aes(x = rich, y = shannon), alpha = 0.3) +
    #geom_smooth(aes(x = rich, y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness", main = "5000 Categories") +
    theme_cowplot()


replica <- sample(1:length(sim.random), 4)

sim.random.5.plot.log1 <- ggplot(sim.random[[replica[1]]]) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", main = "5 Categories (Original)", title = "5 Categories (Original)") +
    theme_cowplot() +
    theme(axis.text.x = element_blank(), #axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          #axis.ticks.x = element_blank()#, axis.ticks.y = element_blank()
          )

sim.random.5.plot.log2 <- ggplot(sim.random[[replica[2]]]) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", main = "5 Categories (Original)") +
    theme_cowplot() +
    theme(axis.text.x = element_blank(), #axis.text.y = element_blank(),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          #axis.ticks.x = element_blank()#, axis.ticks.y = element_blank()
          )

sim.random.5.plot.log3 <- ggplot(sim.random[[replica[3]]]) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", main = "5 Categories (Original)") +
    theme_cowplot()

sim.random.5.plot.log4 <- ggplot(sim.random[[replica[4]]]) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", main = "5 Categories (Original)") +
    theme_cowplot() +
    theme(#axis.text.x = element_blank(), #axis.text.y = element_blank(),
        #axis.title.x = element_blank(),
        axis.title.y = element_blank()#,
          #axis.ticks.x = element_blank()#, axis.ticks.y = element_blank()
          )

sim.random.50.plot.log <- ggplot(sim.random.50) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", title = "50 Categories") +
    theme_cowplot()

sim.random.500.plot.log <- ggplot(sim.random.500) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", title = "500 Categories") +
    theme_cowplot()

sim.random.5000.plot.log <- ggplot(sim.random.5000) +
    geom_point(aes(x = log10(rich), y = shannon), alpha = 0.3) +
    scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = log10(rich), y = shannon), method = "lm", se = FALSE) +
    labs(y = "Mycorrhizal State Diversity Index", x = "Species Richness (log10)", title = "5000 Categories") +
    theme_cowplot()


(sim.random.5.plot + sim.random.50.plot) / (sim.random.500.plot + sim.random.5000.plot)

ggsave("../output/rich_shannon_categories.pdf")

panel1 <- (sim.random.5.plot.log1 + sim.random.5.plot.log2) / (sim.random.5.plot.log3 + sim.random.5.plot.log4)

(panel1 | sim.random.50.plot.log) / (sim.random.500.plot.log + sim.random.5000.plot.log)

ggsave("../output/rich_shannon_categories_log.pdf")




sim.random.5.plot.log.rate <- ggplot(sim.random[[1]]) +
    geom_point(aes(x = shannon, y = r.e09), alpha = 0.3) +
    #scale_y_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = shannon, y = r.e09), method = "lm", se = FALSE) +
    labs(x = "Mycorrhizal State Diversity Index", y = "Net Diversification (e = 0.9)", title = "5 Categories (Original)") +
    theme_cowplot()

sim.random.50.plot.log.rate <- ggplot(sim.random.50) +
    geom_point(aes(x = shannon, y = r.e09), alpha = 0.3) +
    #scale_y_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = shannon, y = r.e09), method = "lm", se = FALSE) +
    labs(x = "Mycorrhizal State Diversity Index", y = "Net Diversification (e = 0.9)", title = "50 Categories") +
    theme_cowplot()

sim.random.500.plot.log.rate <- ggplot(sim.random.500) +
    geom_point(aes(x = shannon, y = r.e09), alpha = 0.3) +
    #scale_y_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = shannon, y = r.e09), method = "lm", se = FALSE) +
    labs(x = "Mycorrhizal State Diversity Index", y = "Net Diversification (e = 0.9)", title = "500 Categories") +
    theme_cowplot()

sim.random.5000.plot.log.rate <- ggplot(sim.random.5000) +
    geom_point(aes(x = shannon, y = r.e09), alpha = 0.3) +
    #scale_y_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    #geom_smooth(aes(x = shannon, y = r.e09), method = "lm", se = FALSE) +
    labs(x = "Mycorrhizal State Diversity Index", y = "Net Diversification (e = 0.9)", title = "5000 Categories") +
    theme_cowplot()


(sim.random.5.plot.log.rate | sim.random.50.plot.log.rate) / (sim.random.500.plot.log.rate + sim.random.5000.plot.log.rate)

ggsave("../output/reply_shannon_rate_log.pdf")


family.data.gen$n.cat <- apply(family.data.gen[, c("AM", "EM", "NM", "OM", "ER")], 1, function(x){sum(x > 0)})


ggplot(family.data.gen) +
    geom_point(aes(x = factor(n.cat), y = log10(rich))) +
    geom_violin(aes(x = factor(n.cat), y = log10(rich)), colour = "lightgrey", fill = NA, size = 2) +
    scale_y_continuous(breaks = c(1, 2, 3, 4), labels = c(10, 100, 1000, 10000)) +
    labs(x = "Number of Mycorrhizal States", y = "Richness (log10)") +
    theme_cowplot()

ggsave("../output/richness_per_number_of_types.pdf")
