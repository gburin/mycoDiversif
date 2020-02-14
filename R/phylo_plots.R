library("ape")
library("phytools")
library("ggtree")
library("ggnewscale")
library("viridis")
library("RColorBrewer")
library("ggforce")

### Tree for the Genus-level analysis


tree <- read.tree("./data/fam_tree_family_full.tre")
age.data <- read.csv("./data/data_all_families.csv", sep = ";")
tree <- drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)])
fulldata <- read.csv("./data/family_data_genus_classif.csv", stringsAsFactors = FALSE)
rownames(fulldata) <- fulldata$family
fulldata$family[fulldata$family == "Compositae"] <- "Asteraceae"
fulldata$family[fulldata$family == "Leguminosae"] <- "Fabaceae"

cleandata <- fulldata[which(!is.na(match(fulldata$family, tree$tip.label))),]
#cleandata <- cleandata[cleandata$MIX.raw.perc != 1,]
cleandata$stem.age <- age.data$stem.age[match(cleandata$family, age.data$familia)]
cleandata$stem.age[cleandata$family == "Cornaceae"] <- tree$edge.length[2 * which(tree$tip.label == "Cornaceae") + 1]
cleandata$stem.age[cleandata$family == "Schlegeliaceae"] <- tree$edge.length[2 * which(tree$tip.label == "Schlegeliaceae") + 1]

cleandata$r.e0 <- geiger::bd.ms(time = cleandata$stem.age, n = cleandata$rich, crown = FALSE, epsilon = 0)
cleandata$r.e09 <- geiger::bd.ms(time = cleandata$stem.age, n = cleandata$rich, crown = FALSE, epsilon = 0.9)

df.rates <- cleandata[, c("r.e0", "r.e09")]
rownames(df.rates)[rownames(df.rates) == "Leguminosae"] <- "Fabaceae"
rownames(df.rates)[rownames(df.rates) == "Compositae"] <- "Asteraceae"
df.myc <- setNames(cleandata[, c("AM.raw.perc", "EM.raw.perc", "NM.raw.perc", "MIX.raw.perc")], c("AM", "EM", "NM", "MIX"))
rownames(df.myc)[rownames(df.myc) == "Leguminosae"] <- "Fabaceae"
rownames(df.myc)[rownames(df.myc) == "Compositae"] <- "Asteraceae"
df.shan <- data.frame(MTDI = as.numeric(vegan::diversity(cleandata[, 3:7])))
rownames(df.shan) <- rownames(df.myc)

df.rates <- df.rates[which(!is.na(match(rownames(df.rates), rownames(df.myc)))),]
df.shan <- data.frame(MTDI = df.shan[which(!is.na(match(rownames(df.shan), rownames(df.myc)))),], row.names = rownames(df.myc))


tree <- drop.tip(tree, tree$tip.label[which(is.na(match(tree$tip.label, rownames(df.rates))))])

p <-
    ggtree(tr = tree, layout = "fan", open.angle = 10) +
    geom_tiplab2(offset = 270, size = 3) +
    geom_treescale(x = 0, width = 50) +
    geom_treescale(x = 100, width = 50) +
    geom_treescale(x = 200, width = 50) +
    geom_treescale(x = 300, width = 50)


p1 <-
    gheatmap(p = p, data = df.rates, offset = 0.8, width = 0.2, colnames_angle = 90, color = "lightgrey", hjust = 0.5) +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[1], "00"), high = paste0(brewer.pal(3, "Dark2")[1], "FF"), name = "Diversification\nRate")

p2 <- p1 + new_scale_fill()

p3 <-
    gheatmap(p = p2, data = df.myc, offset = 90, width = 0.3, colnames_angle = 90, color = "lightgrey") +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[2], "00"), high = paste0(brewer.pal(3, "Dark2")[2], "FF"), name = "Mycorrhizal Type")

p4 <- p3 + new_scale_fill()

p5 <- gheatmap(p = p4, data = df.shan, offset = 205, width = 0.1, colnames_angle = 90, color = "lightgrey") +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[3], "00"), high = paste0(brewer.pal(3, "Dark2")[3], "FF"), name = "Mycorrhizal Type Diversity Index")

ggsave("./output/tree_fig/phylo_shannon_myco_genus.pdf", p5, width = 20, height = 20)
