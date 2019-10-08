library("ape")
library("phytools")
library("ggtree")
library("ggnewscale")
library("viridis")
library("RColorBrewer")
library("ggforce")

tree <- read.tree("family_tree_175.tre")
fulldata <- read.csv("fulldata_mujica.csv")
rownames(fulldata) <- fulldata$familia

cleandata <- fulldata[match(tree$tip.label, rownames(fulldata)),]

df.rates <- cleandata[, c(11, 13)]
names(df.rates) <- c("r0","r0.9")
df.myc <- cleandata[, c(14:16)]
df.myc <- df.myc/apply(df.myc, 1, sum)
df.shan <- data.frame(MTDI = as.numeric(cleandata$Shannon))
rownames(df.shan) <- rownames(df.myc)


p <-
    ggtree(tr = tree, layout = "fan", open.angle = 10) +
    geom_tiplab2(offset = 270) +
    geom_treescale(x = 0, width = 50) +
    geom_treescale(x = 100, width = 50) +
    geom_treescale(x = 200, width = 50) +
    geom_treescale(x = 300, width = 50)


p1 <-
    gheatmap(p = p, data = df.rates, offset = 0.8, width = 0.2, colnames_angle = 90, colnames_offset_y = 178, color = "black", hjust = 0.5) +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[1], "25"), high = paste0(brewer.pal(3, "Dark2")[1], "FF"), name = "Diversification\nRate")

p2 <- p1 + new_scale_fill()

p3 <-
    gheatmap(p = p2, data = df.myc, offset = 90, width = 0.3, colnames_angle = 90, colnames_offset_y = 178, color = "black") +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[2], "25"), high = paste0(brewer.pal(3, "Dark2")[2], "FF"))

p4 <- p3 + new_scale_fill()
    gheatmap(p = p4, data = df.shan, offset = 210, width = 0.1, colnames_angle = 90, colnames_offset_y = 178, color = "black") +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[3], "25"), high = paste0(brewer.pal(3, "Dark2")[3], "FF"))






### Tree for the Genus-level analysis


tree <- read.tree("./data/fam_tree_family_full.tre")
age.data <- read.csv("./data/data_all_families.csv", sep = ";")
tree <- drop.tip(tree, tree$tip.label[duplicated(tree$tip.label)])
fulldata <- read.csv("./data/family_data_genus_classif.csv")
rownames(fulldata) <- fulldata$family

cleandata <- fulldata[which(!is.na(match(rownames(fulldata), tree$tip.label))),]
cleandata$stem.age <- age.data$stem.age[match(cleandata$family, age.data$familia)]
cleandata$stem.age[cleandata$family == "Cornaceae"] <- tree$edge.length[2 * which(tree$tip.label == "Cornaceae") + 1]
cleandata$stem.age[cleandata$family == "Schlegeliaceae"] <- tree$edge.length[2 * which(tree$tip.label == "Schlegeliaceae") + 1]

cleandata$r.e0 <- geiger::bd.ms(time = cleandata$stem.age, n = cleandata$rich, crown = FALSE, epsilon = 0)
cleandata$r.e09 <- geiger::bd.ms(time = cleandata$stem.age, n = cleandata$rich, crown = FALSE, epsilon = 0.9)

df.rates <- cleandata[, c("r.e0", "r.e09")]
df.myc <- setNames(cleandata[, c(22, 23, 26)], c("AM", "EM", "NM"))
df.shan <- data.frame(MTDI = as.numeric(vegan::diversity(cleandata[, c(8, 9, 10, 12, 13)])))
rownames(df.shan) <- rownames(df.myc)

df.myc <- df.myc[-which(is.na(df.myc$AM)),]
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
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[1], "25"), high = paste0(brewer.pal(3, "Dark2")[1], "FF"), name = "Diversification\nRate")

p2 <- p1 + new_scale_fill()

p3 <-
    gheatmap(p = p2, data = df.myc, offset = 90, width = 0.3, colnames_angle = 90, color = "lightgrey") +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[2], "25"), high = paste0(brewer.pal(3, "Dark2")[2], "FF"), name = "Mycorrhizal Type")

p4 <- p3 + new_scale_fill()


p5 <- gheatmap(p = p4, data = df.shan, offset = 210, width = 0.1, colnames_angle = 90, color = "lightgrey") +
    scale_fill_gradient(low = paste0(brewer.pal(3, "Dark2")[3], "25"), high = paste0(brewer.pal(3, "Dark2")[3], "FF"), name = "Mycorrhizal Type Diversity Index")

ggsave("./output/tree_fig/phylo_shannon_myco_genus.pdf", p5, width = 20, height = 20)
