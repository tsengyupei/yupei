#Loading library----
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(vegan)
library(phyloseq)
library(randomcoloR)

#Loading data----
#relative abundance based on traits
bugbase <- read.table("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/predictions.txt")
bugbase$sample_ID <- rownames(bugbase)

#metadata
metadata <- read.table("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/bugbase_metadata.txt", header = T)
colnames(metadata)[1] <- "sample_ID"
metadata <- metadata %>% mutate(sample_ID = sapply(sample_ID, function(x){str_replace(x, "-", "_")})) 
metadata <- metadata %>% mutate(sample_ID = sapply(sample_ID, function(x){str_replace(x, " _", "_")})) %>% as.data.frame()

rownames(metadata) <- metadata$sample_ID
#adding metadata info to abundance data
bugbase <- bugbase %>% cbind(metadata[-69, c(1, 2, 6, 21, 25)]) %>% select(-11)
bugbase_long <- bugbase %>% pivot_longer(Aerobic:Stress_Tolerant, names_to = "Trait", values_to = "Abundance")

#species composition table
spe <- read.table("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/16s_normalized_otus.txt") %>% as.matrix() %>% t() 
rownames(spe) <- sub("X", "", rownames(spe))
spe <- spe[-(66:68), ]

#species with trait information
trait <- read.table("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/contributing_otus.txt")

#species composition greengene
spe_gg <- read.delim("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/greengene_feature_table.txt", header = T, row.names = 1)
tax_gg <- spe_gg[, "taxonomy"] %>% as.data.frame() %>% separate(1, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
rownames(tax_gg) <- rownames(spe_gg)
tax_gg <- tax_gg %>% as.matrix() %>% tax_table()
spe_gg <- spe_gg %>% select(-"taxonomy") 
spe_gg <- spe_gg %>% rename_with(~ mgsub::mgsub(colnames(spe_gg), c("X", "_V6V8", "_"), c("", "", "-"))) %>% otu_table(taxa_are_rows = TRUE)

samdf_gg <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv', header = T) %>% select(-1)
samples.out <- as.data.frame(colnames(spe_gg)) %>% rename("sample_ID" = 1)
samdf_gg <- left_join(samples.out, samdf_gg, by = "sample_ID")
rownames(samdf_gg) <- colnames(spe_gg)


ps <- phyloseq(otu_table(spe_gg, taxa_are_rows = T), sample_data(samdf_gg), tax_table(tax_gg))
ps <- subset_taxa(ps, Kingdom != "Unassigned") %>% 
  subset_taxa(Class != "D_2__Chloroplast") %>% 
  subset_taxa(Family != "D_4__Mitochondria")
ps <- prune_samples(!(sample_names(ps)  %in% c("NC1", "NC2", "NC3", "NC4")), ps) # Remove NC sample

#draw ps rarefaction Curve (raremax = 5235) w/o NA
col <- randomColor(count = 70)
raremax <- min(sample_sums(ps)) #get the minimum point

tab <- otu_table(ps) %>% t()
class(tab) <- "matrix"
rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=F, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_ps <- rarefy_even_depth(ps, rngseed = 1, sample.size = raremax, replace = F)

#rarfied relative abundance/ richness----
trait$OTU <- rownames(trait)
relative_abundance <- function(char = char){
  bugbase_long_gg <- rarefied_ps %>% 
    transform_sample_counts(function(x) {ifelse(x>0, 1, 0)}) %>%
    psmelt() %>% left_join(trait, by = "OTU") %>% select("Sample", char, "Abundance") 
  colnames(bugbase_long_gg)[2] <- "Trait"
  bugbase_long_gg <- bugbase_long_gg %>% group_by(Sample, Trait) %>% summarise(ab = sum(Abundance)) %>% filter(Trait == T) 
  bugbase_long_gg <- bugbase_long_gg %>% as.data.frame( )%>% mutate(Trait = rep(char, length(bugbase_long_gg$Trait)))
  return(bugbase_long_gg)
}
colnames(trait)
relative_abundance_mat <- relative_abundance("Aerobic")
for (i in colnames(trait)[-c(1, 10)]){
  relative_abundance_mat <- rbind(relative_abundance_mat, relative_abundance(i))
}
colnames(relative_abundance_mat)[c(1, 3)] <- c("sample_ID", "Abundance")
relative_abundance_mat <- relative_abundance_mat %>% left_join(samdf, by = "sample_ID")

#Regression
ggscatter(relative_abundance_mat %>% filter(is.na(fern.size) == F & fern.size != "S"), x = 'C_N', y = 'Abundance',
          add = "reg.line", conf.int = TRUE, color = "#DDE5B6", xlab = 'C:N ratio', ylab = 'OTU richness')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size = 7)+
  geom_point(cex = 2, aes( col = layer))+
  scale_colour_manual(values = c("#ADC178", "#DDE5B6", "#997B66", "#F1DCA7"), aesthetics = c("colour", "fill"))+
  theme(text = element_text(size = 15))+
  guides(col=guide_legend(title="Layer"))+
  facet_wrap(~Trait)

#Boxplot with relative abundance depend on layer for each trait----

##saperated plots
# box_m <- function(size = size, trait = trait){
#   stat.test <- bugbase_long %>% filter(fern.size != size & Trait == trait) %>% wilcox_test(Abundance~layer) %>% add_xy_position(x = "layer")
#   p <- ggboxplot(bugbase_long %>% filter(fern.size != size & Trait == trait), x = "layer", y = "Abundance", color = "layer", palette = "jco")+ 
#     stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)+
#     ggtitle(trait)
#   print(p)
# }
# list_lm <- list()
# for(i in colnames(bugbase)[1:9]){
#   list_lm <- append(list, box_m("S", i)) 
# }

box <- function(size = size){
  plotList <- lapply(
    colnames(bugbase)[1:9],
    function(trait) {
      stat.test <- relative_abundance_mat %>% filter(fern.size == size & Trait == trait) %>% wilcox_test(Abundance~layer) %>% add_xy_position(x = "layer")
      p <- ggboxplot(relative_abundance_mat %>% filter(fern.size == size & Trait == trait), x = "layer", y = "Abundance", color = "layer", palette = "jco")+ 
        stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)+
        ggtitle(trait)
      p
    }
  )
  ggarrange(plotlist=plotList, ncol = 3, nrow = 3)
}
box("S")
box("M")
box("L")

#Bar plot for ralative abundance----
bugbase_long %>% 
  filter(is.na(fern.size) == F) %>% 
  filter(Trait == "Aerobic" | Trait == "Anaerobic" | Trait == "Facultatively_Anaerobic") %>%
  ggplot(aes(x = sample_ID, y = Abundance, fill = Trait)) +
  geom_bar(stat="identity")+
  facet_wrap(fern.size ~ layer)

#NMDS----
#transform matrix to relative abundance
spe_lm <- decostand(rarefied_ps@otu_table %>% as.data.frame() %>% t() %>% as.data.frame() %>% filter(colnames(rarefied_ps@otu_table) %in% rownames(samdf_gg)[which(samdf_gg$fern.size == "L" & samdf_gg$layer != "O")]) %>% as.matrix(), method = "total")
NMDS <- metaMDS(spe_lm, noshare = TRUE, autotransform = FALSE, trymax = 500)
# plot(NMDS)
ordisurf(NMDS, samdf_gg[which(samdf_gg$fern.size == "L"& samdf_gg$layer != "O"), "pH"], main="",col="forestgreen")

points(NMDS, display = "site", 
       col = c("#ADC178", "#DDE5B6", "#F1DCA7"), 
       pch = 15, cex = 2)
ordihull(
  NMDS,
  samdf_gg[which(samdf_gg$fern.size == "L" & samdf_gg$layer != "O"), "layer"],
  display = "sites",
  draw = c("polygon"),
  col = NULL,
  border = c("#ADC178", "#DDE5B6", "#F1DCA7"),
  lwd = 2.5
)

points(NMDS, display = "species", 
       col = ifelse(row.names(NMDS$species) %in% rownames(trait)[which(trait$Aerobic == T)], "darkgreen", ifelse(row.names(NMDS$species) %in% rownames(trait)[which(trait$Anaerobic == T)], "skyblue", ifelse(row.names(NMDS$species) %in% rownames(trait)[which(trait$Facultatively_Anaerobic == T)], "red3", "gray"))), 
       pch = "*", cex = 2)

legend(x="topright", legend = unique(samdf_gg[which(samdf_gg$fern.size == "L"& samdf_gg$layer != "O"), "layer"]), col = c("#ADC178", "#DDE5B6", "#F1DCA7"), pch = 15, cex = 0.7, text.font = 20, bg = NULL, box.lty = 0)
legend(x="topleft", legend = c("Aerobic", "Anaerobic", "Facultative aerobic", "Others"), col = c("darkgreen", "skyblue", "red3", "gray"), pch = 15, cex = 0.7, bg = NULL, box.lty = 0, text.font = 20)


# length(which(apply(trait[, c("Aerobic", "Anaerobic","Facultatively_Anaerobic")], 1, FUN = function(x){length(which(x == F)) == 1}) == F))
# length(which(apply(trait[, c("Aerobic", "Anaerobic","Facultatively_Anaerobic")], 1, FUN = function(x){length(which(x == F)) == 2}) == F))
# length(which(apply(trait[, c("Aerobic", "Anaerobic","Facultatively_Anaerobic")], 1, FUN = function(x){length(which(x == F)) == 3}) == T))

rarefied_ps_lm <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'L'| rarefied_ps@sam_data$fern.size == 'M')], rarefied_ps) 
rarefied_ps_lm <- prune_samples(sample_names(rarefied_ps_lm) %in% rownames(rarefied_ps_lm@sam_data)[which(rarefied_ps_lm@sam_data$layer != 'O')], rarefied_ps_lm) 

ph_var <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(var = var(pH, na.rm = T), C_N_var = var(C_N, na.rm = T))
ph_var$fern_ID <- as.character(ph_var$fern_ID)

plotList <- lapply(
  colnames(bugbase)[1:9],
  function(char) {
    otu_trait <- rarefied_ps_lm@otu_table %>% as.data.frame() %>% filter(rownames(rarefied_ps_lm@otu_table) %in% rownames(trait)[which(trait[, char] == T)])
    vec <- c()
    for(i in sort(unique(rarefied_ps_lm@sam_data$fern_ID))){
      dat <- otu_trait[which(rownames(t(otu_trait)) %in% rownames(rarefied_ps_lm@sam_data)[which(rarefied_ps_lm@sam_data$fern_ID== i) ]), ]
      dis_mat <- vegdist(dat, 'bray', na.rm = T)
      vec <- c(vec, i, mean(as.vector(dis_mat)))
    }
    dis <- matrix(vec, ncol = 2, nrow = 16, byrow = T)
    colnames(dis) <- c('fern_ID', 'dis_mean')
    dis[, 1] <- as.character(dis[, 1])
    dis <- as.data.frame(dis)
    
    ph_dis <- left_join(ph_var, dis)
    samdf$fern_ID <- as.character(samdf$fern_ID)
    ph_dis <- unique(left_join(ph_dis, samdf[, c(4, 6)]))
    ph_dis <- ph_dis %>% left_join(read.csv("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi_trait.csv"))
    
    ph_dis$dis_mean <- as.numeric(ph_dis$dis_mean)
    
    
    p <- ggscatter(ph_dis, x = 'C_N_var', y = 'dis_mean',
              add = "reg.line", conf.int = TRUE, color = "#DDE5B6", xlab = 'Variance of C:N ratio within an individual', ylab = 'Mean pairwised disimilarity \n within an individual')+
      stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
               label.y.npc = "top", size =7)+
      geom_point(cex = 5, aes( col = fern.size))+
      scale_colour_manual(values = c("#ADC178", "#997B66"), aesthetics = c("colour", "fill"))+
      theme(text = element_text(size = 15))+
      guides(col=guide_legend(title="Size"))+
      ggtitle(char)
    p
  }
)
ggarrange(plotlist=plotList, ncol = 3, nrow = 3)



