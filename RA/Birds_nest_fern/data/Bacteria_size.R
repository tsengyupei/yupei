#library----
library(stringr)
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(randomcoloR)
library(mgsub)
library(gridExtra)
library(grid)
library(pairwiseAdonis)
library(ggpubr)
library(PerformanceAnalytics)
library(vegan)
library(betapart)
library(dada2)
library(msa)


#Phyloseq data preparation----
#Import data with OTU composition for each sample 
feature_table <- read.csv(file = '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/feature-table.txt', sep = '', header = T, check.names = TRUE, row.names = 1) %>% select(-ID, -taxonomy)
feature_table[feature_table < 6] <- 0
feature_table <- feature_table %>% rename_with(~ mgsub::mgsub(colnames(feature_table), c("X", "_V6V8", "_"), c("", "", "-"))) %>% t() 
feature_table <- feature_table %>% as.data.frame() %>% mutate(fern_ID = as.data.frame(sub("\\-.*", "", rownames(feature_table))))
colnames(feature_table[, 21827]) <- "fern_ID"

feature_table <- 
  feature_table %>% 
  # filter(sub(".*-", "", rownames(feature_table)) != "O") %>% 
  group_by(fern_ID) %>%  summarise_all(sum) 

feature_table <- as.data.frame(feature_table)
rownames(feature_table) <- feature_table[, 1]
feature_table <- feature_table[, -1]


feature_table <- feature_table %>% otu_table(taxa_are_rows = F)

# seqs <- getSequences(t(feature_table_n))
# names(seqs) <- seqs
# mult <- msa(seqs, method="ClustalW", type="dna", order="input")

#Import data with taxonomy information od each OTU 
tax <- read.csv(file='/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/taxonomy.csv', sep = ' ', header = F, check.names = T, row.names = 1)
tax <- tax[-1, ] %>% mutate(all = rownames(tax[-1, ]))
tax[c('ID', 'Kingdom')] <- str_split_fixed(tax$all, ',', 2)
rownames(tax) <- tax[, "ID"]
tax <- tax[, c(9, 1:6)] %>% rename_with(~ c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>% as.matrix() %>% tax_table()

#Import data with sample information 
samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi_trait.csv', header = T) %>% select(-1)
distance <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_distance.csv')
colnames(distance) <- c("tree.ID", "x", "y")
samdf <- left_join(samdf, distance)
rownames(samdf) <- samdf$fern_ID

#construct phyloseq object
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = F), sample_data(samdf), tax_table(tax))
ps <- subset_taxa(ps, Kingdom != "Unassigned") %>% 
  subset_taxa(Class != "D_2__Chloroplast") %>% 
  subset_taxa(Family != "D_4__Mitochondria")
ps <- prune_samples(!(sample_names(ps) %in% c("NC1", "NC2", "NC3", "NC4")), ps) # Remove NC sample

#draw ps rarefaction Curve (raremax = 5235) w/o NA
col <- randomColor(count = 24)
raremax <- min(sample_sums(ps)) #get the minimum point

tab <- otu_table(ps)
class(tab) <- "matrix"
rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=T, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_ps <- rarefy_even_depth(ps, rngseed = 1, sample.size = raremax, replace = F)

##02 alpha diversity----
plot_richness(rarefied_ps, x = "fern.size", measures = c("Observed", "Shannon", "Simpson"), color = "fern.size")+
  theme(text = element_text(size = 30))+
  labs(color = 'Fern size')


#species richness
rar_matrix <- (rarefied_ps@otu_table)

samdf <- samdf[-(25:28), ]
samdf$richness <- apply(rar_matrix>0,1,sum)
chart.Correlation(log10(samdf[, c(14, 24:26, 29)]), histogram=TRUE, pch=19, "spearm")
summary(lm(richness ~ dry_mass + pH_ratio, log10(samdf[, c(14, 24:26, 29)])))
p <- ggscatter(log10(samdf[, c(14, 24:26, 29)]), x = 'dry_mass', y = 'richness', add = "reg.line", conf.int = TRUE, color = "#609067", xlab = 'log(Dry mass)', ylab = 'log(Operational taxonomic unit \n (OTU) richness)')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =10)+
  geom_point(cex = 4, col = "#609067")+
  theme(text = element_text(size = 24))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

#Ordination
ps.prop <- transform_sample_counts(rarefied_ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

ordplot <- plot_ordination(ps.prop, ord.nmds.bray, color = "fern.size", aes(size=5))
ordplot + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  geom_point(cex = 5)+
  labs(colour = "Fern size") +
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#F1DCA7", "#997B66"), label = c("Large", "Medium", "Small"))



bray.dist<-vegdist((ps.prop@otu_table), method='bray')
#PERMANOVA
adonis2(bray.dist~ as.factor(ps.prop@sam_data$fern.size) + ps.prop@sam_data$pH_ratio + ps.prop@sam_data$C_N_ratio, permutations=999)
adonis2(bray.dist~ as.factor(ps.prop@sam_data$fern.size), permutations=999)
permutest(betadisper(bray.dist, group=ps.prop@sam_data$fern.size, type = 'centroid'))



adonis2(dist(disper$distances) ~ as.factor(ps.prop@sam_data$fern.size))
pairwise_size <- pairwise.adonis(bray.dist, as.factor(ps.prop@sam_data$fern.size))

#Distance----
dist.abund <- vegdist(t(prune_samples(rarefied_ps@sam_data$fern.size == "S" | rarefied_ps@sam_data$fern.size == "M", rarefied_ps)@otu_table), method = "bray")
dist.space <- dist(prune_samples(rarefied_ps@sam_data$fern.size == "S" | rarefied_ps@sam_data$fern.size == "M", rarefied_ps)@sam_data[, c(4, 27:28)], method = "euclidean")

mantel(dist.abund, dist.space, method = "spearman", permutations = 9999, na.rm = TRUE)

vec.dist.abun <- as.vector(dist.abund)
vec.dist.space <- as.vector(dist.space)

plot(scale(vec.dist.abun) ~ (scale(vec.dist.space)))
hist(resid(lm((vec.dist.abun) ~ scale(vec.dist.space))))
hist(scale(vec.dist.abun))
hist((vec.dist.space))
summary(lm(scale(vec.dist.abun) ~ (scale(vec.dist.space))))
hist(resid(lm(scale(vec.dist.abun) ~ (scale(vec.dist.space)))))

#all data pair beta and partition into nestedness and turnover----
#transform abundance data to presence absence data
rar_matrix_pa <- rar_matrix
rar_matrix_pa[rar_matrix_pa > 0] <-1
rar_matrix_pa <- as.data.frame(rar_matrix_pa)

#calculate difference of size 
size.dist <- as.vector(dist(scale(samdf[-(25:28), ]$dry_mass), method = "euclidean", diag = FALSE, upper = FALSE))
size.dist <- as.data.frame(size.dist)

#calculate pair jaccard and sorensen
sor_dist <-beta.pair(rar_matrix_pa, index.family="sor")
jac_dist <-beta.pair(rar_matrix_pa, index.family="jac")

#create beta and bind to the difference of size 
size.dist$sor_tun <- as.vector(sor_dist$beta.sim)
size.dist$sor_nes <- as.vector(sor_dist$beta.sne)
size.dist$sor_all <- as.vector(sor_dist$beta.sor)

size.dist$jac_tun <- as.vector(jac_dist$beta.jtu)
size.dist$jac_nes <- as.vector(jac_dist$beta.jne)
size.dist$jac_all <- as.vector(jac_dist$beta.jac)

size.dist$sor_nes_ratio <- size.dist$sor_nes/size.dist$sor_all
size.dist$jac_nes_ratio <- size.dist$jac_nes/size.dist$jac_all

size.dist$sor_tun_ratio <- size.dist$sor_tun/size.dist$sor_all 
size.dist$jac_tun_ratio <- size.dist$jac_tun/size.dist$jac_all 

chart.Correlation(size.dist, histogram=TRUE, pch=19, cex.labels = 1.4, label.pos = 0.8)

#only large and medium beta pair----
rar_matrix_pa_ml <- rar_matrix_pa[-which(samdf$fern.size == 'S'), ]

#calculate difference of size
size.dist_ml <- dist(scale(samdf$nest_volume[-which(samdf$fern.size == 'S')]), method = "euclidean", diag = FALSE, upper = FALSE)
size.dist_ml <- as.data.frame(as.vector(size.dist_ml))

#calculate pair jaccard and sorensen
sor_dist_ml <-beta.pair(rar_matrix_pa_ml, index.family="sor")
jac_dist_ml <-beta.pair(rar_matrix_pa_ml, index.family="jac")

#create beta and bind to the difference of size 
size.dist_ml$sor_tun_ml <- as.vector(sor_dist_ml$beta.sim)
size.dist_ml$sor_nes_ml <- as.vector(sor_dist_ml$beta.sne)
size.dist_ml$sor_all_ml <- as.vector(sor_dist_ml$beta.sor)

size.dist_ml$jac_tun_ml <- as.vector(jac_dist_ml$beta.jtu)
size.dist_ml$jac_nes_ml <- as.vector(jac_dist_ml$beta.jne)
size.dist_ml$jac_all_ml <- as.vector(jac_dist_ml$beta.jac)



chart.Correlation(size.dist_ml, histogram=TRUE, pch=19, cex.labels = 1.4, label.pos = 0.8)

#only small and medium beta pair----
rar_matrix_pa_sm <- rar_matrix_pa[-(25:28), ]
rar_matrix_pa_sm <- rar_matrix_pa[-which(samdf[-(25:28), ]$fern.size == 'L'), ]

#calculate difference of size
size.dist_sm <- dist(scale(samdf[-(25:28), ]$dry_mass[-which(samdf$fern.size == 'L')]), method = "euclidean", diag = FALSE, upper = FALSE)
size.dist_sm <- as.data.frame(as.vector(size.dist_sm))

#calculate pair jaccard and sorensen
sor_dist_sm <-beta.pair(rar_matrix_pa_sm, index.family="sor")
jac_dist_sm <-beta.pair(rar_matrix_pa_sm, index.family="jac")

#create beta and bind to the difference of size 
size.dist_sm$sor_tun_sm <- as.vector(sor_dist_sm$beta.sim)
size.dist_sm$sor_nes_sm <- as.vector(sor_dist_sm$beta.sne)
size.dist_sm$sor_all_sm <- as.vector(sor_dist_sm$beta.sor)

size.dist_sm$jac_tun_sm <- as.vector(jac_dist_sm$beta.jtu)
size.dist_sm$jac_nes_sm <- as.vector(jac_dist_sm$beta.jne)
size.dist_sm$jac_all_sm <- as.vector(jac_dist_sm$beta.jac)

size.dist_sm$sor_nes_ratio <- size.dist_sm$sor_nes_sm/size.dist_sm$sor_all_sm 
size.dist_sm$jac_nes_ratio <- size.dist_sm$jac_nes_sm/size.dist_sm$jac_all_sm 

size.dist_sm$sor_tun_ratio <- size.dist_sm$sor_tun_sm/size.dist_sm$sor_all_sm 
size.dist_sm$jac_tun_ratio <- size.dist_sm$jac_tun_sm/size.dist_sm$jac_all_sm 

chart.Correlation(size.dist_sm, histogram=TRUE, pch=19, cex.labels = 1.4, label.pos = 0.8)


#phylo----

