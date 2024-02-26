# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")
# install.packages('tidyverse')

#BiocManager::install("MicrobiotaProcess")
#library(MicrobiotaProcess)
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(forcats)
# library(MicrobiomeStat) #problem
library(vegan)
# library(ggh4x) #problem
library(randomcoloR)
library(funfuns)
library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library(mobr)
library(tidyr)
library(purrr)
library(broom)
library(ggpubr)
#import data (feature table, metadata, taxanomy)
load('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/seqtab.nochim.RData')
load('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/taxa .RData')
samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv', header = T)
samples.out <- as.data.frame(rownames(seqtab.nochim))

samples.out$sample_ID <- sub('-ITS_R1.fastq.gz.*', '', rownames(seqtab.nochim))
# samples.out <- samples.out[, -1]

distance <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_distance.csv')
colnames(distance) <- c("tree.ID", "x", "y")
samdf <- left_join(samdf, distance)

samdf <- left_join(samples.out, samdf, by = "sample_ID")
rownames(samdf) <- rownames(seqtab.nochim)
samdf <- samdf[, -1]


# cn <- read.csv("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/CN_ratio.csv")
# cn <- cn %>% group_by(ID) %>% summarise(C = mean(C.), N = mean(N.), H = mean(H.)) %>% mutate(C_N = C/N)
# colnames(cn)[1] <- "sample_ID"
# samdf <- left_join(samdf, cn, by = "sample_ID")

boxplot(C_N.x ~ layer, samdf[which(samdf$fern.size == "L"), ])
boxplot(C_N.x ~ layer, samdf[which(samdf$fern.size == "M"), ])
summary(aov(C_N ~ layer, samdf[which(samdf$fern.size == "L"), ]))

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps <- prune_samples(!(sample_names(ps)  %in% c("NC1-ITS_R1.fastq.gz", "NC2-ITS_R1.fastq.gz", "NC3-ITS_R1.fastq.gz", "NC4-ITS_R1.fastq.gz")), ps) # Remove NC sample



ps <- subset_taxa(ps, Kingdom != "Unassigned")
ps <- subset_taxa(ps, Class != "D_2__Chloroplast")
ps <- subset_taxa(ps, Family != "D_4__Mitochondria")
ps

#draw ps rarefaction Curve (raremax = 14231) w/o NA
col <- randomColor(count = 69)
raremax = min(sample_sums(ps)) #get the minimum point

tab <- otu_table(ps)
class(tab) <- "matrix"


rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=F, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_ps = rarefy_even_depth(ps, rngseed = 1, sample.size = raremax, replace = F)

#01_Relative abundance----
pdf(file = "/Users/yupeitseng/Documents/RA/birds_nest_fern/plot/relative_aundance.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 9) # The height of the plot in inches
#Phylum
Stability_phylum<-rarefied_ps %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_phylum[Stability_phylum$Abundance<0.01,]$Phylum="Others"

colourCount = length(unique(Stability_phylum$Phylum))
colourCount
colors_phylum_n <- randomColor(count = 12)
colors_phylum_n <- distinctColorPalette(12)

PlotPhylum<-ggplot(Stability_phylum, aes(x=Sample,y=Abundance, fill=Phylum))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=colors_phylum_n)
  # geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales="free_x", nrow=1)+
  # labs(x="Samples",y="Relative abundance")+
  # scale_fill_manual(values=colors_phylum_n)+
  # guides(fill=guide_legend(ncol=1))+ 
  # theme_minimal(base_size=10)+
  # theme(strip.text=element_text(size=10), axis.text.x=element_text(angle=90, hjust=1), 
  #       legend.text=element_text(size=8), panel.spacing=unit(0.3,"cm"))
PlotPhylum

#Class
Stability_class<-rarefied_ps %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_class[Stability_class$Abundance<0.01,]$Class="Others"

colourCount = length(unique(Stability_class$Class))
colourCount
colors_class_n <- randomColor(count = 26)
colors_class_n <- distinctColorPalette(26)

PlotClass<-ggplot(Stability_class, aes(x=Sample,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=colors_class_n)
PlotClass

#Order
Stability_order<-rarefied_ps %>%
  tax_glom(taxrank = "Order") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_order[Stability_order$Abundance<0.01,]$Order="Others"

colourCount = length(unique(Stability_order$Order))
colourCount
colors_order_n <- randomColor(count = 48)
colors_order_n <- distinctColorPalette(48)

PlotOrder<-ggplot(Stability_order, aes(x=Sample,y=Abundance,fill=Order))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values=colors_order_n)
PlotOrder

#Family
Stability_family<-rarefied_ps %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_family[Stability_family$Abundance<0.01,]$Family="Others"

colourCount = length(unique(Stability_family$Family))
colourCount
colors_family_n <- randomColor(count = 86)
colors_family_n <- distinctColorPalette(86)

PlotFamily<-ggplot(Stability_family, aes(x=Sample,y=Abundance,fill=Family))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_family_n)+
  theme(legend.text=element_text(size=6))
PlotFamily

#Genus
Stability_genus<-rarefied_ps %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_genus
Stability_genus[Stability_genus$Abundance<0.01,]$Genus="Others"

colourCount = length(unique(Stability_genus$Genus))
colourCount
colors_genus_n <- randomColor(count = 112)
colors_genus_n <- distinctColorPalette(112)

PlotGenus<-ggplot(Stability_genus, aes(x=Sample,y=Abundance,fill=Genus))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_genus_n)+
  theme(legend.text=element_text(size=6))
PlotGenus

#Species
Stability_species<-rarefied_ps %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_species
Stability_species[Stability_species$Abundance<0.01,]$Species="Others"

colourCount = length(unique(Stability_species$Species))
colourCount
colors_species_n <- randomColor(count = 120)
colors_species_n <- distinctColorPalette(120)

PlotSpecies<-ggplot(Stability_species, aes(x=Sample,y=Abundance,fill=Species))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_species_n)+
  theme(legend.text=element_text(size=6))
PlotSpecies

dev.off()

##02 alpha diversity----
plot_richness(rarefied_ps, x = "layer", measures = c("Observed", "Shannon", "Simpson"), color = "fern.size")+
  theme(text = element_text(size = 30))+
  labs(color = 'Fern size')

#03 Transform data to proportions as appropriate for Bray-Curtis distances----
ps.prop <- transform_sample_counts(rarefied_ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color = "layer", shape = 'layer', title="Bray NMDS") + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t") + theme_bw()

#medium
ps.prop_m <- prune_samples(sample_names(ps.prop) %in% rownames(ps.prop@sam_data)[which(ps.prop@sam_data$fern.size == 'M')], ps.prop) 
ord.nmds.bray <- ordinate(ps.prop_m, method="NMDS", distance="bray")
plot_ordination(ps.prop_m, ord.nmds.bray, color = 'layer') + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  geom_point(cex = 5)+
  labs(colour = "Layer") +
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#F1DCA7", "#997B66"), label = c("Lower", "Outer", "Upper"))
bray.dist_m<-vegdist((ps.prop_m@otu_table), method='bray')
adonis2(bray.dist_m ~  as.factor(ps.prop_m@sam_data$layer), permutations=9999)
adonis2(bray.dist_m ~  ps.prop_m@sam_data$pH, permutations=9999)
adonis2(bray.dist_m ~  ps.prop_m@sam_data$C_N, permutations=9999)
permutest(betadisper(bray.dist_m, group=ps.prop_m@sam_data$layer))
pairwise.adonis(bray.dist_m, as.factor(ps.prop_m@sam_data$layer))

adonis2(bray.dist_m ~  as.factor(ps.prop_m@sam_data$layer) + ps.prop_m@sam_data$pH + ps.prop_m@sam_data$C_N, permutations=9999)

#large
ps.prop_l <- prune_samples(sample_names(ps.prop) %in% rownames(ps.prop@sam_data)[which(ps.prop@sam_data$fern.size == 'L')], ps.prop)
rarefied_ps_l <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'L')], rarefied_ps)
# ps.prop_l <- prune_samples(sample_names(ps.prop_l) %in% rownames(ps.prop_l@sam_data)[which(!(ps.prop_l@sam_data$sample_ID %in% c("13-L", "13-M", "13-O", "13-U", "19-L", "19-M", "19-O", "19-U")))], ps.prop_l)
ord.nmds.bray <- ordinate(ps.prop_l, method="NMDS", distance="bray")
plot_ordination(ps.prop_l, ord.nmds.bray, color = 'layer') + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  geom_point(cex = 5)+
  labs(colour = "Layer") +
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#DDE5B6", "#F1DCA7", "#997B66"), label = c("Lower", "Middle", "Outer", "Upper"))

bray.dist_l <- vegdist((ps.prop_l@otu_table), method='bray')
adonis2(bray.dist_l ~  as.factor(ps.prop_l@sam_data$layer), permutations=9999)
adonis2(bray.dist_l ~  ps.prop_l@sam_data$pH, permutations=9999, na.action = na.omit)
adonis2(bray.dist_l ~  ps.prop_l@sam_data$C_N, permutations=9999, na.action = na.omit)
permutest(betadisper(bray.dist_l, group=ps.prop_l@sam_data$layer))
pairwise.adonis(bray.dist_l, as.factor(ps.prop_l@sam_data$layer))
adonis2(bray.dist_l ~  as.factor(ps.prop_l@sam_data$layer) + ps.prop_l@sam_data$pH + ps.prop_l@sam_data$C_N, permutations=9999, na.action = na.omit)
#pH range ~ species composition distance (bray)
rarefied_ps_lm <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'L'| rarefied_ps@sam_data$fern.size == 'M')], rarefied_ps)
ps_lm <- prune_samples(sample_names(ps) %in% rownames(ps@sam_data)[which(ps@sam_data$fern.size == 'L'| ps@sam_data$fern.size == 'M')], ps)

ph_range <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(max = max(pH, na.rm = T), min = min(pH, na.rm = T)) %>% mutate(range = max-min)
ph_range <- unique(left_join(ph_range, samdf[, c(5, 7, 18)], by = "fern_ID"))
ph_range$range <- ph_range$max - ph_range$min

#order by pH range with group M and L
level_order <- as.vector(unlist(ph_range[order(ph_range[, 4]), "fern_ID"]))[c(1:2, 4:6, 8:9, 3, 7, 10:16)]

#order by dry mass
level_order <- as.vector(unlist(ph_range[order(ph_range[, 6]), "fern_ID"]))

#individual pH range
ggplot(ph_range) +
  geom_segment(aes(x = factor(fern_ID, level = level_order), xend = factor(fern_ID, level = level_order),
                   y = min, yend = max)) +
  geom_point(aes(x = factor(fern_ID, level = level_order), y = min, col = fern.size), size = 5) +
  geom_point(aes(x = factor(fern_ID, level = level_order), y = max, col = fern.size), size = 5)+
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#997B66"), label = c("Large", "Medium"))+
  labs(colour = "Fern size") +
  # xlab("Dry mass (g)")+
  xlab("Fern ID")+
  ylab("pH")

ggplot(ph_range) +
  geom_segment(aes(x = factor(fern_ID, level = level_order), xend = factor(fern_ID, level = level_order),
                   y = min, yend = max)) +
  geom_point(aes(x = factor(fern_ID, level = level_order), y = min, col = fern.size), size = 5) +
  geom_point(aes(x = factor(fern_ID, level = level_order), y = max, col = fern.size), size = 5)+
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#997B66"), label = c("Large", "Medium"))+
  labs(colour = "Fern size") +
  xlab("Fern ID")+
  ylab("pH")

cn_range <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(max = max(C_N, na.rm = T), min = min(pH, na.rm = T)) %>% mutate(range = max-min)
cn_range <- unique(left_join(cn_range, samdf[, c(5, 7, 18)], by = "fern_ID"))

plot(range ~ dry_mass.y, cn_range)
ggscatter(cn_range, x = 'dry_mass.y', y = 'range',
          add = "reg.line", conf.int = TRUE, color = "#d8b365", xlab = 'Fern size', ylab = 'pH range difference')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =15)+
  geom_point(aes(cex = 8, col = fern.size))+
  theme(text = element_text(size = 38))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )
ggsave('/Users/yupeitseng/Downloads/pH_range.png', p, dpi = 300, height = 12, width = 16)

cn_var <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(var = var(C_N, na.rm = T))
cn_var$fern_ID <- as.character(cn_var$fern_ID)

ph_var <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(var = var(pH, na.rm = T))
ph_var$fern_ID <- as.character(ph_var$fern_ID)
vec <- c()
for(i in sort(unique(rarefied_ps_lm@sam_data$fern_ID))){
  dat <- rarefied_ps_lm@otu_table[which(rownames(rarefied_ps_lm@otu_table) %in% rownames(rarefied_ps_lm@sam_data)[which(rarefied_ps_lm@sam_data$fern_ID== i) ]), ]
  dis_mat <- vegdist(dat, 'bray')
  vec <- c(vec, i, mean(as.vector(dis_mat)))
}
dis <- matrix(vec, ncol = 2, nrow = 16, byrow = T)
colnames(dis) <- c('fern_ID', 'dis_mean')
dis[, 1] <- as.character(dis[, 1])
dis <- as.data.frame(dis)

cn_dis <- left_join(cn_var, dis)
samdf$fern_ID <- as.character(samdf$fern_ID)
cn_dis <- unique(left_join(cn_dis, samdf[, c(5, 7, 18)]))

fit <- lm(dis_mean ~ var, cn_dis)
plot(dis_mean ~ var, cn_dis)
abline(fit)

cn_dis$dis_mean <- as.numeric(cn_dis$dis_mean)
cn_dis$var <- as.numeric(cn_dis$var)
ggscatter(cn_dis, x = 'var', y = 'dis_mean',
          add = "reg.line", conf.int = TRUE, color = "#609067", xlab = 'Variance of pH within an individual', ylab = 'Mean pairwised disimilarity \n within an individual')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =15)+
  geom_point(cex = 8, aes( col = fern.size))+
  theme(text = element_text(size = 38))+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )
ggsave('/Users/yupeitseng/Downloads/pH_dissimilarity.png', p, dpi = 300, height = 12, width = 16)
  

# theme_set(theme_bw())
samdf_lm <- samdf[which(samdf$fern.size == 'L' | samdf$fern.size == 'M'), ]


sam

#03_PcoA #bray
#large
pcoa <- capscale (as.data.frame(ps.prop_l@otu_table)[-(1:4), ] ~ 1, distance = 'bray')
plot(pcoa, type="n", las=1, main = "PCoA (bray)_large BNF")

text(pcoa, display = "sites", 
     labels = ps.prop_l@sam_data$sample_ID)
#only ellipse
# ordiellipse(pcoa, ps.prop_l@sam_data$layer, display="sites", 
#             kind = "se", conf=0.95, label=T, 
#             lwd=2, col=c("red", "blue","orange", "black"))
#only spider
# ordispider(pcoa, ps.prop_l@sam_data$layer, display="sites", 
#            label=T, lwd=2, col=c("red", "blue","orange", "black"))

ordiellipse(pcoa, ps.prop_l@sam_data$layer[-(1:4)], display="sites", 
            kind = "se", conf=0.95, 
            label=F, col=c("red", "blue","orange", "black"), 
            draw = "polygon", alpha = 50)
ordispider(pcoa, ps.prop_l@sam_data$layer[-(1:4)], display="sites", 
           label=T, lwd=2, col=c("red", "blue","orange", "black"))

bray.dist_l<-vegdist(ps.prop_l@otu_table, method='bray')
adonis2(bray.dist_l~as.factor(ps.prop_l@sam_data$layer) + ps.prop_l@sam_data$pH + ps.prop_l@sam_data$C_N, permutations=9999, na.action = na.omit)
permutest(betadisper(bray.dist_l, group=ps.prop_l@sam_data$layer))
a <- pairwise.adonis(bray.dist_l, as.factor(ps.prop_l@sam_data$layer))

write.csv(a, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/pairwise_layer.csv')
bray.dist_m<-vegdist(ps.prop_m@otu_table, method='bray')
adonis2(bray.dist_m~as.factor(ps.prop_m@sam_data$layer)+ ps.prop_m@sam_data$pH + ps.prop_m@sam_data$C_N, permutations=9999, na.action = na.omit)
permutest(betadisper(bray.dist_m, group=ps.prop_m@sam_data$layer))
pairwise.adonis(bray.dist_m, as.factor(ps.prop_m@sam_data$layer))

bray.dist<-vegdist(ps.prop@otu_table, method='bray')
adonis2(bray.dist~as.factor(ps.prop@sam_data$layer), permutations=9999)
permutest(betadisper(bray.dist, group=ps.prop@sam_data$layer))
pairwise.adonis(bray.dist_l, as.factor(ps.prop_l@sam_data$layer))
#rda
spe <- rarefied_ps@otu_table[-3, ]  # rename variables to make them shorter
env <- samdf # select only two explanatory variables
env <- env[-c(3, 66:69), ]

library (vegan)
spe.log <- log1p (spe)
spe.log.hell <- decostand (spe.log, 'hell')
tbRDA <- rda (spe.log.hell ~ pH + layer + fern.size, data = env)
head (summary (tbRDA))  # prints first lines of the summary of tbRDA
ordiplot(tbRDA, type = 'p')
test_axis <- anova(tbRDA, by = 'margin')

test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

write.csv(test_axis.adj, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/rda_layer_margin.csv')

#medium
pcoa <- capscale (as.data.frame(ps.prop_m@otu_table) ~ 1, distance = 'bray')
plot(pcoa, type="n", las=1, main = "PCoA (bray)_medium BNF")

text(pcoa, display = "sites", 
     labels = ps.prop_m@sam_data$sample_ID)
#only ellipse
# ordiellipse(pcoa, ps.prop_l@sam_data$layer, display="sites", 
#             kind = "se", conf=0.95, label=T, 
#             lwd=2, col=c("red", "blue","orange", "black"))
#only spider
# ordispider(pcoa, ps.prop_l@sam_data$layer, display="sites", 
#            label=T, lwd=2, col=c("red", "blue","orange", "black"))

ordiellipse(pcoa, ps.prop_m@sam_data$layer, display="sites", 
            kind = "se", conf=0.95, 
            label=F, col=c("red", "blue","orange"), 
            draw = "polygon", alpha = 50)
ordispider(pcoa, ps.prop_m@sam_data$layer, display="sites", 
           label=T, lwd=2, col=c("red", "blue","orange"))

#04_SAR process (only medium and large)----
#mobr
inv_mob_in <- make_mob_in(ps_lm@otu_table, ps_lm@sam_data)
inv_stats <- get_mob_stats(inv_mob_in, group_var = "fern_ID", index = c("N", "S", "S_n", "S_PIE", "S_asymp"))
plot(inv_stats)

group <- inv_stats[["groups_stats"]] %>% 
  filter(effort == 9825 | is.na(effort) == T) %>% 
  select(-effort) %>% 
  filter(index %in% c('S_n', 'S_PIE', 'S_asymp'))
group$scale <- rep('gamma', nrow(group))
colnames(group)[1] <- 'fern_ID'
group <- left_join(group, samdf[, c(2, 3, 12, 14, 17)])


samples <- inv_stats[["samples_stats"]]
samples <- samples %>% filter(!(index %in% c('N', 'S', 'beta_S', 'S_asymp'))) %>% filter(is.na(index) == F) %>% select(-effort) 
samples <- samples %>% group_by(group, index) %>%  summarise(mean(value))
samples$scale <- rep(c('alpha', 'beta', 'alpha', 'beta'), 16)
colnames(samples)[1] <- 'fern_ID'
samples <- left_join(samples, samdf[, c(2, 3, 12, 14, 17)])
colnames(samples)[3] <- 'value'

dat <- rbind(group, samples)
  
# n_ref <-max(c(max(samples %>% filter(index == "N") %>% select(value)), min(samples %>% filter(index == "N") %>% select(value))))

alpha <- dat[which(dat$scale == 'alpha'), ]
beta<- dat[which(dat$scale == 'beta'), ]

#plotting
png(file="/Users/yupeitseng/Documents/RA/birds_nest_fern/plot/sar_mechanism.png",
    width=2000, height=750)

gamma_p <- ggplot(data = group, aes(x = log10(dry_mass), y = log10(value), colour=index, fill=index,shape=index)) +
  geom_point(size=5)+
  stat_smooth(method="lm",se=TRUE)+
  scale_colour_manual(name="",values = c("S_asymp" = "#D3E4A7","S_n"="#609067", "S_PIE" = "#d8b365"), labels=c("S_asymp"=expression(S["total"]),"S_n"=expression(paste(gamma,"S"["n"])), "S_PIE"=expression(paste(gamma,"S"["PIE"]))))+
  scale_fill_manual(name="",values = c("S_asymp" = "#D3E4A7","S_n"="#609067", "S_PIE" = "#d8b365"), labels=c("S_asymp"=expression(S["total"]),"S_n"=expression(paste(gamma,"S"["n"])), "S_PIE"=expression(paste(gamma,"S"["PIE"]))))+
  scale_shape_manual(name="",values = c("S_asymp" = 8,"S_n"=17, "S_PIE" = 16), labels=c("S_asymp"=expression(S["total"]),"S_n"=expression(paste(gamma,"S"["n"])), "S_PIE"=expression(paste(gamma,"S"["PIE"]))))+
  xlab("lag(Dry mass)") + ylab("log(OTU richness)")+
  theme(text = element_text(size = 40))

beta_p <- ggplot(data = beta, aes(x = log10(dry_mass), y = log10(value), colour=index, fill=index,shape=index)) +
  geom_point(size=5)+
  stat_smooth(method="lm",se=TRUE)+
  scale_colour_manual(name="",values = c("beta_S_n"="#609067", "beta_S_PIE" = "#d8b365"),
                      labels=c("beta_S_n"=expression(paste(beta,"S"["n"])),
                               "beta_S_PIE"=expression(paste(beta,"S"["PIE"]))))+
  scale_fill_manual(name="",values = c("beta_S_n"="#609067", "beta_S_PIE" = "#d8b365"),
                    labels=c("beta_S_n"=expression(paste(beta,"S"["n"])),
                             "beta_S_PIE"=expression(paste(beta,"S"["PIE"]))))+
  scale_shape_manual(name="",values = c("beta_S_n"=17, "beta_S_PIE" = 16),
                     labels=c("beta_S_n"=expression(paste(beta,"S"["n"])),
                              "beta_S_PIE"=expression(paste(beta,"S"["PIE"]))))+
  xlab("lag(Dry mass)") + ylab("log(OTU richness)")+
  theme(text = element_text(size = 40))
 


alpha_p <- ggplot(data = alpha, aes(x = log10(dry_mass), y = log10(value), colour=index, fill=index,shape=index)) +
  geom_point(size=5)+
  stat_smooth(method="lm",se=TRUE)+
  scale_colour_manual(name="",values = c("S_n"="#609067", "S_PIE" = "#d8b365"),
                      labels=c("S_n"=expression(paste(alpha,"S"["n"])),
                               "S_PIE"=expression(paste(alpha,"S"["PIE"]))))+
  scale_fill_manual(name="",values = c("S_n"="#609067", "S_PIE" = "#d8b365"),
                    labels=c("S_n"=expression(paste(alpha,"S"["n"])),
                             "S_PIE"=expression(paste(alpha,"S"["PIE"]))))+
  scale_shape_manual(name="",values = c("S_n"=17, "S_PIE" = 16),
                     labels=c("S_n"=expression(paste(alpha,"S"["n"])),
                              "S_PIE"=expression(paste(alpha,"S"["PIE"]))))+
  xlab("log(Dry mass)") + ylab("log(OTU richness)")+
  theme(text = element_text(size = 40))

ggarrange(gamma_p, alpha_p, beta_p, nrow = 1, ncol = 3)

dev.off()

#
by_index <- dat %>%
  group_by(scale, index) %>%
  nest()

by_index2 <- by_index %>% 
  mutate(model   = map(data, ~ lm(log10(value) ~ log10(dry_mass), data = .)),
         Intercept = map_dbl(model, ~coef(.x)[1]),
         Slope     = map_dbl(model, ~coef(.x)[2]),
         glance    = map(model, glance)
  ) %>% 
  unnest(glance)

names(by_index2)

by_index2a <- by_index2 %>%
  select(index:adj.r.squared, -model, -data, p.value)

write.csv(by_index2a, "/Users/yupeitseng/Documents/RA/birds_nest_fern/plot/sar_mechanism.csv")


#Funguild database----
taxa <- as.data.frame(taxa)
taxa_fun <- do.call(paste, c(taxa[colnames(taxa)], sep = ",")) %>% as.data.frame() 
colnames(taxa_fun) <- "Taxonomy"                              
funguild <- funguild_assign(taxa_fun)
rownames(funguild) <- rownames(taxa)
funguild$OTU <- rownames(funguild)


OTU_funguild <-rarefied_ps %>%
  # tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance != 0)
OTU_funguild <- left_join(OTU_funguild, funguild, by = "OTU")

ggplot(OTU_funguild, aes(x=layer,y=Abundance, fill=trophicMode))+
  geom_bar(stat="identity")+
  facet_grid(fern.size ~ fern_ID)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  xlab("Layer")+
  guides(fill=guide_legend(title="Trophic mode"))



OTU_funguild_group <- OTU_funguild %>% group_by(fern_ID, layer, trophicMode) %>% summarise(richness = n()) %>% left_join(samdf[, c("fern_ID", "dry_mass.y", "litter.and.leaves.weight", "C", "N", "C_N", "layer", "pH", "dry_mass.x")])
stat.test <- OTU_funguild_group %>% group_by(trophicMode) %>% t_test(richness ~ layer)

stat.test <- stat.test %>% add_xy_position(x = "layer")

ggboxplot(OTU_funguild_group, x = "layer", y = "richness", fill = "layer")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)+
  facet_grid(~trophicMode, scales="free")


#function for checking the multiple regression which richness regressed on dry mass, pH and C:N ratio
lm_co <- function(dat = dat, trophic = trophic){
  OTU_funguild_group_fil <- dat %>% filter(trophicMode == trophic) 
  return(summary(lm(log10(richness) ~ layer + log10(dry_mass.x) + scale(pH) + scale(C_N), data = OTU_funguild_group_fil))$coefficients)
}

#function for extrating only the OTU that are in the interesting groups
tophic_group_phy <- function(db = db, keyword = keyword, ps.prop_fil = ps.prop_fil){
  OTU_funguild_fil <- unique(db %>% 
                               # filter(trophicMode %in% grep(keyword, unique(db$trophicMode), value = TRUE)) %>%
                               filter(trophicMode == keyword) %>%
                               select(OTU) %>% unlist %>% as.vector)
  
  allTaxa <- taxa_names(ps.prop_fil)
  allTaxa <- allTaxa[!(allTaxa %in% OTU_funguild_fil)]
  ps.prop_fil <-  prune_taxa(allTaxa, ps.prop_fil)
  return(ps.prop_fil)
}

perm_aov <- function(phy = phy){
  ord.nmds.bray <- ordinate(phy, method="NMDS", distance="bray")
  summary(ord.nmds.bray)
  bray.dist<-vegdist(phy@otu_table, method='bray')
  #PERMANOVA
  return(as.data.frame(adonis2(bray.dist~ (phy@sam_data$layer) + phy@sam_data$pH + phy@sam_data$C_N, permutations=9999, na.action = na.omit)))
}

perm_tab <- function(db = db, ps.prop_fil = ps.prop_fil){
  mat <- perm_aov(tophic_group_phy(db, "Saprotroph", ps.prop_fil))
  for(i in c("Pathotroph", "Symbiotroph")){
    mat <- rbind(mat, perm_aov(tophic_group_phy(db, i, ps.prop_fil)))
  }
  mat <- rbind(mat, perm_aov(ps.prop_fil))
  return(mat)
}

perm_tab_l <- perm_tab(OTU_funguild_l, ps.prop_l)
perm_tab_m <- perm_tab(OTU_funguild_m, ps.prop_m)

write.csv(perm_tab_m, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/perm_tab_m.csv")
write.csv(perm_tab_l, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/perm_tab_l.csv")

#large
OTU_funguild_l <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'L')], rarefied_ps)  %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance != 0)

samdf$fern_ID <- as.numeric(samdf$fern_ID)
OTU_funguild_l <- left_join(OTU_funguild_l, funguild, by = "OTU")
OTU_funguild_group_l <- OTU_funguild_l %>% group_by(fern_ID, layer, trophicMode) %>% summarise(richness = n()) %>% left_join(samdf[, c("fern_ID", "dry_mass.y", "litter.and.leaves.weight", "C", "N", "C_N", "layer", "pH", "dry_mass.x")])

lm_mat_l <- matrix(ncol = 4, nrow = 1)
lm_mat_l <- lm_co(OTU_funguild_group_l, "Pathotroph")
for(i in unique(OTU_funguild_group_l$trophicMode)[-1]){
  lm_mat_l <- rbind(lm_mat_l, lm_co(OTU_funguild_group_l, i))
}

OTU_funguild_all_l <- OTU_funguild_group_l %>% group_by(fern_ID, layer) %>% summarise(richness = sum(richness), pH = mean(pH), C_N = mean(C_N), dry_mass.x = mean(dry_mass.x))
lm_mat_l <- rbind(lm_mat_l, summary(lm(scale(richness) ~ layer + scale(dry_mass.x) + scale(pH) + scale(C_N), data = OTU_funguild_all_l))$coefficients)

write.csv(lm_mat_l, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/lm_mat_l.csv")

#box plot for checking the richness among layers
stat.test <- OTU_funguild_group_l %>% group_by(trophicMode) %>% t_test(richness ~ layer)
stat.test <- stat.test %>% add_xy_position(x = "layer")

ggboxplot(OTU_funguild_group_l, x = "layer", y = "richness", fill = "layer")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)+
  scale_colour_manual(values = c("#ADC178", "#DDE5B6", "#F1DCA7", "#997B66"), aesthetics = c("colour", "fill"))+
  xlab("Layer")+ ylab("OTU richness")+
  guides(fill=guide_legend(title="Layer"))+
  facet_wrap(~trophicMode, scales="free", labeller = label_wrap_gen(width=3))+
  theme(text = element_text(size = 18))

#medium
OTU_funguild_m <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'M')], rarefied_ps)  %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance != 0)
OTU_funguild_m <- left_join(OTU_funguild_m, funguild, by = "OTU")
OTU_funguild_group_m <- OTU_funguild_m %>% group_by(fern_ID, layer, trophicMode) %>% summarise(richness = n()) %>% left_join(samdf[, c("fern_ID", "dry_mass.y", "litter.and.leaves.weight", "C", "N", "C_N", "layer", "pH", "dry_mass.x")])

lm_mat_m <- matrix(ncol = 4, nrow = 1)

lm_mat_m <- lm_co(OTU_funguild_group_m, "Pathotroph")
for(i in unique(OTU_funguild_group_m$trophicMode)[-1]){
  lm_mat_m <- rbind(lm_mat_m, lm_co(OTU_funguild_group_m, i))
}

OTU_funguild_all_m <- OTU_funguild_group_m %>% group_by(fern_ID, layer) %>% summarise(richness = sum(richness), pH = mean(pH), C_N = mean(C_N), dry_mass.x = mean(dry_mass.x))
lm_mat_m <- rbind(lm_mat_m, summary(lm(scale(richness) ~ layer + scale(dry_mass.x) + scale(pH) + scale(C_N), data = OTU_funguild_all_m))$coefficients)

write.csv(lm_mat_m, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/lm_mat_m.csv")

#box plot for checking the richness among layers
stat.test <- OTU_funguild_group_m %>% group_by(trophicMode) %>% t_test(richness ~ layer)
stat.test <- stat.test %>% add_xy_position(x = "layer")

ggboxplot(OTU_funguild_group_m, x = "layer", y = "richness", fill = "layer")+
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01)+
  scale_colour_manual(values = c("#ADC178", "#F1DCA7", "#997B66"), aesthetics = c("colour", "fill"))+
  xlab("Layer")+ ylab("OTU richness")+
  guides(fill=guide_legend(title="Layer"))+
  facet_wrap(~trophicMode, scales="free", labeller = label_wrap_gen(width=3))+
  theme(text = element_text(size = 18))

#large & medium
OTU_funguild_lm <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'M' | rarefied_ps@sam_data$fern.size == 'L')], rarefied_ps)  %>% 
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance != 0)
OTU_funguild_lm <- left_join(OTU_funguild_lm, funguild, by = "OTU")
OTU_funguild_group_lm <- OTU_funguild_lm %>% group_by(fern_ID, layer, trophicMode) %>% summarise(richness = n()) %>% left_join(samdf[, c("fern_ID", "dry_mass.y", "litter.and.leaves.weight", "C", "N", "C_N", "layer", "pH", "dry_mass.x")])

lm_mat_lm <- matrix(ncol = 4, nrow = 1)

lm_mat_lm <- lm_co(OTU_funguild_group_lm, "Pathotroph")
for(i in unique(OTU_funguild_group_lm$trophicMode)[-1]){
  lm_mat_lm <- rbind(lm_mat_lm, lm_co(OTU_funguild_group_lm, i))
}

OTU_funguild_all_lm <- OTU_funguild_group_lm %>% group_by(fern_ID, layer) %>% summarise(richness = sum(richness), pH = mean(pH), C_N = mean(C_N), dry_mass.x = mean(dry_mass.x))
lm_mat_lm <- rbind(lm_mat_lm, summary(lm(scale(richness) ~ layer + scale(dry_mass.x) + scale(pH) + scale(C_N), data = OTU_funguild_all_lm))$coefficients)

write.csv(lm_mat_lm, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/lm_mat_lm.csv")

# OTU_funguild_l %>% filter(trophicMode == "Pathotroph") %>% group_by(guild, layer) %>% summarise(ab = sum(Abundance))
# ggplot(data=OTU_funguild_l %>% filter(trophicMode == "Pathotroph"), aes(x=layer, y=Abundance, fill=guild)) +
#   geom_bar(stat="identity")+
#   facet_grid(~fern_ID)
# OTU_funguild_m %>% filter(trophicMode == "Pathotroph") %>% group_by(guild) %>% summarise(n = n())
# ggplot(data=OTU_funguild_m %>% filter(trophicMode == "Pathotroph"), aes(x=layer, y=Abundance, fill=guild)) +
#   geom_bar(stat="identity")+
#   facet_grid(~fern_ID)



ggscatter(OTU_funguild_group_l, x = 'dry_mass.x', y = 'richness',
          add = "reg.line", conf.int = TRUE, color = "#609067")+
  # geom_point(aes(colour = "#609067"))+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =5)+
  facet_wrap(~trophicMode, scales="free")+
  labs(x = "pH", y = "OTU richness")+
  theme(text = element_text(size = 18))


OTU_funguild_fil <- OTU_funguild %>% filter(trophicMode == "Pathotroph-Saprotroph-Symbiotroph" | trophicMode == "Pathotroph-Symbiotroph") %>% distinct(OTU, .keep_all = T)
unique(OTU_funguild_fil$notes)

# function for calculating pH variance, C:N ratio variance and mean pairwise disimilarity in each individual
ph_CN <- function(keyword = keyword){
  ph_var <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(var = var(pH, na.rm = T), C_N_var = var(C_N), size = unique(fern.size))
  ph_var$fern_ID <- as.character(ph_var$fern_ID)
  if (keyword == "All"){
    rarefied_ps_lm_fil <- rarefied_ps_lm
  }
  else{
    rarefied_ps_lm_fil <- tophic_group_phy(OTU_funguild_lm, keyword, rarefied_ps_lm)
  }
  vec <- c()
  for(i in sort(unique(rarefied_ps_lm_fil@sam_data$fern_ID))){
    dat <- rarefied_ps_lm_fil@otu_table[which(rownames(rarefied_ps_lm_fil@otu_table) %in% rownames(rarefied_ps_lm_fil@sam_data)[which(rarefied_ps_lm_fil@sam_data$fern_ID== i) ]), ]
    dis_mat <- vegdist(dat, 'bray')
    vec <- c(vec, i, mean(as.vector(dis_mat)))
  }
  dis <- matrix(vec, ncol = 2, nrow = 16, byrow = T)
  colnames(dis) <- c('fern_ID', 'dis_mean')
  dis[, 1] <- as.character(dis[, 1])
  dis <- as.data.frame(dis)
  
  ph_dis <- left_join(ph_var, dis)
  samdf$fern_ID <- as.character(samdf$fern_ID)
  ph_dis <- unique(left_join(ph_dis, samdf[, c(5, 18)]))
  
  ph_dis$dis_mean <- as.numeric(ph_dis$dis_mean)
  ph_dis$var <- as.numeric(ph_dis$var)
  return(ph_dis)
  # return(summary(lm(scale(dis_mean) ~  scale(C_N_var), data = ph_dis)))
}
hist(ph_CN("Pathotroph")$var)
ph_CN("Saprotroph")
ph_CN("Symbiotroph")

ggarrange(ph_plot(ph_CN("All")), ph_plot(ph_CN("Pathotroph")), ph_plot(ph_CN("Saprotroph")), ph_plot(ph_CN("Symbiotroph")), nrow = 2, ncol = 2, labels = c("All",  "Pathotroph", "Saprotroph", "Symbiotroph"), common.legend = TRUE, legend = "bottom", vjust = 2)
ph_plot <- function(dat = dat){
  p <- ggscatter(dat, x = 'var', y = 'dis_mean',
            add = "reg.line", conf.int = TRUE, color = "#DDE5B6", xlab = 'Variance of pH within an individual', ylab = 'Mean pairwised disimilarity \n within an individual')+
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
             label.y.npc = "top", size =6)+
    geom_point(cex = 8, aes( col = size))+
    scale_colour_manual(values = c("#ADC178", "#997B66"), aesthetics = c("colour", "fill"))
  return(p)
}

#distance effect with only outer samples
ps.prop_os <- prune_samples(sample_names(ps.prop) %in% rownames(ps.prop@sam_data)[which(ps.prop@sam_data$layer == 'O')], ps.prop)

dist_dis <- function(phy = phy){
  dist.abund <- vegdist(phy@otu_table, method = "bray")
  dist.space <- dist(phy@sam_data[, c(8, 27:28)], method = "euclidean")
  
  mantel(dist.abund, dist.space, method = "spearman", permutations = 9999, na.rm = TRUE)
  vec.dist.abun <- as.vector(dist.abund)
  vec.dist.space <- as.vector(dist.space)
  
  plot(scale(vec.dist.abun) ~ (scale(vec.dist.space)))
  summary(lm(scale(vec.dist.abun) ~ (scale(vec.dist.space))))
}

