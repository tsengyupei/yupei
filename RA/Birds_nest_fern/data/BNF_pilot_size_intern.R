#00 Library----

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")

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
library(PerformanceAnalytics)
library(vegan)
library(igraph)
library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("metagenomeFeatures")
library(metagenomeFeatures)
library(remotes)
install.packages("microbiome")
library(BiocManager)
BiocManager::install("microbiome")
library(microbiome)
#01-1 Data preparation (feature table, metadata, taxanomy)----
load('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/seqtab.nochim.RData') #otu table
load('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/taxa .RData') #taxonomy table

#combine layer-level otu table into size-level otu table
seqtab.nochim <- as.data.frame(seqtab.nochim)
seqtab.nochim$sample <- sub('-.*', '', rownames(seqtab.nochim))
seqtab.nochim <- seqtab.nochim %>% group_by(sample) %>% summarise_all(sum) %>% as.data.frame(seqtab.nochim)
rownames(seqtab.nochim) <- seqtab.nochim[, 1]
seqtab.nochim <- seqtab.nochim[, -1]

#combine layer-level feature table into size-level feature table
samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi.csv', header = T)
samdf <- samdf[, -1]
samdf$fern_ID <- as.character(samdf$fern_ID)
# theme_set(theme_bw())
samples.out <- as.data.frame(rownames(seqtab.nochim))
colnames(samples.out)[1] <- 'fern_ID'
samdf <- left_join(samples.out, samdf, by = "fern_ID")
rownames(samdf) <- rownames(seqtab.nochim)

#01-2 Create phyloseq data----
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#Remove nagetive control
ps <- prune_samples(!(sample_names(ps)  %in% c("NC1", "NC2", "NC3", "NC4")), ps) # Remove mock sample

#Remove invalid taxa groups
ps <- subset_taxa(ps, Kingdom != "Unassigned")
ps <- subset_taxa(ps, Class != "D_2__Chloroplast")
ps <- subset_taxa(ps, Family != "D_4__Mitochondria")
ps

#Draw rarefaction Curve (raremax = 14231) w/o NA
col <- randomColor(count = 24)
raremax = min(sample_sums(ps)) #get the minimum point
tab <- otu_table(ps)
class(tab) <- "matrix"
rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=T, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax

#Create rarefied phyloseq data
rarefied_ps = rarefy_even_depth(ps, rngseed = 1, sample.size = raremax, replace = F)

##02-1 alpha diversity----
plot_richness(rarefied_ps, x = "nest.weight", measures = c("Observed", "Shannon", "Simpson"), color = "fern.size")

#02-2 SAR model----
rar_matrix <- rarefied_ps@otu_table

samdf <- samdf[-(25:28), ]
samdf$richness <- apply(rar_matrix>0,1,sum)
chart.Correlation(scale(samdf[, 4:18]), histogram=TRUE, pch=19)


hist(resid(lm(scale(log(samdf$richness)) ~ scale(log(samdf$leaf_length_mean)))))
chart.Correlation(scale(log(samdf[, c(5, 12:16, 18)])), histogram=TRUE, pch=19)

#Please try the log transformation and simple linear model by yourself----
#Q which morphological measurement and species level will be better to fit the SAR model?
#1) use log() function and scale() function to transform data
#2) use lm() to create simple linear model
#3) extract adjusted r-squared
#4) plotting to see the slope with combination of each species level and morphological measurement
#5) find the best combination and make sure the data fit the model assumption then give the statistic report
phy_family <- rarefied_ps %>% aggregate_taxa(level = 'Family')

#03_Relative abundance----
#Phylum
Stability_phylum<-rarefied_ps %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_phylum[Stability_phylum$Abundance<0.01,]$Phylum="Others"

colourCount = length(unique(Stability_phylum$Phylum))
colourCount
colors_phylum_n <- randomColor(count = 7)
colors_phylum_n <- distinctColorPalette(7)

PlotPhylum<-ggplot(Stability_phylum, aes(x=Sample,y=Abundance, fill=Phylum))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=colors_phylum_n)
PlotPhylum

#Class
Stability_class<-rarefied_ps %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_class[Stability_class$Abundance<0.01,]$Class="Others"

colourCount = length(unique(Stability_class$Class))
colourCount
colors_class_n <- randomColor(count = 16)
colors_class_n <- distinctColorPalette(16)

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
colors_order_n <- randomColor(count = 29)
colors_order_n <- distinctColorPalette(29)

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
colors_family_n <- randomColor(count = 49)
colors_family_n <- distinctColorPalette(49)

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
colors_genus_n <- randomColor(count = 65)
colors_genus_n <- distinctColorPalette(65)

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
colors_species_n <- randomColor(count = 77)
colors_species_n <- distinctColorPalette(77)

PlotSpecies<-ggplot(Stability_species, aes(x=Sample,y=Abundance,fill=Species))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_species_n)+
  theme(legend.text=element_text(size=6))
PlotSpecies


#04 multivariate analysis----
#nMDS (bray-curtis distance)
ps.prop <- transform_sample_counts(rarefied_ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color = "fern.size", title="Bray NMDS")