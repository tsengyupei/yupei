getwd()
setwd("/Users/haijin/Desktop/NTU/02_Bitou/L1_to_L3_silva128/Bitou_Ranalysis_silva128")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("phyloseq"))
#BiocManager::install("MicrobiotaProcess")
#library(MicrobiotaProcess)
install.packages('rlang')
install.packages("ggplot2")
install.packages("vegan")
install.packages("patchwork")
install.packages('Rarefy')
install.packages('randomcoloR')
install.packages('ggh4x')
install.packages('forcats')
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(vegan)
library(ggpubr)
library(dplyr)
library(MicrobiomeStat)
library(Rarefy)
library(randomcoloR)
library(ggh4x)
library(forcats)


#import data (feature table, metadata, taxanomy)
feature_table_L1 <- read.csv(file="feature_table_L1.csv", header = T, check.names = TRUE, row.names = 1)
tax_L1 <- read.csv(file="taxa_table_L1.csv",header=T,check.names = T, row.names = 1)
metadata_L1 <- read.csv(file="metadata_L1.csv", header = T, check.names = T, row.names = 1)
#taxa_2011 <- read.table("taxonomy_2011.tsv", sep = '\t', fill = TRUE)

feature_table_L2 <- read.csv(file="feature_table_L2.csv", header = T, check.names = TRUE, row.names = 1)
tax_L2 <- read.csv(file="taxa_table_L2.csv",header=T,check.names = T, row.names = 1)
metadata_L2 <- read.csv(file="metadata_L2.csv", header = T, check.names = T, row.names = 1)
#taxa_2102 <- read.table("taxonomy_2102.tsv", sep = '\t', fill = TRUE)

feature_table_L3 <- read.csv(file="feature_table_L3.csv", header = T, check.names = TRUE, row.names = 1)
tax_L3 <- read.csv(file="taxa_table_L3.csv",header=T,check.names = T, row.names = 1)
metadata_L3 <- read.csv(file="metadata_L3.csv", header = T, check.names = T, row.names = 1)
#taxa_2102 <- read.table("taxonomy_2102.tsv", sep = '\t', fill = TRUE)

tax_L1 <- as.matrix(tax_L1)
TAX_L1 <- tax_table(tax_L1)
meta_L1 <- as.data.frame(metadata_L1)
META_L1 <- sample_data(meta_L1)
colnames(feature_table_L1) <- row.names(meta_L1)
OTU_L1 <- otu_table(feature_table_L1, taxa_are_rows = TRUE)

physeq_L1 <- phyloseq(OTU_L1, TAX_L1, META_L1)
physeq_L1 <- subset_taxa(physeq_L1, Kingdom != "Unassigned")
physeq_L1 <- subset_taxa(physeq_L1, Class != "D_2__Chloroplast")
physeq_L1 <- subset_taxa(physeq_L1, Class != "D_2__Chloroplast")
physeq_L1 <- subset_taxa(physeq_L1, Family != "D_4__Mitochondria")
physeq_L1

tax_L2 <- as.matrix(tax_L2)
TAX_L2 <- tax_table(tax_L2)
meta_L2 <- as.data.frame(metadata_L2)
META_L2 <- sample_data(meta_L2)
colnames(feature_table_L2) <- row.names(meta_L2)
OTU_L2 <- otu_table(feature_table_L2, taxa_are_rows = TRUE)

physeq_L2 <- phyloseq(OTU_L2, TAX_L2, META_L2)
physeq_L2 <- subset_taxa(physeq_L2, Kingdom != "Unassigned")
physeq_L2 <- subset_taxa(physeq_L2, Class != "D_2__Chloroplast")
physeq_L2 <- subset_taxa(physeq_L2, Class != "D_2__Chloroplast")
physeq_L2 <- subset_taxa(physeq_L2, Family != "D_4__Mitochondria")
physeq_L2

tax_L3 <- as.matrix(tax_L3)
TAX_L3 <- tax_table(tax_L3)
meta_L3 <- as.data.frame(metadata_L3)
META_L3 <- sample_data(meta_L3)
colnames(feature_table_L3) <- row.names(meta_L3)
OTU_L3 <- otu_table(feature_table_L3, taxa_are_rows = TRUE)

physeq_L3 <- phyloseq(OTU_L3, TAX_L3, META_L3)
physeq_L3 <- subset_taxa(physeq_L3, Kingdom != "Unassigned")
physeq_L3 <- subset_taxa(physeq_L3, Class != "D_2__Chloroplast")
physeq_L3 <- subset_taxa(physeq_L3, Class != "D_2__Chloroplast")
physeq_L3 <- subset_taxa(physeq_L3, Family != "D_4__Mitochondria")
physeq_L3

physeq_merge_bitou <- merge_phyloseq(physeq_L1, physeq_L2, physeq_L3)
#result physeq_sample raremax=___ [6828x68]
physeq_Bitou <- subset_samples(physeq_merge_bitou,Type!="Water")


####### Turtle island water #######
getwd()
setwd("/Users/haijin/Desktop/NTU/01_TI/TI_combine_data")
feature_table_2011 <- read.csv(file="feature-table_w_2011.csv", header = T, check.names = TRUE, row.names = 1)
tax_2011<-read.csv(file="taxa_table_2011.csv",header=T,check.names = T, row.names = 1)
metadata_2011 <- read.csv(file = "metadata_2011.csv", header = T, check.names = T, row.names = 1)

tax_2011 <- as.matrix(tax_2011)
TAX_2011 <- tax_table(tax_2011)
meta_2011 <- as.data.frame(metadata_2011)
META_2011 <- sample_data(meta_2011)
colnames(feature_table_2011) <- row.names(meta_2011)
OTU_2011 <- otu_table(feature_table_2011, taxa_are_rows = TRUE)

physeq_2011 <- phyloseq(OTU_2011, TAX_2011, META_2011)
physeq_2011 <- subset_taxa(physeq_2011, Kingdom != "Unassigned")
physeq_2011 <- subset_taxa(physeq_2011, Class != "D_2__Chloroplast")
physeq_2011 <- subset_taxa(physeq_2011, Class != "D_2__Chloroplast")
physeq_2011 <- subset_taxa(physeq_2011, Family != "D_4__Mitochondria")
physeq_2011

physeq_2011_water <- subset_samples(physeq_2011,Type=="Water")
physeq_TI_water <- merge_phyloseq(physeq_2011_water, physeq_L2_water)
#until now should have physeq_2011_water (6), physeq_TI_water (18)

#draw physeq_TI_water rarefaction curve (raremax = 16206)
col <- c("darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
raremax = min(sample_sums(physeq_TI_water)) #get the minimum point
rarecurve(t(otu_table(physeq_TI_water)), step = 50, cex=0.5, xlab = "Sequence sample size", ylab = "Species richness", label = T, col=col)+
  geom_text(title("Water Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_physeq_TI_water = rarefy_even_depth(physeq_TI_water, rngseed = 1, sample.size = raremax, replace = F) 
####### End of Turtle island water #######



#draw physeq_Bitou rarefaction Curve (raremax = 14016)
col <- c("darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
raremax = min(sample_sums(physeq_Bitou)) #get the minimum point

tab <- otu_table(physeq_Bitou)
class(tab) <- "matrix"
tab <- t(tab)
rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=T, col=col)+
  geom_text(title("Bitou Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_physeq_Bitou = rarefy_even_depth(physeq_Bitou, rngseed = 1, sample.size = raremax, replace = F)

#01_Relative abundance
##physeq_Bitou
#Phylum
Stability_phylum<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_phylum[Stability_phylum$Abundance<0.01,]$Phylum="Others"

colourCount = length(unique(Stability_phylum$Phylum))
colourCount
colors_phylum_n <- randomColor(count = 21)
colors_phylum_n <- distinctColorPalette(21)

metadata_L1L2 <- merge(metadata_L1, metadata_L2, by=c("Place","Type","Time","Sample.1"), all=T, sort=F)
metadata_L1toL3 <- merge(metadata_L1L2, metadata_L3, by=c("Place","Type","Time","Sample.1"), all=T, sort=F)

Stability_phylum$Time <- factor(Stability_phylum$Time,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022", "July, 2022", "Sep, 2022"))
x<- c(metadata_L1toL3$Sample.1)
Stability_phylum$Sample.1 = factor(Stability_phylum$Sample.1, levels = x)

PlotPhylum<-ggplot(Stability_phylum, aes(x=Sample,y=Abundance, fill=Phylum))+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
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

#-----subset proteo-----
physeq_CT_NeA_Proteo <- subset_taxa(rarefied_physeq_Bitou, Phylum == "D_1__Proteobacteria")
Stability_class_proteo<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_class_proteo[Stability_class_proteo$Abundance<0.01,]$Class="Others"

colourCount = length(unique(Stability_class_proteo$Class))
colourCount
colors_proteo_n <- randomColor(count = 31)
colors_proteo_n <- distinctColorPalette(31)

Stability_class_proteo$Time <- factor(Stability_class_proteo$Time,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022","July, 2022", "Sep, 2022"))
x<- c(metadata_L1toL3$Sample.1)
Stability_class_proteo$Sample.1 = factor(Stability_class_proteo$Sample.1, levels = x)

PlotClassProteo<-ggplot(Stability_class_proteo, aes(x=Sample,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=colors_proteo_n)
PlotClassProteo

#Class
Stability_class<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_class[Stability_class$Abundance<0.01,]$Class="Others"

colourCount = length(unique(Stability_class$Class))
colourCount
colors_class_n <- randomColor(count = 31)
colors_class_n <- distinctColorPalette(31)

Stability_class$Time <- factor(Stability_class$Time,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022","July, 2022", "Sep, 2022"))
x<- c(metadata_L1toL3$Sample.1)
Stability_class$Sample.1 = factor(Stability_class$Sample.1, levels = x)

PlotClass<-ggplot(Stability_class, aes(x=Sample,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(values=colors_class_n)
PlotClass

#Order
Stability_order<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Order") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_order[Stability_order$Abundance<0.01,]$Order="Others"

colourCount = length(unique(Stability_order$Order))
colourCount
colors_order_n <- randomColor(count = 52)
colors_order_n <- distinctColorPalette(52)

Stability_order$Time <- factor(Stability_order$Time,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022","July, 2022", "Sep, 2022"))
x<- c(metadata_L1toL3$Sample.1)
Stability_order$Sample.1 = factor(Stability_order$Sample.1, levels = x)

PlotOrder<-ggplot(Stability_order, aes(x=Time,y=Abundance,fill=Order))+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=2))+
  scale_fill_manual(values=colors_order_n)
PlotOrder

#Family
Stability_family<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_family[Stability_family$Abundance<0.01,]$Family="Others"

colourCount = length(unique(Stability_family$Family))
colourCount
colors_family_n <- randomColor(count = 73)
colors_family_n <- distinctColorPalette(73)

Stability_family$Time <- factor(Stability_family$Time,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022","July, 2022", "Sep, 2022"))
x<- c(metadata_L1toL3$Sample.1)
Stability_family$Sample.1 = factor(Stability_family$Sample.1, levels = x)

PlotFamily<-ggplot(Stability_family, aes(x=Sample,y=Abundance,fill=Family))+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_family_n)+
  theme(legend.text=element_text(size=6))
PlotFamily

#Genus
Stability_genus<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_genus
Stability_genus[Stability_genus$Abundance<0.01,]$Genus="Others"

colourCount = length(unique(Stability_genus$Genus))
colourCount
colors_genus_n <- randomColor(count = 88)
colors_genus_n <- distinctColorPalette(88)

Stability_genus$Time <- factor(Stability_genus$Time,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022","July, 2022", "Sep, 2022"))
x<- c(metadata_L1toL3$Sample.1)
Stability_genus$Sample.1 = factor(Stability_genus$Sample.1, levels = x)

PlotGenus<-ggplot(Stability_genus, aes(x=Sample,y=Abundance,fill=Genus))+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_genus_n)+
  theme(legend.text=element_text(size=6))
PlotGenus

##pie chart
Stability_genus<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_genus

##### Sep, 2020 #####
sep20 <- subset(Stability_genus, Time=="Sep, 2020")
sep20 <- mutate(sep20, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
sep20_endo_sum <- sum(subset(sep20, new_group=="D_5__Endozoicomonas")$Abundance)
sep20_endo_sum*10 ##71.02938
sep20_others_sum <- sum(subset(sep20, new_group=="Others")$Abundance)
sep20_others_sum*10 ##28.97062

sep20_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(71.03,28.97))

sep20_fig <- 
  ggplot(sep20_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "Sep, 2020")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
#+guides(fill=F)
sep20_fig

##### Jan, 2021 #####
jan21 <- subset(Stability_genus, Time=="Jan, 2021")
jan21 <- mutate(jan21, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
jan21_endo_sum <- sum(subset(jan21, new_group=="D_5__Endozoicomonas")$Abundance)
(jan21_endo_sum/11)*100 ##29.57698
jan21_others_sum <- sum(subset(jan21, new_group=="Others")$Abundance)
(jan21_others_sum/11)*100 ##70.42302

jan21_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(29.58,70.42))

jan21_fig <- 
  ggplot(jan21_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "Jan, 2021")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
jan21_fig

##### May, 2021 #####
may21 <- subset(Stability_genus, Time=="May, 2021")
may21 <- mutate(may21, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
may21_endo_sum <- sum(subset(may21, new_group=="D_5__Endozoicomonas")$Abundance)
(may21_endo_sum/12)*100 ##21.4759
may21_others_sum <- sum(subset(may21, new_group=="Others")$Abundance)
(may21_others_sum/12)*100 ##78.5241

may21_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(21.48,78.52))

may21_fig <- 
  ggplot(may21_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "May, 2021")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
may21_fig

##### Aug, 2021 #####
aug21 <- subset(Stability_genus, Time=="Aug, 2021")
aug21 <- mutate(aug21, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
aug21_endo_sum <- sum(subset(aug21, new_group=="D_5__Endozoicomonas")$Abundance)
(aug21_endo_sum/13)*100 ##12.12153
aug21_others_sum <- sum(subset(aug21, new_group=="Others")$Abundance)
(aug21_others_sum/13)*100 ##87.87847

aug21_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(12.12,87.88))

aug21_fig <- 
  ggplot(aug21_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "Aug, 2021")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
aug21_fig

##### Dec, 2021 #####
dec21 <- subset(Stability_genus, Time=="Dec, 2021")
dec21 <- mutate(dec21, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
dec21_endo_sum <- sum(subset(dec21, new_group=="D_5__Endozoicomonas")$Abundance)
dec21_endo_sum*10 ##11.59061
dec21_others_sum <- sum(subset(dec21, new_group=="Others")$Abundance)
dec21_others_sum*10 ##88.40939

dec21_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(11.6,88.4))

dec21_fig <- 
  ggplot(aug21_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "Dec, 2021")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
dec21_fig

##### May, 2022 #####
may22 <- subset(Stability_genus, Time=="May, 2022")
may22 <- mutate(may22, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
may22_endo_sum <- sum(subset(may22, new_group=="D_5__Endozoicomonas")$Abundance)
may22_endo_sum*10 ##14.91289
may22_others_sum <- sum(subset(may22, new_group=="Others")$Abundance)
may22_others_sum*10 ##85.08711

may22_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(14.91,85.09))

may22_fig <- 
  ggplot(may22_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "May, 2022")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
may22_fig

##### July, 2022 #####
july22 <- subset(Stability_genus, Time=="July, 2022")
july22 <- mutate(july22, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
july22_endo_sum <- sum(subset(july22, new_group=="D_5__Endozoicomonas")$Abundance)
july22_endo_sum*10 ##18.383
july22_others_sum <- sum(subset(july22, new_group=="Others")$Abundance)
july22_others_sum*10 ##81.617

july22_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(18.38,81.62))

july22_fig <- 
  ggplot(july22_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "July, 2022")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
july22_fig

##### Sep, 2022 #####
sep22 <- subset(Stability_genus, Time=="Sep, 2022")
sep22 <- mutate(sep22, new_group = ifelse(Genus!="D_5__Endozoicomonas", "Others", Genus))
sep22_endo_sum <- sum(subset(sep22, new_group=="D_5__Endozoicomonas")$Abundance)
sep22_endo_sum*10 ##26.26431
sep22_others_sum <- sum(subset(sep22, new_group=="Others")$Abundance)
sep22_others_sum*10 ##73.73569

sep22_new = data.frame(
  Genus = c("Endozoicomonas", "Others"),
  Percentage = c(26.26,73.74))

sep22_fig <- 
  ggplot(sep22_new, aes(x="", y=Percentage, fill=Genus))+ 
  geom_bar(stat="identity")+
  coord_polar("y", start=0)+
  scale_fill_manual(values=c("gold1","gray86"))+
  geom_text(aes(label=paste0(Percentage, "%")), position = position_stack(vjust = 0.5))+
  theme_void()+
  labs(title = "Sep, 2022")+
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank(), 
        legend.text = element_text(size = 12))
  #+guides(fill=F)
sep22_fig

ggarrange(sep20_fig, jan21_fig, may21_fig, aug21_fig, dec21_fig, may22_fig, july22_fig, sep22_fig, 
          ncol = 4, nrow = 2, common.legend = T, vjust=0.5)

###Additional codes
#Stability_genus3 <- with(Stability_genus2[Stability_genus2$Genus == "D_5__Endozoicomonas", ], table(Abundance > 0.01, Sample))
#  theme(legend.title=element_blank(), legend.position="right", legend.text = element_text(size = 8))+ 
#  theme(legend.text=element_text(size=6))
#geom_text(aes(label = paste0(x, scales::percent(Stability_genus))))
#sum(sep$new_group=="D_5__Endozoicomonas") ##10
#sum(sep$new_group=="Others") ##11710
#scale_x_discrete(NULL, expand = c(0,0))+
#scale_y_continuous(NULL, expand = c(0,0))+

#Species
Stability_species<-rarefied_physeq_Bitou %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_species[Stability_species$Abundance<0.01,]$Species="Others"
colourCount = length(unique(Stability_species$Species))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colours <- c("#FFA8BB","#F0A3FF", "#993F00","#4C005C","#0075DC","#FFCC99",
             "#00998F","#94FFB5","#8F7C00","#9DCC00","cadetblue","#003380","#FFA405",
             "#2BCE48","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000",
             "cadetblue","#FFA405","#003380","darkgoldenrod1","#426600","#FF0010","#5EF1F2",
             "#FFFF00","#740AFF","#FFA8BB","#990000","#808080", "gold1", "chocolate1");
Stability_species$Month <- factor(Stability_species$Month,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021"))
Stability_species$Sample = factor(Stability_species$Sample, levels = x)

PlotSpecies<-ggplot(Stability_species, aes(x=Sample,y=Abundance,fill=Species))+
  geom_bar(stat="identity")+
  facet_wrap(~ Month, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values = getPalette(colourCount))
  scale_fill_manual(values=colours)+
  guides(fill=guide_legend(ncol=1))
PlotSpecies


##top 20
colors_top20 <- randomColor(count = 21)
colors_top20 <- distinctColorPalette(21)

merge_less_than_top <- function(rarefied_physeq_Bitou, top=20){
  transformed <- transform_sample_counts(rarefied_physeq_Bitou, function(x) x/sum(x)) 
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),] 
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:8] <- "Others"}
  }
  return(merged)
}

ps.gen <- tax_glom(rarefied_physeq_Bitou, "Family..Genus..Species")
ps.gen.top20 <- merge_less_than_top(ps.gen, top=20)

plot_bar(ps.gen.top20, "Sample", fill = "Family..Genus..Species")+
  geom_bar(stat="identity")+
  facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  #facet_nested(~ time + species, scales="free_x", nest_line=element_line(linetype=1))+
  theme(strip.background = element_rect(fill="white"))+
  theme(panel.background = element_rect(fill="white", colour="gray", size=0.5, linetype="solid"))+
  theme(panel.spacing = unit(0.3,"cm"))+
  #theme(strip.text = element_text(face="italic"))+
  labs(x="Samples",y="Relative abundance")+
  scale_fill_manual(values=colors_top20)+
  guides(fill=guide_legend(ncol=1))+
  theme(axis.text.x = element_text(angle=90, hjust=1, size=6))+
  theme(legend.text=element_text(size=8))

#-----subset Endozoicomonas-----
physeq_bitou_endo <- subset_taxa(rarefied_physeq_Bitou, Genus=="D_5__Endozoicomonas")
Stability_species_endo<-physeq_bitou_endo %>%
  tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_species_endo[Stability_species_endo$Abundance<0.01,]$Species="Others"
colourCount = length(unique(Stability_species_endo$Species))
colours <- c("darkgoldenrod1", "cadetblue1", "darkolivegreen4","chocolate1","darkcyan","chartreuse","deeppink2","darkmagenta");
Stability_species$Month <- factor(Stability_species$Month,levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021"))
Stability_species$Sample = factor(Stability_species$Sample, levels = x)

PlotSoeciesEndo<-ggplot(Stability_species_endo, aes(x=Sample,y=Abundance,fill=Species))+
  geom_bar(stat="identity")+
  facet_wrap(~ Month, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_fill_manual(values=colours)
PlotSoeciesEndo

#02_Alpha diversity
sample_data(rarefied_physeq_Bitou)$Month <- factor(sample_data(rarefied_physeq_Bitou)$Month, levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021"))
set.seed(114)
richness <- estimate_richness(rarefied_physeq_Bitou)
head(richness)
rich_all<-plot_richness(rarefied_physeq_Bitou, x="Month", color="Month",
                        measures = c("Observed", "Shannon", "Chao1", "Simpson")) #google it
rich_all+geom_boxplot(stat="boxplot",position="dodge2") + 
  theme(axis.text.x = element_text(angle=90, hjust=1))


#03_PcoA #jaccard
sample_data(rarefied_physeq_Bitou)$Time <- factor(sample_data(rarefied_physeq_Bitou)$Time, levels=c("Sep, 2020","Jan, 2021","May, 2021","Aug, 2021","Dec, 2021","May, 2022","July, 2022","Sep, 2022"))
dist_jaccard_binary = distance(rarefied_physeq_Bitou, method = "jaccard", binary = T)
pcoa_jaccard_binary = ordinate(rarefied_physeq_Bitou, "PCoA", dist_jaccard_binary)
#pcoa_jaccard_binary = ordinate(rarefied_physeq, method = "PCoA", distance = "jaccard")

#R and p value
pcoa_anosim = anosim(dist_jaccard_binary, grouping=sample_data(rarefied_physeq_Bitou)$Time,permutations = 999, distance = "jaccard")
summary(pcoa_anosim) 

col <- c("darkmagenta", "darkcyan", "chocolate1","darkolivegreen4","darkseagreen1","deepskyblue4","gold1","deeppink2")
p_pcoa_jaccard_binary = plot_ordination(rarefied_physeq_Bitou, pcoa_jaccard_binary, color="Time") + scale_color_manual(values = col)
plot(p_pcoa_jaccard_binary)

p_pcoa_jaccard_binary2 = p_pcoa_jaccard_binary +
  geom_point(size=4) +
  #geom_text(aes(label = Month), vjust = 2.5) +
  #scale_color_hue(l=40) +
  coord_fixed(ratio=1) +
  annotate("text",x=-0.3,y=0.6,hjust=0,vjust=0,label="italic(R)==0.4105",parse=T,size=4) +
  annotate("text",x=-0.3,y=0.57,hjust=0,vjust=0,label="italic(p)==0.001",parse=T,size=4) 
#draw circle #+stat_ellipse()
plot(p_pcoa_jaccard_binary2)


##draw rarefied_physeq_TI_water
#Phylum
Stability_phylum<-rarefied_physeq_TI_water %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_phylum[Stability_phylum$Abundance<0.01,]$Phylum="Others"
colourCount = length(unique(Stability_phylum$Phylum))
#getPalette = colorRampPalette(brewer.pal(11, "BrBG"))
colours <- c("darkcyan", "chocolate1","#740AFF","#4C005C",
             "#FFCC99","#8F7C00","deeppink2","#9DCC00","#2BCE48",
             "#808080","cadetblue","#003380","gold1");
Stability_phylum$Month <- factor(Stability_phylum$Month,levels=c("Nov","Feb","July"))
PlotPhylum<-ggplot(Stability_phylum, aes(x=Sample,y=Abundance,fill=Phylum))+
  geom_bar(stat="identity")+
  facet_grid(Month ~ Place, space = "free", scales = "free")+
  #facet_wrap(~ Place, scales = "free_x", nrow = 1)+ #改Type
  #facet_wrap(~ Group + Time, scales = "free_x")
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values = getPalette(colourCount))
  scale_fill_manual(values=colours)
PlotPhylum


#"#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000",
#"#FFFF00","darkgoldenrod1", "#F0A3FF", "#0075DC","darkolivegreen4","darkseagreen1","deepskyblue4",
#"#993F00","darkmagenta"

#-----subset proteo-----
physeq_water_Proteo <- subset_taxa(rarefied_physeq_TI_water, Phylum == "D_1__Proteobacteria")
Stability_class_proteo<-physeq_water_Proteo %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_class_proteo[Stability_class_proteo$Abundance<0.01,]$Class="Others"
colourCount = length(unique(Stability_class_proteo$Class))
colours <- c("darkgoldenrod1", "darkcyan", "chocolate1","darkolivegreen4","deepskyblue4","deeppink2","#FFA8BB","#426600");
Stability_class_proteo$Month <- factor(Stability_class_proteo$Month,levels=c("Nov","Feb","July"))
PlotClassProteo<-ggplot(Stability_class_proteo, aes(x=Sample,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  facet_grid(Month ~ Place, space = "free", scales = "free")+
  #facet_wrap(~ Place, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_fill_manual(values=colours)
PlotClassProteo

#Class
Stability_class<-rarefied_physeq_TI_water %>%
  tax_glom(taxrank = "Class") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_class[Stability_class$Abundance<0.01,]$Class="Others"
colourCount = length(unique(Stability_class$Class))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colours <- c("darkgoldenrod1", "darkcyan", "chocolate1","darkolivegreen4","deepskyblue4","deeppink2",
             "#9DCC00","cadetblue","#003380",
             "#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F",
             "#740AFF","#990000","#FFFF00","darkgoldenrod1", "darkcyan", "chocolate1",
             "darkolivegreen4","darkseagreen1","deepskyblue4","gold1","deeppink2","darkmagenta",
             "#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99",
             "#808080","#94FFB5","#8F7C00","#9DCC00","cadetblue");
Stability_class$Month <- factor(Stability_class$Month,levels=c("Nov","Feb","July"))
PlotClass<-ggplot(Stability_class, aes(x=Sample,y=Abundance,fill=Class))+
  geom_bar(stat="identity")+
  facet_grid(Month ~ Place, space = "free", scales = "free")+
  #facet_wrap(~ Place, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=1))+
  #scale_fill_manual(values = getPalette(colourCount))+
  scale_fill_manual(values=colours)
PlotClass

#Order
Stability_order<-rarefied_physeq_TI_water %>%
  tax_glom(taxrank = "Order") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_order[Stability_order$Abundance<0.01,]$Order="Others"
colourCount = length(unique(Stability_order$Order))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99",
             "#808080","#94FFB5","#8F7C00","#9DCC00","cadetblue","#003380",
             "#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F",
             "#740AFF","#990000","#FFFF00","darkolivegreen4", "deeppink", "chocolate1",
             "darkolivegreen4","darkseagreen1","deepskyblue4","deeppink2","gold1","darkmagenta",
             "aquamarine4","#FFFF00","darkgoldenrod1", "#F0A3FF", "#0075DC","darkolivegreen4",
             "darkseagreen1","deepskyblue4","#993F00");
Stability_order$Month <- factor(Stability_order$Month,levels=c("Nov","Feb","July"))
PlotOrder<-ggplot(Stability_order, aes(x=Sample,y=Abundance,fill=Order))+
  geom_bar(stat="identity")+
  facet_grid(Month ~ Place, space = "free", scales = "free")+
  #facet_wrap(~ Place, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values = getPalette(colourCount))+
  scale_fill_manual(values=colours)+
  guides(fill=guide_legend(ncol=2))
PlotOrder

#Family
Stability_family<-rarefied_physeq_TI_water %>%
  tax_glom(taxrank = "Family") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_family[Stability_family$Abundance<0.01,]$Family="Others"
colourCount = length(unique(Stability_family$Family))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99",
             "#808080","gold1","#8F7C00","#9DCC00","cadetblue","#003380",
             "#740AFF","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F",
             "#FFA405","#990000","#FFFF00","skyblue", "darkcyan", "chocolate1",
             "darkolivegreen4","darkseagreen1","deepskyblue4","#94FFB5","deeppink2","darkmagenta",
             "aquamarine4","darkgoldenrod4","coral2","darksalmon","chartreuse","deeppink2",
             "lightgoldenrod","burlywood","burlywood4","#F0A3FF", "#0075DC",
             "#993F00","#4C005C","#2BCE48","deepskyblue");
Stability_family$Month <- factor(Stability_family$Month,levels=c("Nov","Feb","July"))
PlotFamily<-ggplot(Stability_family, aes(x=Sample,y=Abundance,fill=Family))+
  geom_bar(stat="identity")+
  facet_grid(Month ~ Place, space = "free", scales = "free")+
  #facet_wrap(~ Place, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values = getPalette(colourCount))+
  scale_fill_manual(values=colours)+
  guides(fill=guide_legend(ncol=2))
PlotFamily

#Genus
Stability_genus<-rarefied_physeq_TI_water %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_genus[Stability_genus$Abundance<0.01,]$Genus="Others"
colourCount = length(unique(Stability_genus$Genus))
#getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colours <- c("#F0A3FF", "#0075DC", "#993F00","#4C005C","#2BCE48","#FFCC99",
             "#808080","gold1","#8F7C00","#9DCC00","cadetblue","#003380",
             "#740AFF","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F",
             "#FFA405","#990000","#FFFF00","skyblue", "darkcyan", "chocolate1",
             "darkolivegreen4","darkseagreen1","deepskyblue4","#94FFB5","deeppink2","darkmagenta",
             "aquamarine4","darkgoldenrod4","coral2","darksalmon","chartreuse","deeppink2",
             "lightgoldenrod","burlywood","burlywood4","gold1", "#0075DC",
             "#993F00","#4C005C","#2BCE48","deepskyblue","yellow","tomato","yellowgreen","turquoise");
Stability_genus$Month <- factor(Stability_genus$Month,levels=c("Nov","Feb","July"))
PlotGenus<-ggplot(Stability_genus, aes(x=Sample,y=Abundance,fill=Genus))+
  geom_bar(stat="identity")+
  facet_grid(Month ~ Place, space = "free", scales = "free")+
  #facet_wrap(~ Place, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  #scale_fill_manual(values = getPalette(colourCount))+
  scale_fill_manual(values=colours)+
  guides(fill=guide_legend(ncol=2))
PlotGenus

#Species
colours <- c("lightpink","#990000","#00998F","navy","chocolate1",
             "yellow1","darkseagreen1","deepskyblue4","lightskyblue1","deeppink2",
             "darkmagenta","darkolivegreen4","darkgoldenrod4","darksalmon","chartreuse", 
             "#5EF1F2","antiquewhite","cornflowerblue","mediumorchid","lightgoldenrod1","gray");

merge_less_than_top <- function(rarefied_physeq_TI_water, top=20){
  transformed <- transform_sample_counts(rarefied_physeq_TI_water, function(x) x/sum(x)) 
  otu.table <- as.data.frame(otu_table(transformed))
  otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),] 
  otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
  merged <- merge_taxa(transformed, otu.list, 1)
  for (i in 1:dim(tax_table(merged))[1]){
    if (is.na(tax_table(merged)[i,2])){
      taxa_names(merged)[i] <- "Others"
      tax_table(merged)[i,1:9] <- "Others"}
  }
  return(merged)
}

ps.gen <- tax_glom(rarefied_physeq_TI_water, "Family..Genus..Species")

ps.gen.top20 <- merge_less_than_top(ps.gen, top=20)
sample_data(ps.gen.top20)$Month <- factor(sample_data(ps.gen.top20)$Month, levels=c("Nov","Feb","July"))

plot_bar(ps.gen.top20, "Sample", fill = "Family..Genus..Species")+
  geom_bar(stat="identity")+
  facet_nested( ~ Month, scales = "free_x", nest_line = element_line(linetype = 1))+
  theme(panel.background = element_rect(fill = "white",colour = "lightgray",size = 0.3, linetype = "solid"))+
  theme(panel.spacing = unit(0.1,"cm"))+
  theme(strip.text.x = element_text(size=7))+
  theme(axis.text.x = element_text(size=8))+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  scale_fill_manual(values=colours)+
  guides(fill=guide_legend(ncol=1))

#02_Alpha diversity
set.seed(114)
richness <- estimate_richness(rarefied_physeq_TI_water)
head(richness)
rich_all<-plot_richness(rarefied_physeq_TI_water,x="Place",color = "Place",
                        measures = c("Observed", "Shannon", "Chao1", "Simpson")) #google it
rich_all+geom_boxplot(stat="boxplot",position="dodge2")


#03_PcoA #jaccard
sample_data(rarefied_physeq_TI_water)$Month <- factor(sample_data(rarefied_physeq_TI_water)$Month, levels=c("Nov","Feb","July"))

dist_jaccard_binary = distance(rarefied_physeq_TI_water, method = "bray", binary = T)
pcoa_jaccard_binary = ordinate(rarefied_physeq_TI_water, "PCoA", dist_jaccard_binary)

#R and p value
pcoa_anosim = anosim(dist_jaccard_binary, grouping=sample_data(rarefied_physeq_TI_water)$Month,permutations = 999, distance = "bray")
summary(pcoa_anosim) 

col <- c("darkmagenta", "darkcyan", "chocolate1","darkolivegreen4","gold1","deepskyblue4","gold1","deeppink2")
p_pcoa_jaccard_binary = plot_ordination(rarefied_physeq_TI_water, pcoa_jaccard_binary, color="Month") + scale_color_manual(values = col)
plot(p_pcoa_jaccard_binary)

p_pcoa_jaccard_binary2 = p_pcoa_jaccard_binary +
  geom_point(size=4) +
  #geom_text(aes(label = Month), vjust = 2.5) +
  coord_fixed(ratio=1) +
  annotate("text",x=-0.3,y=0.6,hjust=0,vjust=0,label="italic(R)==0.5181",parse=T,size=4) +
  annotate("text",x=-0.3,y=0.57,hjust=0,vjust=0,label="italic(p)==0.001",parse=T,size=4) 
plot(p_pcoa_jaccard_binary2)



#Bray
pcoa <- ordinate(rarefied_physeq_TI_water, method = "PCoA", distance = "bray")
plot_ordination(rarefied_physeq_TI_water, pcoa, color="Place", shape = 'Group')

#------------------------------
#04_Heatmap
#gpt <- subset_taxa(rarefied_physeq, Phylum == "D_1__Proteobacteria")
gpt <- rarefied_physeq
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:50]), gpt) #top 100 most abundant Proteobacteria taxa
plot_heatmap(gpt, taxa.label="Order")
#plot_heatmap(gpt, "NMDS", "bray",)
