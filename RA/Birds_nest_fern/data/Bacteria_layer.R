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
library(mobr)

#Phyloseq data preparation----
#Import data with OTU composition for each sample 
feature_table <- read.csv(file = '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/feature-table.txt', sep = '', header = T, check.names = TRUE, row.names = 1) %>% select(-ID, -taxonomy)
feature_table <- feature_table %>% rename_with(~ mgsub::mgsub(colnames(feature_table), c("X", "_V6V8", "_"), c("", "", "-"))) %>% otu_table(taxa_are_rows = TRUE)


#Import data with taxonomy information od each OTU 
tax <- read.csv(file='/Users/yupeitseng/Documents/RA/birds_nest_fern/data/Bacteria/taxonomy.csv', sep = ' ', header = F, check.names = T, row.names = 1)
tax <- tax[-1, ] %>% mutate(all = rownames(tax[-1, ]))
tax[c('ID', 'Kingdom')] <- str_split_fixed(tax$all, ',', 2)
rownames(tax) <- tax[, "ID"]
tax <- tax[, c(9, 1:6)] %>% rename_with(~ c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>% as.matrix() %>% tax_table()

#Import data with sample information 
samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv', header = T) %>% select(-1)
samples.out <- as.data.frame(colnames(feature_table)) %>% rename("sample_ID" = 1)
samdf <- left_join(samples.out, samdf, by = "sample_ID")
rownames(samdf) <- colnames(feature_table)

#construct phyloseq object
ps <- phyloseq(otu_table(feature_table, taxa_are_rows = T), sample_data(samdf), tax_table(tax))
ps <- subset_taxa(ps, Kingdom != "Unassigned") %>% 
  subset_taxa(Class != "D_2__Chloroplast") %>% 
  subset_taxa(Family != "D_4__Mitochondria")
ps <- prune_samples(!(sample_names(ps)  %in% c("NC1", "NC2", "NC3", "NC4","9-U")), ps) # Remove NC sample
order(rowSums(ps@otu_table))
order(unique(unlist(as.vector(ps@otu_table))))
#draw ps rarefaction Curve (raremax = 5235) w/o NA
col <- randomColor(count = 69)
raremax <- min(sample_sums(ps)) #get the minimum point

tab <- otu_table(ps) %>% t()
class(tab) <- "matrix"
rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=F, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_ps <- rarefy_even_depth(ps, rngseed = 1, sample.size = raremax, replace = F)
order(unique(rarefied_ps@otu_table))
#01_Relative abundance----
#Kingdom, Phylum, Class, Order, Family
pdf(file = "/Users/yupeitseng/Documents/RA/birds_nest_fern/plot/relative_aundance_bact.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 9) # The height of the plot in inches

for(i in c('Kingdom', 'Phylum', 'Class', 'Order', 'Family')){
  relative_abun_plot(i)
}

relative_abun_plot <- function(tax_level = tax_level){
  Stability <-rarefied_ps %>%
    tax_glom(taxrank = tax_level) %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt()
  Stability[Stability$Abundance<0.01, tax_level] <- "Others"
  colourCount <- Stability %>% select(tax_level) %>% unique() %>% nrow()
  colors_n <- distinctColorPalette(colourCount)
  
  plot_tax <-ggplot(Stability, aes(x=Sample,y=Abundance, fill=get(tax_level)))+
    geom_bar(stat="identity")+
    # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
    labs(x="Samples",y="Relative abundance")+
    theme(axis.text.x = element_text(angle=90, hjust=1))+
    guides(fill=guide_legend(ncol=1))+
    scale_fill_manual(values=colors_n)
  print(plot_tax)
}

#Genus
Stability_genus<-rarefied_ps %>%
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()

Stability_genus$Family_Genus <- paste(Stability_genus$Family, Stability_genus$Genus)
Stability_genus[Stability_genus$Abundance<0.01,]$Family_Genus <- "Others"

colourCount <- Stability_genus %>% select("Family_Genus") %>% unique() %>% nrow()
colors_genus_n <- distinctColorPalette(colourCount)

PlotGenus<-ggplot(Stability_genus, aes(x=Sample,y=Abundance,fill=Family_Genus))+
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
Stability_species$Genus_Species <- paste(Stability_species$Genus, Stability_species$Species)
Stability_species[Stability_species$Abundance<0.01,]$Genus_Species <- "Others"

colourCount <- Stability_species %>% select("Genus_Species") %>% unique() %>% nrow()
colors_species_n <- distinctColorPalette(colourCount)

PlotSpecies<-ggplot(Stability_species, aes(x=Sample,y=Abundance,fill=Genus_Species))+
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
bray.dist_m<-vegdist(t(ps.prop_m@otu_table), method='bray')
adonis2(bray.dist_m ~  as.factor(ps.prop_m@sam_data$layer) + ps.prop_m@sam_data$C_N + ps.prop_m@sam_data$pH, permutations=9999)
permutest(betadisper(bray.dist_m, group=ps.prop_m@sam_data$layer))
pairwise.adonis(bray.dist_m, as.factor(ps.prop_m@sam_data$layer))

adonis2(bray.dist_m ~  as.factor(ps.prop_m@sam_data$layer) + ps.prop_m@sam_data$pH + ps.prop_m@sam_data$C_N, permutations=9999)

#large
ps.prop_l <- prune_samples(sample_names(ps.prop) %in% rownames(ps.prop@sam_data)[which(ps.prop@sam_data$fern.size == 'L')], ps.prop)
ord.nmds.bray <- ordinate(ps.prop_l, method="NMDS", distance="bray")
plot_ordination(ps.prop_l, ord.nmds.bray, color = 'layer') + stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  geom_point(cex = 5)+
  labs(colour = "Layer") +
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#DDE5B6", "#F1DCA7", "#997B66"), label = c("Lower", "Middle", "Outer", "Upper"))

ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color = 'pH') +
  geom_point(cex = 5)+
  labs(colour = "pH") +
  theme_bw()+
  theme(text = element_text(size = 24))+
  scale_colour_continuous(type = "viridis")

bray.dist_l <- vegdist(t(ps.prop_l@otu_table), method='bray')
adonis2(bray.dist_l ~  as.factor(ps.prop_l@sam_data$layer), permutations=9999)
permutest(betadisper(bray.dist_l, group=ps.prop_l@sam_data$layer))
pairwise.adonis(bray.dist_l, as.factor(ps.prop_l@sam_data$layer))

adonis2(bray.dist_l ~  as.factor(ps.prop_l@sam_data$layer) + ps.prop_l@sam_data$pH + ps.prop_l@sam_data$C_N, permutations=9999, na.action = na.omit)

mod <- betadisper(bray.dist_l, group=ps.prop_l@sam_data$layer)
# extract the centroids and the site points in multivariate space.  
centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
vectors<-data.frame(group=mod$group,data.frame(mod$vectors))

# to create the lines from the centroids to each point we will put it in a format that ggplot can handle
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")

# create the convex hulls of the outermost points
grp1.hull<-seg.data[seg.data$group=="L",1:3][chull(seg.data[seg.data$group=="L",2:3]),]
grp2.hull<-seg.data[seg.data$group=="M",1:3][chull(seg.data[seg.data$group=="M",2:3]),]
grp3.hull<-seg.data[seg.data$group=="U",1:3][chull(seg.data[seg.data$group=="U",2:3]),]
grp4.hull<-seg.data[seg.data$group=="O",1:3][chull(seg.data[seg.data$group=="O",2:3]),]
all.hull<-rbind(grp1.hull,grp2.hull,grp3.hull, grp4.hull)

panel.e<-ggplot() + 
  geom_polygon(data=all.hull,aes(x=v.PCoA1,y=v.PCoA2, col = group),alpha=0,linetype="dashed") +
  # geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) + 
  geom_point(data=centroids[,1:3], aes(x=PCoA1,y=PCoA2,col=grps), size=5) + 
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2, col = group),size=2) +
  labs(x="",y="", col = "Layer") +
  coord_cartesian(xlim = c(-0.3,0.45), ylim = c(-0.3,0.4))

#pH range ~ species composition distance (bray)
rarefied_ps_lm <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'L'| rarefied_ps@sam_data$fern.size == 'M')], rarefied_ps) 
rarefied_ps_lm <- prune_samples(sample_names(rarefied_ps_lm) %in% rownames(rarefied_ps_lm@sam_data)[which(rarefied_ps_lm@sam_data$layer != 'O')], rarefied_ps_lm) 


ph_var <- rarefied_ps_lm@sam_data %>% group_by(fern_ID) %>% summarise(var = var(pH, na.rm = T), C_N_var = var(C_N, na.rm = T))
ph_var$fern_ID <- as.character(ph_var$fern_ID)
vec <- c()
for(i in sort(unique(rarefied_ps_lm@sam_data$fern_ID))){
  dat <- rarefied_ps_lm@otu_table[which(rownames(t(rarefied_ps_lm@otu_table)) %in% rownames(rarefied_ps_lm@sam_data)[which(rarefied_ps_lm@sam_data$fern_ID== i) ]), ]
  dis_mat <- vegdist(dat, 'bray')
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


ggscatter(ph_dis, x = 'var', y = 'dis_mean',
               add = "reg.line", conf.int = TRUE, color = "#DDE5B6", xlab = 'Variance of pH within an individual', ylab = 'Mean pairwised disimilarity \n within an individual')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =15)+
  geom_point(cex = 8, aes( col = fern.size))+
  scale_colour_manual(values = c("#ADC178", "#997B66"), aesthetics = c("colour", "fill"))+
  theme(text = element_text(size = 38))+
  guides(col=guide_legend(title="Size"))
  # theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "inches"))+
  # theme(
  #   panel.background = element_rect(fill='transparent'),
  #   plot.background = element_rect(fill='transparent', color=NA),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   legend.background = element_rect(fill='transparent'),
  #   legend.box.background = element_rect(fill='transparent')
  # )

# function for calculating pH variance, C:N ratio variance and mean pairwise disimilarity in each in dividual
ph_CN <- function(fil = fil){
  rarefied_ps_fil <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == "M")], rarefied_ps)
  ph_var <- rarefied_ps_fil@sam_data %>% group_by(fern_ID) %>% summarise(var = var(pH, na.rm = T), C_N_var = var(C_N, na.rm = T), size = unique(fern.size))
  ph_var$fern_ID <- as.character(ph_var$fern_ID)
  
  vec <- c()
  for(i in sort(unique(rarefied_ps_fil@sam_data$fern_ID))){
    dat <- rarefied_ps_fil@otu_table[which(rownames(rarefied_ps_fil@otu_table) %in% rownames(rarefied_ps_fil@sam_data)[which(rarefied_ps_fil@sam_data$fern_ID== 2) ]), ]
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
  # return(ph_dis)
  return(summary(lm(scale(dis_mean) ~  scale(C_N_var) + scale(var), data = ph_dis)))
}
ph_CN("M")


#04_SAR process (only medium and large)----
#mobr
rarefied_ps_lm <- prune_samples(sample_names(rarefied_ps) %in% rownames(rarefied_ps@sam_data)[which(rarefied_ps@sam_data$fern.size == 'L'| rarefied_ps@sam_data$fern.size == 'M')], rarefied_ps)
ps_lm <- prune_samples(sample_names(ps) %in% rownames(ps@sam_data)[which(ps@sam_data$fern.size == 'L'| ps@sam_data$fern.size == 'M')], ps)

inv_mob_in <- make_mob_in(t(rarefied_ps_lm@otu_table), rarefied_ps_lm@sam_data)
inv_stats <- get_mob_stats(inv_mob_in, group_var = "fern_ID", index = c("N", "S", "S_n", "S_PIE", "S_asymp"))
inv_mob_in_ps <- make_mob_in(t(ps_lm@otu_table), ps_lm@sam_data)
inv_stats_ps <- get_mob_stats(inv_mob_in_ps, group_var = "fern_ID", index = c("N", "S", "S_n", "S_PIE", "S_asymp"))
plot(inv_stats)

group <- inv_stats[["groups_stats"]] %>% 
  filter(effort == 5235 | is.na(effort) == T) %>% 
  select(-effort) %>% 
  filter(index %in% c('S_n', 'S_PIE', 'S_asymp'))
group$scale <- rep('gamma', nrow(group))
colnames(group)[1] <- 'fern_ID'
group <- left_join(group, samdf[, c(1, 3, 14)])


samples <- inv_stats[["samples_stats"]]
samples <- samples %>% filter(!(index %in% c('N', 'S', 'beta_S', 'S_asymp'))) %>% filter(is.na(index) == F) %>% select(-effort) 
samples <- samples %>% group_by(group, index) %>%  summarise(mean(value))
samples$scale <- rep(c('alpha', 'beta', 'alpha', 'beta'), 16)
colnames(samples)[1] <- 'fern_ID'
samples <- left_join(samples, samdf[, c(1, 3, 14)])
colnames(samples)[3] <- 'value'

dat <- rbind(group, samples)

# n_ref <-max(c(max(samples %>% filter(index == "N") %>% select(value)), min(samples %>% filter(index == "N") %>% select(value))))

alpha <- dat[which(dat$scale == 'alpha'), ]
beta<- dat[which(dat$scale == 'beta'), ]

#plotting
png(file="/Users/yupeitseng/Documents/RA/birds_nest_fern/plot/sar_mechanism_bact.png",
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