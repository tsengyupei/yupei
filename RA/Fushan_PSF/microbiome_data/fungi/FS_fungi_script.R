#packages----
#library(MicrobiotaProcess)
library(phyloseq)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggrepel)        
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
library(ggnewscale)
library(ape) 
library(spdep) 
library(ade4) 
library(adegraphics) 
library(adespatial) 
library(rstatix)
library(microbiome)
source("/Users/yupeitseng/Documents/RA/FS_PSF/data/fungi/sr.value.R")
source ('http://www.davidzeleny.net/anadat-r/doku.php/en:customized_functions:ordicenter?do=export_code&codeblock=0')

#import data (feature table, metadata, taxanomy)----
load('/Users/yupeitseng/Documents/RA/FS_PSF/data/fungi/track.RData')
load('/Users/yupeitseng/Documents/RA/FS_PSF/data/fungi/seqtab.nochim.RData')
load('/Users/yupeitseng/Documents/RA/FS_PSF/data/fungi/taxa.RData')
samdf <- read.csv('/Users/yupeitseng/Documents/RA/FS_PSF/data/fungi/field_soil_property_data.csv', header = T)

#FS
load("/Users/yupeitseng/Documents/RA/FS_PSF/data/fs.elev.RData")
load("/Users/yupeitseng/Documents/RA/FS_PSF/data/fs.full.RData")

# theme_set(theme_bw())
samples.out <- as.data.frame(rownames(seqtab.nochim))

samples.out$sample_ID <- sub('-ITS_R1.fastq.gz.*', '', rownames(seqtab.nochim))
# samples.out <- samples.out[, -1]

samdf <- left_join(samples.out, samdf, by = "sample_ID")
rownames(samdf) <- samdf[, 2]
samdf <- samdf[, -1]

rownames(seqtab.nochim) <- samples.out$sample_ID

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps <- prune_samples(!(sample_names(ps)  %in% c("M24")), ps)
ps <- prune_samples(!(sample_names(ps)  %in% c("NC1", "NC2", "NC3", "NC4", "NC5", "NC6", "NC7")), ps) # Remove NC sample



ps <- subset_taxa(ps, Kingdom != "Unassigned")
ps <- subset_taxa(ps, Class != "D_2__Chloroplast")
ps <- subset_taxa(ps, Family != "D_4__Mitochondria")
ps

#draw ps rarefaction Curve (raremax = 14231) w/o NA----
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
pdf(file = "/Users/yupeitseng/Documents/RA/FS_PSF/data/plot/relative_aundance.pdf",   # The directory you want to save the file in
    width = 24, # The width of the plot in inches
    height = 18) # The height of the plot in inches
#Phylum
Stability_phylum<-rarefied_ps %>%
  tax_glom(taxrank = "Phylum") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
Stability_phylum[Stability_phylum$Abundance<0.01,]$Phylum="Others"

colourCount = length(unique(Stability_phylum$Phylum))
colors_phylum_n <- randomColor(count = colourCount)
colors_phylum_n <- distinctColorPalette(colourCount)

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
colors_class_n <- randomColor(count = colourCount)
colors_class_n <- distinctColorPalette(colourCount)

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
colors_order_n <- randomColor(count = colourCount)
colors_order_n <- distinctColorPalette(colourCount)

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
colors_family_n <- randomColor(count = colourCount)
colors_family_n <- distinctColorPalette(colourCount)

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
colors_genus_n <- randomColor(count = colourCount)
colors_genus_n <- distinctColorPalette(colourCount)

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
colors_species_n <- randomColor(count = colourCount)
colors_species_n <- distinctColorPalette(colourCount)

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
plot_richness(rarefied_ps, x = "sp", measures = c("Observed", "Shannon", "Simpson"), color = "sp")+
  theme(text = element_text(size = 30))+
  labs(color = 'Species')

sp <- c()
decay_e <- c()
decay_m <- c()
for (i in colnames(ps@tax_table)){
  phy <- rarefied_ps %>% aggregate_taxa(level = i) 
  alpha <- estimate_richness(phy)
  alpha <- cbind(alpha, as.data.frame(phy@sam_data$sp), as.data.frame(phy@sam_data$decay_cat))
  colnames(alpha)[10:11] <- c("sp", "decay_cat")
  
  # Pairwise comparisons
  sp <- c(sp,
         unlist(pairwise_t_test(alpha, Observed ~ sp, p.adjust.method = "bonferroni")[c(8:9)]),
          unlist(pairwise_t_test(alpha, Shannon ~ sp, p.adjust.method = "bonferroni")[c(8:9)]),
          unlist(pairwise_t_test(alpha, Simpson ~ sp, p.adjust.method = "bonferroni")[c(8:9)])
  )
  decay_e <- c(decay_e,
             unlist(pairwise_t_test(alpha[which(alpha$sp == "ENGERO"), ], Observed ~ decay_cat, p.adjust.method = "bonferroni")[, c(8:9)]),
             unlist(pairwise_t_test(alpha[which(alpha$sp == "ENGERO"), ], Shannon ~ decay_cat, p.adjust.method = "bonferroni")[, c(8:9)]),
             unlist(pairwise_t_test(alpha[which(alpha$sp == "ENGERO"), ], Simpson ~ decay_cat, p.adjust.method = "bonferroni")[, c(8:9)])
             )
  decay_m <- c(decay_m,
               unlist(pairwise_t_test(alpha[which(alpha$sp == "MACHZU"), ], Observed ~ decay_cat, p.adjust.method = "bonferroni")[, c(8:9)]),
               unlist(pairwise_t_test(alpha[which(alpha$sp == "MACHZU"), ], Shannon ~ decay_cat, p.adjust.method = "bonferroni")[, c(8:9)]),
               unlist(pairwise_t_test(alpha[which(alpha$sp == "MACHZU"), ], Simpson ~ decay_cat, p.adjust.method = "bonferroni")[, c(8:9)])
  )
}

#Only "Phylum" level, "different species" have significant different "observed" diversity
#Draw the plot to see the pattern
#The difference seems to be 
phy <- rarefied_ps %>% aggregate_taxa(level = 'Phylum')

plot_richness(phy, x = "sp", measures = c("Observed", "Shannon", "Simpson"), color = "sp")+
  theme(text = element_text(size = 30))+
  labs(color = 'Kingdom')



#dbem----
#Species data transformation
ps.prop <- transform_sample_counts(rarefied_ps, function(otu) otu/sum(otu))
spe.log <- log1p(as.data.frame(ps.prop@otu_table))  # species data are in percentage scale which is strongly rightskewed, better to transform them
spe.hell <- decostand (spe.log, 'hell')  # we are planning to do tb-RDA, this is Hellinger pre-transformation
colnames(spe.hell) <- 1:ncol(spe.hell)

#coordination data
xy <- as.matrix(ps.prop@sam_data[, c("gx_new", "gy_new")])


#environment data
env <- as.data.frame(as.matrix(ps.prop@sam_data[, c("sp", "decay_cat", "dbh", "pH")]))
env$dbh <- as.numeric(env$dbh)
env$pH <- as.numeric(env$pH)

decay_cat <- model.matrix(~-1 + env[ ,2])
sp <- model.matrix( ~ env[ ,1])[ ,-1]
env <- cbind(env, sp, decay_cat)
colnames(env)[5:8] <- c("MACHZU", "2009_2014","2014_2019", "Alive")

env_plot <- cbind(env, xy)

env <- env %>% mutate(decay_cat = ifelse(decay_cat == "Alive", "Alive", ifelse(decay_cat == "2014_2019", "Short decay", "Long dacay")))
#boxplot of pH and decay category
ggplot(env, aes(x=decay_cat, y=pH, fill=sp)) + 
  scale_fill_manual(values=c("#A7C084", "#B6A87C"), "Species")+
  labs(x = "Decay category")+
  geom_boxplot()+
  theme(text = element_text(size = 20)) 


env_plot$sp_decay <- paste(env_plot$sp, env_plot$decay_cat)
env_plot$decay_cat_name <- lapply(env_plot$decay_cat, FUN = function(x){
  ifelse (x == "2014_2019", "Short decay", ifelse(x == "Alive", "Alive", "Long decay"))
})
ggplot()+
  geom_contour(fs.elev$col, mapping = aes(x/20, y/20, z = elev), bins = 40, alpha = 0.5,colour = "black")+
  labs(x = "East-West distance (m)", y = "North-South distance (m)", shape = "Decay category")+
  geom_point(data = env_plot, aes(x = gx_new, y = gy_new, color = sp))+
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Species")+
  new_scale_color() +
  geom_point(data = env_plot, aes(x = gx_new, y = gy_new, color = sp_decay, shape = decay_cat, size = 10))+
  geom_text_repel(data = env_plot, aes(x = gx_new, y = gy_new, label = rownames(env_plot)), vjust = 0, nudge_y = 0.5)+
  scale_shape_manual(values=c(15, 16, 17), labels = c("Long decay", "Short decay", "Alive"))+
   scale_color_manual(values=c("#556B37", "#A7C084", "#CAD9B5", "#7D6241", "#B6A87C", "#D3C1AB"))+
  guides(color = FALSE, size = FALSE)+
  theme(text = element_text(size = 20))  

## Step 0.
# Is there a linear trend ?
anova(rda(spe.hell, xy)) # Result: significant trend
# Computation of linearly detrended mite data
spe.hell.det <- resid(lm(as.matrix(spe.hell) ~ ., data = as.data.frame(xy)))

## Step 1. Construct the matrix of dbMEM variables 


dbmem <- as.data.frame(dbmem(ps.prop@sam_data[, c("gx_new", "gy_new")], silent = FALSE)) 

# Truncation distance used above:

thr <- give.thresh(dist(xy))

# Display and count the eigenvalues
attributes(dbmem)$values 
length(attributes(dbmem)$values)

## Step 2. Run the global dbMEM analysis on the detrended Hellinger-transformed mite data
dbmem.rda <- rda(spe.hell.det ~., dbmem) 
anova(dbmem.rda)

## Step 3. Since the R-square is significant, compute the adjusted R2 and run a forward selection of the dbmem variables 
R2a <- RsquareAdj(dbmem.rda)$adj.r.squared 
dbmem.fwd <- forward.sel(spe.hell.det, as.matrix(dbmem), adjR2thresh = R2a)
nb.sig.dbmem <- nrow(dbmem.fwd) # Number of signif. dbMEM
# Identity of the significant dbMEM in increasing order
dbmem.sign <- sort(dbmem.fwd[ ,2])
# Write the significant dbMEM to a new object 
dbmem.red <- dbmem[ ,c(dbmem.sign)]

## Step 4. New dbMEM analysis with 8 significant dbMEM variables ## Adjusted R-square after forward selection: R2adj = 0.02962644
dbmem.rda2 <- rda(spe.hell.det ~ ., data = dbmem.red)
head(summary(dbmem.rda2))
fwd.R2a <- RsquareAdj(dbmem.rda2)$adj.r.squared 
anova(dbmem.rda2)
axes.test <- anova(dbmem.rda2, by = "axis")
# Number of significant axes
nb.ax <- length(which(axes.test[ , ncol(axes.test)] <= 0.05))

## Step 5. Plot the significant canonical axes
rda2.axes <- scores(dbmem.rda2, choices = c(1:nb.ax),
         display = "lc",
         scaling = 1)

par(mfrow = c(1,nb.ax))

for(i in 1:nb.ax){
  sr.value(xy, rda2.axes[ ,i],
           sub = paste("RDA",i), csub = 2)
}

# Interpreting the spatial variation: regression of the significant # canonical axes on the environmental variables, with Shapiro-Wilk
# normality tests of residuals
rda2.axis1.env <- lm(((rda2.axes[ ,1])) ~ ., data = env) 
shapiro.test(resid(rda2.axis1.env)) 
summary(rda2.axis1.env)
rda2.axis2.env <- lm(rda2.axes[ ,2] ~ ., data = env) 
shapiro.test(resid(rda2.axis2.env)) 
summary(rda2.axis2.env)

par(mfrow = c(2, 2))
for(i in 1 : ncol(dbmem.red)){
  sr.value(xy,
           dbmem.red[ ,i],
           sub = paste("dbMEM", dbmem.sign[i]), csub = 2)
}
env_spa <- cbind(env, dbmem.red, xy)
env_spa_nold <- env_spa[, c(4:5, 7:12)]
env_spa_nosd <- env_spa[, c(4:6, 8:12)]
env_spa_noal <- env_spa[, c(4:7, 9:12)]

for (dat in c("env_spa_nold", "env_spa_nosd", "env_spa_noal")){
  rda.env <- rda(spe.hell ~ ., get(dat))
  print(anova(rda.env))
  
  # test_axis <- anova (rda.env, by = 'axis', parallel = 4)
  # #holm's correction (for multiple testing issue)
  # test_axis.adj <- test_axis
  # test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
  # print(test_axis.adj)
  
  test_margin <- anova (rda.env, by = 'margin', parallel = 4)
  test_each <- t(apply (get(dat), 2, FUN = function (x) as.matrix (anova (rda (spe.hell ~ x))[1:4][1,])))
  test_each <- as.data.frame (test_each)
  names (test_each) <- c("Df", "Variance", "F", "Pr(>F)")
  test_each.adj <- test_each
  test_each.adj$`Pr(>F)` <- p.adjust (test_each$`Pr(>F)`, method = 'holm')
  print(test_each.adj[order (test_each.adj$Variance, decreasing = TRUE),])
}

varp <- varpart (spe.hell, env_spa[, 5], env_spa[, 4], env_spa[, 6:7], env_spa[, 9:14])
varp
plot (varp, digits = 2, Xnames = c('species', 'pH', "decay", "spatial"))


adjR2.tbrda <- RsquareAdj (rda.env)$adj.r.squared
sel.fs <- forward.sel (Y = spe.hell, X = env[, 3:7], adjR2thresh = adjR2.tbrda)

pval.adj <- p.adjust (sel.fs$pvalue, method = 'holm')
sel.fs$pval.adj <- pval.adj
sel.fs

#saperate two species
env_plot_e <- env_plot %>% filter(sp == "ENGERO")
spe.hell_e <- spe.hell %>% filter(rownames(spe.hell) %in% rownames(env_plot_e))
rda.e <- rda(spe.hell_e ~ ., env_plot_e[, c(2, 4, 9:10)])
anova(rda.e)
plot(rda.e,
     scaling = 2,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env3 - scaling 1 - lc scores",
     label = T
)

plot(rda.e,
     display = c("sp", "lc", "cn"),
     main = "Triplot RDA spe.hel ~ env3 - scaling 2 - lc scores"
)
spe.sc2 <-
  scores(rda.e, choices = 1:2,
         display = "sp"
  )
arrows(0, 0,
       spe.sc2[, 1] * 0.92,
       spe.sc2[, 2] * 0.92,
       length = 0,
       lty = 1,
       col = "red" )
ordicenter (rda.e, groups = env_plot_e$decay_cat, col = 'orange', cex = 1)

for (i in seq (1, 3)){
  ordispider(rda.e, groups = env_plot_e$decay_cat, show.groups = i, col = i, label = T)
}
ordiplot(rda.e)
points(rda.e, display = "sites", cex = 1, pch = 16, col = "red")
text(rda.e, display = "species", cex = 1, col = "blue")

test_margin <- anova (rda.e, by = 'margin', parallel = 4)
test_each <- t(apply (env_plot_e[, c(4, 6:10)], 2, FUN = function (x) as.matrix (anova (rda (spe.hell_e ~ x))[1:4][1,])))
test_each <- as.data.frame (test_each)
names (test_each) <- c("Df", "Variance", "F", "Pr(>F)")
test_each.adj <- test_each
test_each.adj$`Pr(>F)` <- p.adjust (test_each$`Pr(>F)`, method = 'holm')

ordisurf(rda.e, ps.prop_m@sam_data$pH, add = TRUE)

env_plot_m <- env_plot %>% filter(sp == "MACHZU")
spe.hell_m <- spe.hell %>% filter(rownames(spe.hell) %in% rownames(env_plot_m))



#03 Transform data to proportions as appropriate for Bray-Curtis distances----
ps.prop <- transform_sample_counts(rarefied_ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

plot_ordination(ps.prop, ord.nmds.bray, color = "sp") + 
  geom_point(size = 6)+
  stat_ellipse(type = "norm", linetype = 2) + 
  stat_ellipse(type = "t") + 
  theme_bw() + 
  scale_color_manual(values=c("#A7C084", "#B6A87C"), "Species")+
  theme(text = element_text(size = 40)) 

bray.dist<-vegdist(ps.prop@otu_table, method='bray')
adonis2(bray.dist~as.factor(ps.prop@sam_data$sp), permutations=9999)
permutest(betadisper(bray.dist, group=as.factor(ps.prop@sam_data$sp)))
anosim(ps.prop@otu_table, as.factor(ps.prop@sam_data$sp), distance = "bray", permutations = 9999)

#E
ps.prop_e <- prune_samples(sample_names(ps.prop) %in% rownames(ps.prop@sam_data)[which(ps.prop@sam_data$sp == 'ENGERO')], ps.prop) 
ord.nmds.bray_e <- ordinate(ps.prop_e, method="NMDS", distance="bray")

# scale_color_manual(values=c("#556B37", "#A7C084", "#CAD9B5", "#7D6241", "#B6A87C", "#D3C1AB"))+
  
plot_ordination(ps.prop_e, ord.nmds.bray_e, color = "decay_cat") + 
  geom_point(size = 6)+
  stat_ellipse(type = "norm", linetype = 2) + 
  stat_ellipse(type = "t") + 
  theme_bw() + 
  scale_color_manual(values=c("#556B37", "#A7C084", "#CAD9B5"), "Decay category", labels = c("Long decay", "Short decay", "Alive"))+
  theme(text = element_text(size = 40)) 

bray.dist_e<-vegdist(ps.prop_e@otu_table, method='bray')
adonis2(bray.dist_e ~  as.factor(ps.prop_e@sam_data$decay_cat), permutations=9999)
permutest(betadisper(bray.dist_e, group=ps.prop_e@sam_data$decay_cat))
pairwise.adonis(bray.dist_e, as.factor(ps.prop_e@sam_data$decay_cat))

perm_e <- adonis2(bray.dist_e ~  as.factor(ps.prop_e@sam_data$decay_cat) + ps.prop_e@sam_data$pH + ps.prop_e@sam_data$gx_new +ps.prop_e@sam_data$gy_new, permutations=9999)
write.csv(perm_e, "/Users/yupeitseng/Documents/RA/FS_PSF/data/perm_e_fungi.csv")

#MACHZU
ps.prop_m <- prune_samples(sample_names(ps.prop) %in% rownames(ps.prop@sam_data)[which(ps.prop@sam_data$sp == 'MACHZU')], ps.prop) 
ord.nmds.bray_m <- ordinate(ps.prop_m, method="NMDS", distance="bray")

plot_ordination(ps.prop_m, ord.nmds.bray_m, color = "decay_cat") + 
  geom_point(size = 6)+
  stat_ellipse(type = "norm", linetype = 2) + 
  stat_ellipse(type = "t") + 
  theme_bw() + 
  scale_color_manual(values=c("#7D6241", "#B6A87C", "#D3C1AB"), "Decay category", labels = c("Long decay", "Short decay", "Alive"))+
  theme(text = element_text(size = 40)) 

bray.dist_m<-vegdist(ps.prop_m@otu_table, method='bray')
adonis2(bray.dist_m ~ as.factor(ps.prop_m@sam_data$decay_cat), permutations=9999)
permutest(betadisper(bray.dist_m, group=ps.prop_m@sam_data$zone))
anosim(ps.prop_m@otu_table, as.factor(ps.prop_m@sam_data$decay_cat), permutations=9999)

perm_m <- adonis2(bray.dist_m ~  as.factor(ps.prop_m@sam_data$decay_cat) + ps.prop_m@sam_data$pH + ps.prop_m@sam_data$gx_new + ps.prop_m@sam_data$gy_new, permutations=9999)
write.csv(perm_m, "/Users/yupeitseng/Documents/RA/FS_PSF/data/perm_m_fungi.csv")

plot(ord.nmds.bray_m, type = "n")
points(ord.nmds.bray_m, display = "sites", cex = 1, pch = 16, col = c("blue", "orange", "black") [as.factor(ps.prop_m@sam_data$decay_cat)])
ordiellipse(ord.nmds.bray_m, as.factor(ps.prop_m@sam_data$decay_cat), display="sites", 
            kind = "se", conf=0.95, label=T, 
            lwd=2, col=c("blue","orange", "black"))

text(ord.nmds.bray_m, display = "species", cex = 1, col = "blue")
ordisurf(ord.nmds.bray_m, ps.prop_m@sam_data$pH, add = TRUE)

plot(ord.nmds.bray, type="n", las=1)
points(ord.nmds.bray, display = "sites", cex = 1, pch = 16, col = "red")
ordiellipse(ord.nmds.bray, as.factor(ps.prop_m@sam_data$decay_cat), display="sites", 
            kind = "se", conf=0.95, label=T, 
            lwd=2, col=c("blue","orange", "black"))






spe <- ps.prop_e@otu_table  # rename variables to make them shorter
env <- as.data.frame(as.matrix(ps.prop_e@sam_data))  # select only two explanatory variables

spe.log <- log1p (spe)
spe.log.hell <- decostand (spe.log, 'hell')
tbRDA <- rda (spe.log.hell ~ decay_cat + pH, data = env)
head (summary (tbRDA))  # prints first lines of the summary of tbRDA
test_margin <- anova (tbRDA, by = 'margin', parallel = 4)

test_each <- t(apply (chem, 2, FUN = function (x) as.matrix (anova (rda (vasc.hell ~ x))[1:4][1,])))
test_each <- as.data.frame (test_each)
names (test_each) <- c("Df", "Variance", "F", "Pr(>F)")
test_each.adj <- test_each
test_each.adj$`Pr(>F)` <- p.adjust (test_each$`Pr(>F)`, method = 'holm')
test_each.adj[order (test_each.adj$Variance, decreasing = TRUE),]

#ENGERO
# Transform the data
mite.h <- decostand (ps.prop_e@otu_table, "hellinger")
mite.xy <- ps.prop_e@sam_data[, c('gx_new', 'gy_new')]
mite.xy.c <- scale(mite.xy, center = TRUE, scale = FALSE)
## Univariate spatial correlogram (based on Moran's I)
# Search for neighbours of all points within a radius of 0.7 m # and multiples (i.e., 0 to 0.7 m, 0.7 to 1.4 m and so on). plot.links(mite.xy, thresh = 0.7)
nb1 <- dnearneigh(as.matrix(mite.xy), 0, 0.2)
summary(nb1)
# Correlogram of substrate density
geo_dat <- ps.prop@sam_data[, c('gx_new', 'gy_new')]
ph_dis <- vegdist(ps.prop@sam_data$pH, method="euclidean", diag=T, upper=T) 
spe_dis <- vegdist(ps.prop@otu_table, method="bray", diag=T, upper=T) 
geo_dis <- vegdist(geo_dat, method="euclidean", diag=T, upper=T) 
mantel.rtest(spe_dis, ph_dis, nrepet = 9999)


#rda all
spe <- rarefied_ps@otu_table
env <- samdf # select only two explanatory variables
env <- env[-which(env$sample_ID %in% c("NC1", "NC2", "NC3", "NC4", "NC5", "NC6", "NC7", "M24")), ]

spe.log <- log1p (spe)
spe.log.hell <- decostand (spe.log, 'hell')
tbRDA <- rda (spe.log.hell ~ sp + decay_cat + pH + zone, data = env)
head (summary (tbRDA))  # prints first lines of the summary of tbRDA
ordiplot(tbRDA, type = 'p')
test_axis <- anova(tbRDA, by = 'margin')

test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

#rda ENGERO
env <- samdf # select only two explanatory variables
env <- env[-which(env$sample_ID %in% c("NC1", "NC2", "NC3", "NC4", "NC5", "NC6", "NC7", "M24")), ]
env <- env[which(env$sp == "ENGERO"), ]
spe <- rarefied_ps@otu_table[which(rownames(rarefied_ps@otu_table) %in% env$sample_ID), ]
rownames(spe)

spe.log <- log1p (spe)
spe.log.hell <- decostand (spe.log, 'hell')
tbRDA <- rda (spe.log.hell ~ decay_cat + pH + zone, data = env)
head (summary (tbRDA))  # prints first lines of the summary of tbRDA
ordiplot(tbRDA, type = 'p')
test_axis <- anova(tbRDA, by = 'margin')

test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj

# boxplot(pH ~ zone, samdf)

#rda MACHZU
env <- samdf # select only two explanatory variables
env <- env[-which(env$sample_ID %in% c("NC1", "NC2", "NC3", "NC4", "NC5", "NC6", "NC7", "M24")), ]
env <- env[which(env$sp == "MACHZU"), ]
spe <- rarefied_ps@otu_table[which(rownames(rarefied_ps@otu_table) %in% env$sample_ID), ]
rownames(spe)

spe.log <- log1p (spe)
spe.log.hell <- decostand (spe.log, 'hell')
tbRDA <- rda (spe.log.hell ~ decay_cat + pH, data = env)
head (summary (tbRDA))  # prints first lines of the summary of tbRDA
ordiplot(tbRDA, type = 'p')
test_axis <- anova(tbRDA, by = 'margin')

test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj
