library(PerformanceAnalytics)
library(betapart)
library(vegan)
library(igraph)
library(devtools)
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
library (adespatial)
library(car)
library(dplyr)
library(phyloseq)
library(cooccur)
library(igraph)
library(visNetwork)
library(randomcoloR)
#for funguild database
# devtools::install_github("brendanf/FUNGuildR")
library(FUNGuildR)
library(magrittr)
library(ggplot2)
library(ggpubr)
#function----
avPlots2 <- function(model, terms=~., intercept=FALSE, layout=NULL, ask, 
                     main, xlab, ...){
  terms <- if(is.character(terms)) paste("~",terms) else terms
  vform <- update(formula(model),terms)
  if(any(is.na(match(all.vars(vform), all.vars(formula(model))))))
    stop("Only predictors in the formula can be plotted.")
  terms.model <- attr(attr(model.frame(model), "terms"), "term.labels")
  terms.vform <- attr(terms(vform), "term.labels")
  terms.used <- match(terms.vform, terms.model)
  mm <- model.matrix(model) 
  model.names <- attributes(mm)$dimnames[[2]]
  model.assign <- attributes(mm)$assign
  good <- model.names[!is.na(match(model.assign, terms.used))]
  if (intercept) good <- c("(Intercept)", good)
  nt <- length(good)
  if (nt == 0) stop("No plots specified")
  if (missing(main)) main <- if (nt == 1) paste("Added-Variable Plot:", good) else "Added-Variable Plots"
  if (nt == 0) stop("No plots specified")
  if (nt > 1 & (is.null(layout) || is.numeric(layout))) {
    if(is.null(layout)){
      layout <- switch(min(nt, 9), c(1, 1), c(1, 2), c(2, 2), c(2, 2), 
                       c(3, 2), c(3, 2), c(3, 3), c(3, 3), c(3, 3))
    }
    ask <- if(missing(ask) || is.null(ask)) prod(layout)<nt else ask
    op <- par(mfrow=layout, ask=ask, no.readonly=TRUE, 
              oma=c(0, 0, 1.5, 0), mar=c(5, 4, 1, 2) + .1)
    on.exit(par(op))
  }
  if (missing(xlab)) xlab <- paste(good, "| others")
  if (length(xlab) == 1L) xlab <- rep(xlab, length(good))
  if (length(xlab) > length(good))
    warning("'xlab' not length 1 or the number of model names, truncating")
  res <- as.list(NULL)
  for (i in seq_along(good)) {
    term <- good[[i]]
    res[[term]] <- avPlot(model, term, main="", xlab=xlab[[i]], ...)
  }
  mtext(side=3,outer=TRUE,main, cex=1.2)
  invisible(res)
}


#import data (feature table, metadata, taxanomy)----
load('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/seqtab.nochim.RData')
load('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/taxa .RData')
seqtab.nochim <- as.data.frame(seqtab.nochim)
seqtab.nochim$sample <- sub('-.*', '', rownames(seqtab.nochim))
seqtab.nochim <- seqtab.nochim %>% group_by(sample) %>% summarise_all(sum) %>% as.data.frame(seqtab.nochim)
rownames(seqtab.nochim) <- seqtab.nochim[, 1]
seqtab.nochim <- seqtab.nochim[, -1]

samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi.csv', header = T)
samdf <- samdf[, -1]
samdf$fern_ID <- as.character(samdf$fern_ID)
# theme_set(theme_bw())
samples.out <- as.data.frame(rownames(seqtab.nochim))
colnames(samples.out)[1] <- 'fern_ID'


samdf <- left_join(samples.out, samdf, by = "fern_ID")
rownames(samdf) <- rownames(seqtab.nochim)
samdf$total_weight <- samdf[, 12] + samdf[, 13]

trait <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_trait_clean.csv', header = T)
trait$LA <- 16*(pi*(6.3/2)^2) # mm
trait$SLA <- (trait$Wt.dry*1000)/trait$LA
trait_avg <- trait %>% group_by(fern.ID) %>% summarize(Chl_avg = mean(SPAD_avg),
                                                         LDMC_avg = mean(LDMC),
                                                         Lth_avg = mean(Lth_avg),
                                                         SLA_avg = mean(SLA))
colnames(trait_avg)[1] <- 'fern_ID'
trait_avg$fern_ID <- as.character(trait_avg$fern_ID)
samdf <- left_join(samdf, trait_avg)
rownames(samdf) <- rownames(seqtab.nochim)

samdf$water_content <- (samdf$nest.weight-samdf$dry_mass)/samdf$nest.weight*100


layer_info <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv', header = T)
layer_info <- layer_info %>% mutate(mass_ratio = dry_mass.x/dry_mass.y) %>% mutate(C_N_ratio = C_N*mass_ratio, pH_ratio = pH*mass_ratio) %>% group_by(fern_ID) %>% summarise(C_N_ratio = sum(C_N_ratio, na.rm = T), pH_ratio = sum(pH_ratio, na.rm = T))
layer_info$fern_ID <- as.character(layer_info$fern_ID)

samdf <- left_join(samdf, layer_info, by = "fern_ID")

# write.csv(samdf, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi_trait.csv')

distance <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_distance.csv')
colnames(distance) <- c("tree.ID", "x", "y")
samdf <- left_join(samdf, distance)
rownames(samdf) <- rownames(seqtab.nochim)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

ps <- prune_samples(!(sample_names(ps)  %in% c("NC1", "NC2", "NC3", "NC4")), ps) # Remove NC sample

ps <- subset_taxa(ps, Kingdom != "Unassigned")
ps <- subset_taxa(ps, Class != "D_2__Chloroplast")
ps <- subset_taxa(ps, Family != "D_4__Mitochondria")
ps

#draw ps rarefaction Curve (raremax = 14231) w/o NA
col <- randomColor(count = 24)
raremax = min(sample_sums(ps)) #get the minimum point

tab <- otu_table(ps)
class(tab) <- "matrix"


rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=T, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_ps = rarefy_even_depth(ps, rngseed = 1, sample.size = raremax, replace = F)

##02 SAR----
#alpha diversity
plot_richness(rarefied_ps, x = "nest.weight", measures = c("Observed", "Shannon", "Simpson"), color = "fern.size")

#species richness
rar_matrix <- rarefied_ps@otu_table

samdf <- samdf[-(25:28), ]
samdf$richness <- apply(rar_matrix>0,1,sum)


p <- ggscatter(log10(samdf[, c(14, 29)]), x = 'dry_mass', y = 'richness',
          add = "reg.line", conf.int = TRUE, color = "#609067", xlab = 'log(Dry mass)', ylab = 'log(Operational taxonomic unit \n (OTU) richness)')+
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

ggsave('/Users/yupeitseng/Downloads/SAR.png', p, dpi = 300, height = 12, width = 16)

chart.Correlation(scale(samdf[which(samdf$fern.size == 'L' | samdf$fern.size == 'M'), 4:19]), histogram=TRUE, pch=19)
cor_dat <- samdf[, c(12, 14, 17, 19, 5, 16, 15, 20:23, 24)]
colnames(cor_dat) <- c('Nest fresh mass', 'Nest dry mass', 'Nest volume', 'Total weight', 'Max leaves width', 'Mean leaf length', 'pH', 'Mean chlorophyll\n concentration', 'Mean LDMC', 'Mean leaf\n thickness', 'Mean SLA', 'OTU richness')
chart.Correlation(scale(cor_dat), histogram=TRUE, pch=19, cex.labels = 1.4, label.pos = 0.8)

cor_dat <- samdf[, c(12, 14, 17, 19, 5, 16, 15, 4, 18, 20:23, 24)]
colnames(cor_dat)[1:10] <- c('Nest fresh mass', 'Nest dry mass', 'Nest volume', 'Total weight', 'Max leaves width', 'Mean leaf length', 'pH', 'Height', 'DBH', 'OTU')

cor_dat_size <- samdf[, c(12, 14, 17, 19, 5, 16, 24:31)]
colnames(cor_dat_size)[1:7] <- c('Nest fresh mass', 'Nest dry mass', 'Nest volume', 'Total weight', 'Max leaves width', 'Mean leaf length', 'OTU')
cor_dat_long <- cor_dat_size %>% gather(key="morphological_variables", value="value", c('Nest fresh mass', 'Nest dry mass', 'Nest volume', 'Total weight', 'Max leaves width', 'Mean leaf length')) %>% gather(key="taxa", value="richness", c(colnames(cor_dat_size)[7:14]))
cor_dat_long$value <- (log10(cor_dat_long$value))
cor_dat_long$richness <- (log10(cor_dat_long$richness)) 
cor_dat_long <- cor_dat_long[-which(cor_dat_long$taxa == "Kingdom"), ]
# , levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"))
sp <- ggscatter(cor_dat_long, x = 'value', y = 'richness',
                add = "reg.line", conf.int = TRUE, color = "#609067")+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =5)+
  facet_grid(factor(taxa, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "OTU"))~morphological_variables, scales="free")+
  labs(x = "log(Morphlogical variables)", y = "log(Richness)")+
  theme(text = element_text(size = 16))

# factor(morphological_variables, levels=c('Nest fresh mass', 'Nest dry mass', 'Nest volume', 'Total weight', 'Max leaves width', 'Mean leaf length'))
ggplot(econdatalong, aes(x=log10(value), y=log10(richness)))+
  geom_point()+
  facet_grid(measure~Region, scales="free", space="free_x",  labeller= variable_labeller2)

chart.Correlation((cor_dat), histogram=TRUE, cex.labels = 1.5)
plot(lm(scale(richness) ~ scale(pH) + scale(dry_mass), data = samdf))

rank_names(rarefied_ps)
rarefied_ps_sp <- tax_glom(rarefied_ps, taxrank = 'Species')
rarefied_ps_gn <- tax_glom(rarefied_ps, taxrank = 'Genus')
rarefied_ps_fm <- tax_glom(rarefied_ps, taxrank = 'Family')

richness <- c()
for (i in rank_names(rarefied_ps)){
  rarefied_ps_taxa <- tax_glom(rarefied_ps, taxrank = i)
  richness <- c(richness, apply(rarefied_ps_taxa@otu_table>0,1,sum))
}
richness_mat <- matrix(richness, ncol = 7, nrow = 24, byrow = F)
colnames(richness_mat) <- rank_names(rarefied_ps)
samdf <- cbind(samdf, richness_mat)

lm <- lm(scale(log10(richness)) ~ scale(log10(dry_mass)) + scale(pH) +  scale(Chl_avg), data = as.data.frame(((samdf[, -(1:3)]))))
summary(lm)
performance::check_model(lm)
performance::binned_residuals(lm)
par(bg = NA)
avPlots2(lm, layout = c(1, 3), cex = 3, pch = 16, cex.lab = 2, main = '', cex.axis = 1.8, col.line = '#609067', col = '#609067', ylab = c('log(OTU richness) | others', '', ''), xlab = c('log(Dry mass) | others', 'pH | others', 'Leaf chlorophyll concentration | others'))

# write.csv(samdf, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi_clean.csv')
#03 Transform data to proportions as appropriate for Bray-Curtis distances
#nMDS bry curtis distance
refseq(rarefied_ps)
ps.prop <- transform_sample_counts(rarefied_ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

ordplot <- plot_ordination(ps.prop, ord.nmds.bray, color = "fern.size")
ordplot + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  geom_point(cex = 5)+
  labs(colour = "Fern size") +
  theme_classic()+
  theme(text = element_text(size = 24))+
  scale_colour_manual(values = c("#ADC178", "#F1DCA7", "#997B66"), label = c("Large", "Medium", "Small"))
    

summary(ord.nmds.bray)
bray.dist<-vegdist(ps.prop@otu_table, method='bray')
#PERMANOVA
adonis2(bray.dist~ as.factor(ps.prop_fil@sam_data$fern.size) + ps.prop_fil@sam_data$pH_ratio + ps.prop_fil@sam_data$C_N_ratio, permutations=9999)
disper <-betadisper(bray.dist, group=ps.prop@sam_data$fern.size, type = 'centroid')
anova(disper)
permutest(disper)
plot(disper)

adonis2(dist(disper$distances) ~ as.factor(ps.prop@sam_data$fern.size))
pairwise_size <- pairwise.adonis(bray.dist, as.factor(ps.prop@sam_data$fern.size))
write.csv(pairwise_size, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/pairwise_size.csv')
#rda
spe <- rarefied_ps@otu_table  # rename variables to make them shorter


ggscatter(samdf, x = 'dry_mass', y = 'pH',
          add = "reg.line", conf.int = TRUE, color = "#609067", xlab = 'Dry mass')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =5)+
  theme(text = element_text(size = 16))
ggscatter(samdf, x = 'dry_mass', y = 'Chl_avg',
          add = "reg.line", conf.int = TRUE, color = "#609067", xlab = 'Dry mass', ylab = 'Mean chlorophyll concentration')+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =5)+
  theme(text = element_text(size = 16))

env <- rarefied_ps@sam_data
spe <-  rarefied_ps@otu_table
spe.log <- log1p (spe)
spe.log.hell <- decostand (spe.log, 'hell')
tbRDA <- rda (spe.log.hell ~ dry_mass + water_content, data = samdf[-(25:28),])
head (summary (tbRDA))  # prints first lines of the summary of tbRDA
ordiplot(tbRDA, type = 'p')
par(mfrow = c(1, 1))
anova(tbRDA)
adjR2.tbrda <- RsquareAdj (tbRDA)$adj.r.squared

sel.fs <- forward.sel (Y = spe.log.hell, X = env[,c('Chl_avg', 'dry_mass', 'pH')], adjR2thresh = adjR2.tbrda)

test_axis <- anova (tbRDA, by = 'margin')
test_axis.adj <- test_axis
test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
test_axis.adj



#rda for margin
env <- samdf[1:24, c(14, 15, 20)]  # select only two explanatory variables
test_each <- t(apply (env, 2, FUN = function (x) as.matrix (anova (rda (spe.log.hell ~ x))[1:4][1,])))
test_each <- as.data.frame (test_each)
names (test_each) <- c("Df", "Variance", "F", "Pr(>F)")
test_each.adj <- test_each
test_each.adj$`Pr(>F)` <- p.adjust (test_each$`Pr(>F)`, method = 'holm')
test_each.adj <- test_each.adj[order (test_each.adj$Variance, decreasing = TRUE),]

write.csv(test_axis.adj, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/rda_margin.csv')

boxplot(pH ~ fern.size, data = samdf)
plot(pH ~ dry_mass, data = samdf)
#all data pair beta and partition into nestedness and turnover----
#transform abundance data to presence absence data
rar_matrix_pa <- rar_matrix
rar_matrix_pa[rar_matrix_pa > 0] <-1
rar_matrix_pa <- as.data.frame(rar_matrix_pa)

#calculate difference of size 
size.dist <- as.vector(dist(scale(samdf$nest.weight), method = "euclidean", diag = FALSE, upper = FALSE))
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
rar_matrix_pa_sm <- rar_matrix_pa[-which(samdf$fern.size == 'L'), ]

#calculate difference of size
size.dist_sm <- dist(scale(samdf$nest.weight[-which(samdf$fern.size == 'L')]), method = "euclidean", diag = FALSE, upper = FALSE)
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

size.dist_sm$sample <- rep('S_M', nrow(size.dist_sm))
size.dist_sm_sub <- size.dist_sm[, c(1, 3, 6, 8)]
colnames(size.dist_sm_sub) <- c('Size difference', 'Sorensen', 'Jaccard', 'Sample')
size.dist_sm_sub <- size.dist_sm_sub %>% gather('index', 'nestedness', c(Sorensen, Jaccard))

size.dist$sample <- rep('S_M_L', nrow(size.dist))
size.dist_sub <- size.dist[, c(1, 3, 6, 8)]
colnames(size.dist_sub) <- c('Size difference', 'Sorensen', 'Jaccard', 'Sample')
size.dist_sub <- size.dist_sub %>% gather('index', 'nestedness', c(Sorensen, Jaccard))

combine <- rbind(size.dist_sub, size.dist_sm_sub)
ggscatter(combine, x = 'Size difference', y = 'nestedness',
          add = "reg.line", conf.int = TRUE, color = "#609067")+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =10)+
  facet_grid(Sample~index, scales="free")+
  labs(x = "Size difference", y = "Nestedness-resultant dissimilarity")+
  theme(text = element_text(size = 30))

#mantel test
size.dist_sm_mantel <- dist(scale(samdf$nest.weight[-which(samdf$fern.size == 'L')]), method = "euclidean", diag = T, upper = T)
sor_dist_sm_mantel <- graph.data.frame(as.vector(sor_dist_sm$beta.sne), directed=FALSE)

mantel(sor_dist_sm_mantel$beta.sne, size.dist_sm_mantel, method = "spearman", permutations = 9999, na.rm = TRUE)

#beta multi for different size of BNF
#calculate pair jaccard and sorensen
sor_dist_s <-beta.multi(rar_matrix_pa[which(samdf$fern.size == 'S'), ], index.family="sor")
sor_dist_m <-beta.multi(rar_matrix_pa[which(samdf$fern.size == 'M'), ], index.family="sor")
sor_dist_l <-beta.multi(rar_matrix_pa[which(samdf$fern.size == 'L'), ], index.family="sor")
sor <- matrix(c(sor_dist_s, sor_dist_m, sor_dist_l), ncol = 3, nrow = 3, byrow = F)

jac_dist_s <-beta.multi(rar_matrix_pa[which(samdf$fern.size == 'S'), ], index.family="jac")


barplot(sor[-3, ],
        main = "Beta sor",
        xlab = "Size",
        names.arg = c("S", "M", "L"),
        col = c("grey","white")
)

legend("topleft",
       c("beta nestedness","beta turnover"),
       fill = c("grey","white")
)

###----Distance----
rank_names(rarefied_ps)
rarefied_ps_sp <- tax_glom(rarefied_ps, taxrank = 'Species')
rarefied_ps_gn <- tax_glom(rarefied_ps, taxrank = 'Genus')
rarefied_ps_fm <- tax_glom(rarefied_ps, taxrank = 'Family')
rarefied_ps_or <- tax_glom(rarefied_ps, taxrank = 'Order')
rarefied_ps_cl <- tax_glom(rarefied_ps, taxrank = 'Class')
rarefied_ps_py <- tax_glom(rarefied_ps, taxrank = 'Phylum')

dist.abund <- vegdist(prune_samples(rarefied_ps@sam_data$fern.size == "S" | rarefied_ps@sam_data$fern.size == "M", rarefied_ps)@otu_table, method = "bray")
dist.space <- dist(prune_samples(rarefied_ps@sam_data$fern.size == "S" | rarefied_ps@sam_data$fern.size == "M", rarefied_ps)@sam_data[, c(4, 25:26)], method = "euclidean")

mantel(dist.abund, dist.space, method = "spearman", permutations = 9999, na.rm = TRUE)

dist.abund <- vegdist(rarefied_ps_gn@otu_table, method = "bray")
dist.space <- dist(rarefied_ps_gn@sam_data[,c(4, 25:26)], method = "euclidean")

mantel(dist.abund, dist.space, method = "spearman", permutations = 9999, na.rm = TRUE)
vec.dist.abun <- as.vector(dist.abund)
vec.dist.space <- as.vector(dist.space)

plot(scale(vec.dist.abun) ~ (scale(vec.dist.space)))
hist(resid(lm((vec.dist.abun) ~ scale(vec.dist.space))))
hist(scale(vec.dist.abun))
hist((vec.dist.space))
summary(lm(scale(vec.dist.abun) ~ I(scale(vec.dist.space)^2)))

phy <- rarefied_ps_sp@otu_table %>% as.data.frame()

phy[phy>0] <-1
colnames(phy) <- rarefied_ps_sp@tax_table[which(colnames(phy) %in% rownames(rarefied_ps_sp@tax_table)), 6]

which(colSums(phy) == max(colSums(phy)))
py <- unlist(as.vector(rarefied_ps_sp@tax_table[which(colnames(rarefied_ps_sp@otu_table %>% as.data.frame()) %in% rownames(rarefied_ps_sp@tax_table)), 2]))
cooc_sum <-cooccur(t(phy) , spp_names = T)
cooc_sum <- print(cooc_sum)
plot(cooc_sum)


col <- data.frame(color = rainbow(length(unique(py))), phylum = unique(py))
py_mat <- py %>% as.data.frame() 
colnames(py_mat) <- "phylum"
py_mat <- left_join(py_mat, col, by = "phylum")
nodes <- data.frame(id = 1:ncol(phy), label = colnames(phy), color = py_mat$color, shadow = TRUE) 

edges <- data.frame(from = cooc_sum$sp1, to = cooc_sum$sp2, color = ifelse(cooc_sum$p_lt <= 0.05, "#B0B2C1", "#3C3F51"), dashes = ifelse(cooc_sum$p_lt <= 0.05, TRUE, FALSE))

lnodes <- data.frame(label = col$phylum, color = col$color)


visNetwork(nodes = nodes, edges = edges) %>% visIgraphLayout(layout = "layout_with_kk") %>% visLegend(addNodes = lnodes)

pie(rep(1, length(unique(py))),                                         # Draw colors in pie chart
    col = unique(py_mat$color),
    main = "grDevices Package [rainbow() Function]")

#Funguild database----
richness <- rarefied_ps %>%
  # tax_glom(taxrank = "Species") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt() %>% 
  filter(Abundance != 0)

richness <- richness %>% group_by(fern_ID) %>% summarise(richness = n()) %>% left_join(samdf)

summary(lm(log10(richness$richness) ~ log10(richness$dry_mass) + (richness$pH_ratio)))
plot(log10(richness$dry_mass) ~ log10(richness$pH_ratio))

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

#checking how many proportion of NA OTU is in the whole dataset
na <- OTU_funguild %>% distinct(OTU, .keep_all = T) %>% group_by(trophicMode) %>% summarise(n = n())
na$n[8] / sum(na$n)

#calculating the average proportion of NA OTU in every BNF individuals
OTU_funguild_group %>% group_by(fern_ID) %>% summarise(sum_richness = sum(richness)) %>% left_join(OTU_funguild_group[is.na(OTU_funguild_group$trophicMode), c("fern_ID", "richness")]) %>% mutate(na_prop = richness/sum_richness) %>% select("na_prop") %>% unlist() %>% mean()

ggplot(OTU_funguild, aes(x = as.factor(fern_ID),y = Abundance, fill = trophicMode))+
  geom_bar(stat="identity")+
  scale_colour_manual(values = c("#ADC178", "#DDE5B6", "#BAA587", "#F1DCA7", "#FFCB69", "#E8AC65", "#B58463", "#997B66"), aesthetics = c("colour", "fill"))+
  facet_grid(~fern.size, scales="free_x")+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  xlab("Fern ID")+
  guides(fill=guide_legend(title="Trophic mode"))+
  theme(text = element_text(size = 18))

#SAR in different trophic modes of fungi
OTU_funguild_group <- OTU_funguild %>% group_by(fern_ID, trophicMode) %>% summarise(richness = n()) %>% left_join(samdf) %>% mutate(richness_log = log10(richness), dry_mass_log = log10(dry_mass), C_N_ratio_log = log10(C_N_ratio)) 
OTU_funguild_group %>% group_by(fern_ID) %>% summarise(richness = sum(richness), dry_mass = mean(dry_mass), C_N_ratio = mean(C_N_ratio), pH_ratio = mean(pH_ratio))

lm_richness <- function(keyword = keyword){
  if (keyword == "All"){
    OTU_funguild_group_fil <- OTU_funguild_group %>% group_by(fern_ID) %>% summarise(richness = sum(richness), dry_mass = mean(dry_mass), C_N_ratio = mean(C_N_ratio), pH_ratio = mean(pH_ratio))
  }else{
    OTU_funguild_group_fil <- OTU_funguild_group %>% filter(trophicMode == keyword) 
  }
  return(summary(lm(log10(richness) ~ log10(dry_mass) + scale(pH_ratio) + scale(C_N_ratio), data = OTU_funguild_group_fil))$coefficients)
}
lm_richness_mat <- lm_richness("All")
for(i in unique(OTU_funguild_group$trophicMode)[-8]){
  lm_richness_mat <- rbind(lm_richness_mat, lm_richness(i))
}

write.csv(lm_richness_mat, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/lm_richness_mat.csv") 

p <- ggscatter(OTU_funguild_group, x = 'C_N_ratio', y = 'richness',
          add = "reg.line", conf.int = TRUE, color = "trophicMode")+
  scale_colour_manual(values = c("#ADC178", "#DDE5B6", "#BAA587", "#F1DCA7", "#FFCB69", "#E8AC65", "#B58463", "#997B66"), aesthetics = c("colour", "fill"))+
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
           label.y.npc = "top", size =5)+
  facet_wrap(~trophicMode, ncol = 4)+
  labs(x = "C:N ratio", y = "Richness")+
  theme(text = element_text(size = 18))
 
ggpar(p, legend.title = "Trophic mode")

#
#function for extrating only the OTU that are in the interesting groups
tophic_group_phy <- function(keyword = keyword){
  OTU_funguild_fil <- unique(OTU_funguild %>% 
                               # filter(trophicMode %in% grep(keyword, unique(OTU_funguild$trophicMode), value = TRUE)) %>%
                               filter(trophicMode == i) %>%
                               select(OTU) %>% unlist %>% as.vector)
  
  allTaxa <- taxa_names(ps.prop)
  allTaxa <- allTaxa[!(allTaxa %in% OTU_funguild_fil)]
  ps.prop_fil <-  prune_taxa(allTaxa, ps.prop)
  return(ps.prop_fil)
}

#drawing combine nmds plot with stress value
#function for drawing individual group of nmds plot
ordi <- function(phy = phy){
  ord.nmds.bray <- ordinate(phy, method="NMDS", distance="bray")
  ordplot <- plot_ordination(phy, ord.nmds.bray, color = "fern.size", aes(size=5))
  return(ordplot +
           stat_ellipse(type = "norm", linetype = 2) +
           stat_ellipse(type = "t") +
           geom_point(cex = 5)+
           labs(colour = "Fern size") +
           scale_colour_manual(values = c("#ADC178", "#F1DCA7", "#997B66"), aesthetics = c("colour", "fill"))+
           theme_bw()+
           theme(text = element_text(size = 18))+
           annotate("text", x = 1, y = 2, label = paste("stress = ", as.character(round(ord.nmds.bray$stress, digits = 2))))+
           coord_cartesian(ylim=c(-1.7,1.7),clip="off"))
}

ggarrange(ordi(ps.prop), ordi(tophic_group_phy("Saprotroph")), ordi(tophic_group_phy("Pathotroph")), ordi(tophic_group_phy("Symbiotroph")),
          labels = c("All", "Saprotroph", "Pathotroph", "Symbiotroph"),
          ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom", vjust = 0.1)+
  theme(plot.margin = margin(0.5,0,0.5,0, "cm")) 

perm_aov <- function(phy = phy){
  ord.nmds.bray <- ordinate(phy, method="NMDS", distance="bray")
  summary(ord.nmds.bray)
  bray.dist<-vegdist(phy@otu_table, method='bray')
  #PERMANOVA
  return(as.data.frame(adonis2(bray.dist~ (phy@sam_data$fern.size) + phy@sam_data$pH_ratio + phy@sam_data$C_N_ratio, permutations=9999)))
}

perm_tab <- rbind(perm_aov(ps.prop), perm_aov(tophic_group_phy("Saprotroph")), perm_aov(tophic_group_phy("Pathotroph")), perm_aov(tophic_group_phy("Symbiotroph")))

# write.csv(perm_tab, "/Users/yupeitseng/Documents/RA/birds_nest_fern/data/perm_tab.csv")

for(i in c("Saprotroph", "Pathotroph", "Symbiotroph")){
  #extract the OTU name that belongs to specific trophic mode
  OTU_funguild_fil <- unique(OTU_funguild %>% 
                               filter(trophicMode %in% grep(i, unique(OTU_funguild$trophicMode), value = TRUE)) %>%
                             # filter(trophicMode == i) %>%
                             select(OTU) %>% unlist %>% as.vector)
  
  allTaxa <- taxa_names(ps.prop)
  allTaxa <- allTaxa[!(allTaxa %in% OTU_funguild_fil)]
  ps.prop_fil <-  prune_taxa(allTaxa, ps.prop)
  return(ps.prop_fil)
  
 
  
  ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

  

  summary(ord.nmds.bray)
  bray.dist<-vegdist(ps.prop@otu_table, method='bray')
  #PERMANOVA
  as.data.frame(adonis2(bray.dist~ (ps.prop@sam_data$dry_mass) + ps.prop@sam_data$pH_ratio + ps.prop@sam_data$C_N_ratio, permutations=9999))
  
 
  # checking guilds barplot
  OTU_funguild_fil_dat <-rarefied_ps %>%
    transform_sample_counts(function(x) {x/sum(x)}) %>%
    psmelt() %>%
    filter(Abundance != 0)
  OTU_funguild_fil_dat <- OTU_funguild_fil_dat %>% left_join(funguild, by = "OTU") %>%
    filter(trophicMode == "Pathotroph")
    # filter(trophicMode %in% grep("Symbiotroph", unique(OTU_funguild$trophicMode), value = TRUE))

  ggplot(OTU_funguild_fil_dat, aes(x = as.factor(fern_ID),y = Abundance, fill = guild))+
    geom_bar(stat="identity")+
    labs(x="Samples",y="Relative abundance")+
    theme(axis.text.x = element_text(angle=90, hjust=1))+
    xlab("Fern ID")+
    guides(fill=guide_legend(title="Guild"))+
    theme(text = element_text(size = 10))+
    scale_fill_manual(values = randomColor(count = 92))
  # 
  ps.prop_fil <- tophic_group_phy("Symbiotroph")
  dist.abund <- vegdist(ps.prop_fil@otu_table, method = "bray")
  dist.space <- dist(ps.prop_fil@sam_data[, c(4, 27:28)], method = "euclidean")

  print(mantel(dist.abund, dist.space, method = "spearman", permutations = 9999, na.rm = TRUE))
  lm()
  
spe <- ps.prop@otu_table %>% t() %>% as.data.frame()  %>% mutate(trophic_mode = left_join(colnames(ps.prop@otu_table) %>% as.data.frame() %>% set_colnames("OTU"), unique(OTU_funguild[, c("OTU", "trophicMode")]))$trophicMode) %>% group_by(trophic_mode) %>% summarise_all(sum) %>% as.data.frame()
spe$trophic_mode[8] <- "N"
rownames(spe) <- spe[, 1]
spe <- spe[, -1]
spe <- t(spe)
NMDS <- metaMDS (spe)
ef <- envfit (NMDS, as.data.frame(ps.prop@sam_data)[c(1, 3), -25:28])
ordiplot (NMDS, cex = 1, type = 't')

plot(NMDS, type = "n") #displays empty ordination space
points(NMDS, "sites", pch=19, col="black", cex = 2)
text(NMDS, "species", col="darkgrey", cex = 2)
ordiellipse(NMDS, unlist(as.data.frame(ps.prop@sam_data)[-(25:28), 3]), kind="se", conf=0.95, lwd=2, col=c("#ADC178", "#F1DCA7", "#997B66"), label=T, cex = 3)



summary(ord.nmds.bray)
bray.dist<-vegdist(spe, method='bray')
#PERMANOVA
as.data.frame(adonis2(bray.dist~ (ps.prop@sam_data$dry_mass) + ps.prop@sam_data$pH_ratio + ps.prop@sam_data$C_N_ratio, permutations=9999))


 #RDA 
  # vec.dist.abun <- as.vector(dist.abund)
  # vec.dist.space <- as.vector(dist.space)
  # combine <- as.data.frame(cbind(vec.dist.abun, vec.dist.space))
  # hist(scale(vec.dist.abun))
  # hist(log10(vec.dist.space))
  # print(plot(scale(vec.dist.abun) ~ (log10(vec.dist.space))))
  # 
  # ggscatter(combine, x = 'vec.dist.space', y = 'vec.dist.abun',
  #           add = "reg.line", conf.int = TRUE, color = "#609067")+
  #   stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, label.x.npc = "left",
  #            label.y.npc = "top", size =10)+
  #   labs(x = "Distance", y = "Bray-curtis distance of OTU composition\n of trophic mode include symbiotroph")+
  #   theme(text = element_text(size = 30))
  # 
  # allTaxa <- taxa_names(rarefied_ps)
  # allTaxa <- allTaxa[!(allTaxa %in% OTU_funguild_fil)]
  # ps.prop_fil <-  prune_taxa(allTaxa, rarefied_ps)
  # env <- rarefied_ps@sam_data
  # spe <-  rarefied_ps@otu_table
  # spe.log <- log1p (spe)
  # spe.log.hell <- decostand (spe.log, 'hell')
  # tbRDA <- rda (spe.log.hell ~ dry_mass + pH_ratio + C_N_ratio, data = samdf[-(25:28),c(3, 14, 25, 26)])
  # head (summary (tbRDA))  # prints first lines of the summary of tbRDA
  # ordiplot(tbRDA, type = 'p')
  # par(mfrow = c(1, 1))
  # anova(tbRDA)
  # adjR2.tbrda <- RsquareAdj (tbRDA)$adj.r.squared
  # 
  # sel.fs <- forward.sel (Y = spe.log.hell, X = env[,c('dry_mass', 'pH_ratio', 'C_N_ratio')], adjR2thresh = adjR2.tbrda)
  # 
  # test_axis <- anova (tbRDA, by = 'terms')
  # test_axis.adj <- test_axis
  # test_axis.adj$`Pr(>F)` <- p.adjust (test_axis$`Pr(>F)`, method = 'holm')
  # test_axis.adj
  # 
  # 
  # 
  # #rda for margin
  # env <- samdf[1:24, c(14, 15, 20)]  # select only two explanatory variables
  # test_each <- t(apply (env, 2, FUN = function (x) as.matrix (anova (rda (spe.log.hell ~ x))[1:4][1,])))
  # test_each <- as.data.frame (test_each)
  # names (test_each) <- c("Df", "Variance", "F", "Pr(>F)")
  # test_each.adj <- test_each
  # test_each.adj$`Pr(>F)` <- p.adjust (test_each$`Pr(>F)`, method = 'holm')
  # test_each.adj <- test_each.adj[order (test_each.adj$Variance, decreasing = TRUE),]
}

chart.Correlation(samdf[, c("dry_mass", "pH_ratio", "C_N_ratio")])

disper <-betadisper(bray.dist, group=ps.prop@sam_data$fern.size, type = 'centroid')
anova(disper)
permutest(disper)
plot(disper)

adonis2(dist(disper$distances) ~ as.factor(ps.prop@sam_data$fern.size))
pairwise_size <- pairwise.adonis(bray.dist, as.factor(ps.prop@sam_data$fern.size))


PlotSpecies<-ggplot(Stability_species, aes(x=Sample,y=Abundance,fill=trophicMode))+
  geom_bar(stat="identity")+
  # facet_wrap(~ Time, scales = "free_x", nrow = 1)+
  labs(x="Samples",y="Relative abundance")+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  guides(fill=guide_legend(ncol=3))+
  scale_fill_manual(values=colors_species_n)+
  theme(legend.text=element_text(size=6))
PlotSpecies


funguild_query("endophyte", "guild", db = Stability_species)$fern.size
