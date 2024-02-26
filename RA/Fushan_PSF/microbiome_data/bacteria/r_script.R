#import data (feature table, metadata, taxanomy)
feature_table_L1 <- read.csv(file = '/home/po/Documents/FS_bac_all/dada2_output_232-230/feature-table.txt', sep = '', header = T, check.names = TRUE, row.names = 1)
feature_table_L1 <- feature_table_L1[, -8]
tax_L1 <- read.csv(file='/home/po/Documents/FS_bac_all/taxa/taxonomy.csv', sep = '', header = F, check.names = T, row.names = 1)
tax_L1 <- tax_L1[-1, ]
tax_L1 <- tax_L1[, -8]

tax_L1 <- as.matrix(tax_L1)
TAX_L1 <- tax_table(tax_L1)
colnames(TAX_L1) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

metadata_L1 <- read.csv(file = '/home/po/Documents/field_soil_property_data.csv', header = T, check.names = T, row.names = 1) %>% filter(sample_ID %in% colnames(feature_table_L1))
meta_L1 <- as.data.frame(metadata_L1)
META_L1 <- sample_data(meta_L1)
rownames(META_L1) <- META_L1$sample_ID
OTU_L1 <- otu_table(feature_table_L1, taxa_are_rows = TRUE)

physeq_L1 <- phyloseq(OTU_L1, TAX_L1, META_L1)
physeq_L1 <- subset_taxa(physeq_L1, Kingdom != "Unassigned")
physeq_L1 <- subset_taxa(physeq_L1, Class != "D_2__Chloroplast")
physeq_L1 <- subset_taxa(physeq_L1, Family != "D_4__Mitochondria")

#draw physeq_L1 rarefaction Curve (raremax = 7877)
col <- c("darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
raremax = min(sample_sums(physeq_L1)) #get the minimum point

tab <- otu_table(physeq_L1)
class(tab) <- "matrix" 

tab <- t(tab)



rarecurve(tab, step=50, cex=0.5, xlab="Sequence sample size", ylab="Species richness", label=T, col=col)+
  geom_text(title("Rarefaction Curve"))+
  abline(v=raremax) #add vertical line with the minimum point
raremax
rarefied_physeq_L1 = rarefy_even_depth(physeq_L1, rngseed = 1, sample.size = raremax, replace = F)

##02 alpha diversity
plot_richness(rarefied_physeq_L2, x = "sample_ID", measures = c("Shannon", "Simpson"), color = "decay_cat", shape = 'soil')

#03 Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(rarefied_physeq_L1, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="PCoA", distance="jaccard")

plot_ordination(ps.prop, ord.nmds.bray, color = "pH", shape = 'sp',  title="Bray NMDS")
erie_bray <- phyloseq::distance(ps.prop, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(rarefied_physeq_L1))

# Adonis test
adonis2(erie_bray ~ decay_cat+sp, data = sampledf)
beta <- betadisper(erie_bray, sampledf$pH)
permutest(beta)

