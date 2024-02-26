library(stringr)
library(multcompView)
library(hrbrthemes)
indi <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_individual.csv', header = T)
layer <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_layer.csv', header = T)
tree <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_tree.csv')

#indi_info
layer$dry_mass <- rowSums(layer[, 2:4], na.rm = T) 
layer <- layer[, -(2:4)]
layer$fern_ID <- sub('-.*', '', layer$sample_ID)
pH <- layer
layer <- layer %>% group_by(fern_ID) %>% summarise(sum(dry_mass))
pH <- pH %>% group_by(fern_ID) %>% summarise(mean(pH, na.rm = T))
layer <- left_join(layer, pH)
colnames(layer)[2:3] <- c('dry_mass', 'pH')
indi <- indi[, -3]
colnames(indi)[2] <- 'fern_ID'
indi$fern_ID <- as.character(indi$fern_ID)
indi <- left_join(indi, layer)


indi$leaf_length_mean <- rowSums(indi[, 6:8])/3
indi$nest_volume <- indi$nest.length * indi$nest.height * indi$nest.width
indi <- left_join(indi, tree)

write.csv(indi, file = '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi.csv')

#BNF_info (layer info)----
indi <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info_indi.csv', header = T)
layer <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_layer.csv', header = T)
layer$dry_mass <- rowSums(layer[, 2:4], na.rm = T) 
layer <- layer[, -(2:4)]
layer$fern_ID <- sub('-.*', '', layer$sample_ID)
colnames(indi)[4] <- 'fern_ID' 
indi$fern_ID <- as.character(indi$fern_ID)
layer <- left_join(layer, indi, by = 'fern_ID')
layer <- layer[, -5]
layer$layer <- sub('.*-', '', layer$sample_ID)
layer[which(!(layer$layer %in% c('U', 'M', 'L', 'O'))), 21] <- 'S'
layer[66:69, c(3, 4, 21)] <- NA

write.csv(layer, file = '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv')

#pH
samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv', header = T)
samdf <- samdf[, -1]
# theme_set(theme_bw())
samdf <- samdf[-23, ]

samdf_l <- samdf[which(samdf$fern.size == 'L'), ]
samdf_m <- samdf[which(samdf$fern.size == 'M'), ]

samdf_l$fern_ID <- as.character(samdf_l$fern_ID)
samdf_m$fern_ID <- as.character(samdf_m$fern_ID)

ggplot(samdf_l, aes(x=factor(fern_ID, level=c('12', '13', '19', '22', '25', '10', '7', '9', '11')), y=pH, color=factor(layer, level = c('U', 'M', 'L', 'O')))) + 
  geom_point(size=6) +
  geom_line(aes(group=layer))+
  theme(text = element_text(size = 30))+
  labs(x = "Fern ID", y = "pH", colour = 'Layer')+
  scale_color_manual(values = c("#D1C092",
                                         "#ACA08B",
                                         "#615117",
                                         "lightgrey"))
                                         
ggplot(samdf_m, aes(x=factor(fern_ID, level=c('2', '17', '20', '23', '8', '24', '5')), y=pH, color=factor(layer, level = c('U', 'L', 'O')))) + 
  geom_point(size=6) +
  geom_line(aes(group=layer))+
  theme(text = element_text(size = 30))+
  labs(x = "Fern ID", y = "pH", colour = 'Layer')+
  scale_color_manual(values = c("#D1C092",
                                         "#615117",
                                         "lightgrey"))

boxplot(pH ~ fern_ID, data = samdf_l)

box_plot <- function(dat = dat, main = main){
  aov_l <- aov(pH ~ layer , data = dat)
  summary(aov_l)
  tukey <- TukeyHSD(aov_l)
  return(tukey)
  cld <- multcompLetters4(aov_l, tukey)
  Tk <- group_by(dat, layer) %>%
    summarise(mean=mean(pH), quant = quantile(pH, probs = 0.75)) %>%
    arrange(desc(mean))
  cld <- as.data.frame.list(cld$layer)
  Tk$cld <- cld$Letters
  
  ggplot(dat, aes(layer, pH)) + 
    geom_boxplot() +
    labs(x="Layer", y="pH", title = main) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_text(data = Tk, aes(x = layer, y = quant, label = cld), size = 10, vjust=-1, hjust =-1) +
    theme(text = element_text(size = 30))
}
dat = samdf_m
main = 'Medium'
box_plot(samdf_l, 'Large')
box_plot(samdf_m, 'Medium')

#CN ratio----
samdf <- read.csv('/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv', header = T)
cn <- read.csv("/Users/yupeitseng/Documents/RA/birds_nest_fern/data/CN_ratio.csv")
cn <- cn %>% group_by(ID) %>% summarise(C = mean(C.), N = mean(N.), H = mean(H.)) %>% mutate(C_N = C/N)
colnames(cn)[1] <- "sample_ID" 
samdf <- left_join(samdf, cn, by = "sample_ID")
samdf <- samdf[, -1]

write.csv(samdf, '/Users/yupeitseng/Documents/RA/birds_nest_fern/data/BNF_info.csv')
