barcode_list_generator <- function(sample = sample, barcode = barcode, output = output){
  barcode_F <- barcode[seq(2, 70, 2), ]
  final_F_1 <- c() 
  final_F_1[seq(1, 70, 2)] <- sample
  final_F_1[seq(2, 70, 2)] <- barcode_F
  write.csv(final_F_1 , file = output)
}

sample <- read.csv('/Users/yupeitseng/Documents/RA/FS_PSF/Wet_lab/PSF_bac_sample.csv', header = F)
sample <- apply(sample, 2, FUN = function(x){paste('>', x, sep = '')})
barcode_F <- read.csv('/Users/yupeitseng/Documents/RA/FS_PSF/bio_info/barcode/B2_barcode_F.txt', header = F)
barcode_R <- read.csv('/Users/yupeitseng/Documents/RA/FS_PSF/bio_info/barcode/B2_barcode_R.txt', header = F)
barcode_list_generator(sample[, 1], barcode_F, '/Users/yupeitseng/Downloads/L1_barcode_F.csv')
barcode_list_generator(sample[, 1], barcode_R, '/Users/yupeitseng/Downloads/L1_barcode_R.csv')
barcode_list_generator(sample[, 2], barcode_F, '/Users/yupeitseng/Downloads/L2_barcode_F.csv')
barcode_list_generator(sample[, 2], barcode_R, '/Users/yupeitseng/Downloads/L2_barcode_R.csv')
