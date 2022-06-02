#imaris analysis script:
#postn fiber quantification
# for the 20211202:
#read in the min dist transform to postn fibers
b16_s9_1600_600_z_2nd_Intensity_Min_Ch_10_Img_1 <- read_csv("D:/Krummel Lab/Microscopy/SP8/20211207_stain8_9_trem2_postn_mc38_v_b16/b16_stain98_apc/imaris_analysis/b16_s9-_1600_600_z_2nd_Statistics/b16_s9-_1600_600_z_2nd_Intensity_Min_Ch=10_Img=1.csv")
mc38_s9_1600_600_z_Intensity_Min_Ch_5_Img_1 <- read_csv("D:/Krummel Lab/Microscopy/SP8/20211207_stain8_9_trem2_postn_mc38_v_b16/mc38_stain9_apc_4/mc38_s9_1600_600_z_Statistics/mc38_s9_1600_600_z_Intensity_Min_Ch=5_Img=1.csv")

b16 = b16_s9_1600_600_z_2nd_Intensity_Min_Ch_10_Img_1$`Intensity Min`
mc38 = mc38_s9_1600_600_z_Intensity_Min_Ch_5_Img_1$`Intensity Min`
b16 = b16_s8_z_Intensity_Min_Ch_5_Img_1$`Intensity Min`
mc38 = mc38_s9_z1_Intensity_Min_Ch_5_Img_1$`Intensity Min`
mc38 = mc38_s9_1600_600_zoom2_Detailed$`Intensity Min`

dist = c(mc38,b16)
type = c(rep('MC38',length(mc38)),rep('B16F10',length(b16)))
df = data.frame(dist,type)
library(ggplot2)
df$type = factor(df$type,levels = c('MC38','B16F10'))
ggplot(data = df,aes(x=dist,fill=type))+geom_histogram(aes(y=..density..), color="white",alpha=0.5,position = 'identity',binwidth = 2)+
  geom_density(alpha=0.3,position = 'identity',bw=2)+xlab('Dist. to nearest POSTN surface')+theme_bw()+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=16),legend.position = 'None')
#for KS testing
ks.out = ks.test(mc38,b16)

#dist of cd11b+ mhcii hi lo to cd31 vessels histogram
# for the 20211202:
#read in the min dist transform to postn fibers
#take anytuing with y value<600 micron
pos_lo = read_csv("E:/Krummel Lab/Microscopy/SP8/20211116_cleared_stain2/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_lo_stats_whole_Statistics/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_lo_stats_whole_Position.csv")
pos_hi <- read_csv("E:/Krummel Lab/Microscopy/SP8/20211116_cleared_stain2/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_hi_stats_whole_Statistics_Statistics/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_hi_pos.csv")
dist_lo <- read_csv("E:/Krummel Lab/Microscopy/SP8/20211116_cleared_stain2/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_lo_stats_whole_Statistics/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_lo_stats_dist_vessel.csv")
dist_hi <- read_csv("E:/Krummel Lab/Microscopy/SP8/20211116_cleared_stain2/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_hi_stats_whole_Statistics_Statistics/mc38_stain2_512_2_10_TileScan 1 Merged_cd11b_mhcii_hi_stats_dist_vessel.csv")
#crop to specified region of interest
dist_lo = dist_lo[pos_lo$`Position Y`<600,]
dist_hi = dist_hi[pos_hi$`Position Y`<600,]
dist_lo = dist_lo[pos_lo$`Position X`<1000,]
dist_hi = dist_hi[pos_hi$`Position X`<1000,]
lo = dist_lo$`Intensity Min`
hi = dist_hi$`Intensity Min`
dist = c(lo,hi)
type = c(rep('MHCII Lo',length(lo)),rep('MHCII Hi',length(hi)))
df = data.frame(dist,type)
library(ggplot2)
#df$type = factor(df$type,levels = c('MC38','B16F10'))
ggplot(data = df,aes(x=dist,fill=type))+geom_histogram(aes(y=..density..), color="white",alpha=0.5,position = 'identity',binwidth = 5)+
  geom_density(alpha=0.3,position = 'identity',bw=2)+xlab('Dist. to nearest CD31 surface')+theme_bw()+
  theme(axis.text = element_text(size=16),axis.title = element_text(size=16))+scale_fill_manual(values = c('chartreuse4','magenta3'))
