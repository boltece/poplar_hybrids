######################################## RDA ########################################

#straight from forester paper
#https://popgen.nescent.org/2018-03-27_RDA_GEA.html
install.packages("devtools")
library(devtools)
install.packages("data.table")
library(data.table)
install.packages("tidyverse")
library(tidyverse)
install.packages("vegan")
library(vegan)
install.packages("psych")
library(psych)
library(viridis)
install.packages("ggsci")
library(ggsci)
install.packages("ggcorrplot")
library(ggcorrplot)
require(ggrepel)
library(lme4)
library(readr)

source('~/Desktop/Imports.R')


############################### ENV data for RDA ##################################
setwd("./RDA/")


df_mod <- read.csv("./Modified_transect_assignments_with_ClimGeo_data.csv")
df_mod <- df_mod[,-1]

#env correlations

env_all <- df_mod[,3:30]
corr <- cor(env_all, use="complete.obs", method="pearson")


corr_NA <- round(corr,2)
write.csv(corr_NA,"./Correlations_env_variables_546_samples.csv")
corr_NA[which(abs(corr) < .75)] <- NA
abs(corr_NA)
write.csv(corr_NA, "./DF_correlations_env_546samples_NA_is_under0.75.csv")



head(df_mod)

MAT <- df_mod[,6]
TD <- df_mod[,9]
MWMT <- df_mod[,7]
MAP <- df_mod[,10]
RH <- df_mod[,28]
PAS <- df_mod[,22]
CMD <- df_modn[,27]

climate_7var <- cbind(MAT,TD, MAP, MWMT, PAS, CMD, RH)
clim7_scaled <- apply(climate_7var, 2, scale)
clim_df <-as.data.frame(clim7_scaled)


lat <- df_mod$Latitude
long <- df_mod$Longitude
elev <- df_mod$Elevation
geo_df <- cbind(long, lat, elev)
geo_scaled <- apply(geo_df,2, scale)
geo_scaled <- as.data.frame(geo_scaled)

geo_clim <- cbind(geo_scaled, clim_df)

################################ GENETIC data for RDA ##################################


df012<- fread("./poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012",sep="\t",data.table = F)
df012 <- df012[,-1]
dim(df012)

formated<-apply(df012,2,function(df) as.numeric(df))
dim(formated)


### writing genetic data to a text file 
df012_txt <- df012
df012_indv <- read.table("./poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012.indv")

rownames(df012_txt)<- df012_indv$V1
View(df012_txt[1:100,1:100])

df012_pos <- read.table("./poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012.pos")
df012_pos_v2 <- paste(df012_pos$V1, df012_pos$V2, sep="_")
df012_pos_v2<- as.data.frame(df012_pos_v2)

colnames(df012_pos_v2) <- "SNP_position"
df012_pos_v2_t <- t(df012_pos_v2)
colnames(df012_txt) <- df012_pos_v2_t

write.table(df012_txt, "./poplar_546trees_mac2_missingNONE_LDpruned.txt", sep="\t", col.names = T, row.names=T, quote=F)


############################### SCALE genetic data ####################################

colmean<-apply(formated,2,mean,na.rm=TRUE) 
## all missing data was filtered out using vcftools prior to analysis so na.rm=TRUE is irrelevant for this df

normalize<-matrix(nrow = nrow(formated),ncol=ncol(formated))

af<-colmean/2

for (m in 1:length(af)){
  nr<-formated[ ,m]-colmean[m]
  dn<-sqrt(af[m]*(1-af[m]))
  normalize[ ,m]<-nr/dn
}


scaled_df012 <- as.data.frame(normalize)
View(scaled_df012[1:20,1:20])







########################### partitioning of variance ################################



vp<- varpart(scaled_df012, ~ as.matrix(clim_df), ~ as.matrix(geo_scaled))
vp




#### RDA of geo and clim ####

n <- rda(formula = scaled_df012 ~.,scale=FALSE, center=TRUE, data = geo_clim)


anova(n)



saveRDS(n,'./RDA_poplar_546trees_clim_geo_scaled.RDS')

RsquareAdj(n)

summary(eigenvals(n,model='constrained'))

RDA1_gen <- sort(abs(summary(n)$biplot[,1]),decreasing=TRUE)
RDA1_gen

RDA2_gen <- sort(abs(summary(n)$biplot[,2]),decreasing=TRUE)
RDA2_gen

sort(abs(summary(n)$biplot[,1]) + abs(summary(n)$biplot[,2]),decreasing=TRUE)

RDA_df <- as.data.frame(summary(n)$biplot)
write.csv(RDA_df,'RDA_poplar_546trees_clim_geo_scaled.csv',row.names = F)


anova(n,by="axis")


###################################### ggplot ############################################


n <- readRDS('./RDA_poplar_546trees_clim_geo_scaled.RDS')


### file prep 

n_sum_pve <- summary(n)

RDA1_pve <- paste("RDA1 (",round((n_sum_pve$concont$importance[2,1]*100),digits=2),"%)", sep="")
RDA2_pve <- paste("RDA2 (",round((n_sum_pve$concont$importance[2,2]*100),digits=2),"%)", sep="")
RDA3_pve <- paste("RDA3 (",round((n_sum_pve$concont$importance[2,3]*100),digits=2),"%)", sep="")
n_sum <- summary(n)
n_sum <- scores(n, display=c("sp","sites", "bp")) 
rda_snp <- as.data.frame(n_sum$species)
rda_indv <- as.data.frame(n_sum$sites)


pca_df <- read.csv("./pca_df_final_546.csv")

rda_df2 <- data.frame(ID=as.character(df_mod$ID),ANC=as.character(pca_df$Pop),
                      TRANSECT=as.character(df_mod$Transect_SL)) 
### some samples collected to increase the sample size of balsamifera were labeled with a region
### We reassigned these samples to the closest transect
              
rda_indv <- cbind(rda_df2,rda_indv)
write.csv(rda_indv, "./RDA_indv_assignment_loadings_clim_geo_scaled_modified_transects.csv")
rda_indv <-read.csv("./RDA_indv_assignment_loadings_clim_geo_scaled_modified_transects.csv")
rda_biplot <- as.data.frame(n_sum$biplot)
rda_biplot$var <- row.names(rda_biplot)

basplot <- plot(n)
mult <- attributes(basplot$biplot)$arrow.mul




### label by ancestral coefficient


rda_plot <- ggplot(data=rda_indv, aes(x=RDA1, y=RDA2)) + 
  geom_point(aes(colour=ANC)) +
  scale_color_gradient(low="deepskyblue", high="blue4") +
  xlim(-15, 15) +
  ylim(-15,15) +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult * RDA1 * 0.5,y = 0, yend = mult * RDA2 * 0.5),
               arrow = arrow(length = unit(0.15, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult + mult/5) * RDA1, y = (mult + mult/5) * RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 4,fontface = "bold") + 
  xlab(RDA1_pve) + ylab(RDA2_pve) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

rda_plot  

ggsave('RDA_label_clim_geo_scaled.pdf',rda_plot,height=6,width=7,units='in')




### label by transect

col_trans <- c("#FFDB6D", "#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

rda_plot_transect <- ggplot(data=rda_indv, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv, aes(x=RDA1, y=-1*RDA2,fill=TRANSECT),pch=21,col='black',size=3) +
  xlim(-15, 15) +
  ylim(-15,15) +
  scale_fill_manual(values=col_trans) +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult * RDA1 * 0.95,y = 0, yend = mult * -1*RDA2 * 0.95), linewidth=1.0,
               arrow = arrow(length = unit(0.35, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult + mult/7) * RDA1, y = (mult + mult/7) * -1*RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 5,fontface = "bold") + 
  xlab(RDA1_pve) + ylab(RDA2_pve) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=12), 
    axis.title = element_text(size = 14, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

rda_plot_transect 

ggsave('./RDA_Transects_label_V2_clim_geo_scaled.pdf',rda_plot_transect,height=6,width=9,units='in')



