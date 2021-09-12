setwd("E:/UK doc/MRes Project 2 AstraZeneca/A Thesis")
df <- read.csv(file = 'new_ori_for_impute.csv', row.names = 1, header= TRUE)

df <- df[,-c(1,2,3,4,5,6,7,8)] # remove unwanted cols for imputation

library(mice)
md.pattern(df, rotate.names = TRUE) # visualize missing value pattern

library(VIM)
aggr_plot <- aggr(df, col=c('navyblue','red'), 
                  numbers=TRUE, sortVars=TRUE, 
                  labels=names(data), cex.axis=1, gap=3, 
                  ylab=c("Histogram of missing data","Pattern")) # another way to look at the missing value pattern


# imputation by mice package
imputed <- mice(df, 
                m = 1,   # default multiple imputation
                method="pmm",        # predictive mean matching
                maxit = 50,
                seed = 2021)

summary(imputed)


densityplot(imputed) # plot predicted value and original/observed value
densityplot(imputed, subset=.imp==1) # plot predicted value only

completedData <- complete(imputed,1) # impute with the 1st imputation
write.csv(completedData,'imputation.csv')


##### PCA
library(ropls)
df_o <- read.csv(file = 'new_264.csv', row.names = 1, header= TRUE) # mouse data set
df_o <- read.csv(file = 'new_110.csv', row.names = 1, header= TRUE) # rat data set

df_o <- read.csv(file = 'new_250.csv', row.names = 1, header= TRUE) # mouse(E3b ligase) data set
df_o <- read.csv(file = 'new_44.csv', row.names = 1, header= TRUE) # rat (E3b ligase data set)

df <- df_o[,-c(1,2,3,4,5,6,7)]
df <- scale(log(df)) # scaling and log-transformation

df = subset(df, select = -c(Rat_Bioav.F..,Mouse_Bioav.F..) ) # remove Bioav (Y) here, identify outliers solely on X-descriptors
Bioav.pca <- opls(df,
                  predI = 4, # PC components
                  algoC = c("svd"),  # for pca without missing values
                  scaleC = c("none") # no scaling needed
)

# coloured PCA by categorical variables
E3 <- df_o[, "E3_Ligase"]
plot(Bioav.pca, typeVc = c("x-score"), parAsColFcVn = E3, parEllipsesL = F,parPaletteVc = c("green4",'magenta'))
plot(Bioav.pca, typeVc = c("x-score"), parAsColFcVn = E3, parEllipsesL = F,parPaletteVc = c("black", "green4",'magenta'))
plot(Bioav.pca, typeVc = c("x-score"), parAsColFcVn = E3, parEllipsesL = F,parPaletteVc = c('magenta'))

ion <- df_o[, "Ion_Class"]
plot(Bioav.pca, typeVc = "x-score", parCompVi = c(1, 2), parAsColFcVn = ion, parEllipsesL = F,parPaletteVc = c("blue", "black",'orange'))

F_cate <- df_o[, "Mouse_F_cate"]
plot(Bioav.pca, typeVc = "x-score", parCompVi = c(1, 2), parAsColFcVn = F_cate, parEllipsesL = T,parPaletteVc = c("orange",'green','blue','red'))

F_cate <- df_o[, "Rat_F_cate"]
plot(Bioav.pca, typeVc = "x-score", parCompVi = c(1, 2), parAsColFcVn = F_cate, parEllipsesL = T ,parPaletteVc = c("orange",'green','blue','red'))


plot(Bioav.pca, typeVc = c("overview"))
plot(Bioav.pca, typeVc = c("x-loading"),parCompVi = c(1, 2))

Bioav.pca@summaryDF
Bioav.pca@descriptionMC

write.csv(Bioav.pca@modelDF, "264 PCA variance explained.csv")
write.csv(Bioav.pca@loadingMN, "264 PCA loading.csv")

write.csv(Bioav.pca@modelDF, "110 PCA variance explained.csv")
write.csv(Bioav.pca@loadingMN, "110 PCA loading.csv")

write.csv(Bioav.pca@modelDF, "250 PCA variance explained.csv")
write.csv(Bioav.pca@loadingMN, "250 PCA loading.csv")

write.csv(Bioav.pca@modelDF, "44 PCA variance explained.csv")
write.csv(Bioav.pca@loadingMN, "44 PCA loading.csv")


##### OPLS
library(dplyr)
seed = 2021
df_o <- read.csv(file = 'new_250.csv', row.names = 1, header= TRUE)
df_o <- read.csv(file = 'new_44.csv', row.names = 1, header= TRUE)

row_names_df_to_remove<-c('AZ5220','AZ2903','AZ5291') # for 250
row_names_df_to_remove<-c('AZ3435') # for 44
df_o <- df_o[!(row.names(df_o) %in% row_names_df_to_remove),]

df <- df_o[,-c(1,2,3,4,5,6)] # for mouse
df <- df_o[,-c(1,2,3,4,5,7)] # for rat
df <- scale(log(df))

df_x <- as.matrix(subset(df, select = -c(Mouse_Bioav.F..) ))
df_y <- as.matrix(subset(df, select = c(Mouse_Bioav.F..) ))

df_x <- as.matrix(subset(df, select = -c(Rat_Bioav.F..) ))
df_y <- as.matrix(subset(df, select = c(Rat_Bioav.F..) ))


n <- as.integer(0.8*nrow(df_x)) # prepare for 80% 20% splitting
Bioav.opls <- opls(df_x, # x-variable
                   df_y, # y-response
                   predI = 1, # no. of predictive component in OPLS
                   orthoI = 3, # no. of orthogonal component in OPLS
                   algoC = c("nipals"),
                   crossvalI = 20,
                   permI = 0, 
                   scalC = c("none"),
                   subset = 1:n # train/test partition selected
)


plot(Bioav.opls, typeVc = c("overview"))
plot(Bioav.opls, typeVc = c("x-score"))
plot(Bioav.opls, typeVc = c("x-loading"))
plot(Bioav.opls, typeVc = c("predict-train"))
plot(Bioav.opls, typeVc = c("predict-test"))

write.csv(Bioav.opls@modelDF, "247 OPLS components.csv")
write.csv(Bioav.opls@vipVn , "247 OPLS VIPs.csv")
write.csv(Bioav.opls@coefficientMN, "247 OPLS coefficients.csv")
write.csv(Bioav.opls@weightMN, "247 OPLS weights.csv")
Bioav.opls@descriptionMC


write.csv(Bioav.opls@modelDF, "43 OPLS components.csv")
write.csv(Bioav.opls@vipVn , "43 OPLS VIPs.csv")
write.csv(Bioav.opls@coefficientMN, "43 OPLS coefficients.csv")
write.csv(Bioav.opls@weightMN, "43 OPLS weights.csv")
Bioav.opls@descriptionMC


# plot of VIP scores in descending order
library(ggplot2)
Vip <- read.csv(file = '247 OPLS VIPs.csv', header= TRUE)
Vip <- read.csv(file = '43 OPLS VIPs.csv', header= TRUE)

colnames(Vip)
names(Vip)[names(Vip) == "X"] <- "Descriptors"
names(Vip)[names(Vip) == "x"] <- "VIP"

p<-ggplot(data=Vip, aes(x=reorder(Descriptors, -VIP), y=VIP)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(size = 9, angle=45, hjust=1))+
  labs(title = 'VIP scores in the mouse OPLS model') +
  theme(plot.title = element_text(color = "black", size = 15, face = "bold", hjust = 0.5),
        plot.caption = element_text(color = "dark grey", size = 10, face = "italic", hjust = 1)) +
  geom_hline(yintercept=1, linetype="dashed", color = "red", size=2)
p


p<-ggplot(data=Vip, aes(x=reorder(Descriptors, -VIP), y=VIP)) +
  geom_bar(stat="identity") +
  theme(axis.text.x=element_text(size = 9, angle=45, hjust=1))+
  labs(title = 'VIP scores of descriptors in the OPLS model (rat)') +
  theme(plot.title = element_text(color = "blue4", size = 13, face = "bold", hjust = 0.5),
        plot.caption = element_text(color = "dark grey", size = 10, face = "italic", hjust = 1)) +
  geom_hline(yintercept=1, linetype="dashed", color = "blue4", size=2)
p



# plot comparison between 'before log-transformation' and 'after log-transformation'
library(ggplot2)
library(cowplot)
fig1 <- qplot(x= Mouse_Bioav.F..,data = df_train,bins = 50) + 
  ggtitle("Distribution of mouse F% before log-transformation") + xlab("Mouse_Bioav_F%") + ylab("Count")+
  labs(caption = ' ') +
  theme(plot.title = element_text(color = "black", size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = "black", size = 10, face = "italic", hjust = 0.7)) 

fig2 <- qplot(x= log(Mouse_Bioav.F..),data = df_train, bins = 50)+ 
  ggtitle("Distribution of mouse F% after log-transformation") + xlab("log(Mouse_Bioav_F%)") + ylab("Count")+ 
  labs(caption = 'Data: 80% training set from E3b-PROTACs') +
  theme(plot.title = element_text(color = "red", size = 13, face = "bold", hjust = 0.5),
        plot.caption = element_text(color = "dark grey", size = 10, face = "italic", hjust = 1)) 
plot_grid(fig1, fig2, labels = "AUTO")



















