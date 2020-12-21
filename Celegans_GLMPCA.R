#Set working environment 
setwd("~/Box/Comp_bio/Final project/Biol_536_final_project_code/")

#Dependencies - many of these are not on CRAN, but on Bioconductor
library(SingleCellExperiment)
library(ggplot2); theme_set(theme_bw())
library(DuoClustering2018)
library(scry)
library(BiocGenerics)
library(Biobase)
library(Matrix)
library(ggExtra)
library(glmpca)
library(Seurat)
library(mclust)
#Load in some scripts from Townes et al. 2019 - https://github.com/willtownes/scrna2019 
source("./util/functions.R")
source("./util/functions_genefilter.R")
source("./util/clustering.R")


#Load in the data - the data is contained in an EXpression set object
cur_eset <- readRDS("./data/eset_ABprp.rds")
c_elegans <- SingleCellExperiment(assays=list(counts=cur_eset@assayData$exprs, logcounts=cur_eset@assayData$norm_exprs), 
                                  colData = list(nUMI=pData(cur_eset)$n.umi, cell.lin=pData(cur_eset)$lineage, embryo.time=pData(cur_eset)$embryo.time))

#Check some of the predictions of a multinomial distribution on the C. elegans dataset 
#1. Test whether the fraction of zeros in a droplet decreases with the total UMIs per droplet 
pd<-data.frame(sz=colData(c_elegans)$nUMI,pz=colMeans(counts(c_elegans)==0))
plt <- ggplot(pd,aes(x=sz,y=pz))+
  geom_point()+
  theme_bw()+
  xlab("total UMI per droplet")+
  ylab("fraction of zeros per droplet")+
  scale_x_log10()+
  theme(text = element_text(size=20))
#Add some histograms to get a nice spread of the data 
ggExtra::ggMarginal(plt,type="histogram",fill="white",bins=100)
#It appears that it does 

#prediction 2: The number of times that a zero is observed is a decreasing function of its expression across cells. 
Yds<-Down_Sample_Matrix(counts(c_elegans))
Yds<-Yds[rowSums(Yds)>0,]
#variance=mean, suggests poisson
m<-rowMeans(Yds); v<-apply(Yds,1,var)
summary(v/m)
plot(log(m),log(v),xlab="log(mean)",ylab="log(var)")
abline(0,1,col="blue") #poi
curve(x+log1p(exp(x)/5),from=-8,to=3,add=TRUE,lty=2,col="red") #nb

N<-median(colSums(Yds))
predict_zeros_binom<-function(x){(1-exp(x)/N)^N} #binomial
predict_zeros_poi<-function(x){exp(-exp(x))}
predict_zeros_nb<-function(x,phi=2){
  exp(-phi*log1p(exp(x-log(phi))))
}
pd<-data.frame(log_mean=log(m),frac_zero=rowMeans(Yds==0))
xlo<-min(pd$log_mean)
xhi<-max(pd$log_mean)
xcv<-data.frame(x=c(xlo,xhi))
ggplot(xcv)+geom_point(data=pd,aes(x=log_mean,y=frac_zero),alpha=.5) +stat_function(aes(x,color="bin"),fun=predict_zeros_binom) +stat_function(aes(x,color="poi"),fun=predict_zeros_poi) +stat_function(aes(x,color="nb"),fun=predict_zeros_nb) #+scale_color_manual("model",breaks=c("bin","poi","nb"),values=c("blue","green","red"))

# png(fp(pth,"logmean_pzero_monocytes.png"),width=6,height=4)
# #same plot but using base plot
# with(pd,plot(log_mean,frac_zero,xlab="log of mean expression",ylab="fraction of zero droplets",cex=1.5, cex.lab = 1.5, cex.axis = 1.5))
# curve(predict_zeros_binom,from=xlo,to=xhi,col="blue",lwd=4,add=TRUE)
# curve(predict_zeros_poi,from=xlo,to=xhi,col="green",lwd=3,lty=2,add=TRUE)
# curve(predict_zeros_nb(x,phi=4),from=xlo,to=xhi,col="red",lwd=3,lty=3,add=TRUE)
# legend("bottomleft",c("Multinomial","Poisson","Negative Binomial"),lty=c(1,2,3),lwd=c(4,3,3),col=c("blue","green","red"), cex = 1.5)
# dev.off()

#Ok, that also matches our predictions! 

#Next up is to perform feature selection. We'll be using the deviance criteria - but this may not be the best measure to capture cell to cell 
# variation when there are so many cell types and perhaps just differentiation would be best. 
cl <- scran::quickCluster(c_elegans)
c_elegans<-scran::computeSumFactors(c_elegans,clusters=cl)
c_elegans<-scater::logNormCounts(c_elegans)
#This script from Townes et al. performs gene ranking by constant expression, differential expression (mean to variance ratio) and deviance
gm<-rank_all_genes(c_elegans)
rowData(c_elegans)<-gm


#Lets just visualize this: 
pd<-meanvar_plotdata(c_elegans,G=1000)
ggplot(pd,aes(x=m,y=vmr,colour=criteria))+
  geom_point(alpha=.9)+
  xlab("average normalized expression")+
  ylab("variance to mean ratio")+
  theme(legend.position=c(0.2,.8))+
  scale_color_manual(values=c("orange","red","blue","gray"))+
  scale_x_log10()+
  scale_y_log10()+
  theme(text = element_text(size=16))


#Now that the featers have been picked out, we can run a glm pca
#We'll use the poison liklihood for now. 
# Still need to figure out exactly what chaning the dimension parameter (L) does

L<-10
Y <- counts(c_elegans)
#get rid of the zeros 
Y <- Y[rowSums(counts(c_elegans))!=0,]
dev_ind <- rownames(rowData(c_elegans))[order(rowData(c_elegans)$dev)][1:500]
Y <- Y[rownames(Y) %in% dev_ind,]
system.time(res<-glmpca(Y,L,fam="poi",verbose=TRUE)) #25 seconds - 113 iterations 
plot(res$dev,type="l",log="y")
factors<-res$factors

#Visualize.
sz<-colMeans(Y)
pd<-cbind(factors,celltype=colData(c_elegans)$cell.lin,pz=colMeans(Y==0),z=log10(sz))
ggplot(pd,aes(x=dim1,y=dim2,colour=celltype))+geom_point(show.legend=TRUE)
#if(sp) ggsave(fp(pth,"zheng8eq_glmpca12.pdf"),width=6,height=4)
ggplot(pd,aes(x=dim2,y=dim3,colour=celltype))+geom_point(show.legend=T)
ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=T)

#3d visualization of latent dimensions - to see if we can determine any further structure in the 
# emebedding
fig <- plot_ly(pd, x = ~dim1, y = ~dim2, z = ~dim3, color = ~celltype, 
               marker= list(size=10))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'dim1'),
                                   yaxis = list(title = 'dim2'),
                                   zaxis = list(title = 'dim3')))
fig

#Run clustering algorthms and an ARI test. 
seurat_clust <- seurat_cluster(pd[,1:10], res=2)
# res_0.1 <- subset(seurat_clust, resolution == 0.1)
# res_0.5 <- subset(seurat_clust, resolution == 0.5)
# res_1 <- subset(seurat_clust, resolution == 1)

#We then compare those clusters to the true cell identies using the Rand index
ari <- adjustedRandIndex(seurat_clust$cluster, c_elegans$cell.lin)
# [1] 0.2941941

# Lets look at embryo time. 
sz<-colMeans(Y)
pd<-cbind(factors,embryo.time=colData(c_elegans)$embryo.time,pz=colMeans(Y==0),z=log10(sz))
ggplot(pd,aes(x=dim1,y=dim2,colour=embryo.time))+geom_point(show.legend=TRUE)+
  scale_color_gradientn(colours = rainbow(5))
#if(sp) ggsave(fp(pth,"zheng8eq_glmpca12.pdf"),width=6,height=4)
ggplot(pd,aes(x=dim2,y=dim3,colour=celltype))+geom_point(show.legend=T)
ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=T)



### OLD STUFF 

#NWe can try running the GLM-PCA on the genes with the most variable expression.  
# L<-5
# Y <- counts(c_elegans)
# #get rid of the zeros 
# Y <- Y[rowSums(counts(c_elegans))!=0,]
# dev_ind <- rownames(rowData(c_elegans))[order(rowData(c_elegans)$expr)][1:1000]
# Y <- Y[rownames(Y) %in% dev_ind,]
# system.time(res<-glmpca(Y,L,fam="poi",verbose=TRUE)) #25 seconds - 113 iterations 
# plot(res$dev,type="l",log="y")
# factors<-res$factors
# 
# sz<-colMeans(Y)
# pd<-cbind(factors,celltype=colData(c_elegans)$cell.lin,pz=colMeans(Y==0),z=log10(sz))
# ggplot(pd,aes(x=dim1,y=dim2,colour=celltype))+geom_point(show.legend=F)
# ggplot(pd,aes(x=dim2,y=dim3,colour=celltype))+geom_point(show.legend=F)
# ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=F)
# 
# library(plotly)
# 
# 
# 
# 
# While it looks there is some structure to the data, we are not getting 
# specific clusters by cell type. - which makes some sense given that we 
# collecting data on a time series. 
# To evaluate the "chops" of this method - lets just try filtering out 
# only the differentiated cell types

#Subset by the "differentiated" cell types - this is lineages that have 10 characters - each 
# step appears to increase by one character. 
# 
# dif_cells <- names(table(colData(c_elegans)$cell.lin))[nchar(names(table(colData(c_elegans)$cell.lin))) == 10]
# 
# c_elegans_dif <- c_elegans[,colData(c_elegans)$cell.lin %in% dif_cells]
# #First we can try it on the "globally" ranked genes - meaning genes ranked for their deviance, etc. on 
# # all cells. We can try reranking on just these cells if that fails. 
# L<-10
# Y <- counts(c_elegans_dif)
# #get rid of the zeros 
# Y <- Y[rowSums(counts(c_elegans_dif))!=0,]
# dev_ind <- rownames(rowData(c_elegans_dif))[order(rowData(c_elegans_dif)$dev)][1:500]
# Y <- Y[rownames(Y) %in% dev_ind,]
# system.time(res<-glmpca(Y,L,fam="poi",verbose=TRUE)) #25 seconds - 113 iterations 
# plot(res$dev,type="l",log="y")
# factors<-res$factors
# 
# sz<-colMeans(Y)
# pd<-cbind(factors,celltype=colData(c_elegans_dif)$cell.lin,pz=colMeans(Y==0),z=log10(sz))
# ggplot(pd,aes(x=dim1,y=dim2,colour=celltype))+geom_point(show.legend=F)
# ggplot(pd,aes(x=dim2,y=dim3,colour=celltype))+geom_point(show.legend=F)
# ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=F)
# 
# fig <- plot_ly(pd, x = ~dim1, y = ~dim2, z = ~dim3, color = ~celltype, 
#                marker= list(size=10))
# fig <- fig %>% add_markers()
# fig <- fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
#                                    yaxis = list(title = 'Gross horsepower'),
#                                    zaxis = list(title = '1/4 mile time')))
# fig
# 
# #There is decent separation. Lets try redoing variable gene analysis and see if that gives
# # more separation. 
# 
# cl <- scran::quickCluster(c_elegans_dif)
# c_elegans_dif<-scran::computeSumFactors(c_elegans_dif,clusters=cl)
# c_elegans_dif<-scater::logNormCounts(c_elegans_dif)
# #This script from Townes et al. performs gene ranking by constant expression, differential expression (mean to variance ratio) and deviance
# gm<-rank_all_genes(c_elegans_dif)
# rowData(c_elegans_dif)<-gm
# 
# L<-10
# Y <- counts(c_elegans_dif)
# #get rid of the zeros 
# Y <- Y[rowSums(counts(c_elegans_dif))!=0,]
# dev_ind <- rownames(rowData(c_elegans_dif))[order(rowData(c_elegans_dif)$dev)][1:1500]
# Y <- Y[rownames(Y) %in% dev_ind,]
# system.time(res<-glmpca(Y,L,fam="poi",verbose=TRUE)) #25 seconds - 113 iterations 
# plot(res$dev,type="l",log="y")
# factors<-res$factors
# 
# sz<-colMeans(Y)
# pd<-cbind(factors,celltype=colData(c_elegans_dif)$cell.lin,pz=colMeans(Y==0),z=log10(sz))
# ggplot(pd,aes(x=dim1,y=dim2,colour=celltype))+geom_point(show.legend=F)
# ggplot(pd,aes(x=dim2,y=dim3,colour=celltype))+geom_point(show.legend=F)
# ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=F)
# 
# fig <- plot_ly(pd, x = ~dim1, y = ~dim2, z = ~dim3, color = ~celltype, 
#                marker= list(size=10))
# fig <- fig %>% add_markers()
# fig <- fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
#                                    yaxis = list(title = 'Gross horsepower'),
#                                    zaxis = list(title = '1/4 mile time')))
# fig
# 
# 
# #These appear maybe marginally better. z
# #What happens if we apple a clustering algorthm and a discriminant function. 
# #Use the seurat clustering (Which implementes the Louvain algorthm) function to begin adjusting for our performance
# seurat_clust <- seurat_cluster(pd[,1:10], res=2)
# # res_0.1 <- subset(seurat_clust, resolution == 0.1)
# # res_0.5 <- subset(seurat_clust, resolution == 0.5)
# # res_1 <- subset(seurat_clust, resolution == 1)
# 
# #We then compare those clusters to the true cell identies using the Rand index
# ari <- adjustedRandIndex(seurat_clust$cluster, c_elegans_dif$cell.lin)

