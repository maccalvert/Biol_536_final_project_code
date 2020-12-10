#Set working environment 
setwd("~/Box/Comp_bio/Final project/")
#load in the feature identification scripts from Townes et al. 2019 - this uses the multinomial deviance estimation
# better than constant experions / highly variable expression 
#source("scrna2019-master/util/functions_genefilter.R")

#Basic workflow for performing a GLM-PCA dimension reduction using the scry package developed throug Townes et al. 2019

#First, we'll make sure that the data follows a normal distrubtion: below is the code reproduced from Townes 

#Make sure you follow the data downloading steps in 01_data_loading.Rmd

# Data loading 
library(SingleCellExperiment)
source("./util/functions.R") #get_10x_readcounts function
fp<-file.path
bp<-"./real/zheng_2017_monocytes"

sce<-get_10x_readcounts(fp(bp,"data/hg19"),fp(bp,"data/cd14_monocytes_molecule_info.h5"))
saveRDS(sce,fp(bp,"data/01_sce_all_genes_all_cells.rds"))




suppressPackageStartupMessages(library(SingleCellExperiment))
library(ggplot2); theme_set(theme_bw())
library(DuoClustering2018)

require(scry)

sce<-sce_full_Zhengmix4eq()

celegans <- readRDS("~/Box/Comp_bio/Final project/lineage_tree_reconstruction-main/data/eset_ABprp.rds")

#Lets see the distribution of our UMIs / reads. 
library(BiocGenerics)
library(Biobase)
library(Matrix)
library(ggExtra)
library(ggplot2)
cur_eset <- readRDS("~/Box/Comp_bio/Final project/lineage_tree_reconstruction-main/data/eset_ABprp.rds")
cur_tree <- readRDS("~/Box/Comp_bio/Final project/lineage_tree_reconstruction-main/data/subtree_ABprp.rds")
cur_eset@assayData$exprs[1:5, 1:5] # The exprs slot contains the raw count matrix in sparse matrix format
cur_eset@assayData$norm_exprs[1:5, 1:5] # The norm_exprs slot contains the log2 transformed normalized count matrix in sparse matrix format

#Make a really simple version of the single cell experiment objects so that we can 
# use the feature identification functions from Townes et al. 2019. 
#Counts 
mat <- as.matrix(cur_eset@assayData$exprs)
log_mat <- log(mat)
#Lets see if it follows some of the predictions of a multinomial distribution: 
#Prediction 1: Fraction of zeros is inversely related to the total UMIs in a sample
#First get all of the UMI counts for each droplet 
umi_droplet<-colSums(mat)
#Next get the fraction of zeros 
zeroFrac_droplet <- colMeans(mat==0)
#Lets plot the relationship now. 
pd<-data.frame(sz=umi_droplet,pz=zeroFrac_droplet)
plt <- ggplot(pd,aes(x=sz,y=pz))+
  geom_point()+
  theme_bw()+
  xlab("total UMI per droplet")+
  ylab("fraction of zeros per droplet")+
  scale_x_log10()+
  theme(text = element_text(size=20))
#Add some histograms to demonstrate the data. 
ggExtra::ggMarginal(plt,type="histogram",fill="white",bins=100)
#Ok, it mostly matches the prediction and looks very similar to the data we saw 
# in Townes

#Next we'll test the predction that probability of a gene having a zero count is a 
# decreasing function of its mean expression. 

#downsample to normalize droplet size (total UMI)
Yds<-Down_Sample_Matrix(mat)
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
#ggs("logmean_pzero_binom_monocytes.pdf")

png(fp(pth,"logmean_pzero_monocytes.png"),width=6,height=4)
#same plot but using base plot
with(pd,plot(log_mean,frac_zero,xlab="log of mean expression",ylab="fraction of zero droplets",cex=1.5, cex.lab = 1.5, cex.axis = 1.5))
curve(predict_zeros_binom,from=xlo,to=xhi,col="blue",lwd=4,add=TRUE)
curve(predict_zeros_poi,from=xlo,to=xhi,col="green",lwd=3,lty=2,add=TRUE)
curve(predict_zeros_nb(x,phi=4),from=xlo,to=xhi,col="red",lwd=3,lty=3,add=TRUE)
legend("bottomleft",c("Multinomial","Poisson","Negative Binomial"),lty=c(1,2,3),lwd=c(4,3,3),col=c("blue","green","red"), cex = 1.5)
dev.off()

#These also look great! 

#Ok next step is to look at the deviance for feature selection. 
#Lets get ranks for the deviance, highly expressed, and highly variable genes 
#Need to make a singlecellExpression object 
umi_droplet <- as.data.frame(umi_droplet)
colnames(umi_droplet)<-"nUMI"
#add cell linage data
cell_lin <- data.frame("lineage" = cur_eset@phenoData@data$lineage, row.names = rownames(cur_eset@phenoData@data))
#add time points
time_points <- data.frame("time_points" = pData(cur_eset)$time.point, row.names = row.names(pData(cur_eset)))
c_elegans <- SingleCellExperiment(assays=list(counts=mat, logcounts=log_mat), colData = list(umi_droplet, cell_lin))
cl <- scran::quickCluster(c_elegans)
c_elegans<-scran::computeSumFactors(c_elegans,clusters=cl)
c_elegans<-scater::logNormCounts(c_elegans)
gm<-rank_all_genes(c_elegans)
rowData(c_elegans)<-gm

#Plot the data from these 
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
  
#Ok, no that we have our features picked out, time for a GLM-PCA?
#First we'll run it on all of the data, and then we'll run it on a feature selected set


L<-10 #number of dimensions 
K<-8 #number of known clusters
Y <- counts(c_elegans)
#get rid of the zeros 
Y <- Y[rowSums(counts(c_elegans))!=0,]
system.time(res<-glmpca(Y,L,fam="poi",verbose=TRUE)) #16 min - 173 iterations 
plot(res$dev,type="l",log="y")
factors<-res$factors

sz<-colMeans(Y)
pd<-cbind(factors,celltype=cell_lin$lineage,pz=colMeans(Y==0),z=log10(sz))
ggplot(pd,aes(x=dim1,y=dim2,colour=celltype))+geom_point(show.legend=TRUE)
#if(sp) ggsave(fp(pth,"zheng8eq_glmpca12.pdf"),width=6,height=4)
ggplot(pd,aes(x=dim3,y=dim4,colour=celltype))+geom_point(show.legend=FALSE)
ggplot(pd,aes(x=pz,y=dim1,colour=celltype))+geom_point(show.legend=FALSE)
#if(sp) ggsave(fp(pth,"zheng8eq_pz_glmpca1.pdf"),width=4,height=4)
# kmeans_res<-kmeans(factors,K,nstart=25)
# cl<-kmeans_res$cluster
# ari(cl,pd$celltype) #.74
# mcl_res<-Mclust(factors,K)
# ari(mcl_res$classification,pd$celltype) #.74

#Ok, those looked terrible, just a point cloud with nothing really meaningful differentiating

#Lets try with the most "deviant" 500 genes - they are already ranked in the rowData field
#so we'll just get the rownames / gene names for the top 5. 
L<-5
Y <- counts(c_elegans)
#get rid of the zeros 
Y <- Y[rowSums(counts(c_elegans))!=0,]
dev_ind <- rownames(rowData(c_elegans))[order(rowData(c_elegans)$dev)][1:500]
Y <- Y[rownames(Y) %in% dev_ind,]
system.time(res<-glmpca(Y,L,fam="poi",verbose=TRUE)) #25 seconds - 113 iterations 
plot(res$dev,type="l",log="y")
factors<-res$factors

sz<-colMeans(Y)
pd<-cbind(factors,celltype=cell_lin$lineage,pz=colMeans(Y==0),z=log10(sz))
ggplot(pd,aes(x=dim1,y=dim2,colour=celltype))+geom_point(show.legend=TRUE)
#if(sp) ggsave(fp(pth,"zheng8eq_glmpca12.pdf"),width=6,height=4)
ggplot(pd,aes(x=dim2,y=dim3,colour=celltype))+geom_point(show.legend=T)
ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=T)

#Ok, if we try for time points, what do we get. 
pd<-cbind(factors,time.point=time_points,pz=colMeans(Y==0),z=log10(sz))
ggplot(pd,aes(x=dim1,y=dim2,colour=time_points))+geom_point(show.legend=TRUE)


#We can try the neg binomial too. 
system.time(res<-glmpca(Y,L,fam="nb",verbose=TRUE)) #25 seconds - 113 iterations 
factors<-res$factors

sz<-colMeans(Y)
#pd<-cbind(factors,celltype=cell_lin$lineage,pz=colMeans(Y==0),z=log10(sz))
pd<-cbind(factors,time.point=time_points,pz=colMeans(Y==0),z=log10(sz))
pdf("test.pdf")
ggplot(data = pd, aes(x=dim1, y=dim2, colour=time_points))+
  geom_point()
dev.off()
ggplot(pd,aes(x=dim1,y=dim2,colour=time_points))+geom_point(show.legend=F)
#if(sp) ggsave(fp(pth,"zheng8eq_glmpca12.pdf"),width=6,height=4)
ggplot(pd,aes(x=dim2,y=dim3,colour=time_points))+geom_point(show.legend=T)
ggplot(pd,aes(x=dim1,y=dim3,colour=celltype))+geom_point(show.legend=T)

ggplot(pd,aes(x=dim3,y=dim4,colour=celltype))+geom_point(show.legend=FALSE)
ggplot(pd,aes(x=pz,y=dim1,colour=celltype))+geom_point(show.legend=FALSE)


#We can see whether these point clouds correspond at all to the "true tree" - this will be 
# on the other data type
pData(cur_eset)$nb_pca1 <- pd$dim1
pData(cur_eset)$nb_pca2 <- pd$dim2
pData(cur_eset)$nb_pca3 <- pd$dim3
#tree plust scatter
res <- make_tree_dimr(proj=pData(cur_eset)[, c("nb_pca2", "nb_pca3", "lineage")], left_tree = NULL, right_tree = NULL, top_tree = top_tree, shared.col = "lineage", colorBy = "lineage", tree.color = lineage_color, label.time.cut = 120,
                      plot.link = "top", shift.y.scale = 1/20, 
                      return_coords = F) 
res

pd <- data.frame("dim1"=pData(cur_eset)$nb_pca1,"dim2"=pData(cur_eset)$nb_pca2, "time.point"=pData(cur_eset)$embryo.time)
ggplot(pd,aes(x=dim1,y=dim2,colour=time.point))+geom_point(show.legend=TRUE)+
  scale_color_gradientn(colours = rainbow(5))








res<-list()
res$expr<-filterExpr(c_elegans)
res$hvg<-filterHVG(c_elegans,total_umi="total_counts")
res$devb<-filterDev(sce,dev="binomial")
res$devp<-filterDev(sce,dev="poisson")
f<-function(sn){
  s<-res[[sn]]
  x<-data.frame(gene=rownames(s),rank=seq_len(nrow(s)))
  colnames(x)[2]<-sn
  x
}
res2<-lapply(names(res),f)
rk<-Reduce(function(x,y){merge(x,y,by="gene",all=TRUE)},res2)
rownames(rk)<-rk$gene
rk$gene<-NULL