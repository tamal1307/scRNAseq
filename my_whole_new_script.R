setwd("A:/try data/my analysis")
memory.limit(size=30000)
#load("A:/try data/my analysis/my_analysis_environment.RData")
library(ggplot2)
library(ggbiplot)
library(gplots)
library(ggrepel)

#reading files##########
primary= read.table('cancer_cell_non-met_scRNA_HNSCC.txt',
                    header=T, row.names = 1, sep= '\t')
metastasis= read.table('cancer_cell_met_scRNA_HNSCC.txt', 
                       header = T, row.names= 1, sep= '\t')

#creating lables####
m= rep('metastasis', dim(metastasis)[2])
p= rep('primary', dim(primary)[2])
cell_type= c(p,m)

#transposing the variables ####
prima1= t(primary)
meta1= t(metastasis)

#attaching the data #######
mydata= rbind(prima1, meta1)
mydata.labled= data.frame(cell_type, mydata, stringsAsFactors = T)

#trimming genes with zero values for all columns
zero.mean= which(apply(mydata, 2, mean)==0)
mydata.a= mydata[,-zero.mean] 

prima1.nz= prima1[,-zero.mean]
meta1.nz= meta1[,-zero.mean]

#running for loop for doing wilcox test for each gene 
dim(prima1.nz)
b= ncol(meta1.nz)#b stores number of columns for *.nz data  

#running a for loop for wilcox test####
genes= c()
for (i in 1:b){
  genes[i]= wilcox.test(prima1.nz[,i], meta1.nz[,i],
                        paired=F, alternative= 'two.sided')$p.value
}
length(genes)
g.names= colnames(prima1.nz)
g.list= data.frame( g.names, genes, stringsAsFactors = T)
head(g.list)
colnames(g.list)= c('gene_id', 'pval')

#creating log2fold change list, primary tumor cells control
lgf= c()
for (i in 1:b){
  lgf[i]= log2(mean(meta1.nz[,i])) - log2(mean(prima1.nz[,i]))
}

#A complete list with gene id, pvalue, logfold change####
g.list= data.frame( g.list, lgf,stringsAsFactors = T)
colnames(g.list)= c('gene_id', 'pval', 'log2_fold')
class(g.list)
head(g.list)

#multiple testing corection#######
g.list$pval= p.adjust(g.list$pval)

#crearing log2 fold gene list with FDR values####
which(g.list$log2_fold==1)
inf_pos= (1/0)
inf_neg= -(1/0)
up_in_primary= subset(g.list,(g.list[,3] > 2) 
                      & (g.list[,3]<inf_pos)
                      & (g.list[,2] <0.05))
down_in_primary= subset(g.list,(g.list[,3] < (-2))
                        & (g.list[,3]> inf_neg)
                        & (g.list[,2] <0.05))

#subsetting sorted gene list from whole data#######
#all the genes with zero values in samples from both groups
#are removed
mydata.nz= rbind(prima1.nz, meta1.nz)
my_gene_list= rbind(up_in_primary, down_in_primary)

#saving the gene list
write.table(my_gene_list, 'genes_fold_chng.txt', sep= '\t')

#subsetting####
pq= my_gene_list$gene_id
pq=as.character(pq)
x=mydata.nz[,(colnames(mydata.nz) %in% pq)]
dim(mydata.nz)
dim(x)
m1= rep('metastasis', dim(meta1.nz)[1])
p1= rep('primary', dim(prima1.nz)[1])
cell_type= c(p1,m1)
sorted_data= data.frame(cell_type, x, stringsAsFactors = T)

#PCA####
pc= prcomp(sorted_data[,-1],
           scale.= T,
           center = T)

g <- ggbiplot(pc,
              choices = 1:2,
                  obs.scale = 1,
                  var.scale = 1,
                  groups= sorted_data$cell_type,
                  var.axes = F,
                  )

print(g)

#volcano plotting using g.list####
g.list$difexpressed <- 'NO'
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
g.list$diffexpressed[g.list$log2_fold > 0.6 & g.list$pval < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
g.list$diffexpressed[g.list$log2_fold < -0.6 & g.list$pval < 0.05] <- "DOWN"
g.list= g.list[,-4]
# Re-plot but this time color the points with "diffexpressed"
plt <- ggplot(data=g.list, aes(x=log2_fold, y=-log10(pval), col=diffexpressed)) + geom_point() + theme_minimal()
print(plt)
