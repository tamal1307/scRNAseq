#setting working directory and loading data
setwd("A:/try data/my analysis/")
memory.limit(size=30000)
load("A:/try data/my analysis/tsne_env.RData")
library(Seurat)
library(ReactomeGSA)
#save.image("A:/try data/my analysis/tsne_env.RData")

#installing t-sne package####
#if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
#BiocManager::install("ReactomeGSA", force = T)

#library(Rtsne)

#running t-sne ####
my.tsne= Rtsne(sorted_data[,-1], check_duplicates= F)
plot(my.tsne$Y,col=sorted_data$cell_type, asp=1, pch=20,
     xlab= 'Dimension1', ylab = 'Dimension2', main= 't-sne plot')

#plotting heatmap####
library(dplyr)
library(ggplot2)
library(reshape)

#melting the data using reshape package
my.dat=data.frame(rownames(sorted_data), sorted_data[,-1],
                  stringsAsFactors = T)   
colnames(my.dat)
g3<-melt(my.dat)
colnames(g3)=c('cell_id', 'variable', 'value')
colnames(g3)
#Initial Heatmap 
plot1<-ggplot(g3, aes(variable, cell_id, fill= value)) + 
  geom_tile()

#plot1
plot2= plot1+ scale_fill_distiller(palette = 'Accent')
plot2

a= t(sorted_data[,-1])
head(a)
which(rownames(a)== 'cell_type')

#working with seurat####

x= CreateSeuratObject(a)
x1= FindVariableFeatures(x)
x1= ScaleData(x1)
x1= RunPCA(x1)
x1= FindNeighbors(x1)
x1= FindClusters(x1)
x1= RunTSNE(x1, check_duplicates=F)
#x1= RunUMAP(x1, dims = 1:5) #dims can be changed upto max PC
#x1-> x2

#tsne and PCA plotting####
DimPlot(x1, reduction= 'tsne', 
        label = T, pt.size = 0.75) + NoLegend()
#DimPlot(x1, reduction= 'pca')
#DimPlot(x1, reduction= 'umap')

#rename cluster ids####
#new.cluster.ids <- c('clust0','clust1','clust2', 'clust3',
                     'clust4', 'clust4', 'clust5', 'clust6',
                     'clust7','clust8', 'clust9', 'clust10',
                     'clust11', 'clust12', 'clust13',
                     'clust14', 'clust15', 'clust16',
                     'clust17')
#names(new.cluster.ids) <- levels(x1)
#x1 <- RenameIdents(x1, new.cluster.ids)
#levels(x1)
#x1-> x3
#creat heatmap###

#DimHeatmap(x1, dims=1:2, balanced=T) 

#dims-> corresponds to principal component number

head(Embeddings(x1, reduction = 'pca')[,1:2])

#creating cells types####
#Idents(x1, cells= 1:3114) <- 'primary'
#Idents(x1, 3115: 3902) <- 'metastasis'

head(Idents(x1))

#Cluster biomarkers####
clst0.markr= FindMarkers(x1, ident.1 = 0,
                              min.pct = 0.25 )
clst1.markr= FindMarkers(x1, ident.1 = 1,
                         min.pct = 0.25 )
clst2.markr= FindMarkers(x1, ident.1 = 2,
                         min.pct = 0.25 )
clst3.markr= FindMarkers(x1, ident.1 = 3,
                         min.pct = 0.25 )
clst4.markr= FindMarkers(x1, ident.1 = 4,
                         min.pct = 0.25 )
clst5.markr= FindMarkers(x1, ident.1 = 5,
                         min.pct = 0.25 )
clst6.markr= FindMarkers(x1, ident.1 = 6,
                         min.pct = 0.25 )
clst7.markr= FindMarkers(x1, ident.1 = 7,
                         min.pct = 0.25 )
clst8.markr= FindMarkers(x1, ident.1 = 8,
                         min.pct = 0.25 )
clst9.markr= FindMarkers(x1, ident.1 = 9,
                         min.pct = 0.25 )
clst10.markr= FindMarkers(x1, ident.1 = 10,
                         min.pct = 0.25 )
clst11.markr= FindMarkers(x1, ident.1 = 11,
                         min.pct = 0.25 )
clst12.markr= FindMarkers(x1, ident.1 = 12,
                         min.pct = 0.25 )
clst13.markr= FindMarkers(x1, ident.1 = 13,
                         min.pct = 0.25 )

#CLUSTER SPECIFIC GENE EXPRESSOIN

VlnPlot(x1, features = "SOX2", sort = T) + NoLegend()
FeaturePlot(x1, 'SOX2')

#geneset enrichment using ReactomeGSA####
#https://bioconductor.org/packages/release/bioc/manuals/ReactomeGSA/man/ReactomeGSA.pdf
file.gsa= analyse_sc_clusters(
  x1,
  use_interactors = TRUE,
  include_disease_pathways = T,
  create_reactome_visualization = T,
  create_reports = FALSE,
  report_email = NULL,
  verbose = FALSE,
)

#file.gsa
#open_reactome(file.gsa)
my.paths=pathways(file.gsa)

#heatmap of pathways
plot_gsva_heatmap(
  file.gsa,
  pathway_ids = NULL,
  max_pathways = 50,
  truncate_names = TRUE)

#generating plot for individual pathway####
plot_gsva_pathway(file.gsa, 'R-HSA-2892245')


