

#setwd("~/Documents/EPHE_nm/Analyses_TrKOSpg11_2021")
setwd("/Users/daniel/Documents/EPHE_nm/Analyses_RNAseq_Gio/Synthese")
#load("geneinfo4.rda")
dir()
#load("out1_table_FDR_25421x14.rda")
#load("out2_table_logFC_25421x14.rda")
load("gout.rda")
contrast_names =c("Cer_KO_06_vs_Cer_WT_06", "Cer_KO_16_vs_Cer_WT_16", "Cer_KO_32_vs_Cer_WT_32",
                  "Cor_KO_06_vs_Cor_WT_06", "Cor_KO_16_vs_Cor_WT_16", "Cor_KO_32_vs_Cor_WT_32",
                  "Hip_KO_06_vs_Hip_WT_06", "Hip_KO_16_vs_Hip_WT_16", "Hip_KO_32_vs_Hip_WT_32",
                  "Spi_KO_06_vs_Spi_WT_06")
# Creation des fichiers de sortie de l'analyse differentielle
table_FDR = gout[[1]]$table[,c("Chr","gene_name","ensembl_gene_id","entrezgene_id")] # starts by taking the columns with info genes
for(i in 1:10) table_FDR=cbind(table_FDR,gout[[i]]$table[,"FDR"]) #add "logFC" column for each element of outi_list
colnames(table_FDR) = c("Chr","gene_name","ensembl_gene_id","entrezID",paste0("FDR",names(contrast_names)))
paste0("FDR",names(contrast_names))
table_logFC = gout[[1]]$table[,c("Chr","gene_name","ensembl_gene_id","entrezgene_id")] # starts by taking the columns with info genes
for(i in 1:10) table_logFC=cbind(table_logFC,gout[[i]]$table[,"logFC"]) #add "logFC" column for each element of outi_list
colnames(table_logFC) = c("Chr","gene_name","ensembl_gene_id","entrezID",paste0("logFC",names(contrast_names)))

out1 = table_FDR
out2 = table_logFC

# colonne calculant la plus petite FDR
out1$FDRmin=apply(out1[,5:14],1,min)
sum(out1$FDRmin < 0.001)
o = order(out1$FDRmin,decreasing = FALSE)
out1 = out1[o,]
out2 = out2[o,]

row_num = 1000
neworder = rev(1:row_num)
mat12 = cbind(as.matrix(out1[neworder,5:13]),as.matrix(out2[neworder,5:13]))
ma = match(rownames(mat12),rownames(out1))
rownames(mat12) = out1$gene_name[ma]
mat12_norm=scale(mat12)

require(Rtsne)
set.seed(5)
tsne <- Rtsne(mat12_norm, check_duplicates = FALSE, pca = FALSE, perplexity=30, theta=0.1, dims=2)
dat = data.frame(tsne1 = tsne$Y[,1],
                 tsne2 = tsne$Y[,2],
                 genename = rownames(mat12_norm),
                 color = rgb(0,0,0,0.5))
library(ggplot2)
library(ggrepel)

mycolors = ifelse(out1[neworder,"Chr"]=="chr2",rgb(1,0,0,0.5),rgb(0,0,0,0.2))


p <- ggplot(dat, aes(x=tsne1, y=tsne2, label = genename)) +
  xlim(c(-50,40))+
  ylim(c(-40,40))+
  geom_point(col = rgb(0,0,0,0.5))+
  theme_bw()
p

pdf(file = paste0("tsne_all_v3.pdf"), width = 20, height = 20)
print(p)
dev.off()

p2 <- p + geom_text_repel(col = rgb(0,0,0,0.5), max.overlaps=1000)
pdf(file = paste0("tsne_all_names_v3.pdf"), width = 20, height = 20)
print(p2)
dev.off()

plot_tsne=function(the_comparison = 1,
                   threshold = 1e-4,
                   minvalue = -3,
                   maxvalue = 3){
  
  title = paste(colnames(mat12)[the_comparison],"<",threshold,"\n",
                "minvalue in red = ",minvalue," and maxvalue in green =", maxvalue)
  
  myvalues = mat12[,the_comparison+9]
  myvalues_norm = (myvalues - minvalue)/(maxvalue-minvalue)
  
  myvalues_norm[myvalues_norm<0]=0
  myvalues_norm[myvalues_norm>1]=1
  
  rgb_table = colorRamp(c("blue","red"))(myvalues_norm)
  mycolors = rgb(red= rgb_table[,1]/255, green = rgb_table[,2]/255, blue = rgb_table[,3]/255,alpha=1)
  
  mycolors[mat12[,the_comparison]>threshold]=rgb(0.7,0.7,0.7,0.1)
  
  p <- ggplot(dat, aes(x=tsne1, y=tsne2, label = genename)) +
    xlim(c(-50,40))+
    ylim(c(-40,40))+
    geom_point(col = mycolors)+
    ggtitle(title) +
    theme_bw()
  p2 <- p + geom_text_repel(col = mycolors, max.overlaps=1000)
  ?geom_text_repel

  #png(file = paste0("tsne_",colnames(mat12)[the_comparison],"_v3.png"), width = 600, height = 600)
  pdf(file = paste0("tsne_",colnames(mat12)[the_comparison],"_v3.pdf"), width = 20, height = 20)
  print(p2)
  dev.off()
  print(paste0("tsne_",colnames(mat12)[the_comparison]))
}


for(i in 1:9) plot_tsne(the_comparison = i, threshold = 1e-4,  minvalue = -2, maxvalue = 2)


