
load("table_FDR.rda")
load("table_logFC.rda")

# additional column to calculate the minimum FDR
table_FDR$FDRmin=apply(table_FDR[,5:14],1,min)
sum(table_FDR$FDRmin < 0.001)
o = order(table_FDR$FDRmin,decreasing = FALSE)
table_FDR = table_FDR[o,]
table_logFC = table_logFC[o,]

row_num = 1000
neworder = rev(1:row_num)
mat12 = cbind(as.matrix(table_FDR[neworder,5:13]),as.matrix(table_logFC[neworder,5:13]))
ma = match(rownames(mat12),rownames(table_FDR))
rownames(mat12) = table_FDR$gene_name[ma]
mat12_norm=scale(mat12)

library(Rtsne)
set.seed(5)

tsne <- Rtsne(mat12_norm, check_duplicates = FALSE, pca = FALSE, perplexity=30, theta=0.1, dims=2)
dat = data.frame(tsne1 = tsne$Y[,1],
                 tsne2 = tsne$Y[,2],
                 genename = rownames(mat12_norm),
                 color = rgb(0,0,0,0.5))

library(ggplot2)
library(ggrepel)

p <- ggplot(dat, aes(x=tsne1, y=tsne2, label = genename)) +
  xlim(c(-50,40))+
  ylim(c(-40,40))+
  geom_point(col = rgb(0,0,0,0.5))+
  theme_bw()
p

pdf(file = paste0("tsne_all.pdf"), width = 20, height = 20)
print(p)
dev.off()

p2 <- p + geom_text_repel(col = rgb(0,0,0,0.5), max.overlaps=1000)
pdf(file = paste0("tsne_all_names.pdf"), width = 20, height = 20)

print(p2)
dev.off()
p2


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


#for(i in 1:9) plot_tsne(the_comparison = i, threshold = 1e-4,  minvalue = -2, maxvalue = 2)


