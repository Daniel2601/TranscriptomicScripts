
# generation of a more complete geneinfo
# Addition of 5 columns for the geninfo table:
# the 8th column gives the Ensemble name without the dot and sub number
# the 9th column gives a short Wiki description of the gene
# the 10th column gives the Wiki name of the gene
# the 11th column gives the entrez gene identifier
load("geneinfo.rda")
#BiocManager::install("biomaRt")
library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
genecleaned = unlist(lapply(geneinfo$Geneid,  sub, pattern = "\\.\\d+$", replacement = ""))
full_names <- getBM(attributes=c('ensembl_gene_id','wikigene_description','wikigene_name',
                                 'entrezgene_id'),
                    filters = 'ensembl_gene_id', values = genecleaned, mart = ensembl)
geneinfo$ensembl_no_version = unlist(lapply(geneinfo$Geneid,  sub, pattern = "\\.\\d+$", replacement = ""))
ma = match(geneinfo$ensembl_no_version, full_names$ensembl_gene_id)
geneinfo = cbind(geneinfo, full_names[ma,])
geneinfo$Chr = unlist(lapply(strsplit(geneinfo$Chr,";"),function(x)x[1]))
geneinfo$Start = unlist(lapply(strsplit(geneinfo$Start,";"),function(x)x[1]))
geneinfo$End = unlist(lapply(strsplit(geneinfo$End,";"),function(x)x[1]))
geneinfo$Strand = unlist(lapply(strsplit(geneinfo$Strand,";"),function(x)x[1]))
geneinfo = subset(geneinfo, select=-ensembl_no_version)

rm(full_names)
rm(ma)
rm(ensembl)
rm(genecleaned)
save(geneinfo,file = "geneinfo4.rda")
