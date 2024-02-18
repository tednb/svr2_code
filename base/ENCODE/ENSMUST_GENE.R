load("Diestrus_ENSMUST.rda")
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(biomaRt)
tranpt$TID <- sub("\\..*$", "", tranpt$TID)
ensembl <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl",mirror = "useast")
gene_info <- getBM(attributes=c('ensembl_transcript_id', 'mgi_symbol'),
                   filters='ensembl_transcript_id', 
                   values=tranpt$TID, 
                   mart=ensembl)

df_with_gene_names <- merge(tranpt, gene_info, by.x="TID", by.y='ensembl_transcript_id')
idx<- which(df_with_gene_names$mgi_symbol=="")
df_with_gene_names <- df_with_gene_names[-idx,]
gene_tpm <- aggregate(tpm ~ mgi_symbol, data=df_with_gene_names, FUN=sum)
uterus_RNA <- data.table(tpm = log2(gene_tpm$tpm+1),symbol = gene_tpm$mgi_symbol)
