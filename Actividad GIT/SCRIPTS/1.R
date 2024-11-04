#Secuencias 
library("ggmsa")
library(Biostrings)
library(BiocManager)
library(msa)
library(ape)
library(phangorn)
library(phytools)
library(seqinr)

globulinas<-readDNAStringSet("Data/DivergentGlobins.fasta")
globulinas
insulinas<-readDNAStringSet("Data/Insulinas.fasta")
insulinas

alignments <- msa(globulinas, method = "ClustalW")
alignments1 <- msa(insulinas, method = "Muscle")
alignments1

ggmsa(globulinas, 0, 10, front = "DroidSansMono", color = "Zappo_AA", char_width= 1, seq_name = TRUE )
ggmsa(insulinas, 0, 10, font = "DroidSansMono", color = "Zappo_AA", char_width= 1, seq_name = TRUE )

#GLOBULINAS
pdf("GRAFICOS/GLOBULINAS.pdf")
arbol<- msaConvert(alignments, type = "seqinr::alignment")
matriz<-dist.alignment(arbol, "identity")
matriz
as.matrix(matriz)[5:116, "Homo_sapiens", drop=FALSE]

tree<-nj(matriz)
plot(tree)
dev.off()

#INSULINAS
pdf("GRAFICOS/INSULINAS.pdf")
arbol1<- msaConvert(alignments1, type = "seqinr::alignment")
matriz1<-dist.alignment(arbol1, "identity")
matriz1
as.matrix(matriz1)[14:72, "Insulinas", drop=FALSE]

tree1<-nj(matriz1)
plot(tree1)
dev.off()
