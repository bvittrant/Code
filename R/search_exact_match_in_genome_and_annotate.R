
###############################################################################

### IMPORTANT A LIRE

# Faire attention aux version de genome. SI on commence avec hg19 il faut 
# toujours garder hg19 sinon gros problèmes

###############################################################################

# DL package and library

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

# BSgenome
biocLite("BSgenome")
library(BSgenome.Hsapiens.UCSC.hg19)

#citation https://bioconductor.org/packages/release/bioc/html/BSgenome.html

# CHIPseeker
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(clusterProfiler)
library(org.Hs.eg.db)

#citation https://www.bioconductor.org/packages/release/bioc/html/ChIPpeakAnno.html

###############################################################################


# On met le genome humain dans un objet
genome = BSgenome.Hsapiens.UCSC.hg19
seqnames = seqnames(genome)

# On crée la liste avec les séquence
seq = c("TAGTGCTAAGCTGGGA","TCCTGAGGTCTAACCT","TGCGCAGGGAGGCGCC") # 16 bp
seq = as.character(seq)

# On crée une data frame vide pour recevoir les résultats
res = as.data.frame(matrix(ncol=15))
colnames(res) = c("seqnames", "start","end","width",
	"strand","annotation","geneChr","geneStart",
	"geneEnd","geneLength","geneStrand","geneId",
	"transcriptId","distanceToTSS",'seq')

# On lance les boucles pour chaque sequence et chaque Chromosomes
for(i in seq){
	print(i)
	for(seqname in seqnames){
		print(seqname)
		match = matchPattern(i, genome[[seqname]])
		if(dim(as.data.frame(match, row.names = NULL, optional = FALSE))[1]==0){next}
		print('matched')
		tmp2 = GRanges(seqname, IRanges(start(match), end(match)))

		# Ici on choisit tssRegion=c(-10,10) donc on assigne l'objet le plus proche 
		# a possible -10 ou +10 base de la sequence
		peakAnno <- annotatePeak(tmp2,tssRegion=c(-10,10),TxDb = txdb, level = "transcript",
                           genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"))
		df_tmp = as.data.frame(peakAnno, row.names = NULL, optional = FALSE)
		
		# On rajoute le nom de la sequence a l'arrache ... 
		df_tmp$seq = i

		# On met tout ça dans notre DF
		res = rbind(res, df_tmp)
	}
}

# On enlève la ligne de NA
res = res[-1,]
res

write.table(res, file='res.tsv', row.names=F, quote=F,sep='\t')




