dMMR = c("MSH2","MSH6","MLH1","PMS2")
DDR = c("FANCC", "FANCD2", "PALB2", "BRCA1", "POLE", "RB1", "TP53", "ATM", "ERCC2")
list_g = c(dMMR, DDR)

library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id", values = list_g,
                 mart = mart)
