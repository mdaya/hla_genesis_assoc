library(stringr)

#Set arguments
args <- commandArgs(trailingOnly = TRUE)
plink.prefix <- args[1]
hibag.files <- args[-1]

#Merge all the HIBAG files
hibag.frame <- data.frame()
for (file in hibag.files) {
  gene.frame <- read.delim(file, stringsAsFactors = F)
  gene.frame <- gene.frame[gene.frame$prob >= 0.5,]
  gene.start.pos <- str_locate(file, "_HLA_")[,2] + 1
  gene.end.pos <- str_locate(file, "_value")[,1] - 1
  gene <- substring(file, gene.start.pos, gene.end.pos)
  gene.frame$allele1 <- paste0(gene, "*", gene.frame$allele1)
  gene.frame$allele2 <- paste0(gene, "*", gene.frame$allele2)
  gene.frame$gene <- gene
  hibag.frame <- rbind(hibag.frame, gene.frame)
}

#Get list of sample IDs
sample.ids <- unique(hibag.frame$sample.id)

#Create initial ped frame with all values set to missing
called.alleles <- unique(c(hibag.frame$allele1, hibag.frame$allele2))
called.alleles <- called.alleles[order(called.alleles)]
ped.frame <- data.frame(FID=sample.ids, IID=sample.ids,
                        FATHER=rep(0, length(sample.ids)), MOTHER=rep(0, length(sample.ids)),
                        SEX = rep(-9, length(sample.ids)), PHENO = rep(-9, length(sample.ids)))
for (called.allele in called.alleles) {
  ped.frame$ALLELE <- "0 0"
  names(ped.frame)[dim(ped.frame)[2]] <- called.allele
}

#For each sample and possible allele, set value to "1 1" for non-missing non-carriers, and "1 2" for carriers
for (sample.id in sample.ids) {
  sample.allele.frame <- hibag.frame[hibag.frame == sample.id,]
  for (called.allele in called.alleles) {
    gene <- unlist(str_split(called.allele, "\\*"))[1]
    sample.allele.row <- sample.allele.frame[(sample.allele.frame$gene == gene),]
    if (nrow(sample.allele.row) == 1) {
      if ((called.allele %in% sample.allele.row$allele1) | (called.allele %in% sample.allele.row$allele2)) {
        ped.frame[ped.frame$IID == sample.id,names(ped.frame) == called.allele] <- "1 2"
      } else {
        ped.frame[ped.frame$IID == sample.id,names(ped.frame) == called.allele] <- "1 1"
      } 
    }
  }
}

#Create map file
nr.alleles <- length(called.alleles)
map.frame <- data.frame( CHR=rep(6, nr.alleles),
                         SNP=called.alleles,
                         CM=rep(0, nr.alleles),
                         BP=1:nr.alleles)

#Write the output
write.table(ped.frame, paste0(plink.prefix, ".ped"),  sep=" ", quote=F, row.names=F, col.names=F)
write.table(map.frame, paste0(plink.prefix, ".map"),  sep=" ", quote=F, row.names=F, col.names=F)