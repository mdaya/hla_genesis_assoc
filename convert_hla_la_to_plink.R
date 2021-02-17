###################
# Setup parameters
###################
#Note: it is assumed that the HLA-LAoutput files have the formaat <sample-id>_*.txt; 
#the sample-id is extracted as the first "_" seperated filed in the filename
args <- commandArgs(trailingOnly = TRUE)
hla.la.out.dir <- paste0(args[1], "/")
min.Q1 <- args[2]
min.avg.cov <- args[3]
incl.genes.str <- args[4]
incl.genes <- unlist(strsplit(gsub(" ", "", incl.genes.str), ","))
out.prefix <- args[5]

############
# Functions
############

#Function to extract first two alleles
extractAntigenAllele <- function(a) {
  if (is.na(a)) {
    return (NA)
  }
  colon.pos <- unlist(gregexpr(":", a))
  if (length(colon.pos) == 1) {
    if (nchar(a) == 5) {
      return (a)
    } else {
      return(a)
    }
  } else {
    return (substring(a, 1, colon.pos[2]-1))
  }
}

#Function to do a first parse of HLA allele output
parseFile <- function(hla.la.out.file) {
  sample <- unlist(strsplit(hla.la.out.file, "_"))[1]
  print(paste("Parsing", sample))
  calls.out.file <- paste0(out.prefix, "_HLA_LA.txt")
  hla.la.frame <- read.delim(paste0(hla.la.out.dir, hla.la.out.file), stringsAsFactors = F, head=T)
  hla.la.frame <- hla.la.frame[hla.la.frame$Locus %in% incl.genes,]
  hla.la.frame$Sample <- sample
  failed.call.rows <- (hla.la.frame$Q1 < min.Q1) | (hla.la.frame$perfectG != 1) | (hla.la.frame$AverageCoverage < min.Q1)
  hla.la.frame$Allele[failed.call.rows] <- NA
  hla.la.frame$Allele <- sapply(hla.la.frame$Allele, extractAntigenAllele)
  hla.la.frame <- hla.la.frame[,c("Sample", "Locus", "Chromosome", "Allele", "Q1","AverageCoverage")]
  write.table(hla.la.frame, file=calls.out.file,  
              sep="\t", quote=F, row.names=F, 
              col.names=!file.exists(calls.out.file), append=file.exists(calls.out.file))

}

##############
# Run parsing
##############

for (file in list.files(hla.la.out.dir)) {
  parseFile(file)
}

#########################
# Convert to PLINK format
#########################

#Read in the parsed file
assoc.calls <- read.delim(paste0(out.prefix, "_HLA_LA.txt"))

#Create initial ped frame with all values set to missing
called.alleles <- unique(assoc.calls$Allele[!is.na(assoc.calls$Allele)])
called.alleles <- called.alleles[order(called.alleles)] 
sample.ids <- unique(assoc.calls$Sample)
ped.frame <- data.frame(FID=sample.ids, IID=sample.ids, 
                        FATHER=rep(0, length(sample.ids)), MOTHER=rep(0, length(sample.ids)),
                        SEX = rep(-9, length(sample.ids)), PHENO = rep(-9, length(sample.ids)))
for (called.allele in called.alleles) {
  ped.frame$ALLELE <- "0 0"
  names(ped.frame)[dim(ped.frame)[2]] <- called.allele
}

#For each sample and possible allele, set value to "1 1" for non-missing non-carriers, and "1 2" for carriers
for (sample.id in sample.ids) {
  sample.allele.frame <- assoc.calls[assoc.calls$Sample == sample.id,]
  for (called.allele in called.alleles) {
    locus <- unlist(strsplit(called.allele, "\\*"))[1]
    sample.allele.row <- sample.allele.frame[sample.allele.frame$Locus == locus,]
    if (sum(is.na(sample.allele.row$Allele)) == 0) {
      if (called.allele %in% sample.allele.row$Allele) {
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
write.table(ped.frame, paste0("ped_", out.prefix, ".ped"),  sep=" ", quote=F, row.names=F, col.names=F)
write.table(map.frame, paste0("ped_", out.prefix, ".map"),  sep=" ", quote=F, row.names=F, col.names=F)



