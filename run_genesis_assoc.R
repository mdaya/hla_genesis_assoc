#Load required libraries
library("GENESIS")
library("GWASTools")
library("SNPRelate")

#Set parameters
args <- commandArgs(trailingOnly = TRUE)
null.rdata.file <- args[1]
bed.fn <- args[2]
bim.fn <- args[3]
fam.fn <- args[4]
out.file.name <- args[5]

#Load null model
load(null.rdata.file)

#Convert PLINK file to GDS
gdsfile <- tempfile()
snpgdsBED2GDS(bed.fn=bed.fn, 
              bim.fn=bim.fn,
              fam.fn=fam.fn,
              gdsfile)
gds <- GdsGenotypeReader(gdsfile)
genoData <- GenotypeData(gds)
iterator <- GenotypeBlockIterator(genoData)

#Run association model
assoc <- assocTestSingle(gdsobj=iterator, null.model=nullmod, test="Score")

#Write association output
write.table(assoc,file=out.file.name, sep="\t",row.names=F,col.names=T,quote=F)


