## When plink1.9 is used to generate the freq file
#Sys.setenv(TMPDIR="/wibicomfs/STBS/carene/tmp")
#
#maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)
#pdf("MAF_distribution.pdf")
#hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
#dev.off()

# When plink1.9 is used to generate the freq file

maf_freq <- read.table("MAF_check.afreq", header =TRUE, as.is=T)
pdf("MAF_distribution.pdf")
hist(maf_freq[,6],main = "MAF distribution", xlab = "MAF")
dev.off()


