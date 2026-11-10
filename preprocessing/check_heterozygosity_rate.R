# When esxecuted by plink2

het <- read.table("R_check.het", header=TRUE, comment.char = "")
print(colnames(het))

het$HET_RATE = (het$OBS_CT - het$O.HOM.) / het$OBS_CT

pdf("heterozygosity.pdf")
hist(het$HET_RATE,
     xlab="Heterozygosity Rate",
     ylab="Frequency",
     main="Heterozygosity Rate")
dev.off()

