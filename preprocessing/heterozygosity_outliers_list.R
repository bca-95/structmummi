# plink2
het <- read.table("R_check.het", header=TRUE, comment.char = "")
print(colnames(het))

het$HET_RATE = (het$OBS_CT - het$O.HOM.) / het$OBS_CT

het_fail = subset(
  het,
  (HET_RATE < mean(HET_RATE) - 3*sd(HET_RATE)) |
  (HET_RATE > mean(HET_RATE) + 3*sd(HET_RATE))
)

het_fail$HET_DST = (het_fail$HET_RATE - mean(het$HET_RATE)) / sd(het$HET_RATE)

write.table(het_fail, "fail-het-qc.txt", row.names=FALSE)
