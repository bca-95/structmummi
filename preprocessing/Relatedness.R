pdf("relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
par(pch=16, cex=1)
with(relatedness,plot(Z0,Z1, xlim=c(0,1), ylim=c(0,1), type="n"))
with(subset(relatedness,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness,RT=="UN") , points(Z0,Z1,col=3))
with(subset(relatedness,RT=="OT") , points(Z0,Z1,col=5))
#legend(1,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))
if (!is.null(relatedness$RT) && length(unique(relatedness$RT)) > 0) {
  legend(1, 1, xjust = 1, yjust = 1, legend = unique(relatedness$RT))
}
dev.off()

pdf("zoom_relatedness.pdf")
relatedness_zoom = read.table("zoom_pihat.genome", header=T)
par(pch=16, cex=1)
with(relatedness_zoom,plot(Z0,Z1, xlim=c(0,0.02), ylim=c(0.98,1), type="n"))
with(subset(relatedness_zoom,RT=="PO") , points(Z0,Z1,col=4))
with(subset(relatedness_zoom,RT=="UN") , points(Z0,Z1,col=3))
with(subset(relatedness_zoom,RT=="OT") , points(Z0,Z1,col=5))
#legend(0.02,1, xjust=1, yjust=1, legend=levels(relatedness$RT), pch=16, col=c(4,3))
if (!is.null(relatedness_zoom$RT) && length(unique(relatedness_zoom$RT)) > 0) {
  legend(0.02, 1, xjust=1, yjust=1, legend=unique(relatedness_zoom$RT), pch=16, col=c(4,3))
}
dev.off()

pdf("hist_relatedness.pdf")
relatedness = read.table("pihat_min0.2.genome", header=T)
hist(relatedness[,10],main="Histogram relatedness", xlab= "Pihat")  
dev.off() 

