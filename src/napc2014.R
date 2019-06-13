require(ape)
require(geiger)

source("~/Dropbox/code/R/common_src/phy_dateTree.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R")
source("~/Dropbox/code/R/artioLimbs/src/ratesFromDiscreteCharacters.R")
source("~/Dropbox/code/R/common_src/ratesOverIntervals.R")

oneRep<-function(rep, tree, dates, intervals, thisVec, thisBreaks, do.parallel=FALSE) {
	# thisVec <- array(thisMat[,chr], dimnames=list(rownames(thisMat)))
	thisTree <- ladderize(multi2di(tree), right=FALSE)

	thisVec <- thisVec[match(thisTree$tip.label,names(thisVec))]
	thisVec <- thisVec[is.finite(thisVec)]
	thisTree<-drop.tip(thisTree, thisTree$tip.label[!(thisTree$tip.label %in% names(thisVec))])

	dates <- dates[match(thisTree$tip.label, row.names(dates)),]   #reorders dat to match the taxon order of thisTree
	thisTree <- dateTreeWithRanges(thisTree, dates, within.error=TRUE, random=TRUE, extra=0.5, rootAdj=20.0, include.ranges=TRUE)
	# thisTree<-dateTreeWithList(thisTree, dates, offset=0.05, extendNodes=TRUE, rootAdj=20.0, include.ranges=FALSE)
	# plot(thisTree, cex=0.3, no.margin=TRUE)

	testIntervalRates(thisTree, x= thisVec, intervals, do.parallel=do.parallel)

}

# quartz(width=10, height=6)
# plot(tree, edge.color="white", tip.color="white", edge.width=1.5, cex=0.3, no.margin=TRUE, label.offset=1)
do.parallel<-FALSE

tree<-read.nexus("~/Dropbox/code/R/amandaTeeth/dat/ungulates.tre")
dates<-read.csv("~/Dropbox/code/R/amandaTeeth/dat/ranges_20130918_total.csv", row.names=1)
rownames(dates)<-gsub(" ", "_", row.names(dates))

tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% row.names(dates))])
dates<-dates[(row.names(dates) %in% tree$tip.label),]	#cull taxa in dates list that aren't in the tree
rownames(thisMat)<-gsub(" ", "_", rownames(thisMat))
thisMat<-thisMat[rownames(thisMat) %in% tree$tip.label,]

chr<-"bodyMass"
thisBreaks <- bmBreaks
reps=10

# rez <- oneRep(1, tree=tree, intervals=intervals, dates=dates, thisVec=array(thisMat[,chr], dimnames=list(rownames(thisMat))), thisBreaks=thisBreaks, do.parallel=TRUE)
	thisVec <- array(thisMat[,chr], dimnames=list(rownames(thisMat)))
	thisTree <- ladderize(multi2di(tree), right=FALSE)

	thisVec <- thisVec[match(thisTree$tip.label,names(thisVec))]
	thisVec <- thisVec[is.finite(thisVec)]
	thisTree<-drop.tip(thisTree, thisTree$tip.label[!(thisTree$tip.label %in% names(thisVec))])

	dates <- dates[match(thisTree$tip.label, row.names(dates)),]   #reorders dat to match the taxon order of thisTree
	thisTree <- dateTreeWithRanges(thisTree, dates, within.error=TRUE, random=TRUE, extra=0.5, rootAdj=20.0, include.ranges=TRUE)
	# thisTree<-dateTreeWithList(thisTree, dates, offset=0.05, extendNodes=TRUE, rootAdj=20.0, include.ranges=FALSE)
	# plot(thisTree, cex=0.3, no.margin=TRUE)

	rez <- testIntervalRates(thisTree, x= thisVec, intervals, do.parallel=TRUE)


# repList<-list()
# if (do.parallel) {
	# require(parallel)
	# repList<-mclapply(seq_len(reps), oneRep, tree=tree, dates=dates, intervals=intervals, thisVec=array(thisMat[,chr], dimnames=list(rownames(thisMat))), thisBreaks=thisBreaks, mc.cores=detectCores()-2)
# } else repList<-lapply(seq_len(reps), oneRep, tree=tree, intervals=intervals, dates=dates, thisVec=array(thisMat[,chr], dimnames=list(rownames(thisMat))), thisBreaks=thisBreaks)

# # # save(repList, file = "~/Desktop/replist.R")
# repList<-load("~/Desktop/replist.R")

# histMat <- sapply(repList, getHistFromOneRepList, simplify="array")

# # bigList <- repList[[3]]
# countList <- lapply(repList, getCountsFromOneRepList, intervals, thisBreaks)
# countSums <- sapply(countList, function(x) apply(x, 3, rowSums), simplify="array")

# # thisCounts<-getCountsFromOneRepList(repList[[1]], intervals, thisBreaks)

# ifonly <- array(NA, dim=c(nrow(intervals)-1, 4, reps))
# for (i in seq_len(dim(countSums)[3])) {
	# for (intv in 1:(nrow(intervals)-1)) {
		# ifonly[intv,,i]<-(countSums[,intv,i]/(sum(histMat[intv+1,,i])+(2*countSums[,intv,i])))[2:5]
	# }
# }
# medIfOnly <- apply(ifonly, c(1,2), quantile, probs=c(0.025, 0.50, 0.975))
# quartz(width=10, height=6)
# plot(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,1], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0,max(medIfOnly)), xlab="Time (Ma)", ylab="Contribution in Isolation", fg="gray75", col.axis="gray75", col.lab="gray75")
# overlayCzTimescale()
# polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medIfOnly[1,,1], medIfOnly[3,(nrow(intervals)-1):1,1]), col=rgb(t(col2rgb("dodgerblue4"))/255, alpha=0.25), border=NA)
# polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medIfOnly[1,,4], medIfOnly[3,(nrow(intervals)-1):1,4]), col=rgb(t(col2rgb("firebrick4"))/255, alpha=0.25), border=NA)
# polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(2*medIfOnly[1,,2], 2*medIfOnly[3,(nrow(intervals)-1):1,2]), col=rgb(t(col2rgb("green4"))/255, alpha=0.25), border=NA)
# lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,1], lwd=4, col="dodgerblue4")
# lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,4], lwd=4, col="firebrick4")
# lines(rowMeans(intervals)[1:(nrow(intervals)-1)], 2*medIfOnly[2,,2], lwd=4, col="green4")
# # polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medIfOnly[1,,3], medIfOnly[3,(nrow(intervals)-1):1,3]), col="goldenrod1", border=NA)
# # lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,3], lwd=3, col="goldenrod4")

# medDist <- apply(apply(histMat, 3, getDistanceBetweenSubsequentIntervals), 1, quantile, probs=c(0.025, 0.5, 0.975))
# quartz(width=10, height=6.5)
# par(mar=c(5,6,4,2))
# plot(rowMeans(intervals)[1:(nrow(intervals)-1)], medDist[2,], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0,max(medDist)), xlab="Time (Ma)", ylab="Change in Body Size Distribution\n(Bray-Curtis Distance)", fg="gray75", col.axis="gray75", col.lab="gray75")
# overlayCzTimescale()
# polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medDist[1,], medDist[3,(nrow(intervals)-1):1]), col=rgb(t(col2rgb("plum4"))/255, alpha=0.25), border=NA)
# lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medDist[2,], lwd=3, col="plum4")







# # # bigList<-repList[[3]]
# # t(apply(countSums, c(1,2), median))

# colSums(this[1:3,])-colSums(this[4:5,])
# thisHist[1:2,]
# thisHist[2,] + colSums(this[2:3,]) - colSums(this[4:5,])


# sum(this[,2:3]) - sum(this[4,5]) 



# # rateMat<-foreach(i=seq_len(reps)) %do% oneRep(i,tree,dates,thisMat)

# # tmo<-transformPhylo.ML(as.matrix(thisVec), thisTree, model="tm1")

# thisMedians<-t(apply(rateMat, 1, quantile, probs=c(0.0275, 0.5, 0.975), na.rm=TRUE))

# thisMax<-max(thisMedians[,2], na.rm=TRUE)
# par(mar=c(5,6,4,1), fg="gray75")
# plot(rowMeans(intervals), thisMedians[,3], xlim=c(max(intervals), min(intervals)), ylim=c(0,thisMax), type="n", xlab="Time (Mybp)", ylab="Rate of Evolution", cex.lab=1.0, col.axis="gray75", col.lab="gray75")
# overlayCzTimescale()
# text(rowMeans(intervals), thisMedians[,2], paste(rowMeans(intervals), "Ma", sep=""), pos=3, col="gray50", cex=0.5)
# polygon(c(rowMeans(intervals)[1:nrow(intervals)], rowMeans(intervals)[nrow(intervals):1]), c(thisMedians[1:nrow(intervals),1], thisMedians[nrow(intervals):1,3]), col=rgb(0.6,0.96,1,0.5), border="cadetblue")
# lines(rowMeans(intervals), thisMedians[,2], lwd=2, col="cadetblue4")

# # fdBM<- evBM[1:(nrow(evBM)-1),2]-evBM[2:(nrow(evBM)),2]



# # as.matrix(vegdist(bigHist))[cbind(1:(nrow(intervals)-1), 2:nrow(intervals))]