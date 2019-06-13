require(ape)
require(geiger)

source("~/Dropbox/code/R/common_src/phy_dateTree.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R")
source("~/Dropbox/code/R/artioLimbs/src/ratesFromDiscreteCharacters.R")

ancDescOneBranch<-function(edge, thisTree, thisVec, asr) {
	if (edge[1] %in% names(asr)) ad <- asr[names(asr)==edge[1]] else ad <- thisVec[thisTree$tip.label[edge[2]]]
	if (edge[2] %in% names(asr)) ad <- c(ad, asr[names(asr)==edge[2]]) else ad <- c(ad, thisVec[thisTree$tip.label[edge[2]]])
	ad
}


getClassFromBreaks <- function(x, breaks) {
	for (i in 1:(length(breaks)-1)) {
		if (x > breaks[i] & x[1] <= breaks[i+1]) return (i)
	}
	i
}

getADClasses <- function(x, thisBreaks) {
	if (any(!is.finite(x))) return(c(NA, NA))
	c(getClassFromBreaks (x[1], breaks=thisBreaks), getClassFromBreaks (x[2], breaks=thisBreaks))
}		

getBranchClasses <- function(thisTree, thisVec, asr, thisBreaks) {
	# branchDates <- t(apply(thisTree$edge, 1, getBranchAges, thisTree, dates, use.ranges=FALSE))
	branchProp <- getProportionsOfBranchesWithinIntervals(thisTree, intervals)

	ad <- t(apply(thisTree$edge, 1, ancDescOneBranch, thisTree, thisVec, asr))	# ancestor and descendant states for each branch
	branchChange <- ad[,2] - ad[,1] 		# total (linear) change along branch?

	branchPropChange <- branchChange * branchProp

	adClass <- t(apply(ad, 1, getADClasses, thisBreaks))	# ancestor and descendant chr classes for each branch
	branchClassChanges <- adClass[,1] != adClass[,2]		# does the lineage change character classes along branch?
	
	adMat<-array(NA, dim=c(nrow(thisTree$edge), nrow(intervals), 2), dimnames=list(NULL, rownames(intervals), c("start", "end")))
	
	for (i in seq_len(nrow(ad))) {
		for (j in nrow(intervals):1) {
			if (branchProp[i,j] > 0) {
				if (j == nrow(intervals)) { adMat[i, j, "start"] <- ad[i,1] + ((1 - sum(branchProp[i,]))*branchChange[i]) # 1 - sum(branchProp[i,]) is the proportion of the branch before the window of observation (determined by intervals)
				} else if (branchProp[i,j+1] == 0) { adMat[i, j, "start"] <- ad[i, 1]	# if this is the FA, the start of this branch is the ancState
				} else adMat[i, j, "start"] <- adMat[i, j+1, "end"]				# if this is not the FA, the start of this branch is the end from the last interval
			}
			adMat[i, j, "end"] <- adMat[i, j, "start"] + branchPropChange[i,j]	# the start of this branch is the end from the last interval
		}
	}
	
	apply(adMat, c(1,2), getADClasses, thisBreaks)
	
	# if completely unsampled, add branch as taxon
		# what if it spans > 1 character class?
	# if terminal has FA in same interval, then DO NOT add taxon
	
	# problem is when single branch spans > 1 chr class...
	
}

getIntervalForDate <- function(x, intervals) {
	if (x < min(intervals)) { return(0) 
	} else if (x > max(intervals)) { return(nrow(intervals)+1)
	} else for (i in 1:nrow(intervals)) if (x >= intervals[i,"ageTop"] & x < intervals[i,"ageBase"]) return (i)
	return(NA)
}

getClassCounts <- function(x) {
	sapply(1:5, function(y) sum(x==y, na.rm=TRUE))	#why are there branches with NAs?
}

oneRep<-function(rep, tree, dates, intervals, thisVec, thisBreaks) {
	# thisVec <- array(thisMat[,chr], dimnames=list(rownames(thisMat)))
	thisTree <- ladderize(multi2di(tree), right=FALSE)

	thisVec <- thisVec[match(thisTree$tip.label,names(thisVec))]
	thisVec <- thisVec[is.finite(thisVec)]
	thisTree<-drop.tip(thisTree, thisTree$tip.label[!(thisTree$tip.label %in% names(thisVec))])

	dates <- dates[match(thisTree$tip.label, row.names(dates)),]   #reorders dat to match the taxon order of thisTree
	thisTree <- dateTreeWithRanges(thisTree, dates, within.error=TRUE, random=TRUE, extra=0.5, rootAdj=20.0, include.ranges=TRUE)
	# thisTree<-dateTreeWithList(thisTree, dates, offset=0.05, extendNodes=TRUE, rootAdj=20.0, include.ranges=FALSE)
	# plot(thisTree, cex=0.3, no.margin=TRUE)
	
	asr <- getAncStates(thisVec, thisTree)	# geiger function
	
	branchClasses <- getBranchClasses(thisTree, thisVec, asr, thisBreaks)
	# branchClassChanges <- apply(branchClasses, c(2,3), function (x) if (any(!is.finite(x))) return(FALSE) else x[1] != x[2])
	# thisIndex<-which(branchClassChanges, arr.ind=TRUE)
	# branchClasses[1:2,thisIndex[,1], thisIndex[,2]]
	
	allClasses <- c(sapply(thisVec, getClassFromBreaks, breaks=thisBreaks), sapply(asr, getClassFromBreaks, breaks=thisBreaks))

	# spn = nodes
	# ext/emm = LA
	# chrin/out = branchClass changes
	indates <- thisTree$node.date[(length(thisTree$tip.label)+1):length(thisTree$node.date)]
	inIntv <- sapply(indates, getIntervalForDate, intervals)
	allIntv <- sapply(thisTree$node.date, getIntervalForDate, intervals)
	taxIntv <- apply(thisTree$strat.range, 2, function(x) sapply(x, getIntervalForDate, intervals))
	# taxIntv[taxIntv[,1]>0 & is.na(taxIntv[,2]),2]<-0
	# taxIntv[taxIntv[,2]<=nrow(intervals) & is.na(taxIntv[,1]),1] <- nrow(intervals)+1
	# lapply(1:nrow(intervals), function(x) thisTree$edge[thisTree$edge[,1] %in% (length(thisTree$tip.label)+(which(inIntv==x))),])

	bigList<-list()

	for (intv in 1:nrow(intervals)) {
		intBranches <- which(is.finite(branchClasses[1,,intv]))	# row of edge matrix in this interval
		terminals <- intBranches[thisTree$edge[intBranches,2] <= length(thisTree$tip.label)] #row of edge matrix corresponding to terminal branches
		# newNodeBranches <- which(thisTree$edge[,1] %in% newNodes) #row of edge matrix with a new node
		# thisTree$edge[newNodeBranches,]
		
		# FL <- intBranches[allIntv[thisTree$edge[intBranches,1]] == intv & allIntv[thisTree$edge[intBranches,2]] == intv]
		# bt <- intBranches[allIntv[thisTree$edge[intBranches,1]] > intv & allIntv[thisTree$edge[intBranches,2]] < intv]
		# Ft <- intBranches[allIntv[thisTree$edge[intBranches,1]] == intv & allIntv[thisTree$edge[intBranches,2]] < intv]
		# bL <- intBranches[allIntv[thisTree$edge[intBranches,1]] > intv & allIntv[thisTree$edge[intBranches,2]] == intv]

		newNodes <- which(is.finite(inIntv) & inIntv==intv) + length (thisTree$tip.label) # new nodes in this interval
		spn<- allClasses[newNodes]
				
		changedBranches <- which(apply(branchClasses[,,intv], 2, function(x) all(is.finite(x)) & diff(x) != 0))
		chrIn <- branchClasses[2,changedBranches,intv]
		chrOut <- branchClasses[1,changedBranches,intv]

		b <- intBranches[allIntv[thisTree$edge[intBranches,1]] > intv & allIntv[thisTree$edge[intBranches,2]] <= intv]	# all branches that enter this interval
		b <- b[!b%in%changedBranches]
		bTax <- which(taxIntv[,1] > intv & taxIntv[,2] <= intv)	# all taxa that enter this interval
		same <- c(allClasses[bTax], branchClasses[1,b,intv])		# branches/taxa in this interval that were also present in the previous interval
		
		ext <- allClasses[which(apply(taxIntv, 1, function(x) all(is.finite(x))) & taxIntv[,2]==(intv+1))] # taxa whose interval of last appearance is in the previous interval
		
		bigList[[intv]] <- list(same=same, spn=spn, chrIn=chrIn, chrOut=chrOut, ext=ext)

	}

	return(bigList)

	# newTaxa<-list()
	# oldTaxa<-list()
	# sameTaxa<-list()s
	# for (i in 1:(length(intSp)-1)) {
		# array(NA, dim=5, dimnames=c("ext", "chrOut", "spn", "imm", "chrIn"))
		# newTaxa[[i]] <- intSp[[i]][!(intSp[[i]] %in% intSp[[i+1]])]		#taxa that first appeared in this interval
		# oldTaxa[[i]] <- intSp[[i+1]][!(intSp[[i+1]] %in% intSp[[i]])]	#taxa from previous interval that didnt' make it into this one
		# sameTaxa[[i]] <- intSp[[i]][intSp[[i]] %in% intSp[[i+1]]]		#taxa in this and previous intervalstree
	# }
	# sapply(intSp, length)
	
	# intBr <- apply(branchClasses[1,,], 2, function(x) which(is.finite(x)))
	# sapply(intBr, length)
	# # clears out branch segments in the same interval as their terminal taxon 
	# for (i in seq_len(nrow(intervals))) {
		# this <- array(TRUE, dim = length(intBr[[i]])) 
		# for (j in seq_along(intBr[[i]])) {
			# if (thisTree$edge[intBr[[i]][j], 2] < length(thisTree$tip.label)) {
				# if (i == 1) { this[j] <- FALSE
				# } else if (!(intBr[[i]][j] %in% intBr[[(i-1)]])) this[j] <- FALSE
			# }
		# }
		# intBr[[i]] <- intBr[[i]][this]
	# }
	# sapply(intBr, length)
	
	# newBr<-list()
	# oldBr<-list()
	# sameBr<-list()
	# for (i in 1:(nrow(intervals)-1)) {
		# newBr[[i]] <- intBr[[i]][!(intBr[[i]] %in% intBr[[i+1]])]		#taxa that first appeared in this interval
		# oldBr[[i]] <- intBr[[i+1]][!(intBr[[i+1]] %in% intBr[[i]])]		#taxa from previous interval that didnt' make it into this one
		# sameBr[[i]] <- intBr[[i]][intBr[[i]] %in% intBr[[i+1]]]			#taxa in this and previous intervals
	# }
		
	# # lapply(newTaxa, function(x) rownames(thisMat)[x])
	# # lapply(sameTaxa, function(x) rownames(thisMat)[x])
	# # lapply(intSp, function(x) rownames(thisMat)[x])
	# newMeans<-lapply(newTaxa, function(x) thisVec[x])
	# oldMeans<-lapply(oldTaxa, function(x) thisVec[x])

	# nhist<-lapply(newTaxa, function(x) hist(thisVec[x], breaks=thisBreaks, plot=FALSE))
	# ohist<-lapply(oldTaxa, function(x) hist(thisVec[x], breaks=thisBreaks, plot=FALSE))
	# shist<-lapply(sameTaxa, function(x) hist(thisVec[x], breaks=thisBreaks, plot=FALSE))
	# mwList<-lapply(1:(length(intSp)-1), function(x) wilcox.test(thisVec[newTaxa[[x]]], thisVec[oldTaxa[[x]]]))
	# ksList<-lapply(1:(length(intSp)-1), function(x) ks.test(thisVec[newTaxa[[x]]], thisVec[oldTaxa[[x]]]))
	# cbind(intervals[1:(length(intSp)-1),], sapply(mwList, function(x) x$p.value), sapply(ksList, function(x) x$p.value))
	
	# par(mfrow=c((nrow(intervals)-1),2), mar=c(0.5,0.5,0.5,0.5), cex=0.5)
	# par(mfrow=c((nrow(intervals)-1),1), mar=c(0.5,0.5,0.5,0.5), cex=0.5)
	# for (i in 1:(length(intSp)-1)) {
		# # barplot(nhist[[i]]$counts+ohist[[i]]$counts+shist[[i]]$counts, col="firebrick", main=rownames(intervals)[i], cex.main=0.75, ylim=c(0,max(sapply(1:(length(intSp)-1), function(x) max(nhist[[x]]$counts+ohist[[x]]$counts+shist[[x]]$counts)))), border=NA)
		# barplot(thisHist[i,], col="firebrick", main=rownames(intervals)[i], cex.main=0.75, ylim=c(0,max(thisHist)), border=NA)
		# # if (mwList[[i]]$p.value<0.10) box(col="yellow", lwd=3)
		# # if (ksList[[i]]$p.value<0.10) box(col="red", lwd=1)
		# # barplot(nhist[[i]]$counts+shist[[i]]$counts, col="dodgerblue", add=TRUE, border=NA)
		# # barplot(shist[[i]]$counts, col=rgb(0.5,0.5,0.5), add=TRUE, border=NA)
	
		# # boxplot(x=list(extinct=thisVec[oldTaxa[[i]]], same=thisVec[sameTaxa[[i]]], new=thisVec[newTaxa[[i]]]), col=c("firebrick", "gray50", "dodgerblue"), horizontal=TRUE)
	# }

	# nodeValues <- rbind(as.matrix(thisVec),as.matrix(asr))
	# changeMat <- packContinuousChangeIntoIntervals(tree=thisTree, nodeValues, intervals)
	# lmaMat <- colSums(t(apply(apply(thisTree$edge, 1, getBranchAges, thisTree, dates, use.ranges=TRUE), 2, getLMAForSingleBranch, intervals)))
	# changeMat / lmaMat

}

getCountsFromOneRepList <- function(bigList, intervals, thisBreaks) {
	counts<-array(NA, dim=c(5, (length(thisBreaks)-1), (nrow(intervals)-1)), dimnames=list(c("same", "spn", "chrIn", "chrOut", "ext"), NULL, paste(intervals[1:(nrow(intervals)-1), 2], "Ma", sep="")))
	for (i in 1:(nrow(intervals)-1)) {
		counts[,,i] <- t(sapply(bigList[[i]], function(x) getClassCounts(x)))
	}
	counts
}

getHistFromOneRepList <- function(bigList) {
	t(sapply(bigList, function(x) getClassCounts(unlist(x[c("same", "spn", "chrIn")])))) 
}

# getDistanceBetweenSubsequentIntervals <- function(thisVec, intSp, thisBreaks) {
	# shortVecs <- t(sapply(intSp, function(x) hist(thisVec[x], breaks=thisBreaks)$counts, plot=FALSE))
	# # allTaxa <- unique(unlist(lapply(intSp, function(x) names(thisVec)[x])))
	# dVec <- vector()
	# for (i in 1:(nrow(shortVecs)-1)) dVec[i]<-vegdist(shortVecs[i:(i+1),])
	# dVec
	# # plot(rowMeans(intervals)[1:(nrow(shortVecs)-1)], dVec, type="o")
# }

getDistanceBetweenSubsequentIntervals <- function(x) {
	# return(dim(x))
	distVec<-vector()
	for (intv in 1:(nrow(intervals)-1)) {
		distVec[intv] <- vegdist(x[intv:(intv+1),])
	}
	distVec
}


quartz(width=10, height=6)
plot(tree, edge.color="white", tip.color="white", edge.width=1.5, cex=0.3, no.margin=TRUE, label.offset=1)
do.parallel<-FALSE

tree<-read.nexus("~/Dropbox/code/R/amandaTeeth/dat/ungulates.tre")
dates<-read.csv("~/Dropbox/code/R/common_dat/ranges(6Jan2012).csv", row.names=1)
rownames(dates)<-gsub(" ", "_", row.names(dates))

tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% row.names(dates))])
dates<-dates[(row.names(dates) %in% tree$tip.label),]	#cull taxa in dates list that aren't in the tree
rownames(thisMat)<-gsub(" ", "_", rownames(thisMat))
thisMat<-thisMat[rownames(thisMat) %in% tree$tip.label,]

chr<-"bodyMass"
thisBreaks <- bmBreaks
reps=100

repList<-list()
if (do.parallel) {
	require(parallel)
	repList<-mclapply(seq_len(reps), oneRep, tree=tree, dates=dates, intervals=intervals, thisVec=array(thisMat[,chr], dimnames=list(rownames(thisMat))), thisBreaks=thisBreaks, mc.cores=detectCores()-2)
} else repList<-lapply(seq_len(reps), oneRep, tree=tree, intervals=intervals, dates=dates, thisVec=array(thisMat[,chr], dimnames=list(rownames(thisMat))), thisBreaks=thisBreaks)

# save(repList, file = "~/Desktop/replist.R")
repList<-load("~/Desktop/replist.R")

histMat <- sapply(repList, getHistFromOneRepList, simplify="array")

# bigList <- repList[[3]]
countList <- lapply(repList, getCountsFromOneRepList, intervals, thisBreaks)
countSums <- sapply(countList, function(x) apply(x, 3, rowSums), simplify="array")

# thisCounts<-getCountsFromOneRepList(repList[[1]], intervals, thisBreaks)

ifonly <- array(NA, dim=c(nrow(intervals)-1, 4, reps))
for (i in seq_len(dim(countSums)[3])) {
	for (intv in 1:(nrow(intervals)-1)) {
		ifonly[intv,,i]<-(countSums[,intv,i]/(sum(histMat[intv+1,,i])+(2*countSums[,intv,i])))[2:5]
	}
}
medIfOnly <- apply(ifonly, c(1,2), quantile, probs=c(0.025, 0.50, 0.975))
quartz(width=10, height=6)
plot(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,1], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0,max(medIfOnly)), xlab="Time (Ma)", ylab="Contribution in Isolation", fg="gray75", col.axis="gray75", col.lab="gray75")
overlayCzTimescale()
polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medIfOnly[1,,1], medIfOnly[3,(nrow(intervals)-1):1,1]), col=rgb(t(col2rgb("dodgerblue4"))/255, alpha=0.25), border=NA)
polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medIfOnly[1,,4], medIfOnly[3,(nrow(intervals)-1):1,4]), col=rgb(t(col2rgb("firebrick4"))/255, alpha=0.25), border=NA)
polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(2*medIfOnly[1,,2], 2*medIfOnly[3,(nrow(intervals)-1):1,2]), col=rgb(t(col2rgb("green4"))/255, alpha=0.25), border=NA)
lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,1], lwd=4, col="dodgerblue4")
lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,4], lwd=4, col="firebrick4")
lines(rowMeans(intervals)[1:(nrow(intervals)-1)], 2*medIfOnly[2,,2], lwd=4, col="green4")
# polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medIfOnly[1,,3], medIfOnly[3,(nrow(intervals)-1):1,3]), col="goldenrod1", border=NA)
# lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medIfOnly[2,,3], lwd=3, col="goldenrod4")

medDist <- apply(apply(histMat, 3, getDistanceBetweenSubsequentIntervals), 1, quantile, probs=c(0.025, 0.5, 0.975))
quartz(width=10, height=6.5)
par(mar=c(5,6,4,2))
plot(rowMeans(intervals)[1:(nrow(intervals)-1)], medDist[2,], type="n", xlim=c(max(intervals), min(intervals)), ylim=c(0,max(medDist)), xlab="Time (Ma)", ylab="Change in Body Size Distribution\n(Bray-Curtis Distance)", fg="gray75", col.axis="gray75", col.lab="gray75")
overlayCzTimescale()
polygon(c(rowMeans(intervals)[1:(nrow(intervals)-1)], rowMeans(intervals)[(nrow(intervals)-1):1]), c(medDist[1,], medDist[3,(nrow(intervals)-1):1]), col=rgb(t(col2rgb("plum4"))/255, alpha=0.25), border=NA)
lines(rowMeans(intervals)[1:(nrow(intervals)-1)], medDist[2,], lwd=3, col="plum4")







# # bigList<-repList[[3]]
# t(apply(countSums, c(1,2), median))

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