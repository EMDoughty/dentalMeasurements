require(picante)
source("~/Dropbox/code/R/common_src/sampling.R")

makeCommFromIntervals<-function(thisMat, intervals) {
	thisArray<-array(0, dim=c(nrow(intervals), nrow(thisMat)), dimnames=list(rownames(intervals), rownames(thisMat)))
	for (int in seq_len(nrow(intervals))) thisArray[int,]<-as.numeric(is.finite(thisMat$FO) & is.finite(thisMat$LO) & thisMat$FO>intervals$ageTop[int] & thisMat$LO<intervals$ageBase[int])
	as.matrix(thisArray)
}

makeCommFromCollections<-function(occs, min.occs=10) {
	colList<-unique(occs$collection_no[order(occs$collections.ma_rand)])
	comm<-array(0, dim=c(length(colList), length(unique(occs$taxon))), dimnames=list(colList, gsub("[[:space:]]", "_", unique(occs$taxon))))
	for (thisCol in seq_along(colList)) if (dim(occs[occs$collection_no==colList[thisCol],])[1]>min.occs) comm[thisCol,]<-unique(occs$taxon) %in% occs$taxon[occs$collection_no==colList[thisCol]]
	comm<-as.matrix(comm[rowSums(comm)>0,])
	comm<-comm[sapply(as.numeric(rownames(comm)), function(x) { (unique(occs$collections.ma_rand[occs$collection_no==x])>min(intervals) & unique(occs$collections.ma_rand[occs$collection_no==x])<max(intervals)) }),]
	comm
}

makeCommFromFormations<-function(occs, min.taxa=0) {
	cols<-read.csv("~/Dropbox/code/R/common_dat/allNAMammalsFlat_collections_(2012Feb22).csv", stringsAsFactors = FALSE)
	cols<-cols[cols$collection_no%in%occs$collection_no,]
	# members<-dim(apply(cbind(cols$collections.formation, cols$collections.member), 1, function(x) cat(x[1], x[2])))
	fms<-sort(unique(gsub("[[:punct:]] ", "", gsub("\\\"", "", cols$collections.formation))))
	fms<-fms[-1]
	# colList<-unique(occs$collection_no[order(occs$collections.ma_rand)])
	comm<-array(0, dim=c(length(fms), length(unique(occs$taxon))), dimnames=list(fms, gsub("[[:space:]]", "_", unique(occs$taxon))))
	for (thisFm in seq_along(fms)) comm[thisFm,]<-unique(occs$taxon) %in% occs$taxon[occs$collection_no%in%cols$collection_no[cols$collections.formation==fms[thisFm]]]
	comm<-as.matrix(comm[rowSums(comm)>min.taxa,])
	# comm<-comm[sapply(as.numeric(rownames(comm)), function(x) { (unique(occs$collections.ma_rand[occs$collection_no==x])>min(intervals) & unique(occs$collections.ma_rand[occs$collection_no==x])<max(intervals)) }),]
	comm
}

makePhyListInt<-function(phy, intervals) {
	phyList<-list()
	for (int in seq_len(nrow(intervals))) {
		phyList[[int]]<-drop.tip(phy, phy$tip.label[is.finite(phy$strat.range[,1]) & is.finite(phy$strat.range[,2]) & !(phy$strat.range[,1]>intervals$ageTop[int] & phy$strat.range[,1]<intervals$ageBase[int])])
	}
	phyList
}

makePhyListCol<-function(phy, colList) {
	phyList<-list()
	for (thisCol in seq_len(nrow(colList))) {
		phyList[[thisCol]]<-drop.tip(phy, phy$tip.label[is.finite(phy$strat.range[,1]) & is.finite(phy$strat.range[,2]) & !(phy$strat.range[,1]>=colList[thisCol, 2] & phy$strat.range[,2]<=colList[thisCol, 2])])
	}
	names(phyList)<-colList[,1]
	phyList
}

makePhyListFm<-function(phy, fms) {
	phyList<-list()
	for (thisFm in seq_len(nrow(fms))) {
		phyList[[thisFm]]<-drop.tip(phy, phy$tip.label[is.finite(phy$strat.range[,1]) & is.finite(phy$strat.range[,2]) & (phy$strat.range[,1]<fms[thisFm, 1] | phy$strat.range[,2]>fms[thisFm, 2])])
	}
	names(phyList)<-rownames(fms)
	phyList
}

# poolCols<-function(cols, occs) {
	# i=1
	# while(i<nrow(cols)) {
		# for (j in (i+1):nrow(cols)) {
			# **** if (length(c(agrep("Barstow", "b"), agrep("b", "Barstow")))>0 & distKm(cols$collections.latdec[i], cols$collections.lngdec[i], cols$collections.latdec[j], cols$collections.lngdec[j])<10)
		# }
	# }
# }

occs<-appendMissingOrders(read.csv("~/Dropbox/code/R/common_dat/allNAMammalsFlat-occs-(6Jan2012).csv", stringsAsFactors = FALSE))
occs<-occs[occs$occurrences.order_name%in%c("Artiodactyla", "Perissodactyla"),]
# occs<-occs[occs$occurrences.species_name=="sp.",]
occs<-appendTaxonNames(occs, taxonomic.level="species")
occs<-cbind(occs, collections.ma_rand=getCollectionDatesFromOccs(occs, age.determination="random"))

tree<-read.nexus("~/Dropbox/code/R/amandaTeeth/dat/ungulates.tre")
dates<-read.csv("~/Dropbox/code/R/common_dat/ranges(6Jan2012).csv", row.names=1)
rownames(dates)<-gsub(" ", "_", row.names(dates))

# dates<-t(apply(ranges, 1, function(x) c(FO=runif(1, x[2], x[1]), LO=runif(1, x[4], x[3]))))
tree<-drop.tip(tree, tree$tip.label[!(tree$tip.label %in% row.names(dates))])
dates<-dates[(row.names(dates) %in% tree$tip.label),]	#cull taxa in dates list that aren't in the tree
rownames(thisMat)<-gsub(" ", "_", rownames(thisMat))
thisMat<-thisMat[rownames(thisMat) %in% tree$tip.label,]
thisTree<-ladderize(multi2di(tree), right=FALSE)
dates<-dates[match(thisTree$tip.label, row.names(dates)),]   #reorders dat to match the taxon order of thisTree
thisTree<-dateTreeWithList(thisTree, dates, random=TRUE, offset=0.05, extendNodes=TRUE, rootAdj=20.0, include.ranges=TRUE)

# comm<-makeComm(thisMat, intervals)
comm<-makeCommFromFormations(occs, min.taxa=10)
phy<-prune.sample(comm, thisTree)
# phyList<-makePhyListInt(phy, intervals)
comm<-comm[,match(phy$tip.label, colnames(comm))]
# colList<-unique(occs[occs$collection_no%in%as.numeric(rownames(comm)), c("collection_no",  "collections.ma_rand")])
# colList<-colList[order(colList[,2]),]
# phyList<-makePhyListCol(phy, colList)
fms<-t(sapply(rownames(comm), function(x) range(occs$collections.ma_rand[occs$collection_no%in%cols$collection_no[cols$collections.formation==x]])))
phyList<-makePhyListFm(phy, fms)

commList<-list()
for (i in seq_along(phyList)) commList[[i]]<-matrix(comm[names(phyList)[i],phyList[[i]]$tip.label], nrow=1, dimnames=list(NULL, phyList[[i]]$tip.label))
reps=999
rez<-mapply(ses.mpd, commList, lapply(phyList, cophenetic), MoreArgs=list(null.model = "taxa.labels", abundance.weighted = FALSE, runs = reps), SIMPLIFY=FALSE)
rezStats<-data.frame(collection=rownames(comm), p.value=sapply(rez, function(x) sum(x$mpd.obs.rank>1)/reps), sig=sapply(rez, function(x) sum(x$mpd.obs.rank>1)/reps)<0.05)
# ses.mpd(commList[[1]], cophenetic(phyList[[1]]), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)

# par(mfrow = c(ceiling(sqrt(nrow(colList))), ceiling(sqrt(nrow(colList)))))
par(mfrow = c(ceiling(sqrt(nrow(comm))), ceiling(sqrt(nrow(comm)))))
for (i in row.names(comm)) {
	plot(phyList[[i]], show.tip.label = FALSE, direction="up", no.margin = TRUE, main = i)
	if (rezStats$sig[rezStats$collection==i]) rect(par()$usr[1], par()$usr[3], par()$usr[2], par()$usr[4], border="red", col=rgb(1,0,0,0.25), lwd=3)
	tiplabels(tip = match(colnames(comm)[comm[i,] > 0], phyList[[i]]$tip.label), pch = 19, cex = 2, col=rainbow(sum(comm[i, ] > 0), alpha=0.5))
	text(0,0, labels=i, pos=4, cex=0.75)
}
cbind(colList, rezStats)
cbind(fms, rezStats)

# pd.result<-pd(comm, phy)
# plot(rowMeans(intervals), pd.result$PD, type="l", xlim<-c(max(intervals), min(intervals)))

# # commList<-list()
# for (i in seq_len(nrow(comm))) commList[[i]]<-matrix(comm[i,phyList[[i]]$tip.label], nrow=1, dimnames=list(NULL, phyList[[i]]$tip.label))
# reps=999
# rez<-mapply(ses.mpd, commList, lapply(phyList, cophenetic), MoreArgs=list(null.model = "taxa.labels", abundance.weighted = FALSE, runs = reps), SIMPLIFY=FALSE)
# data.frame(collection=as.numeric(rownames(comm)), p.value=sapply(rez, function(x) sum(x$mpd.obs.rank>1)/reps), sig=sapply(rez, function(x) sum(x$mpd.obs.rank>1)/reps)<0.05)
# # ses.mpd(commList[[1]], cophenetic(phyList[[1]]), null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)

# phydist <- cophenetic(phy)
# ses.mpd.result <- ses.mpd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
# ses.mpd.result <- ses.mpd(comm[1,], phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
# ses.mpd.result

# ses.mntd.result <- ses.mntd(comm, phydist, null.model = "taxa.labels", abundance.weighted = FALSE, runs = 999)
# ses.mntd.result

# plot(rowMeans(intervals), ses.mpd.result$mpd.obs, type="l", xlim<-c(max(intervals), min(intervals)), ylim=c(0, max(ses.mpd.result$mpd.obs)))
# plot(rowMeans(intervals), ses.mntd.result$mntd.obs, type="l", xlim<-c(max(intervals), min(intervals)), ylim=c(0, max(ses.mntd.result$mntd.obs)))
# plot(ses.mntd.result$ntaxa, ses.mntd.result$mntd.obs)
# lines(rowMeans(intervals), ses.mntd.result$mntd.obs)
