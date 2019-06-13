require(ape)
require(geiger)
require(phytools)

source("~/Dropbox/code/R/NAPC2014/src/NAPC2014_src.R")
source("~/Dropbox/code/R/common_src/phy_dateTree.R")
source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/amandaTeeth/src/amandaSynonymize.R")
source("~/Dropbox/code/R/amandaTeeth/src/amandaSrc.R")
source("~/Dropbox/code/R/blasto/src/blasto_Birlenbach.R")
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R")

ranges <- read.csv("~/Dropbox/code/R/amandaTeeth/dat/ranges_20130918_total.csv", row.names=1)

specimenMat <- getSpecimenMatFromMeasurements()
# specimenMat <- merge(specimenMat, getBlastoSpecimenMat(), all=TRUE)
specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(), all=TRUE)
specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(), all=TRUE)
specimenMat <- specimenMat[!is.na(specimenMat$species)& specimenMat$species!="",]

specimenMat$species <- synonymize(as.character(specimenMat$species))
rownames(specimenMat) <- make.unique(specimenMat$species)
specimenMat[,sapply(specimenMat, is.numeric)] <- specimenMat[,sapply(specimenMat, is.numeric)]/10  #converts mm measurements to cm for compatibility with Janis regressions
specimenMat <- transform(specimenMat, p4_a=p4_l*p4_w, m1_a=m1_l*m1_w, m2_a=m2_l*m2_w, m3_a=m3_l*m3_w, M2_A=M2_L*M2_W)
# specimenMat <- specimenMat[-(which(specimenMat$specimen=="UNSM 54805")),]

# thisMat <- aggregate(specimenMat, by=list(species=specimenMat$species), mean, na.rm=TRUE)
thisMat <- aggregate(specimenMat[,-(c(1:2,37:46))], by=list(species=specimenMat$species), median, na.rm=TRUE)
thisMat[sapply(thisMat, is.nan)] <- NA
# rownames(specimenMat) <- make.unique(paste(specimenMat$species, specimenMat$specimen))
thisMat <- appendMissingPaleoDBSpecies(thisMat, ranges)		# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
# thisMat[,"bodyMass"] <- getBodyMassVectorFromThisMat(thisMat)
# thisMat$bodyMass <- fillMissingBodyMasses(thisMat)				# this fills taxa missing their body mass with the average body mass of its cogeners

rownames(thisMat) <- gsub(pattern="[[:space:]]", replacement="_", x=rownames(thisMat))
rownames(ranges) <- gsub(pattern="[[:space:]]", replacement="_", x=rownames(ranges))

# do.parallel=FALSE
# if (do.parallel) require(parallel)
# reps=10
tree <- read.nexus("~/Dropbox/code/R/NAPC2014/dat/analysis.tre")
# plot.phylo(phy, direction="upwards", no.margin=TRUE, cex=0.25)
shortTree <- getNAPCTree(tree, ranges, resolve=TRUE)

upLabels<-c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
loLabels <- casefold(upLabels)
shortMat <- thisMat[shortTree$tip.label,loLabels]
# shortMat <- thisMat[shortTree$tip.label,upLabels]
# shortMat <- thisMat[shortTree$tip.label,c(upLabels, loLabels)]
shortMat <- log(shortMat[complete.cases(shortMat),])

shortTree <- drop.tip(shortTree, tip=shortTree$tip.label[!shortTree$tip.label %in% rownames(shortMat)])

ppca <- phyl.pca(tree=shortTree, Y=shortMat)

	pc_h<- 2
	pc_v<- 3
	# par(mfrow=c(3,1), mar=c(4, 4, 1, 1))
	famList <- read.csv("~/Dropbox/code/common_dat/taxonomy.csv", stringsAsFactors=FALSE)

   #by family
   famList$taxon <- gsub(pattern="[[:space:]]", replacement="_", x=famList$taxon)
	bigList <- famList[famList$taxon %in% unique(rownames(shortMat)),c("taxon", "occurrences.family_name", "occurrences.parent_name")]
	bigList[bigList[,2] =="",2] <- bigList[bigList[,2] =="",3]

   #by parent
	# bigList <- cbind(as.character(famList$taxon[famList$taxon %in% unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]), as.character(famList$occurrences.parent_name[famList$taxon%in%unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]))
	# bigList[bigList[,2] =="",2] <- as.character(famList$occurrences.family_name[famList$taxon %in% unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))])[bigList[,2] ==""]

	shortFam <- sort(unique(bigList[,2]))
	famColors <- rainbow(length(shortFam))
	symbolVec <- array(data=c(15:17, 21:22,24), dim=length(shortFam))

	# plot(ppca$S[,pc_h], ppca$S[,pc_v], xlim=c(min(ppca$S[,pc_h]),1.5*max(ppca$S[,pc_h])), xlab=paste("PC",pc_h," (",round(100*(pcaAll$sdev[pc_h]^2/sum(pcaAll$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaAll$sdev[pc_v]^2/sum(pcaAll$sdev^2)), digits=1),"%)", sep=""), type="n", main="Upper and Lower P3-M3")
	plot(ppca$S[,pc_h], ppca$S[,pc_v], xlim=c(min(ppca$S[,pc_h]),1.5*max(ppca$S[,pc_h])), type="n", main="Upper and Lower P3-M3")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(ppca$S[,pc_h], ppca$S[,pc_v], labels=rownames(ppca$S), cex=0.5, col= famColors[match(bigList[match(rownames(ppca$S), bigList[,1]),2], shortFam)])
	text(ppca$S[,pc_h], ppca$S[,pc_v], labels=rownames(ppca$S), pos=4, cex=0.3, col=famColors[match(bigList[match(rownames(ppca$S), bigList[,1]),2], shortFam)])
	points(ppca$S[,pc_h], ppca$S[,pc_v], cex=1.0, pch=symbolVec[match(bigList[match(rownames(ppca$S), bigList[,1]),2], shortFam)], col=famColors[match(bigList[match(rownames(ppca$S), as.character(bigList[,"taxon"])),2], shortFam)])
	legend("bottomright", legend=shortFam, pch=symbolVec, col=famColors, box.col="gray50", bg="white", cex=0.55)
