# source("~/Dropbox/code/R/amandaTeeth/src/amandaSynonymize.R")
source("~/Dropbox/code/R/amandaTeeth/src/amandaSrc.R")
source("~/Dropbox/code/R/blasto/src/blasto_Birlenbach.R")
# source("http://dl.dropbox.com/s/oe6pzhsuopy32mc/amandaSynonymize.R")

specimenMat <- getSpecimenMatFromMeasurements()
# specimenMat <- merge(specimenMat, getBlastoSpecimenMat(), all=TRUE)
specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(), all=TRUE)
specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(), all=TRUE)
# specimenMat$species <- synonymize(as.character(specimenMat$species))
specimenMat <- specimenMat[!is.na(specimenMat$species),]
rownames(specimenMat) <- make.unique(specimenMat$species)
# rownames(specimenMat) <- make.unique(paste(specimenMat$species, specimenMat$specimen))

oneSpeciesMat <- makeOneSpeciesMatFromSpecimenMat(specimenMat)

upLabels<-c("P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
loLabels <- casefold(upLabels)

#PCA by specimen
	# pcaUp <- prcomp(log(specimenMat[complete.cases(specimenMat[,upLabels]),upLabels]))
	# pcaLo <- prcomp(log(specimenMat[complete.cases(specimenMat[,loLabels]),loLabels]))
	# pcaAll <- prcomp(log(specimenMat[complete.cases(specimenMat[,c(upLabels, loLabels)]),c(upLabels, loLabels)]))
	
#PCA by species
	pcaUp <- prcomp(log(oneSpeciesMat[complete.cases(oneSpeciesMat[,upLabels]),upLabels]))
	pcaLo <- prcomp(log(oneSpeciesMat[complete.cases(oneSpeciesMat[,loLabels]),loLabels]))
	pcaAll <- prcomp(log(oneSpeciesMat[complete.cases(oneSpeciesMat[,c(upLabels, loLabels)]),c(upLabels, loLabels)]))

	# # # # # quartz()
	pc_h<- 2
	pc_v<- 3
	# par(mfrow=c(3,1), mar=c(4, 4, 1, 1))
	famList <- read.csv("~/Dropbox/code/common_dat/taxonomy.csv")

   #by family
	bigList <- cbind(as.character(famList$taxon[famList$taxon%in%unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]), as.character(famList$occurrences.family_name[famList$taxon%in%unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]))
	bigList[bigList[,2] =="",2] <- as.character(famList$occurrences.parent_name[famList$taxon%in%unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))])[bigList[,2] ==""]

   #by parent
	# bigList <- cbind(as.character(famList$taxon[famList$taxon %in% unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]), as.character(famList$occurrences.parent_name[famList$taxon%in%unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]))
	# bigList[bigList[,2] =="",2] <- as.character(famList$occurrences.family_name[famList$taxon %in% unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))])[bigList[,2] ==""]

	shortFam <- sort(unique(bigList[,2]))
	famColors <- rainbow(length(shortFam))
	symbolVec <- array(data=c(15:17, 21:22,24), dim=length(shortFam))

	plot(pcaAll$x[,pc_h], pcaAll$x[,pc_v], xlim=c(min(pcaAll$x[,pc_h]),1.5*max(pcaAll$x[,pc_h])), xlab=paste("PC",pc_h," (",round(100*(pcaAll$sdev[pc_h]^2/sum(pcaAll$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaAll$sdev[pc_v]^2/sum(pcaAll$sdev^2)), digits=1),"%)", sep=""), type="n", main="Upper and Lower P3-M3")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(pcaAll$x[,pc_h], pcaAll$x[,pc_v], labels=rownames(pcaAll$x), cex=0.5, col= famColors[match(bigList[match(rownames(pcaAll$x), bigList[,1]),2], shortFam)])
	text(pcaAll$x[,pc_h], pcaAll$x[,pc_v], labels=rownames(pcaAll$x), pos=4, cex=0.3, col=famColors[match(bigList[match(rownames(pcaAll$x), bigList[,1]),2], shortFam)])
	points(pcaAll$x[,pc_h], pcaAll$x[,pc_v], cex=1.0, pch=symbolVec[match(bigList[match(rownames(pcaAll$x), bigList[,1]),2], shortFam)], col=famColors[match(bigList[match(rownames(pcaAll$x), bigList[,1]),2], shortFam)])
	legend("bottomright", legend=shortFam[shortFam%in%bigList[bigList[,1]%in%rownames(pcaAll$x),2]], pch=symbolVec[shortFam%in%bigList[bigList[,1]%in%rownames(pcaAll$x),2]], col=famColors[shortFam%in%bigList[bigList[,1]%in%rownames(pcaAll$x),2]], box.col="gray50", bg="white", cex=0.55)
	
# # # # # pcaAll$rotation

	plot(pcaUp$x[,pc_h], pcaUp$x[,pc_v], xlim=c(min(pcaUp$x[,pc_h]),1.5*max(pcaUp$x[,pc_h])), xlab=paste("PC",pc_h, " (",round(100*(pcaUp$sdev[pc_h]^2/sum(pcaUp$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaUp$sdev[pc_v]^2/sum(pcaUp$sdev^2)), digits=1),"%)", sep=""), type="n", main="Upper P3-M3", bg="gray90")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(pcaUp$x[,pc_h], pcaUp$x[,pc_v], labels=rownames(pcaUp$x), cex=0.33, col=famColors[match(bigList[match(rownames(pcaUp$x), bigList[,1]),2], shortFam)])
	# text(pcaUp$x[,pc_h], pcaUp$x[,pc_v], labels=gsub("[[:space:]]", "\n", rownames(pcaUp$x)), cex=0.3, col=famColors[match(bigList[match(sapply(strsplit(rownames(pcaUp$x),"[[:punct:]]"), function(x) x[1]), bigList[,1]),2], shortFam)])
	text(pcaUp$x[,pc_h], pcaUp$x[,pc_v], labels=rownames(pcaUp$x), pos=4, cex=0.3, col=famColors[match(bigList[match(sapply(strsplit(rownames(pcaUp$x),"[[:punct:]]"), function(x) x[1]), bigList[,1]),2], shortFam)])
	points(pcaUp$x[,pc_h], pcaUp$x[,pc_v], cex=1.0, pch=symbolVec[match(bigList[match(rownames(pcaUp$x), bigList[,1]),pc_h], shortFam)], col=famColors[match(bigList[match(rownames(pcaUp$x), bigList[,1]),2], shortFam)])
	legend("bottomright", legend=shortFam[shortFam%in%bigList[bigList[,1]%in%rownames(pcaUp$x),2]], pch=symbolVec[shortFam%in%bigList[bigList[,1]%in%rownames(pcaUp$x),2]], col=famColors[shortFam%in%bigList[bigList[,1]%in%rownames(pcaUp$x),2]], box.col="gray50", bg="white", cex=0.55)

 # # # pcaUp$rotation

	plot(pcaLo$x[,pc_h], pcaLo$x[,pc_v], xlim=c(min(pcaLo$x[,pc_h]),1.5*max(pcaLo$x[,pc_h])), xlab=paste("PC",pc_h," (",round(100*(pcaLo$sdev[pc_h]^2/sum(pcaLo$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaLo$sdev[pc_v]^2/sum(pcaLo$sdev^2)), digits=1),"%)", sep=""), type="n", main="Lower P3-M3")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(pcaLo$x[,pc_h], pcaLo$x[,pc_v], labels=gsub("[[:space:]]", "\n", rownames(pcaLo$x)), cex=0.3, col=famColors[match(bigList[match(rownames(pcaLo$x), bigList[,1]),2], shortFam)])
	text(pcaLo$x[,pc_h], pcaLo$x[,pc_v], labels=rownames(pcaLo$x), pos=4, cex=0.3, col=famColors[match(bigList[match(sapply(strsplit(rownames(pcaLo$x),"[[:punct:]]"), function(x) x[1]), bigList[,1]),2], shortFam)])
	points(pcaLo$x[,pc_h], pcaLo$x[,pc_v], cex=0.5, pch=symbolVec[match(bigList[match(rownames(pcaLo$x), bigList[,1]),2], shortFam)], col=famColors[match(bigList[match(rownames(pcaLo$x), bigList[,1]),2], shortFam)])
	legend("bottomright", legend=shortFam[shortFam%in%bigList[bigList[,1]%in%rownames(pcaLo$x),2]], pch=symbolVec[shortFam%in%bigList[bigList[,1]%in%rownames(pcaLo$x),2]], col=famColors[shortFam%in%bigList[bigList[,1]%in%rownames(pcaLo$x),2]], box.col="gray50", bg="white", cex=0.55)

# # pcaLo$rotation

