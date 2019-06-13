
source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/amandaTeeth/src/amandaSrc.R")
source("~/Dropbox/code/R/blasto/src/blasto_Birlenbach.R")
source("~/Dropbox/code/R/common_src/sampling.R")
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R")

specimenMat <- getSpecimenMatFromMeasurements()
specimenMat <- merge(specimenMat, getBirlenbachBlastoSpecimens(), all=TRUE)
specimenMat <- merge(specimenMat, getLiteratureSpecimenMat(), all=TRUE)

specimenMat[sapply(specimenMat, is.nan)] <- NA
# namesMat <- cbind(sort(unique(specimenMat$species)), as.character(getCurrentNames(sort(unique(specimenMat$species)))))
# namesMat[apply(namesMat, 1, function(x) x[1] != x[2]),]
specimenMat$species <- getCurrentTaxa(specimenMat$species)

upLabels<-c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
loLabels <- casefold(upLabels)

# thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), mean, na.rm=TRUE)	#aggregate by MEAN
thisMat <- aggregate(specimenMat[,c(upLabels, loLabels)], by=list(species=specimenMat$species), median, na.rm=TRUE)	#aggregate by MEDIAN

thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species)
rownames(thisMat) <- thisMat$species
thisMat[,sapply(thisMat, is.numeric)] <- thisMat[,sapply(thisMat, is.numeric)] / 10  #converts mm measurements to cm for compatibility with Janis regressions
thisMat <- transform(thisMat, p4_a=p4_l*p4_w, m1_a=m1_l*m1_w, m2_a=m2_l*m2_w, m3_a=m3_l*m3_w, M2_A=M2_L*M2_W)
thisMat[sapply(thisMat, is.nan)] <- NA


# occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=attr,class,genus,ident,coords,loc,paleoloc,stratext,lithext,geo,rem,ref&limit=all")
# occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
# occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Artiodactyla,Perissodactyla&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
# occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)

occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
# occs <- appendTaxonNames1.2(occs, taxonomic.level="species", keep.indet=FALSE)
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)

ranges <- getTaxonRangesFromOccs(occs=occs[occs$order %in% c("Artiodactyla","Perissodactyla"),], random=TRUE)
rownames(ranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(ranges))
# ranges[order(ranges[,1]),]
# as.matrix(sort(unique(rownames(thisMat)[!rownames(thisMat) %in% rownames(ranges)])))
# as.matrix(sort(unique(as.character(occs$taxon)[!as.character(occs$taxon) %in% as.character(thisMat$species)])))
# as.matrix(sort(unique(rownames(ranges)[!rownames(ranges) %in% rownames(thisMat)])))
# as.matrix(unique(specimenMat$species[!specimenMat$species %in% rownames(ranges)]), ncol=1) #taxa with measures, but no ranges

uniqTax <- getTaxonomyForTaxa(gsub("_", " ", rownames(thisMat)))
uniqTax$species <- gsub(" ", "_", as.character(uniqTax$species))
family.names <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$species)]
regMat <- read.csv("~/Dropbox/ungulate_RA/BodyMassRegressionAssignment/regressionLabelsJDM.csv")
reg.vec <- regMat$reg[!is.na(regMat$family)][match(family.names, regMat$family[!is.na(regMat$family)])]  # this is the regression "labels" for the species from thisMat in the correct order, based on family name
names(reg.vec) <- rownames(thisMat)    

genus.names <- uniqTax$genus[match(x=names(reg.vec)[is.na(reg.vec)], table=uniqTax$species)]
reg.vec[is.na(reg.vec)] <- regMat$reg[match(genus.names, regMat$genus)]   # this is the regression "labels" for the species from thisMat in the correct order, based on genus name
    
    # species.names <- uniqTax$accepted_name[match(x=rownames(thisMat), table=uniqTax$accepted_name)]
    # reg.vec[is.na(reg.vec)] <- regMat$reg[match(species.names, regMat$taxon)][is.na(reg.vec)]   # this is the regression "labels" for the species from thisMat in the correct order, based on genus name
    
    #match_genus<- ut$genus[match(rownames(thisMat),ut$species)]
    #vec[is.na(vec)] <- regMat$reg[match(match_genus[is.na(vec)], regMat$taxon)]
    #species.names <- uniqTax$accepted_name[match(x=rownames(thisMat),table=uniqTax$accepted_name)]
    #reg.vec[is.na(reg.vec)] <- regMat$reg[match(species.names,regMat$accepted_names)][is.na(reg.vec)]
thisMat <- cbind(thisMat, reg.vec)
thisMat[,"bodyMass"] <- getBodyMassVectorFromThisMatAllMeasures(thisMat, linked.files=TRUE)
# thisMat <- appendMissingPaleoDBSpecies(thisMat, ranges)		# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
thisMat$bodyMass <- fillMissingBodyMasses(thisMat)			# this fills taxa missing their body mass with the average body mass of its cogeners
thisMat[!sapply(thisMat, is.finite)] <- NA

intervals <- makeIntervals(1, 56.5, 2)
intList <- listifyMatrixByRow(intervals)

	richness <- matrix(0, nrow=nrow(intervals), ncol=4)
	for (intv in seq_len((nrow(intervals)))) {
		thisInt <- rownames(ranges)[ranges[,"FO"] >= intervals$ageTop[intv] & ranges[, "LO"] < intervals$ageBase[intv]]
		richness[intv,1] <- length(thisInt)
		richness[intv,2] <- sum(thisInt %in% rownames(thisMat)[apply(thisMat, 1, function(x) any(is.finite(as.numeric(x))))])
		richness[intv,3] <- sum(thisInt %in% rownames(thisMat)[is.finite(thisMat$bodyMass)])
		# print(thisInt[!thisInt %in% rownames(thisMat)])	
	}
	# quartz("Taxonomic Richness", width=3.3, height=1.65)
	# par(mfrow=c(2,1), mar=c(3,3,0.5,0.5), cex=0.5)

	# plot raw richnesses
		# plot(rowMeans(intervals), richness[,1], type="n", xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)), ylim=c(0, max(richness, na.rm=TRUE)), main="", xlab="Time (Ma)", ylab="Number of Species")
		# overlayCzTimescale()
		# lines(rowMeans(intervals), richness[,1], type="l", lty=3)
		# # lines(rowMeans(intervals), richness[,2], type="l", lwd=1.5,col="firebrick")
		# lines(rowMeans(intervals), richness[,3], type="l", lwd=1.5,col="darkolivegreen")
		
	# plot proportional richnesses
		par(mar=c(3,4,0.5,0.5), cex=0.66)
		plot(rowMeans(intervals), richness[,3]/richness[,1], type="n", xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)), ylim=c(0,1), main="", xlab="Time (Ma)", ylab="% sampled")
		overlayCzTimescale()
		lines(rowMeans(intervals), richness[,2]/richness[,1], lwd=1.5, type="o", pch=15, col="firebrick")
		lines(rowMeans(intervals), richness[,3]/richness[,1], lwd=1, type="o", pch=21, col="green4", bg="green1")

		richness<-NULL

# thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species) 

# pcVec <- pcaLo$x[,2]
# pcVec <- c(pcVec,-(pcaUp$x[!rownames(pcaUp$x)%in%rownames(pcaLo$x),2])) # gathers PC2 scores from pcaUp (upper PCA)
# thisMat[,"PC2"] <- pcVec[match(thisMat$species, names(pcVec))]
# thisMat[,"PC3"] <- pcaLo$x[match(thisMat$species, rownames(pcaLo$x)),3]

# famList <- read.csv("~/Dropbox/code/common_dat/taxonomy.csv")
# write.csv(cbind(rownames(thisMat), famList[match(rownames(thisMat), famList$taxon),]), file="~/Desktop/famList.csv")

# write.csv(as.matrix(thisMat[is.finite(thisMat$bodyMass),]$bodyMass, ncol=1), file="~/Dropbox/code/R/amandaTeeth/results/bodyMassEstimates.csv")
# ranges[!(rownames(ranges) %in% rownames(thisMat[is.finite(thisMat$bodyMass),])),] # taxa with ranges, but no body mass estaimate (yet)
# ranges[!(rownames(ranges) %in% rownames(thisMat[is.finite(thisMat$PC3),])),] # taxa with ranges, but no PCA (yet)
# # thisMat<-thisMat[rownames(thisMat)%in%occs$taxon,]

bigList <- unique(occs[occs$order %in% c("Artiodactyla","Perissodactyla"),c("accepted_name","family","order")])
# bigList[order(bigList$family, bigList$accepted_name),]
# thisMat <- thisMat[rownames(thisMat) %in% bigList$accepted_name[bigList$order=="Artiodactyla"],]

reps <- 10
quota <- 0.4
bootstrapSpecimens <- FALSE
bootstrapSpecies <- FALSE
bootstrapSpeciesWithinIntervals <- FALSE
plotHist <- FALSE
do.range.through <- TRUE
do.heuristic <- FALSE
do.subsample <- TRUE
do.disparity <- FALSE

do.parallel <- FALSE
if (do.parallel) require(parallel)

if (bootstrapSpecies) holderMat <- thisMat

if (plotHist) {
	quartz("Guild Histograms")
	par(mfrow=c((nrow(intervals)), 3), mar=c(0,0,0.75,0), cex.axis=0.5, cex.main=0.75)
}

# bmBreaks <- seq(from=min(thisMat$bodyMass, na.rm=TRUE), to=max(thisMat$bodyMass, na.rm=TRUE), length.out=11)
# pc2Breaks <- seq(from=min(pcVec), to=max(pcVec), length.out=11)
# pc3Breaks <- seq(from=min(pcaLo$x[,3]), to=max(pcaLo$x[,3]), length.out=11)
bmBreaks <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf) #Janis 2000  max(thisMat$bodyMass, na.rm=TRUE)
# bmBreaks <- c(-Inf, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, Inf) #Badgely and Fox 2000
# bmBreaks <- hist(thisMat$bodyMass, plot=FALSE)$breaks
# pc2Breaks <- hist(thisMat$PC2, plot=FALSE)$breaks
# pc3Breaks <- hist(thisMat$PC3, plot=FALSE)$breaks
# bmBox <- array(NA, dim=c(nrow(intervals), length(bmBreaks)-1, reps), dimnames=list(rownames(intervals), (bmBreaks[1:(length(bmBreaks)-1)]+bmBreaks[2:(length(bmBreaks))])/2, paste("Rep", seq_len(reps), sep="")))
# pc2Box <- array(NA, dim=c(nrow(intervals), length(pc2Breaks)-1, reps), dimnames=list(rownames(intervals), (pc2Breaks[1:(length(pc2Breaks)-1)]+pc2Breaks[2:(length(pc2Breaks))])/2, paste("Rep", seq_len(reps), sep="")))
# pc3Box <- array(NA, dim=c(nrow(intervals), length(pc3Breaks)-1, reps), dimnames=list(rownames(intervals), (pc3Breaks[1:(length(pc3Breaks)-1)]+pc3Breaks[2:(length(pc3Breaks))])/2, paste("Rep", seq_len(reps), sep="")))
# bstatBox <- array(NA, dim=c(nrow(intervals), 6, reps), dimnames=list(rownames(intervals), c("mean", "median", "variance", "skewness", "range_min", "range_max"), paste("Rep", seq_len(reps), sep="")))
# pc2statBox <- array(NA, dim=c(nrow(intervals), 6, reps), dimnames=list(rownames(intervals), c("mean", "median", "variance", "skewness", "range_min", "range_max"), paste("Rep", seq_len(reps), sep="")))
# pc3statBox <- array(NA, dim=c(nrow(intervals), 6, reps), dimnames=list(rownames(intervals), c("mean", "median", "variance", "skewness", "range_min", "range_max"), paste("Rep", seq_len(reps), sep="")))
# if(do.disparity) dispMat <- matrix(NA, nrow=reps, ncol=nrow(intervals), dimnames=list(paste("Rep", seq_len(reps), sep=""), rownames(intervals))) else dispMat <- NULL

# par(mfrow=c(1,2))
# hist(thisMat$bodyMass)
# hist(thisMat$bodyMass, breaks=c(min(thisMat$bodyMass, na.rm=TRUE), bmBreaks[2:5], max(thisMat$bodyMass, na.rm=TRUE)))

repIntSp <- list()
thisCol <- thisMat[,"bodyMass"]
breaks <- bmBreaks
# thisCol <- which(colnames(thisMat)=="PC3")
# breaks <- pc3Breaks
# intList <- listifyMatrixByRow(intervals)

print("Beginning sampling routine...")

for (rep in seq_len(reps)) {
	cat("Beginning Rep", rep, "of", reps, "...\r")
	if (bootstrapSpecimens) {
		thisMat <- specimenMat[sample.int(nrow(specimenMat), size=nrow(specimenMat), replace=TRUE),]
		thisMat <- aggregate(thisMat, by=list(species=specimenMat$species), mean, na.rm=TRUE)
		# thisMat <- thisMat[,apply(!sapply(thisMat, is.na), 2, any)]
		rownames(thisMat) <- thisMat$species
		thisMat[sapply(thisMat, is.nan)] <- NA
		# thisMat<-cbind(thisMat, cbind(FO=vector(length=nrow(thisMat), mode="numeric"), LO=vector(length=nrow(thisMat), mode="numeric")))
		thisMat[,"reg"] <- as.character(famList$reg[match(thisMat$species,famList$taxon)])
		thisMat[,"bodyMass"] <- makeBodyMasses(thisMat, regList, best.only=TRUE)
		thisMat[,"PC2"] <- pcVec[match(thisMat$species, names(pcVec))]
		thisMat[,"PC3"] <- pcaLo$x[match(thisMat$species, rownames(pcaLo$x)),3]
	}
	if (bootstrapSpecies) thisMat <- holderMat[sample.int(n=nrow(thisMat), size=nrow(thisMat), replace=TRUE),]

	col.dates <- getCollectionAgesFromOccs(occs=occs[, c("collection_no", "max_ma", "min_ma")], random=TRUE)
	occDates <- col.dates$collection_age[match(occs$collection_no, col.dates$collection_no)]
	intOccs <- apply(intervals, 1, function(thisIntv) occs$occurrence_no[occDates > thisIntv[1] & occDates <= thisIntv[2]])
	# intTaxa <- sapply(intOccs, function(x) unique(occs$accepted_name[occs$occurrence_no %in% x]))
	intSp <- sapply(intOccs, function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$accepted_rank =="species" & occs$occurrence_no %in% x]))), rownames(thisMat)))
	
	if (do.subsample) { 
		nOccs <- sapply(intOccs, length)
	 	nTaxa <- sapply(intSp, length)
	 	nTaxa <- 0
		quota <- max(c(max(nTaxa), min(nOccs)))
		cat("Subsampling quota set to", quota, "occurrences")

		# intSp <- sapply(intSp, function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$occurrence_no %in% sample(x=x, size=quota)]))), rownames(thisMat)))
		# ssOccs <- lapply(intOccs, sample, size=quota)
		intSp <- sapply(lapply(intOccs, sample, size=quota), function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$accepted_rank =="species" & occs$occurrence_no %in% x]))), rownames(thisMat)))
	}
	if (do.range.through) intSp <- makeRangeThroughOneRep(intSp)
	repIntSp[[rep]] <- intSp
}

# x <- repIntSp[[1]]
# sum(sapply(x, length))

# varAdjSp <- function(spIndex, traitVec) {
	# intVec <- sort(traitVec[spIndex])
	# null_var <- c(min(intVec, na.rm=TRUE), sort(sample(traitVec[traitVec >= min(intVec, na.rm=TRUE) & traitVec < max(intVec, na.rm=TRUE)], size=length(intVec)-2)), max(intVec, na.rm=TRUE))
	# null_unif <- c(min(intVec, na.rm=TRUE), sort(runif(n=length(intVec)-2, min=min(intVec), max=max(intVec))), max(intVec, na.rm=TRUE))
	# null_equip <- c(min(traitVec, na.rm=TRUE), sort(sample(traitVec[is.finite(traitVec)], size=length(intVec))), max(traitVec, na.rm=TRUE))
	# if (length(intVec) > 1) c(variance=var(intVec), null_var=var(null_var), varAdj=var(intVec[2:length(intVec)] - intVec[1:(length(intVec)-1)]), null_unif=var(null_unif[2:length(null_unif)] - null_unif[1:(length(null_unif)-1)]), null_equip=var(null_equip[2:length(null_equip)] - null_equip[1:(length(null_equip)-1)])) else return(NA)
# }

# varCube <- sapply(repIntSp, function(x) sapply(x, varAdjSp, thisMat$bodyMass), simplify="array")
# varBox <- apply(varCube, c(1,2), quantile, probs=c(0.025, 0.50, 0.975))
# par(mfrow=c(3,1))
# plot(rowMeans(intervals), varBox[3,"varAdj",], xlim=c(55,0), ylim=c(0, max(varBox[3,"varAdj",])), type="n")
# overlayCzTimescale()
# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(varBox[1,"null_unif",], rev(varBox[3,"null_unif",])), col=alphaColor("red1", 0.25), border="red1")
# lines(rowMeans(intervals), varBox[2,"varAdj",], col="green4")
# # lines(rowMeans(intervals), varBox[2,"null_unif",], col="black")
# plot(rowMeans(intervals), varBox[3,"varAdj",], xlim=c(55,0), ylim=c(0, max(varBox[3,"varAdj",])), type="n")
# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(varBox[1,"null_equip",], rev(varBox[3,"null_equip",])), col=alphaColor("steelblue1", 0.25), border="steelblue1")
# # lines(rowMeans(intervals), varBox[2,"null_equip",], col="steelblue4")
# # polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(varBox[1,2,], rev(varBox[3,2,])), col=alphaColor("green1", 0.25), border="green1")
# lines(rowMeans(intervals), varBox[2,"varAdj",], col="green4")

# plot(rowMeans(intervals), varBox[3,"variance",], xlim=c(55,0), ylim=c(min(varBox[,1,]), max(varBox[,1,])), type="n")
# overlayCzTimescale()
# # polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(varBox[1,1,], rev(varBox[3,1,])), col=alphaColor("green1", 0.25), border="green1")
# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(varBox[1,"null_var",], rev(varBox[3,"null_var",])), col=alphaColor("green1", 0.25), border="green1")
# lines(rowMeans(intervals), varBox[2,"variance",])

# load('~/Desktop/fullSearch100Reps.R')
# load('~/Desktop/mostRecentResults.R')

####################################################################################################################################
### Handley analysis of taxonomic distributions
print("Beginning taxonomic Handley analysis...")
shortFam <- sort(unique(bigList$family))

taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
med.n <- median(sapply(repIntSp, function(x) length(unique(unlist(sapply(x, function(y) rownames(thisMat)[y]))))))
# optList_tax_median <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, do.parallel=do.parallel)	# based on means
optList_tax_median_heuristic <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on means
optList_tax_median_full <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=FALSE, do.parallel=do.parallel)	# based on means
optList_tax_median <- optList_tax_median_full
# optList_tax_allReps <- list()
optList_tax_allReps_heuristic <- list()
optList_tax_allReps_full <- list()
for (this.rep in seq_len(reps)) {
	taxCube <- sapply(repIntSp, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% rownames(thisMat)[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
	this.n <- length(unique(unlist(sapply(repIntSp [[this.rep]], function(x) rownames(thisMat)[x]))))
	# optList_tax_allReps[[this.rep]] <- doHandleyTest(taxCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, do.parallel=do.parallel)	# based on means
	optList_tax_allReps_heuristic[[this.rep]] <- doHandleyTest(taxCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on means
	optList_tax_allReps_full[[this.rep]] <- doHandleyTest(taxCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=FALSE, do.parallel=do.parallel)	# based on means
}

####################################################################################################################################
### Handley analysis of body mass distributions
print("Beginning body mass Handley analysis...")

countCube <- sapply(repIntSp, function(y) sapply(y, function(x) hist(thisCol[x], breaks=breaks, plot=FALSE)$counts), simplify = "array")
countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 

# optList_bm_median <- doHandleyTest(countBox[2,,], n=nrow(thisMat), do.heuristic=do.heuristic)
optList_bm_median_heuristic <- doHandleyTest(countBox[2,,], n=nrow(thisMat), do.heuristic=TRUE)
optList_bm_median_full <- doHandleyTest(countBox[2,,], n=nrow(thisMat), do.heuristic=FALSE)
optList_bm_median <- optList_bm_median_full
	# # quartz(width=3.3, height=9.8)
	# # par(mfrow=c(ncol(countBox[2,,])/2,2), mar=c(0,3,0.5, 0.5), mfg=c(2,1))
	# # for (i in seq(from=1, to=ncol(countBox[2,,]), by=1)) barplot(countBox[2,,][,i], width=c(0.68897, 0.68897, 0.778151, 0.522879, 1.18786),space=0, cex.axis=0.5, ylim=c(0,30))
	
	# quartz()
	# thisTab <- table(unlist(sapply (optList_bm_allReps, function(x) x[[(length(x) - 1)]]$optBreaks)))
	# thisTab <- array(thisTab[match(seq_len(nrow(intervals)), names(thisTab))], dimnames=list(rownames(intervals)))
	# barplot(rev(thisTab)/reps, cex.names=0.5, ylim=c(0,1))
	# abline(h=c(0.95, 0.75, 0.5), lty=3, col="gray50")

# load(file="~/Desktop/optList_bm_allReps.R")

# optList_bm_allReps <- list()
optList_bm_allReps_heuristic <- list()
optList_bm_allReps_full <- list()
for (this.rep in seq_len(reps)) {
	this.n <- length(unique(unlist(sapply(repIntSp [[this.rep]], function(x) rownames(thisMat)[x]))))
	# optList_bm_allReps[[this.rep]] <- doHandleyTest(countCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, do.parallel=do.parallel)	# based on means
	optList_bm_allReps_heuristic[[this.rep]] <- doHandleyTest(countCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on means
	optList_bm_allReps_full[[this.rep]] <- doHandleyTest(countCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=FALSE, do.parallel=do.parallel)	# based on means
}
# sapply(optList_bm_allReps, function(x) length(x) - 1)
# sapply(optList_bm_allReps_heuristic, function(x) length(x) - 1)
# sapply(optList_bm_allReps_full, function(x) length(x) - 1)
save(repIntSp, optList_tax_median_heuristic, optList_tax_median_full, optList_tax_allReps_heuristic, optList_tax_allReps_full, optList_bm_median_heuristic, optList_bm_median_full, optList_bm_allReps_heuristic, optList_bm_allReps_full, file="~/Desktop/mostRecentResults.R")

	### Handley-block histogram series
	this.rep <- 1
	old.shift <- nrow(intervals)
	shift.ints <- rev(which(intervals$ageBase %in% optList_bm_allReps[[this.rep]][[length(optList_bm_allReps[[this.rep]]) - 1]]$optBreaks))
	par(mfrow=c(1,length(shift.ints)+1), mar=c(4,2,0.5,0.5), col.axis="gray50", col.lab="gray50", fg="gray50")
	for (this.shift in c(shift.ints, 0)) {
		hist(thisMat$bodyMass[unique(unlist(repIntSp[[this.rep]][seq(from=old.shift, to=(this.shift+1))]))], freq=FALSE, breaks=c(min(thisMat$bodyMass, na.rm=TRUE),bmBreaks[2:6],max(thisMat$bodyMass, na.rm=TRUE)), col=rainbow(n=length(shift.ints)+1), main="", xlab="log Body Mass", ylab="", ylim=c(0,1))
	}

par(mfrow=c(2,1), mar=c(5, 4, 4, 1))
hist(sapply(optList_tax_allReps, function(x) length(x) - 2), breaks=seq(0.5, 10.5, 1.0), col="orchid4", main="Number of taxonomic shifts in each rep", xlab="Number of Shifts")
hist(sapply(optList_bm_allReps, function(x) length(x) - 2), breaks=seq(0.5, 10.5, 1.0), col="firebrick4", main="Number of body mass shifts in each rep", xlab="Number of Shifts")
	
quants <- apply(sapply(repIntSp, function(y) sapply(y, function(x) quantile(thisMat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
# quants <- quants[,-56]


####################################################################################################################################
	### number of Replicates with a shift in that interval
	quartz()
	par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
	optList_tax_allReps <- optList_tax_allReps_heuristic
	optList_tax_allReps <- optList_tax_allReps_full
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_tax_allReps, function(x) x[[length(x) - 1]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), ylim=c(0,reps), main="Number of Replicates with a Taxonomic Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="orchid4", border="orchid1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=c(55, 0), ylim=c(0,reps), add=TRUE)

	par(mfrow=c(2,1), mar=c(3.5, 3.5, 1, 1))
	optList_bm_allReps <- optList_bm_allReps_heuristic
	optList_bm_allReps <- optList_bm_allReps_full
	### number of Replicates with a shift in that interval
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_bm_allReps, function(x) unique(x[[length(x) - 1]]$optBreaks)))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), ylim=c(0,reps), main="Number of Replicates with a Body Mass Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="firebrick4", border="firebrick1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=c(55, 0), ylim=c(0,reps), add=TRUE)


####################################################################################################################################
	quartz(width=6.89)
		par(mfrow=c(3,1), mar=c(0,4,0.5,0.5), mgp=c(2, 1,0))
		
		### isotope panel
		source("~/Dropbox/code/R/common_src/isotopes.R")
		optList_topes <- doTopesRateAnalysis(intervals)
		plotTopesRateAnalysis(optList_topes, intervals, x.axis=FALSE) #
		box(lwd=3)
		# getAlroyStatistics(intervals)
		
	### taxonomy panel
		par(mar=c(0,4,0,0.5))
		# for (intv in seq_len(nrow(intervals))) rawSp[[intv]] <- which(thisRanges[,"FO"] > intervals$ageTop[intv] & thisRanges[,"LO"] < intervals$ageBase[intv])
		
		# thisRich <- colSums(apply(countCube, c(1,2), mean, na.rm=TRUE))
		# thisRich <- colSums(apply(countCube, c(1,2), mean, na.rm=TRUE))
		# plot(rowMeans(intervals), thisRich, xlim=c(max(intervals), min(intervals)), type="n", ylab="Number of Species", xaxt="n", ylim=c(0,max(colSums(countBox[3,,])))) #sapply(rawSp, length)
		# overlayCzTimescale(do.subepochs=TRUE)
		# # lines(rowMeans(intervals), sapply(rawSp, length), col=alphaColor("black", 1.0), lwd=1, lty=2)
		# polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(colSums(countBox[1,,]), rev(colSums(countBox[3,,]))), col=alphaColor("green1", 0.25), border="green1")
		# lines(rowMeans(intervals), thisRich, col=alphaColor("green4", 1.0), lwd=3)
		
		prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
		dimnames(prop) <- list(rownames(intervals), shortFam)
		source("~/Dropbox/code/R/common_src/taxonomicEv.R")
		proportionalRichness(occs[occs$order %in% c("Artiodactyla","Perissodactyla"),], intervals, classCol="bigList", prop, age.determination="midpoint", range.through=TRUE, do.prop=FALSE, do.log=FALSE, overlay.labels=FALSE, legend=FALSE, xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)), do.parallel=FALSE)
		
		abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median) - 1]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median) - 1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median) - 1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median) - 1]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
		box(lwd=3)
	
	### body mass panel
		# quartz(width=10, height=6)
		# par(mar=c(4,4,2,2))
		par(mar=c(3,4,0,0.5))
		thisRanges <- getTaxonRangesFromOccs(occs=occs[occs$order %in% c("Artiodactyla","Perissodactyla"),], random=FALSE)
		thisMat[,c("FO","LO")] <- thisRanges[match(rownames(thisMat), rownames(thisRanges)),]
		# plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(55,5,10), xlab="Time (Ma)")
		plot(thisMat$FO, thisMat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(55,5,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75")
		rect(-10e6, -10e6, 10e6, 10e6, col="white")
		overlayCzTimescale(do.subepochs=TRUE)
		
		famColors <- rainbow(length(shortFam))
		colorList <- famColors[match(bigList$family[as.character(bigList$accepted_name) %in% names(thisCol)], shortFam)]
		colorList[is.na(colorList)] <- "gray25"
		
		orderColors <- array(NA, dim=nrow(thisMat))
		orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)]=="Perissodactyla"] <- "dodgerblue4"
		orderColors[bigList$order[match(rownames(thisMat), bigList$accepted_name)] =="Artiodactyla"] <- "deeppink4"
		for (i in seq_len(nrow(thisMat))) {
			# lines(x=c(this["FO"], x["LO"]), y=c(x["bodyMass"], x["bodyMass"]), lwd=3, pch=21, col=famColors[match(bigList[match(rownames(thisMat), bigList[,1]),2], shortFam)])
			# lines(x=c(thisMat$FO[i], thisMat$LO[i]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor(colorList[i], 0.75))
			# lines(x=c(thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"FO"], thisRanges[match(rownames(thisMat)[i], rownames(thisRanges)),"LO"]), y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor("gray0", 0.75)) #
			if (is.finite(thisMat$FO[i]) & is.finite(thisMat$LO[i]) & thisMat$FO[i] != thisMat$LO[i]) lines(x=thisMat[i,c("FO","LO")], y=c(thisMat$bodyMass[i], thisMat$bodyMass[i]), lwd=0.75, pch=21, col=alphaColor("gray0", 0.5)) #alphaColor(orderColors[i], 0.5)
		}
		points(thisMat[complete.cases(thisMat[ ,c("FO","LO")]) & thisMat$FO == thisMat$LO, c("FO", "bodyMass")], pch=21, col=alphaColor("gray0", 0.5), cex=0.25)
		
		# optList_bm_median <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, do.parallel=do.parallel)	# based on means
		# optList_bm_median <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=do.heuristic, do.parallel=do.parallel)	# based on median
		abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks,2], range(intervals))), lwd=1.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks,2])), "Ma"), adj=c(0,0),cex=0.5, col="firebrick4")
		
		polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
		polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
		lines(rowMeans(intervals), quants[3,], col=alphaColor("goldenrod1", 0.5), lwd=5)
		lines(rowMeans(intervals), quants[3,], col=alphaColor("darkorange4", 1.0), lwd=3)
		points(rowMeans(intervals), quants[3,], col=alphaColor("darkorange1", 0.5), cex=0.5)
		box(lwd=3)




##############

# # bm_stdev <- apply(sapply(repIntSp, function(y) sapply(y, function(x) sd(thisMat[x,"bodyMass"], na.rm=TRUE)), simplify = "array"), 1, median, na.rm=TRUE)
# topes <- getAlroyStatistics(intervals)
# topes$FD <- c(topes[1:(nrow(topes)-1),"midpoint"] - topes[2:nrow(topes),"midpoint"], NA)
# quartz(width=3.445, height=3.445)
# par(mfrow=c(2,1), mar=c(3,3,0,0.5), mgp=c(1.25, 0.5, 0), cex=0.75, cex.axis=0.5, cex.lab=0.75)
# # plot(rowMeans(intervals), -topes$midpoint, type="o", lwd=2, pch=21, col="dodgerblue4", bg="dodgerblue1", xlim=c(max(intervals), min(intervals)), cex.axis=0.75, xaxp=c(55,5,10), cex=0.75)
# # plot(rowMeans(intervals), abs(-topes$FD), type="o", lwd=2, pch=21, col="dodgerblue4", bg="dodgerblue1", xlim=c(max(intervals), min(intervals)), cex.axis=0.75, xaxp=c(55,5,10), cex=0.75)
# # abline(h=0.0, lty=3, col="gray50")

# # raw correlation
# # plot(topes$midpoint, quants[3,], pch=21, col="black", bg="gray50")
# # abline(lm(formula=quants[3,] ~ topes$midpoint), col="black")
# # cor.test(topes$midpoint, quants[3,], method="spearman")

# # first-differences
# FD_median <- c(quants[3,1:(ncol(quants) - 1)] - quants[3,2:ncol(quants)], NA)
# # plot(-topes$FD, FD_median, pch=21, col="black", bg="gray50", xlim=c(max(intervals), min(intervals)))
# # plot(topes$FD, FD_median, pch=21, col="black", bg="gray50", xlab=expression(paste(plain("first-differences "), delta^18, plain(O))), ylab="first-differences\nmedian body mass",)
# plot(topes$FD, FD_median, pch=21, col="dodgerblue4", bg="dodgerblue1", xlab=expression(paste(plain("first-differences "), delta^18, plain(O))), ylab="first-differences\nmedian body mass", fg="gray75", col.axis="gray75", col.lab="gray75")
# abline(lm(FD_median ~ topes$FD), col="black")
# cor.test(topes$FD, FD_median, method="spearman")

# # first differences bodyMass to that interval with stdev of that interval
# # plot(topes$stdev, FD_median, pch=21, col="black", bg="gray50", xlab=expression(paste(plain("Std. Dev. "), delta^18, plain(O))), ylab="first-differences\nmedian body mass")
# plot(topes$stdev, FD_median, pch=21, col="firebrick4", bg="firebrick1", xlab=expression(paste(plain("Std. Dev. "), delta^18, plain(O))), ylab="first-differences\nmedian body mass", fg="gray75", col.axis="gray75", col.lab="gray75")
# # plot(topes$stdev[2:nrow(topes)], FD_median[1:(length(FD_median)-1)], pch=21, col="black", bg="gray50", xlab=expression(paste(plain("Std. Dev. "), delta^18, plain(O))), ylab="first-differences\nmedian body mass")
# abline(lm(FD_median ~ topes$stdev), col="black")
# cor.test(topes$stdev, FD_median, method="spearman")

# # first differences bodyMass to that interval with stdev of that interval
# plot(topes$stdev, bm_stdev, pch=21, col="black", bg="gray50", xlab=expression(paste(plain("Std. Dev. "), delta^18, plain(O))), ylab="first-differences\nmedian body mass")
# # plot(topes$stdev[2:nrow(topes)], FD_median[1:(length(FD_median)-1)], pch=21, col="black", bg="gray50", xlab=expression(paste(plain("Std. Dev. "), delta^18, plain(O))), ylab="first-differences\nmedian body mass")
# abline(lm(FD_median ~ bm_stdev), col="black")
# # cor.test(topes$stdev, FD_median, method="spearman")

# # cbind(bigList[match(rownames(thisMat), bigList[,1]),2], rownames(thisMat))

# # quartz("Univariate Statistics", width=7, height=5.5)
# # bstatBox <- bstatBox[,colnames(pc2statBox)!="skewness",]
# # pc2statBox <- pc2statBox[,colnames(pc2statBox)!="skewness",]
# # pc3statBox <- pc3statBox[,colnames(pc3statBox)!="skewness",]
# # if (do.disparity) par(mfrow=c(3,dim(bstatBox)[2])) else par(mfrow=c(3,dim(bstatBox)[2]-1), mar=c(3,4,2,1)) 
# # plotUnivStats(bstatBox, intervals=intervals, dispMat=dispMat, new.window=FALSE, thisLab="log Body Mass (kg)")
# # plotUnivStats(pc2statBox, intervals=intervals, dispMat=dispMat, new.window=FALSE, thisLab="PC2 score")
# # plotUnivStats(pc3statBox, intervals=intervals, dispMat=dispMat, new.window=FALSE, thisLab="PC3 score")

# # # # thisCol <- which(colnames(thisMat)=="bodyMass")
# # # breaks <- bmBreaks
# # # thisCol <- which(colnames(thisMat)=="PC3")
# # # breaks <- pc3Breaks
# # optList_bm_allReps <- doHandleyTest_heuristic(dframe=thisMat, thisCol=thisCol, breaks=breaks, intList=listifyMatrixByRow(intervals), sig=0.01, do.parallel=do.parallel)
# # optList_bm_allReps

# #######################

require(vegan)
quartz(width=6.89, height=6)
bmBox <- sapply(repIntSp, function(thisRep) t(sapply(thisRep, function(taxa) hist(thisMat$bodyMass[taxa], freq=TRUE, xlim=c(min(bmBreaks), max(bmBreaks)), ylim=c(0,22), breaks=bmBreaks, col="firebrick1", ylab="", plot=plotHist)$counts)), simplify="array") #, main=paste(intervals[intv, "ageTop"], "-", intervals[intv, "ageBase"], "Ma (n = ", length(intSp[[intv]]), ")")
par(mfrow=c(2,1), mar=c(3,3,0.5,0.5), mgp=c(1.25, 0.5, 0))
tax <- plotNMDS(thisBox=t(apply(taxCube, c(1,2), mean, na.rm=TRUE)), intervals, polygon.ints=optList_tax_median[[(length(optList_tax_median) - 1)]]$optBreaks, scaler=5, title="taxonomy") #, filename="~/Desktop/NMDS_bodyMass.pdf"
bm <- plotNMDS(thisBox=bmBox, intervals, polygon.ints= optList_bm_median[[(length(optList_bm_median) - 1)]]$optBreaks, scaler=5, title="Body Mass") #, filename="~/Desktop/NMDS_bodyMass.pdf"
# # pc2<-plotNMDS(thisBox=pc2Box, intervals, scaler=5, title="PC2") #, filename="~/Desktop/NMDS_pc2.pdf"
# # pc3<-plotNMDS(thisBox=pc3Box, intervals, scaler=5, title="PC3") #, filename="~/Desktop/NMDS_pc3.pdf"


# plot(vegdist(t(apply(taxCube,c(1,2),mean,na.rm=TRUE))), vegdist(apply(bmBox,c(1,2),mean,na.rm=TRUE)))
# abline(a=0, b=1, lty=3)
# bm<-NMDSDist(bmBox)
# pc2<-NMDSDist(pc2Box)
# pc3<-NMDSDist(pc3Box)

	# quartz(width=13, height=4)
	# # pdf(file="NMDSDist.pdf", width=11, height=3)
	# par(mfrow=c(1,3))
	# plot(rowMeans(bm[,c("ageTop", "ageBase")]), bm[,5], xlim=c(max(intervals), min(intervals)), ylim=c(min(bm[,3], na.rm=TRUE), max(bm[,5], na.rm=TRUE)), type="n", main="Body Mass")
	# overlayCzTimescale()
	# polygon(x=c(rowMeans(bm[,c("ageTop", "ageBase")])[1:(nrow(bm)-1)], rowMeans(bm[,c("ageTop", "ageBase")])[(nrow(bm)-1):1]), y=c(bm[1:(nrow(bm)-1),3],bm[(nrow(bm)-1):1,5]), border=NA, col=rgb(t(col2rgb("firebrick")), alpha=64, maxColorValue=255))
	# lines(rowMeans(bm[,c("ageTop", "ageBase")]), bm[,4], col="firebrick", lwd=2, )
	 
	# plot(rowMeans(pc2[,c("ageTop", "ageBase")]), pc2[,5], xlim=c(max(intervals), min(intervals)), ylim=c(min(pc2[,3], na.rm=TRUE), max(pc2[,5], na.rm=TRUE)), type="n", main="PC2")
	# overlayCzTimescale()
	# polygon(x=c(rowMeans(pc2[,c("ageTop", "ageBase")])[1:(nrow(pc2)-1)], rowMeans(pc2[,c("ageTop", "ageBase")])[(nrow(pc2)-1):1]), y=c(pc2[1:(nrow(pc2)-1),3],pc2[(nrow(pc2)-1):1,5]), border=NA, col=rgb(t(col2rgb("forestgreen")), alpha=64, maxColorValue=255))
	# lines(rowMeans(pc2[,c("ageTop", "ageBase")]), pc2[,4], col="forestgreen", lwd=2, )
	
	# plot(rowMeans(pc3[,c("ageTop", "ageBase")]), pc3[,5], xlim=c(max(intervals), min(intervals)), ylim=c(min(pc3[,3], na.rm=TRUE), max(pc3[,5], na.rm=TRUE)), type="n", main="PC3")
	# overlayCzTimescale()
	# polygon(x=c(rowMeans(pc3[,c("ageTop", "ageBase")])[1:(nrow(pc3)-1)], rowMeans(pc3[,c("ageTop", "ageBase")])[(nrow(pc3)-1):1]), y=c(pc3[1:(nrow(pc3)-1),3],pc3[(nrow(pc3)-1):1,5]), border=NA, col=rgb(t(col2rgb("dodgerblue")), alpha=64, maxColorValue=255))
	# lines(rowMeans(pc3[,c("ageTop", "ageBase")]), pc3[,4], col="dodgerblue", lwd=2, )
	# # dev.off()
	
# quartz(width=14, height=11)
# par(mfrow=c(2,3), mar=c(4,4,2,0.5))
# box<-bmBox
# box<-pc2Box
# box<-pc3Box

# ksbm<-pairwiseKSTestsSubsequent(bmBox)
# kspc2<-pairwiseKSTestsSubsequent(pc2Box)
# kspc3<-pairwiseKSTestsSubsequent(pc3Box)

# x<-list()
# for (intv in seq_len(nrow(intervals))) {
	# x[[intv]]<-kruskal.test(list(thisMat$bodyMass[!is.na(thisMat$FO) & thisMat$FO<intervals$ageBase[intv]], thisMat$bodyMass[!is.na(thisMat$FO) & !is.na(thisMat$LO) & thisMat$LO	<intervals$ageBase[intv]]))
	# # x[[intv]]<-ks.test(thisMat$bodyMass[!is.na(thisMat[,"FO"]) & !is.na(thisMat[, "LO"]) & thisMat[,"FO"]>intervals[intv, "ageBase"]], thisMat$bodyMass[!is.na(thisMat[,"FO"]) & !is.na(thisMat[, "LO"]) & thisMat[,"LO"]<intervals[intv, "ageBase"]])
	# # x[[intv-1]]<-ks.test(thisMat$bodyMass[!is.na(thisMat[,"FO"]) & !is.na(thisMat[, "LO"]) & thisMat[,"FO"]>intervals[intv, "ageTop"] & thisMat[, "LO"]<intervals[intv, "ageBase"]], thisMat$bodyMass[!is.na(thisMat[,"FO"]) & !is.na(thisMat[, "LO"]) & thisMat[,"FO"]>intervals[intv-1, "ageTop"] & thisMat[, "LO"]<intervals[intv-1, "ageBase"]])
# }
# cbind(intervals, sapply(x, function (xx) { xx$p.value } ))

# # allHist<-hist(thisMat$bodyMass, plot=FALSE)
# # thisECDF<-ecdf(allHist$counts)

# # require(MASS)
# # thisDist<-list()
# # thisDist<-c(thisDist, normal=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"normal"))
# # thisDist<-c(thisDist, cauchy=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"cauchy"))
# # thisDist<-c(thisDist, gamma1=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"gamma"))
# # thisDist<-c(thisDist, gamma2=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),dgamma,list(shape=1,rate=0.1)))
# # thisDist<-c(thisDist, sn=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),VGAM::dsnorm,list(location=mean(thisMat$bodyMass], na.rm=TRUE), scale=var(thisMat$bodyMass], na.rm=TRUE), shape=0)))
# # thisDist<-c(thisDist, Poisson=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"Poisson"))
# # thisDist<-c(thisDist, geometric=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"geometric"))
# # thisDist<-c(thisDist, logistic=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"logistic"))
# # thisDist<-c(thisDist, weibull=fitdistr(scale(thisMat$bodyMass[is.finite(thisMat$bodyMass)]),"weibull"))

# taxList<-rownames(pcaLo$x)[rownames(pcaLo$x)%in%rownames(pcaUp$x)]
# x<-pcaLo$x[taxList,1]
# y<-pcaUp$x[taxList,1]
# x <- thisMat$bodyMass
# # pc1Vec <- pcaAll$x[,1]
# # pc1Vec <- c(pc1Vec,pcaLo$x[!rownames(pcaLo$x) %in% names(pc1Vec),1])
# # pc1Vec <- c(pc1Vec,pcaUp$x[!rownames(pcaUp$x) %in% names(pc1Vec),1])
# # y<- -pc1Vec[match(rownames(x), names(pc1Vec))]
# # cbind(x,y)
# # plot(x, y)

# taxList <- rownames(pcaLo$x)[rownames(pcaLo$x) %in% rownames(pcaUp$x)]
# x <- pcaLo$x[taxList,3]
# y <- pcaUp$x[taxList,3]
# plot(x, y)

# plot(x, -pcaAll$x[match(rownames(x), rownames(pcaAll$x)),1], col="gray50", xlab="thisMat$bodyMass", ylab="PC1")
# points(x, -pcaLo$x[match(rownames(x), rownames(pcaLo$x)),1], col="blue")
# points(x, -pcaUp$x[match(rownames(x), rownames(pcaUp$x)),1], col="red")
# abline(h=0, lty=3)
# abline(v=mean(x, na.rm=TRUE), lty=3)
# # abline(a=0, b=1, lty=2)
# xylm <- lm(-pcaUp$x[match(rownames(x), rownames(pcaUp$x)),1] ~ x)
# abline(lm(-pcaAll$x[match(rownames(x), rownames(pcaAll$x)),1] ~ thisMat$bodyMass))
# abline(lm(-pcaUp$x[match(rownames(x), rownames(pcaUp$x)),1] ~ thisMat$bodyMass), col="blue")
# abline(lm(-pcaLo$x[match(rownames(x), rownames(pcaLo$x)),1] ~ thisMat$bodyMass), col="red")
# legend("topleft", fill=c("gray50", "blue", "red"), legend=c("both", "lower", "upper"))
# cor.test(x, y)
# # hist(xylm$residuals)
# thisMat[names(xylm$residuals[xylm$residuals< -1]),]

# # plot(sapply(names(pres), function(x) mean(thisMat$bodyMass[grep(x, rownames(thisMat))], na.rm=TRUE)), pres)
# # cor.test(sapply(names(pres), function(x) mean(thisMat$bodyMass[grep(x, rownames(thisMat))], na.rm=TRUE)), pres, method="spearman")
# xylm$residuals[order(xylm$residuals)]
