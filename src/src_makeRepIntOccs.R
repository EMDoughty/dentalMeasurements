################################################################################################################################################
#### get occurrences within intervals
################################################################################################################################################

getRepIntOccs <- function(settings, 
						intervals, 
						return.box=TRUE, 
						save.to.file=TRUE, 
						file.path="~/Desktop/EcologyResults/", 
						file.name="coreData") {
	
	if (bootstrapSpecies) holderMat <- measure.mat
	
	repIntOccs <- list()

	for (rep in seq_len(n.reps)) {
		cat("Beginning Rep", rep, "of", n.reps, "...\r")
		##################################################We need to update this bootstrap section
		if (bootstrapSpecimens) {
			measure.mat <- specimenMat[sample.int(nrow(specimenMat), size=nrow(specimenMat), replace=TRUE),]
			measure.mat <- aggregate(measure.mat, by=list(taxon=specimenMat$taxon), mean, na.rm=TRUE)
			rownames(measure.mat) <- measure.mat$taxon
			measure.mat[sapply(measure.mat, is.nan)] <- NA
			measure.mat[,"reg"] <- as.character(famList$reg[match(measure.mat$taxon,famList$taxon)])
			measure.mat[,"bodyMass"] <- makeBodyMasses(measure.mat, regList, best.only=TRUE)
		}
		if (bootstrapSpecies) measure.mat <- holderMat[sample.int(n=nrow(measure.mat), size=nrow(measure.mat), replace=TRUE),]
	
		col.dates <- getCollectionAgesFromOccs(occs=occs[, c("collection_no", "max_ma", "min_ma")], random=TRUE)
		occDates <- col.dates$collection_age[match(occs$collection_no, col.dates$collection_no)]
		intOccs <- apply(intervals, 1, function(thisIntv) occs$occurrence_no[occDates > thisIntv[1] & occDates <= thisIntv[2]]) # greater than ageTop, less than or equal to ageBase
		# intTaxa <- sapply(intOccs, function(x) unique(occs$accepted_name[occs$occurrence_no %in% x]))
		# x <- intOccs
		# intSp <- sapply(intOccs, function(x) match(sort(unique(gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name[occs$accepted_rank  == "species" & occs$occurrence_no %in% x]))), measure.mat$taxon))
		# intSp <- sapply(intOccs, function(x) match(sort(unique(occs$accepted_name[occs$accepted_rank  %in% c("genus", "species") & occs$occurrence_no %in% x]))), measure.mat$taxon))
		# intSp <- sapply(intOccs, function(x) sort(unique(as.character(occs$accepted_name[occs$occurrence_no %in% x]))))
	
		#which(occs$occurrence_no %in% x == TRUE) # none are being returned as TRUE
		
		if (do.subsample) { 
			nOccs <- sapply(intOccs, length)
		 	# nTaxa <- sapply(intSp, length)			### if you want to set the quota no lower than the maximum number of SIB taxa; intSp is required for this to work, so has to be done above
		 	nTaxa <- 0									### set to zero to simply set the quota to the minimum number of occurrences
			quota <- max(c(max(nTaxa), min(nOccs)))		### quota is either the maximum number of observed taxa, or the minimum number of occurrences
			cat("Subsampling quota set to", quota, "occurrences")
	
			intOccs <- lapply(X=intOccs, FUN=sample, size=quota)
		}
		repIntOccs[[rep]] <- intOccs 
	}
	
###################################################################################################################################
	repIntTaxa <- getRepIntTaxaFromRepIntOccs(repIntOccs, this.rank=this.rank, do.rangethrough=settings$do.rangethrough)
	
	# all.taxa <- all.taxa[all.taxa %in% as.character(bigList$accepted_name)]
	
	taxRangeCube <- sapply(repIntTaxa, getIntRangesOneRep, all.taxa, simplify="array")
	taxRangeBox <- array(data=FALSE, dim=c(nrow(taxRangeCube), nrow(intervals), settings$n.reps), dimnames=list(rownames(taxRangeCube), rownames(intervals)))
	for(this.rep in seq_len(settings$n.reps)) {
		for (this.taxon in seq_len(nrow(taxRangeCube))) {
			if (all(is.finite(taxRangeCube[this.taxon,1,this.rep]) & is.finite(taxRangeCube[this.taxon,2,this.rep]))) taxRangeBox[this.taxon, taxRangeCube[this.taxon,1,this.rep]:taxRangeCube[this.taxon,2,this.rep], this.rep] <- TRUE
		}
	}

	if (save.to.file) {
		if(Sys.info()["sysname"] == "Darwin"){
			save(settings, occs, measure.mat, taxRangeBox, file=paste0(file.path, file.name, "_standardized=", do.subsample, timestamp(),".Rdata"))
			#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
		} 
		# else if(Sys.info()["sysname"] == "Windows"){
			# save(repIntOccs, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/repIntOccs_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
			# # load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
		# }
		
		print("**** repIntOccs saved to file...")
	}
	
	if (return.box) return(repIntOccs)
}

