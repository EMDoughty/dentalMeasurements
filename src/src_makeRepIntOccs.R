source("~/Dropbox/code/R/common_src/strat.R")
### occs, interval must be global

################################################################################################################################################
#### get occurrences within intervals
################################################################################################################################################

dateOccsOneRep <- function(mclapply.dummy=0, draw.random.dates=TRUE) {
	col.dates <- getCollectionAgesFromOccs(occs=unique(occs[, c("collection_no", "max_ma", "min_ma")]), random=draw.random.dates)
	occ.dates <- col.dates$collection_age[match(occs$collection_no, col.dates$collection_no)]
	occ.dates
}

getRepIntOccsOneRep <- function(settings, occ.dates, draw.random.dates=TRUE) {
	intOccs <- apply(intervals, 1, function(thisIntv) occs$occurrence_no[occ.dates > thisIntv[1] & occ.dates <= thisIntv[2]]) # greater than ageTop, less than or equal to ageBase
	return(intOccs)
}

getRepIntOccsAllReps <- function(settings,
						# occs, 
						# intervals, 
						occ.dates.from.file=NULL,
						return.box=TRUE, 
						save.to.file=TRUE, 
						file.path="~/Desktop/EcologyResults/", 
						file.name="repIntOccs") {
	repIntOccs <- list()
	
	#### occs.dates is now required before calling getRepIntOccsOneRep. 
	#### occs.dates.from.file should be a character string of the name of a Rdata file including an object, repOccDates
	#### if occs.dates.from.file is missing, will generate a new repIntOccs

	if (is.null(occ.dates.from.file)) { 
		print("Generating random occurrence dates...")
		if (settings$do.parallel) { repOccDates <- simplify2array(mclapply(X=seq_len(settings$n.reps), FUN=dateOccsOneRep, mc.cores=detectCores()-2))
		} else repOccDates <- replicate(n=settings$n.reps, expr=dateOccsOneRep())
		if (save.to.file) save(repOccDates, file=paste0(file.path, "repOccDates_", gsub(pattern=" ", replacement="_", x=date()),".Rdata"))
		print("Done...")
	} else {
		paste("Reading occurrence dates from file:", occ.dates.from.file, "...\n")
		load(occ.dates.from.file)
	}

	for (rep in seq_len(settings$n.reps)) {
		cat("Beginning rep", rep, "of", settings$n.reps, "...\r")
		################################################## We need to update this bootstrap section
		repIntOccs[[rep]] <- getRepIntOccsOneRep(settings, occ.dates=repOccDates[,rep])
	} 
	
	if (save.to.file) {
		if(Sys.info()["sysname"] == "Darwin"){
			save(settings, repIntOccs, file=paste0(file.path, file.name, "_", gsub(pattern=" ", replacement="_", x=date()),".Rdata"))
		} 
		print("**** repIntOccs saved to file...")
	}
	print("Done...")
	
	if (return.box) return(repIntOccs)
}

subsampleRepIntOccsOneRep <- function(intOccs) {
		n.occs <- sapply(intOccs, length)
	 	# nTaxa <- sapply(intSp, length)			### if you want to set the quota no lower than the maximum number of SIB taxa; intSp is required for this to work, so has to be done above
	 	nTaxa <- 0									### set to zero to simply set the quota to the minimum number of occurrences
		quota <- max(c(max(nTaxa), min(n.occs)))		### quota is either the maximum number of observed taxa, or the minimum number of occurrences
		cat("Subsampling quota set to", quota, "occurrences\r")

		intOccs <- lapply(X=intOccs, FUN=sample, size=quota, replace=FALSE)
}

subsampleRepIntOccs <- function(repIntOccs) {
	lapply(repIntOccs, subsampleRepIntOccsOneRep)
}

##################################################################################################################################
################################################################################################################
#### Marcot 2019 06 18
################################################################################################################
makeRangeThroughOneRep <- function(intTaxa) {
	sampled.taxa <- sort(unique(unlist(intTaxa)))					## vector of all sampled taxa in this rep
	# sampled.taxa <- sampled.taxa[grep(pattern="_", x=sampled.taxa)]		## only those at the species level - do not want to range through occurrences at higher taxonomic ranks
	first.last <- t(sapply(sampled.taxa, function(x) array(data=range(which(sapply(intTaxa, function(y, x) x %in% y, x=x))), dimnames=list(c("LI", "FI")))))
	first.last <- data.frame(taxon=rownames(first.last), first.last, stringsAsFactors=FALSE)
	for (this.taxon in seq_along(sampled.taxa)) {
		for (this.intv in seq(from=first.last$LI[this.taxon], to=first.last$FI[this.taxon])) {								# for all intervals in which this.taxon should range-through
			if (!first.last$taxon[this.taxon] %in% intTaxa[[this.intv]]) intTaxa[[this.intv]] <- c(intTaxa[[this.intv]], first.last$taxon[this.taxon])	# if it is not currently in the interval, add it
		}
	}
	intTaxa
}

getIntTaxaFromOneRepIntOccs <- function(this.repIntOccs, settings) {
	if (settings$this.rank=="species") {
		intTaxa <- lapply(this.repIntOccs, function(x) sort(unique(as.character(occs$accepted_name[occs$occurrence_no %in% x & occs$accepted_rank=="species"]))))
	} else if (settings$this.rank=="genus") {
		intTaxa <- lapply(this.repIntOccs, function(x) sort(unique(as.character(occs$genus[occs$occurrence_no %in% x & occs$accepted_rank %in% c("genus", "species")]))))
	}

	if (settings$do.rangethrough) intTaxa <- makeRangeThroughOneRep(intTaxa)
	intTaxa
}

getRepIntTaxaFromRepIntOccs <- function(repIntOccs, 
										settings,
										force.rangethrough=FALSE,
										return.box=TRUE, 
										save.to.file=TRUE, 
										file.path="~/Desktop/EcologyResults/", 
										file.name="repIntTaxa") {

	if (force.rangethrough) settings$do.rangethrough <- TRUE
	if (settings$do.parallel) { repIntTaxa <- mclapply(repIntOccs, getIntTaxaFromOneRepIntOccs, settings, mc.cores=detectCores()/2) 
	} else repIntTaxa <- lapply(repIntOccs, getIntTaxaFromOneRepIntOccs, settings)
	
	if (save.to.file) {
		if(Sys.info()["sysname"] == "Darwin"){
			save(settings, repIntTaxa, file=paste0(file.path, file.name, "_standardized=", settings$do.subsample, "_rangethrough=", settings$do.rangethrough, "_", gsub(pattern=" ", replacement="_", x=date()),".Rdata"))
		} 
		print("**** repIntTaxa saved to file...")		
	}
	if (return.box) return(repIntTaxa)
}

################################################################################################################

	
