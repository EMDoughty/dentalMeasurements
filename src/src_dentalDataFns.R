require(compiler)
# require(MASS)
source("~/Dropbox/code/R/common_src/occFns.R", chdir = TRUE)

getSpecimenMatFromAmandaMeasurements <- function(filename = "~/Dropbox/code/R/dentalMeasurements/dat/amanda_specimens.csv") {
	require(abind)
	dat <- read.csv(filename, strip.white = TRUE)

	measure.labels <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W", "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w", "upP", "upM", "loP", "loM")
	# specimen.mat <- array(dim = c(0, length(measure.labels)), dimnames = list(NULL, measure.labels))
	specimen.vec <- unique(dat[,c("species", "Specimen.no.")])
	specimen.mat <- t(sapply(X=specimen.vec$Specimen.no., FUN=function(x) colMeans(dat[dat$Specimen.no.==x, measure.labels], na.rm=TRUE)))
	specimen.mat[!is.finite(specimen.mat)] <- NA
	specimen.mat <- data.frame(identified.name=specimen.vec$species, specimen=specimen.vec$Specimen.no., specimen.mat)

	# cube <- abind::abind(dat[seq(from = 1, to = nrow(dat), by = 3), measure.labels], dat[seq(from = 2, to = nrow(dat), by = 3), measure.labels], along = 3) 
	# cube <- abind::abind(cube, dat[seq(from = 3, to = nrow(dat), by = 3), measure.labels], along = 3)
	# replicatemeasure.mat <- apply(cube, c(1, 2), mean, na.rm = TRUE)
	# rownames(replicatemeasure.mat) <- dat$species[seq(from = 1, to = nrow(dat), by = 3)]
	# replicatemeasure.mat[!is.finite(replicatemeasure.mat)] <- NA
	# # replicatemeasure.mat<-data.frame(species=dat$species[seq(from=1, to=nrow(dat), by=3)], specimen=dat$Specimen.no.[seq(from=1, to=nrow(dat), by=3)], replicatemeasure.mat, stringsAsFactors=FALSE, row.names=NULL)
	
	# theseSpecimenNos <- dat$Specimen.no.[seq(from = 1, to = nrow(dat), by = 3)]
	# specimenNos <- unique(theseSpecimenNos)
	# specimen.mat <- matrix(nrow = 0, ncol = ncol(replicatemeasure.mat))
	# for (i in seq_len(length(specimenNos))) {
		# index <- which(theseSpecimenNos == specimenNos[i])
		# thisSpecimen <- replicatemeasure.mat[index, ]
		# if (!is.null(dim(thisSpecimen))) 
			# thisSpecimen <- apply(thisSpecimen, 2, mean, na.rm = TRUE)
		# specimen.mat <- rbind(specimen.mat, thisSpecimen)
		# rownames(specimen.mat)[nrow(specimen.mat)] <- rownames(replicatemeasure.mat)[index[1]]
	# }
	# specimen.mat <- data.frame(species = rownames(specimen.mat), specimen = specimenNos, specimen.mat, stringsAsFactors = FALSE, row.names = NULL)
	# specimen.mat
}

getBlastoMeasuresOneSpecimen <- function(this.block) {
	this.dim <- unique(this.block$dim)
	this.m <- sapply(this.dim, function(x) mean(this.block$length[this.block$dim == x], na.rm = TRUE))
	data.frame(matrix(this.m, nrow = 1, dimnames = list(NULL, this.dim)))
}

getSpecimenMatFromBlastoMeasurements <- function(filename="~/Dropbox/code/R/dentalMeasurements/dat/blasto_Birlenbach20140207.csv") {
	dat <- read.csv(filename, strip.white = TRUE)
	specimen.vec <- sort(unique(dat$specimen))
	this.dim <- unique(dat$dim)
	measure.labels <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W", "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w", "dp3_l", "dp3_w", "dp4_l", "dp4_w")
	specimen.mat <- array(dim = c(0, length(measure.labels)), dimnames = list(NULL, measure.labels))
	for (i in seq_along(specimen.vec)) specimen.mat <- merge(getBlastoMeasuresOneSpecimen(dat[dat$specimen == specimen.vec[i], ]), specimen.mat, all = TRUE, sort=FALSE)
	specimen.mat <- data.matrix(specimen.mat)
	specimen.mat[!is.finite(specimen.mat)] <- NA
	specimen.mat <- data.frame(specimen=specimen.vec, specimen.mat)

	specimen.mat <- merge(read.csv("~/Dropbox/code/R/dentalMeasurements/dat/blasto_info2.csv", strip.white = TRUE), specimen.mat, by = "specimen", all = TRUE, sort=FALSE)
	colnames(specimen.mat)[colnames(specimen.mat) == "sp_current"] <- "identified.name"
	specimen.mat
}

getSpecimenMatFromLiteratureMeasurements <- function(filename = "~/Dropbox/code/R/dentalMeasurements/dat/ungulate_literature.csv") {
	dat <- read.csv(filename, strip.white = TRUE)
	# dat <- dat[apply(is.finite(data.matrix(dat[,3:ncol(dat)])), 1, any),] # removes taxa with no measurements
	dat
}

# makeOneSpeciesMatFromSpecimenMat <- function(specimen.mat) {
	# # species <-unique(specimen.mat$species)
	# oneSpeciesMat <- aggregate(specimen.mat, by = list(specimen.mat$species), mean, na.rm = TRUE)
	# # oneSpeciesMat <- aggregate(specimen.mat, by=list(specimen.mat$species), median, na.rm=TRUE)
	# oneSpeciesMat <- data.frame(oneSpeciesMat, row.names = oneSpeciesMat[, 1])
	# oneSpeciesMat[sapply(oneSpeciesMat, is.nan)] <- NA
	# colnames(oneSpeciesMat)[1] <- "taxon"
	# oneSpeciesMat[, -(c(2:3, which(colnames(oneSpeciesMat) %in% c("Published.name", "n", "Reference", "PaleoDB.ref", "Notes", "X", "X.1", "X.2", "X.3", 
		# "X.4"))))]
# }

makeOneGenusMatFromSpecimenMat <- function(measure.mat) {
	#need to make it so only the genus is being used to aggregate as some entries have a subgenus included within the genus field in pbdb. 10/6/2023 It seems that this is not an issue with measure.mat but in occs and repIntTaxa
  oneGenusMat <- aggregate(measure.mat, by = list(measure.mat$genus), mean, na.rm = TRUE)		# the column "genus" is appended with the reg.vec, so might need a better way of getting this...
	# oneGenusMat <- aggregate(measure.mat, by=list(measure.mat$genus), median, na.rm=TRUE)
	oneGenusMat <- data.frame(oneGenusMat, row.names = oneGenusMat[, 1])
	oneGenusMat[sapply(oneGenusMat, is.nan)] <- NA
	colnames(oneGenusMat)[1] <- "taxon"
	oneGenusMat[,-2]
}

getToothRowLengths <- function(taxon) {
	thisRow <- vector(mode = "numeric", length = 0)
	if (!is.nan(taxon$upP)) {
		thisRow <- c(thisRow, upP = taxon$upP)
	} else {
		thisList <- c(taxon$P2_L, taxon$P3_L, taxon$P4_L)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, upP = sum(thisList))
		else thisRow <- c(thisRow, upP = NA)
	}

	if (!is.nan(taxon$upM)) {
		thisRow <- c(thisRow, upM = taxon$upM)
	} else {
		thisList <- c(taxon$M1_L, taxon$M2_L, taxon$M3_L)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, upM = sum(thisList))
		else thisRow <- c(thisRow, upM = NA)
	}

	if (!is.nan(taxon$loP)) {
		thisRow <- c(thisRow, loP = taxon$loP)
	} else {
		thisList <- c(taxon$p2_l, taxon$p3_l, taxon$p4_l)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, loP = sum(thisList))
		else thisRow <- c(thisRow, loP = NA)
	}

	if (!is.nan(taxon$loM)) {
		thisRow <- c(thisRow, loM = taxon$loM)
	} else {
		thisList <- c(taxon$m1_l, taxon$m2_l, taxon$m3_l)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, loM = sum(thisList))
		else thisRow <- c(thisRow, loM = NA)
	}
	return(thisRow)
}

getExantMat <- function() {
	lab <- c("p2", "p3", "p4", "m1", "m2", "m3")
#	widths <- read.csv("~/Dropbox/code/R/cooperExtantDental/dat/widths_cooper.csv", stringsAsFactors=FALSE, strip.white=TRUE) original line
	widths <- read.csv("~/Dropbox/code/R/dentalMeasurements/dat/widths_cooper.csv", stringsAsFactors=FALSE, strip.white=TRUE) # changed this line so it matches with directory on github
#	widths <- merge(widths, read.csv("~/Dropbox/code/R/cooperExtantDental/dat/widths_roberto.csv", stringsAsFactors=FALSE, strip.white=TRUE), all=TRUE, sort=FALSE) original line
	widths <- merge(widths, read.csv("~/Dropbox/code/R/dentalMeasurements/dat/widths_roberto.csv", stringsAsFactors=FALSE, strip.white=TRUE), all=TRUE, sort=FALSE) # changed this line so it matches with directory on github
	w_means <-sapply(lab, function(x) rowMeans(widths[,grepl(x, colnames(widths))], na.rm=TRUE))
	w_means[!is.finite(data.matrix(w_means))] <- NA
	w_means <- data.frame(specimen=widths[,1], w_means)
	w_means <- aggregate(w_means, by=list(w_means$specimen), FUN=mean, na.rm=TRUE)
	w_means <- w_means[,-2]
	colnames(w_means) <- sapply(colnames(w_means), function(x) paste(x, "_w", sep=""))
	colnames(w_means)[1] <- "specimen"
	
#	lengths <- read.csv("~/Dropbox/code/R/cooperExtantDental/dat/lengths_cooper.csv", stringsAsFactors=FALSE, strip.white=TRUE) original line
	lengths <- read.csv("~/Dropbox/Code/R/dentalMeasurements/dat/lengths_cooper.csv", stringsAsFactors=FALSE, strip.white=TRUE) # changed this line so it matches with directory on github
#	lengths <- merge(lengths, read.csv("~/Dropbox/code/R/cooperExtantDental/dat/lengths_roberto.csv", stringsAsFactors=FALSE, strip.white=TRUE), all=TRUE, sort=FALSE) original line	
	lengths <- merge(lengths, read.csv("~/Dropbox/code/R/dentalMeasurements/dat/lengths_roberto.csv", stringsAsFactors=FALSE, strip.white=TRUE), all=TRUE, sort=FALSE) # changed this line so it matches with directory on github
	l_means <- data.frame(sapply(lab, function(x) rowMeans(lengths[,grepl(x, colnames(lengths))], na.rm=TRUE)))
	l_means <- cbind(specimen=lengths[,1], l_means)
	w_means[is.nan(data.matrix(w_means))] <- NA
	l_means <- aggregate(l_means, by=list(l_means[,1]), FUN=mean, na.rm=TRUE)
	l_means <- l_means[,-2]
	colnames(l_means) <- sapply(colnames(l_means), function(x) paste(x, "_l", sep=""))
	colnames(l_means)[1] <- "specimen"
	
	extant.mat <- merge(w_means, l_means, by="specimen", all=TRUE, sort=FALSE)
	# extant.mat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = widths$species[match(extant.mat$specimen, widths$specimen)]) 
	extant.mat$identified.name <- widths$species[match(extant.mat$specimen, widths$specimen)] 

	extant.mat
}

getArchaicMat <- function() { 
  #changed 2023_9_6 to added if statement as I (Evan) added identified.name as a field when cleaning the archaic dataset
  #2023_9_26 reverted to original code (Evan) since using the identified "verbatim" name will fail to keep the taxonomy of partial synonyms correct.  Must use the accepted name in its place.
	dat <- read.csv("~/Dropbox/code/R/dentalMeasurements/dat/ArchaicUngulate_UploadFile_Master.csv", stringsAsFactors=TRUE, strip.white=TRUE)
	
	dat$verbatim.name <- dat$identified.name
	if(any(colnames(dat) %in% "accepted_name"))
	{
	 dat$identified.name <- gsub("_", " ", dat$accepted_name)
	} else {
	 dat$identified.name  <- apply(X=dat, MARGIN=1, FUN=function(x) if (x["Accepted.Species"]=="") x["Accepted.Genus"] else paste(x["Accepted.Genus"], x["Accepted.Species"], sep=" "))
	}
	
	names(dat)[names(dat)=="Catalog.Number"] <- "specimen"
	dat
}

############################################################################################################################################

#function is meant to be run prior to analysis to bring all measurement datasets together into a single entity
getSingleSpeciesMatrix <- function(append.archaic=TRUE, append.extant=TRUE) {
	#compile and lable dental measurments for specimens
	specimen.mat <- getSpecimenMatFromAmandaMeasurements()
	specimen.mat <- merge(x = specimen.mat, y = getSpecimenMatFromBlastoMeasurements(), all = TRUE, sort=FALSE)
	specimen.mat <- merge(x = specimen.mat, y = getSpecimenMatFromLiteratureMeasurements(), all = TRUE, sort=FALSE)
	specimen.mat[sapply(specimen.mat, is.nan)] <- NA
		
	if (append.archaic) specimen.mat <- merge(specimen.mat, getArchaicMat(), all=TRUE, sort=FALSE)
	if (append.extant) specimen.mat <- merge(specimen.mat, getExantMat(), all=TRUE, sort=FALSE)

	options(timeout=300)
	specimen.mat$species <- getCurrentTaxa(tax.vec = specimen.mat$identified.name)
	specimen.mat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = specimen.mat$species)

	############################
	#### note that specimens are aggregated by their medians, so as to minimize the effect of outlier measurements
	############################
	
	upLabels <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W", "upM") #\P2_L\",\"P2_W\","
	loLabels <- casefold(upLabels)
	loLabels[loLabels=="upm"] <- "loM"

	# measure.mat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(taxon = specimen.mat$species), mean, na.rm = TRUE)
	measure.mat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(taxon = specimen.mat$species), median, na.rm = TRUE)

	measure.mat[, sapply(measure.mat, is.numeric)] <- measure.mat[, sapply(measure.mat, is.numeric)]/10 #converts mm measurements to cm for compatibility with Janis & Damuth (1990) regressions
	# measure.mat <- transform(measure.mat, p4_a = p4_l * p4_w, m1_a = m1_l * m1_w, m2_a = m2_l * m2_w, m3_a = m3_l * m3_w, M2_A = M2_L * M2_W)
	measure.mat[sapply(measure.mat, is.nan)] <- NA

	rownames(measure.mat) <- measure.mat$taxon

	return(measure.mat)
}

appendMissingPaleoDBSpecies <- function(measure.mat, tax.vec) {
	# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
	tax.vec <- tax.vec[!tax.vec %in% measure.mat$taxon]
	tax.frame <- data.frame(array(NA, dim = c(length(tax.vec), ncol(measure.mat)), dimnames = list(tax.vec, colnames(measure.mat))))
	tax.frame$taxon <- tax.vec
	if (any(tax.frame$taxon %in% measure.mat$taxon)) {
		# measure.mat[match(rangesNotmeasure.mat$taxon, measure.mat$taxon),c("FO", "LO")] <- ranges[,c("FO", "LO")]
		} else measure.mat <- merge(measure.mat, tax.frame, all = TRUE, sort = FALSE)
	rownames(measure.mat) <- measure.mat$taxon
	measure.mat
}

getMeasureMatWithBodyMasses <- function(settings, append.archaic=TRUE, append.extant=TRUE) {
	print("Building measurement matrix...")
	measure.mat <- getSingleSpeciesMatrix(append.archaic = append.archaic, append.extant = append.extant)

	tax.vec <- sort(unique(occs$accepted_name[occs$accepted_rank %in% c("genus", "species")]))
	tax.vec <- tax.vec[!tax.vec %in% rownames(measure.mat)]
	measure.mat <- appendMissingPaleoDBSpecies(measure.mat, tax.vec) # species missing from the measurements are appended, and will receive a body mass estimate based on their cogeners

	source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
	measure.mat <- appendRegressionCategories(measure.mat = measure.mat, regMat=read.csv(file="~/Dropbox/code/R/dentalMeasurements/dat/regressionLabelsJDM.csv"), focal.tax=settings$focal.tax)
	measure.mat <- approxBodyMass(measure.mat = measure.mat)
	measure.mat <- measure.mat[is.finite(measure.mat$bodyMass),]		#### clears taxa without body mass estimate

	return(measure.mat)
}
