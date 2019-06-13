require(compiler)
# require(MASS)
source("~/Dropbox/code/R/common_src/occFns.R", chdir = TRUE)

getSpecimenMatFromAmandaMeasurements <- function(filename = "~/Dropbox/code/R/dentalMeasurements/dat/amanda_specimens.csv") {
	require(abind)
	datMeasures <- read.csv(filename, strip.white = TRUE)

	m <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W", "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", 
		"p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w", "upP", "upM", "loP", "loM")
	cube <- abind::abind(datMeasures[seq(from = 1, to = nrow(datMeasures), by = 3), m], datMeasures[seq(from = 2, to = nrow(datMeasures), by = 3), 
		m], along = 3)
	cube <- abind::abind(cube, datMeasures[seq(from = 3, to = nrow(datMeasures), by = 3), m], along = 3)
	replicateMeasureMat <- apply(cube, c(1, 2), mean, na.rm = TRUE)
	rownames(replicateMeasureMat) <- datMeasures$species[seq(from = 1, to = nrow(datMeasures), by = 3)]
	replicateMeasureMat[!is.finite(replicateMeasureMat)] <- NA
	# replicateMeasureMat<-data.frame(species=datMeasures$species[seq(from=1, to=nrow(datMeasures), by=3)], specimen=datMeasures$Specimen.no.[seq(from=1, to=nrow(datMeasures), by=3)], replicateMeasureMat, stringsAsFactors=FALSE, row.names=NULL)
	
	theseSpecimenNos <- datMeasures$Specimen.no.[seq(from = 1, to = nrow(datMeasures), by = 3)]
	specimenNos <- unique(theseSpecimenNos)
	specimen.mat <- matrix(nrow = 0, ncol = ncol(replicateMeasureMat))
	for (i in seq_len(length(specimenNos))) {
		index <- which(theseSpecimenNos == specimenNos[i])
		thisSpecimen <- replicateMeasureMat[index, ]
		if (!is.null(dim(thisSpecimen))) 
			thisSpecimen <- apply(thisSpecimen, 2, mean, na.rm = TRUE)
		specimen.mat <- rbind(specimen.mat, thisSpecimen)
		rownames(specimen.mat)[nrow(specimen.mat)] <- rownames(replicateMeasureMat)[index[1]]
	}
	specimen.mat <- data.frame(species = rownames(specimen.mat), specimen = specimenNos, specimen.mat, stringsAsFactors = FALSE, row.names = NULL)
	specimen.mat
}

getBlastoMeasuresOneSpecimen <- function(this.block) {
	this.dim <- unique(this.block$dim)
	this.m <- sapply(this.dim, function(x) mean(this.block$length[this.block$dim == x], na.rm = TRUE))
	data.frame(matrix(this.m, nrow = 1, dimnames = list(NULL, this.dim)))
}

getSpecimenMatFromBlastoMeasurements <- function() {
	dat <- read.csv("~/Dropbox/code/R/dentalMeasurements/dat/blasto_Birlenbach20140207.csv", strip.white = TRUE)
	specimen.vec <- sort(unique(dat$specimen))
	this.dim <- unique(dat$dim)
	dim.names <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W", "p2_l", "p2_w", "p3_l", "p3_w", 
		"p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w", "dp3_l", "dp3_w", "dp4_l", "dp4_w")
	specimen.mat <- array(dim = c(0, length(dim.names)), dimnames = list(NULL, dim.names))
	for (i in seq_along(specimen.vec)) specimen.mat <- merge(getBlastoMeasuresOneSpecimen(dat[dat$specimen == specimen.vec[i], ]), specimen.mat, all = TRUE)
	specimen.mat <- data.frame(specimen = specimen.vec, specimen.mat)

	specimen.mat[!sapply(specimen.mat, is.finite)] <- NA
	specimen.mat <- merge(read.csv("~/Dropbox/code/R/dentalMeasurements/dat/blasto_info2.csv", strip.white = TRUE), specimen.mat, by = "specimen", 
		all = TRUE)
	colnames(specimen.mat)[colnames(specimen.mat) == "sp_current"] <- "species"
	specimen.mat
}

getSpecimenMatFromLiteratureMeasurements <- function(filename = "~/Dropbox/code/R/dentalMeasurements/dat/ungulate_literature.csv") {
	dat <- read.csv(filename, strip.white = TRUE)
	# dat <- dat[apply(is.finite(data.matrix(dat[,3:ncol(dat)])), 1, any),] # removes taxa with no measurements
	dat
}

makeOneSpeciesMatFromSpecimenMat <- function(specimen.mat) {
	# species<-unique(specimen.mat$species)
	oneSpeciesMat <- aggregate(specimen.mat, by = list(specimen.mat$species), mean, na.rm = TRUE)
	# oneSpeciesMat <- aggregate(specimen.mat, by=list(specimen.mat$species), median, na.rm=TRUE)
	oneSpeciesMat <- data.frame(oneSpeciesMat, row.names = oneSpeciesMat[, 1])
	oneSpeciesMat[sapply(oneSpeciesMat, is.nan)] <- NA
	colnames(oneSpeciesMat)[1] <- "species"
	oneSpeciesMat[, -(c(2:3, which(colnames(oneSpeciesMat) %in% c("Published.name", "n", "Reference", "PaleoDB.ref", "Notes", "X", "X.1", "X.2", "X.3", 
		"X.4"))))]
}

getToothRowLengths <- function(species) {
	thisRow <- vector(mode = "numeric", length = 0)
	if (!is.nan(species$upP)) {
		thisRow <- c(thisRow, upP = species$upP)
	} else {
		thisList <- c(species$P2_L, species$P3_L, species$P4_L)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, upP = sum(thisList))
		else thisRow <- c(thisRow, upP = NA)
	}

	if (!is.nan(species$upM)) {
		thisRow <- c(thisRow, upM = species$upM)
	} else {
		thisList <- c(species$M1_L, species$M2_L, species$M3_L)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, upM = sum(thisList))
		else thisRow <- c(thisRow, upM = NA)
	}

	if (!is.nan(species$loP)) {
		thisRow <- c(thisRow, loP = species$loP)
	} else {
		thisList <- c(species$p2_l, species$p3_l, species$p4_l)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, loP = sum(thisList))
		else thisRow <- c(thisRow, loP = NA)
	}

	if (!is.nan(species$loM)) {
		thisRow <- c(thisRow, loM = species$loM)
	} else {
		thisList <- c(species$m1_l, species$m2_l, species$m3_l)
		if (!any(is.nan(thisList))) 
			thisRow <- c(thisRow, loM = sum(thisList))
		else thisRow <- c(thisRow, loM = NA)
	}
	return(thisRow)
}

############################################################################################################################################

#function is meant to be run prior to analysis to bring all measurement datasets together into a single entity
getSingleSpeciesMatrix <- function() {
	#compile and lable dental measurments for specimens
	specimen.mat <- getSpecimenMatFromAmandaMeasurements()
	specimen.mat <- merge(x = specimen.mat, y = getSpecimenMatFromBlastoMeasurements(), all = TRUE)
	specimen.mat <- merge(x = specimen.mat, y = getSpecimenMatFromLiteratureMeasurements(), all = TRUE)

	specimen.mat[sapply(specimen.mat, is.nan)] <- NA
	specimen.mat$species <- getCurrentTaxa(tax.vec = specimen.mat$species)

	upLabels <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W") #\P2_L\",\"P2_W\","
	loLabels <- casefold(upLabels)

	### node that specimens are aggregated by their medians, so as to minimize the effect of outlier measurements
	thisMat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(species = specimen.mat$species), median, na.rm = TRUE)

	thisMat[, sapply(thisMat, is.numeric)] <- thisMat[, sapply(thisMat, is.numeric)]/10 #converts mm measurements to cm for compatibility with Janis regressions
	thisMat <- transform(thisMat, p4_a = p4_l * p4_w, m1_a = m1_l * m1_w, m2_a = m2_l * m2_w, m3_a = m3_l * m3_w, M2_A = M2_L * M2_W)
	thisMat[sapply(thisMat, is.nan)] <- NA

	#thisMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = thisMat$species)
	rownames(thisMat) <- thisMat$species

	return(thisMat)
}

appendMissingPaleoDBSpecies <- function(thisMat, tax.vec) {
	# this adds taxa that are in PaleoDB (i.e., occurrence data), but not in the measurement files
	tax.vec <- tax.vec[!tax.vec %in% thisMat$species]
	tax.frame <- data.frame(array(NA, dim = c(length(tax.vec), ncol(thisMat)), dimnames = list(tax.vec, colnames(thisMat))))
	tax.frame$species <- tax.vec
	if (any(tax.frame$species %in% thisMat$species)) {
		# thisMat[match(rangesNotThisMat$species, thisMat$species),c("FO", "LO")] <- ranges[,c("FO", "LO")]
		} else thisMat <- merge(thisMat, tax.frame, all = TRUE, sort = FALSE)
	rownames(thisMat) <- thisMat$species
	thisMat
}

getMeasureMat <- function() {
	print("Building measurement matrix...")
	measureMat <- getSingleSpeciesMatrix()
	# print("getSingleSpMat completed")
	rownames(measureMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = measureMat$species)

	measureMat <- appendRegressionCategories(thisMat = measureMat, regMat = read.csv(file="~/Dropbox/code/R/dentalMeasurements/dat/regressionLabelsJDM.csv"))
	# print("regression append completed")
	#head(measureMat)
	measureMat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = measureMat$species)
	#Approxiate body masses
	measureMat <- approxBodyMass(measureMat = measureMat)
	# print("bodymass approximated")
	measureMat <- measureMat[is.finite(measureMat $bodyMass),]		#### clears taxa without body mass estimate

	return(measureMat)
}
