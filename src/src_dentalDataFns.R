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
	replicatemeasure.mat <- apply(cube, c(1, 2), mean, na.rm = TRUE)
	rownames(replicatemeasure.mat) <- datMeasures$species[seq(from = 1, to = nrow(datMeasures), by = 3)]
	replicatemeasure.mat[!is.finite(replicatemeasure.mat)] <- NA
	# replicatemeasure.mat<-data.frame(species=datMeasures$species[seq(from=1, to=nrow(datMeasures), by=3)], specimen=datMeasures$Specimen.no.[seq(from=1, to=nrow(datMeasures), by=3)], replicatemeasure.mat, stringsAsFactors=FALSE, row.names=NULL)
	
	theseSpecimenNos <- datMeasures$Specimen.no.[seq(from = 1, to = nrow(datMeasures), by = 3)]
	specimenNos <- unique(theseSpecimenNos)
	specimen.mat <- matrix(nrow = 0, ncol = ncol(replicatemeasure.mat))
	for (i in seq_len(length(specimenNos))) {
		index <- which(theseSpecimenNos == specimenNos[i])
		thisSpecimen <- replicatemeasure.mat[index, ]
		if (!is.null(dim(thisSpecimen))) 
			thisSpecimen <- apply(thisSpecimen, 2, mean, na.rm = TRUE)
		specimen.mat <- rbind(specimen.mat, thisSpecimen)
		rownames(specimen.mat)[nrow(specimen.mat)] <- rownames(replicatemeasure.mat)[index[1]]
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
	# species <-unique(specimen.mat$species)
	oneSpeciesMat <- aggregate(specimen.mat, by = list(specimen.mat$species), mean, na.rm = TRUE)
	# oneSpeciesMat <- aggregate(specimen.mat, by=list(specimen.mat$species), median, na.rm=TRUE)
	oneSpeciesMat <- data.frame(oneSpeciesMat, row.names = oneSpeciesMat[, 1])
	oneSpeciesMat[sapply(oneSpeciesMat, is.nan)] <- NA
	colnames(oneSpeciesMat)[1] <- "taxon"
	oneSpeciesMat[, -(c(2:3, which(colnames(oneSpeciesMat) %in% c("Published.name", "n", "Reference", "PaleoDB.ref", "Notes", "X", "X.1", "X.2", "X.3", 
		"X.4"))))]
}

makeOneGenusMatFromSpecimenMat <- function(measure.mat) {
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

############################################################################################################################################

#function is meant to be run prior to analysis to bring all measurement datasets together into a single entity
getSingleSpeciesMatrix <- function() {
	#compile and lable dental measurments for specimens
	specimen.mat <- getSpecimenMatFromAmandaMeasurements()
	specimen.mat <- merge(x = specimen.mat, y = getSpecimenMatFromBlastoMeasurements(), all = TRUE)
	specimen.mat <- merge(x = specimen.mat, y = getSpecimenMatFromLiteratureMeasurements(), all = TRUE)

	specimen.mat[sapply(specimen.mat, is.nan)] <- NA
	specimen.mat$species <- getCurrentTaxa(tax.vec = specimen.mat$species)
	specimen.mat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = specimen.mat$species)

	upLabels <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W") #\P2_L\",\"P2_W\","
	loLabels <- casefold(upLabels)

	############################
	#### note that specimens are aggregated by their medians, so as to minimize the effect of outlier measurements
	############################
	
	# measure.mat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(taxon = specimen.mat$species), mean, na.rm = TRUE)
	measure.mat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(taxon = specimen.mat$species), median, na.rm = TRUE)

	measure.mat[, sapply(measure.mat, is.numeric)] <- measure.mat[, sapply(measure.mat, is.numeric)]/10 #converts mm measurements to cm for compatibility with Janis regressions
	measure.mat <- transform(measure.mat, p4_a = p4_l * p4_w, m1_a = m1_l * m1_w, m2_a = m2_l * m2_w, m3_a = m3_l * m3_w, M2_A = M2_L * M2_W)
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

getMeasureMatWithBodyMasses <- function() {
	print("Building measurement matrix...")
	measure.mat <- getSingleSpeciesMatrix()

	tax.vec <- sort(unique(occs$accepted_name[occs$accepted_rank %in% c("genus", "species")]))
	tax.vec <- tax.vec[!tax.vec %in% rownames(measure.mat)]
	measure.mat <- appendMissingPaleoDBSpecies(measure.mat, tax.vec)

	measure.mat <- appendRegressionCategories(measure.mat = measure.mat, regMat = read.csv(file="~/Dropbox/code/R/dentalMeasurements/dat/regressionLabelsJDM.csv"))
	measure.mat <- approxBodyMass(measure.mat = measure.mat)
	measure.mat <- measure.mat[is.finite(measure.mat$bodyMass),]		#### clears taxa without body mass estimate

	return(measure.mat)
}
