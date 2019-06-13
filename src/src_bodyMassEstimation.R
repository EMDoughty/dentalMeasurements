#body mass estimation functions 

getBMEstimateFromMeasurementAndRegression <- function(m, regRow) {
	10^((log10(m) * regRow$slope) + regRow$intercept)
}
getBMEstimateFromMeasurementAndRegression_compiled <- cmpfun(getBMEstimateFromMeasurementAndRegression)

getBodyMassFromMeasurement <- function(thisCol, regRow) {
	sapply(thisCol, getBMEstimateFromMeasurementAndRegression_compiled, regRow)
}

getMendozaBodyMasses <- function (this) {
	sevenOne <- (2.543 * this$m1_w) + (1.827 * this$loM) + (1.402 * this$m2_w) + 1.280 #7.1
	sevenTwo <- (1.850 * this$m1_w) + (1.883 * this$loM) + (1.073 * this$m2_w) + (0.308 * this$p4_w) + 1.216 #7.2
	sevenThree <- (1.460 * this$m1_w) + (1.363 * this$loM) + (1.182 * this$m2_w) + (0.387 * this$p4_w) + (0.955 * this$m2_l) + 1.578 #7.3
	sevenFour <- (1.226 * this$m1_w) + (1.313 * this$loM) + (1.090 * this$m2_w) + (0.320 * this$p4_w) + (1.095 * this$m2_l) + (0.150 * this$loP) + 1.361 #7.4
	sevenFive <- (1.355 * this$m1_w) + (1.427 * this$loM) + (1.322 * this$m2_w) + (0.457 * this$p4_w) + (1.177 * this$m2_l) + (0.234 * this$loP) + (0.340 * this$p4_l) + 1.168 #7.5
	sevenSix <- (2.045 * this$m1_l) + (1.073 * this$loM) + 1.507	#7.6*
	sevenSeven <- (1.213 * this$m1_l) + (1.421 * this$loM) + (0.422 * this$p4_w) + 1.380	#7.7*
	c(sevenOne, sevenTwo, sevenThree, sevenFour, sevenFive, sevenSix, sevenSeven)
}

appendStDevToReg <- function(thisReg) {
	cbind(thisReg, array((log10(100 + thisReg["see"]) - 2), dimnames=list(rownames(thisReg), "stdev")))
	# cbind(thisReg, array((log10(100 + thisReg["see"]) - 2) * sqrt(thisReg["n"]), dimnames=list(rownames(thisReg), "stdev")))
	# thisReg["see"]=(2+see)	
}

getLikelihoodOneBodyMassFromMeasurementAndReg <- function(bm, this, thisReg) {
	lnl <- vector()
	for (i in seq_along(this)) lnl[i] <- dnorm(bm, mean = (log10(this[i]) * thisReg$slope[thisReg$m == colnames(this)[i]]) + thisReg$intercept[thisReg$m == colnames(this)[i]], sd=thisReg$stdev[thisReg$m == colnames(this)[i]], log=TRUE)
	sum(-lnl)
}

getMLbodyMassForOneSpecies <- function(this, thisReg, best.only=FALSE) {
	if (best.only) thisReg <- thisReg[thisReg$j %in% c("FLML","FLMA","SLML","SLMA","SUML","SUMA","TLMA","LMRL") ,]		#Janis's "best" variables Table 13.3, p.281
	this <- this[sapply(this, is.finite) & names(this) %in% thisReg$m]
	if (ncol(this)==0) return (NA)
	# bmVec <- vector()
	# for (i in seq_along(this)) bmVec <- c(bmVec, getBMEstimateFromMeasurementAndRegression_compiled(this[i], thisReg[match(colnames(this)[i], thisReg$m),]))
	# # sapply(this, getBMEstimateFromMeasurementAndRegression, thisReg[match(names(this[i]), as.character(thisReg$m)),])
	# mean(bmVec, na.rm=TRUE)
	optim(2, fn=getLikelihoodOneBodyMassFromMeasurementAndReg, this=data.matrix(this), thisReg=thisReg, method="L-BFGS-B", lower=-Inf, upper=Inf)$par
}
getMLbodyMassForOneSpecies_compiled <- cmpfun(getMLbodyMassForOneSpecies)

getMLBodyMasses <- function(thisMat, regList, best.only=FALSE) {
	bmVec <- array(NA, dim=c(nrow(thisMat), 1), dimnames=list(rownames(thisMat), "bodyMass"))
	for (i in seq_len(nrow(thisMat))) {
		if (!is.na(thisMat$reg[i])) {
			bmVec[i] <- getMLbodyMassForOneSpecies_compiled(this=thisMat[i,], thisReg=regList[[which(names(regList) == thisMat$reg[i])]])
		} else bmVec[i] <- NA
	} 
	bmVec
} 
getMLBodyMasses_compiled <- cmpfun(getMLBodyMasses)

appendRegTypeToThisMat <- function(thisMat) {
		famList <- unique(read.csv("~/Dropbox/code/common_dat/taxonomy.csv"), strip.white=TRUE)
		thisMat[,"reg"] <- famList$reg[match(x=thisMat$species, famList$taxon)]	# now assuming reg will already be a part of thisMat
		thisMat
}

getBodyMassVectorFromMeasureMatAllMeasures <- function(thisMat, linked.files=FALSE) {
	#######################################################################################################################################
	##### read Janis 1990 regression parameters from file, and append standard deviations
	#######################################################################################################################################
	if (linked.files) {
		# regList <- list(ruminantia=read.csv("https://dl.dropbox.com/s/dcd0bs1x5v9e7lh/regRuminantia.csv"), perissodactyla=read.csv("https://dl.dropbox.com/s/04k387q7yh4wp9u/regPerissodactyla.csv"), ungulate=read.csv("https://dl.dropbox.com/s/310ayur1s1dc8sl/regAllUngulates.csv"))
	} else {
		regList <- list(ruminantia=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regRuminantia.csv"), perissodactyla=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regPerissodactyla.csv"), ungulate=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regAllUngulates.csv"))
	}
	regList <- lapply(regList, appendStDevToReg)

	#######################################################################################################################################
	##### get body mass for only those taxa that have all (i.e., are not missing any) of the published measurements (about 325 species)
	#######################################################################################################################################
	thisMat$reg[thisMat$reg==""] <- NA
	theseColumns <- c(as.character(regList[[1]]$m)[-which(as.character(regList[[1]]$m) %in% c("loP", "loM"))], "reg") # theseColumns is the names of columns that are also in reg - published measuremnts
	bm <- getMLBodyMasses_compiled(thisMat[complete.cases(thisMat[,theseColumns]), theseColumns], regList, best.only=FALSE)
	
	#######################################################################################################################################
	##### get regression parameters for measurements not in published regression
	#######################################################################################################################################
	other_m <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M3_L", "M3_W")
# 	other_m <- colnames(thisMat[,sapply(thisMat, is.numeric)])[!colnames(thisMat[,sapply(thisMat, is.numeric)]) %in% theseColumns]

	otherReg <- lapply(X=names(regList), FUN=function(this.group) {
		shortMat <- thisMat[rownames(thisMat) %in% rownames(bm) & thisMat$reg==this.group, other_m]
		short.bm <- bm[rownames(bm) %in% rownames(shortMat),]
		apply(log10(shortMat), MARGIN=2, FUN=function(x, bm) {
			lm(bm ~ x) } , bm=short.bm ) 
		} )

	names(otherReg) <- names(regList)

	otherList <- lapply(otherReg, function(thisReg) t(sapply(thisReg, function(x) c(coef(x), summary(x)$sigma))))
	
	#######################################################################################################################################
	##### merge regression parameters for unpublished measurements (otherList) with those from published (regList)
	#######################################################################################################################################
	for (i in seq_along(regList)) {
		colnames(otherList[[i]]) <- c("intercept", "slope", "stdev")
		otherList[[i]] <- data.frame(m=rownames(otherList[[i]]), otherList[[i]], stringsAsFactors=FALSE)
		regList[[i]] <- merge(regList[[i]], otherList[[i]], all=TRUE, sort=FALSE)
	}

	#######################################################################################################################################
	##### recalculate body masses of all taxa with all (merged) measurements
	#######################################################################################################################################
	bm <-  getMLBodyMasses_compiled(thisMat, regList, best.only=FALSE)
	bm[match(thisMat$species, rownames(bm))]	
}

fillMissingBodyMasses <- function(thisMat) {
	require(stringr)
	noMass <- rownames(thisMat)[!is.finite(thisMat$bodyMass)]
	noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) mean(thisMat[grep(str_split(x, pattern=" ")[[1]][1], rownames(thisMat)),"bodyMass"], na.rm=TRUE)))
	thisMat$bodyMass[match(rownames(noMass), rownames(thisMat))] <- noMass[,"bodyMass"]
	array(thisMat$bodyMass, dimnames=list(rownames(thisMat)))
}

############################################################################################################################################

#Compile and format matrix of all measurements from multiple sources
appendRegressionCategories <- function(thisMat, regMat) {
	
	uniqTax <- lapply(c("Artiodactyla", "Perissodactyla"), FUN=getTaxonomyForOneBaseTaxon)
	uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]])
	uniqTax$taxon_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$taxon_name)
	
	#delete unique taxa that lack family and genus
	thisMat$family <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	# thisMat$family <- sapply(thisMat$family, as.character)
	thisMat$family[thisMat$family == ""] <- NA
	
	thisMat$genus <- uniqTax$genus[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	# thisMat$genus <- sapply(thisMat$genus, as.character)
	thisMat$genus[thisMat$genus == ""] <- NA
	
	#### 
	#append regression catagories to each species
	####
	
	# reg.vec <- array(dim=nrow(thisMat))
	# apply(thisMat, 1, function(x) {
		# if (!is.na(x["family"]) & x["family"] %in% regMat$family) 
	
	family.names <- uniqTax$family[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	reg.vec <- regMat$reg[!is.na(regMat$family)][match(family.names, regMat$family[!is.na(regMat$family)])]  # this is the regression "labels" for the species from measure.mat in the correct order, based on family name
	
	genus.names <- uniqTax$genus[match(x=rownames(thisMat), table=uniqTax$taxon_name)]
	reg.vec[is.na(reg.vec)] <- regMat$reg[!is.na(regMat$genus)][match(genus.names, regMat$genus[!is.na(regMat$genus)])][is.na(reg.vec)]   #this is the regression "labels" for the species from measure.mat in the correct order, based on genus name; it appears that having regMat$genus[!is.na(regMat$genus)] will cause the index to improperly assign regressions
	
	thisMat$reg.vec <- reg.vec

	#check for taxa that are not recieving a regression
	#missingReg <- measure.mat[is.na(measure.mat $reg.vec),]
	#missingReg <- missingReg[!grepl("sp.",missingReg$species),]
	#missingReg <- missingReg[!grepl("indet.",missingReg$species),]
	#missingReg <- missingReg[!grepl("cf.",missingReg$species),]
	#write.csv(missingReg, '/Users/evandoughty/Dropbox/ungulate_RA/RCode/JonCode/2017_2_27_missingReg.csv')
	
	# nrow(thisMat) #803 species remain after final removal 
	
	# rownames(thisMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisMat))
	
	return(thisMat)
}

############################################################################################################################################

approxBodyMass <- function(measureMat) {
	#Approximate body mass
	print("Building body mass estimates...")
	measureMat[!is.na(measureMat $family) & measureMat $family=="Entelodontidae", c("P2_L", "P3_L", "p2_w", "m2_w", "m3_w")] <- NA	### entelodont tooth widths were generating >5 ton body masses, so dropped here.
	measureMat[,"bodyMass"] <- getBodyMassVectorFromMeasureMatAllMeasures(measureMat, linked.files=FALSE)
	measureMat $bodyMass <- fillMissingBodyMasses(measureMat)	# this fills taxa missing their body mass with the average body mass of its cogeners
	measureMat[!sapply(measureMat, is.finite)] <- NA
	#rownames(measureMat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(measureMat))
	measureMat $species <- rownames(measureMat)
	return(measureMat)
}

############################################################################################################################################

#### Function for checking for, isolating, and appending species entries that were not present in the initial regression catagorizing matrix.
#### Entries are merely located and appended. The reg.vec column in the regMat file will retain all NA values for the taxa missing those entries
#### unless the user changes or removes them via manual or automatic means.

checkMissingReg<- function(measReg, uniqTax) {
	missingReg <- measure.matReg[is.na(measure.matReg$reg.vec),]
	#merge/match with occs fiel to get order, genus, and family columns
	occurrences.order_name <- uniqTax$order[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.order_name)
	
	occurrences.family_name <- uniqTax$family[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.family_name)
	
	occurrences.genus_name <- uniqTax$genus[match(x=rownames(missingReg), table=uniqTax$accepted_name)]
	missingReg <- cbind(missingReg,occurrences.genus_name)
	
	missingReg <- missingReg[,c("occurrences.order_name","occurrences.family_name","occurrences.genus_name","species","reg.vec")]
	head(missingReg)
	colnames(missingReg)[colnames(missingReg) == "species"] <- "taxon"
	colnames(missingReg)[colnames(missingReg) == "reg.vec"] <- "reg"
	
	regMat <- rbind(regMat, missingReg)
	
	return(regMat) 
}
