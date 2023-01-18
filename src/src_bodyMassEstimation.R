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
	thisReg.StDev <- cbind(thisReg, (log10(100 + thisReg["see"]) - 2))
  colnames(thisReg.StDev) <- c(colnames(thisReg),"stdev")
  
 # rownames(thisReg) <- thisReg$m
	
#  cbind(thisReg, array((log10(100 + thisReg["see"]) - 2), dimnames=list("stdev")))   #can't get this to cbind correctly dwspite fixing the dim[2] issues of lkine below
#  cbind(thisReg, array((log10(100 + thisReg["see"]) - 2), dimnames=list(rownames(thisReg), "stdev"))) #11/30/2022 this is currently causing an error as dimanmes is anot proper length.  rownames are currently 1-n rather than any string.
  # cbind(thisReg, array((log10(100 + thisReg["see"]) - 2) * sqrt(thisReg["n"]), dimnames=list(rownames(thisReg), "stdev")))
	# thisReg["see"]=(2+see)
  
  return(thisReg.StDev)
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

getMLBodyMasses <- function(measure.mat, regList, best.only=FALSE) {
	bmVec <- array(NA, dim=c(nrow(measure.mat), 1), dimnames=list(rownames(measure.mat), "bodyMass"))
	for (i in seq_len(nrow(measure.mat))) {
		if (!is.na(measure.mat$reg[i])) {
			bmVec[i] <- getMLbodyMassForOneSpecies_compiled(this=measure.mat[i,], thisReg=regList[[which(names(regList) == measure.mat$reg[i])]])
		} else bmVec[i] <- NA
	} 
	bmVec
} 
getMLBodyMasses_compiled <- cmpfun(getMLBodyMasses)

appendRegTypeTomeasure.mat <- function(measure.mat) {
		famList <- unique(read.csv("~/Dropbox/code/common_dat/taxonomy.csv"), strip.white=TRUE)
		measure.mat[,"reg"] <- famList$reg[match(x=measure.mat$taxon, famList$taxon)]	# now assuming reg will already be a part of measure.mat
		measure.mat
}

getBodyMassVectorFromMeasureMatAllMeasures <- function(measure.mat, linked.files=FALSE) {
	#######################################################################################################################################
	##### read Janis 1990 regression parameters from file, and append standard deviations
	#######################################################################################################################################
	if (linked.files) {
		# regList <- list(ruminantia=read.csv("https://dl.dropbox.com/s/dcd0bs1x5v9e7lh/regRuminantia.csv"), perissodactyla=read.csv("https://dl.dropbox.com/s/04k387q7yh4wp9u/regPerissodactyla.csv"), ungulate=read.csv("https://dl.dropbox.com/s/310ayur1s1dc8sl/regAllUngulates.csv"))
	} else {
		regList <- list(ruminantia=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regRuminantia.csv"), 
		                perissodactyla=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regPerissodactyla.csv"), 
		                ungulate=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regAllUngulates.csv"), 
		                ArchaicAllSelenodonts=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicAllSelenodonts.csv"), 
		                ArchaicNonselenodonts=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicNonselenodonts.csv"), 
		                ArchaicSelenodontNonBrowsers=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicSelenodontNonBrowsers.csv"), 
		                ArchaicSelenodontBrowsers=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicSelenodontBrowsers.csv"), 
		                DasmuthAllUngulates=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regAllArchaicUngulates.csv"))
	}
	regList <- lapply(regList, appendStDevToReg)

	#######################################################################################################################################
	##### get body mass for only those taxa that have all (i.e., are not missing any) of the published measurements (about 325 species)
	#######################################################################################################################################
	measure.mat$reg[measure.mat$reg==""] <- NA
	theseColumns <- c(as.character(regList[[1]]$m)[-which(as.character(regList[[1]]$m) %in% c("loP", "loM"))], "reg") # theseColumns is the names of columns that are also in reg - published measuremnts
	bm <- getMLBodyMasses_compiled(measure.mat = measure.mat[complete.cases(measure.mat[,theseColumns]), theseColumns], regList = regList, best.only=FALSE)
	
	#######################################################################################################################################
	##### get regression parameters for measurements not in published regression
	#######################################################################################################################################
#	other_m <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M3_L", "M3_W")
# 	other_m <- colnames(measure.mat[,sapply(measure.mat, is.numeric)])[!colnames(measure.mat[,sapply(measure.mat, is.numeric)]) %in% theseColumns]

#	otherReg <- lapply(X=names(regList), FUN=function(this.group) {
#		shortMat <- measure.mat[rownames(measure.mat) %in% rownames(bm) & measure.mat$reg==this.group, other_m]
#		short.bm <- bm[rownames(bm) %in% rownames(shortMat),]
#		apply(log10(shortMat), MARGIN=2, FUN=function(x, bm) {
#			lm(bm ~ x) } , bm=short.bm ) 
#		} )

#	names(otherReg) <- names(regList)

#	otherList <- lapply(otherReg, function(thisReg) t(sapply(thisReg, function(x) c(coef(x), summary(x)$sigma))))
	
	#######################################################################################################################################
	##### merge regression parameters for unpublished measurements (otherList) with those from published (regList)
	#######################################################################################################################################
#	for (i in seq_along(regList)) {
#		colnames(otherList[[i]]) <- c("intercept", "slope", "stdev")
#		otherList[[i]] <- data.frame(m=rownames(otherList[[i]]), otherList[[i]], stringsAsFactors=FALSE)
#		regList[[i]] <- merge(regList[[i]], otherList[[i]], all=TRUE, sort=FALSE)
#	}

	#######################################################################################################################################
	##### recalculate body masses of all taxa with all (merged) measurements
	#######################################################################################################################################
	bm <-  getMLBodyMasses_compiled(measure.mat, regList, best.only=FALSE)
	bm[match(measure.mat$taxon, rownames(bm))]	
}

fillMissingBodyMasses <- function(measure.mat) {
	# require(stringr)
	noMass <- measure.mat$taxon[!is.finite(measure.mat$bodyMass)]

		#### mean of congeners
		noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) mean(measure.mat[grep(pattern=strsplit(x=x, split="_")[[1]][1], x=measure.mat$taxon),"bodyMass"], na.rm=TRUE)))
		# #### median of congeners
		# noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) median(measure.mat[grep(pattern= strsplit(x=x, split="_")[[1]][1], x=measure.mat$taxon),"bodyMass"], na.rm=TRUE)))
		# #### randomly select a congener
		# noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) sample(x=measure.mat[grep(pattern= strsplit(x=x, split="_")[[1]][1], x=measure.mat$taxon),"bodyMass"], size=1)))

	measure.mat$bodyMass[match(noMass$taxon, measure.mat$taxon)] <- noMass[,"bodyMass"]
	measure.mat
}

############################################################################################################################################

#Compile and format matrix of all measurements from multiple sources
appendRegressionCategories <- function(measure.mat, regMat, regAuthor = c("Janis")) {
	
	uniqTax <- lapply(c("Artiodactyla", "Perissodactyla", "Arctocyonidae", "Hyopsodontidae", "Periptychidae", "Phenacodontidae"), FUN=getTaxonomyForOneBaseTaxon)
	uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]], uniqTax[[3]], uniqTax[[4]], uniqTax[[5]], uniqTax[[6]])
	uniqTax$taxon_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$taxon_name)
	

	#delete unique taxa that lack family and genus
	measure.mat$family <- uniqTax$family[match(x=rownames(measure.mat), table=uniqTax$taxon_name)]
	# measure.mat$family <- sapply(measure.mat$family, as.character)
	measure.mat$family[measure.mat$family == ""] <- NA
	
	measure.mat$genus <- uniqTax$genus[match(x=rownames(measure.mat), table=uniqTax$taxon_name)]
	# measure.mat$genus <- sapply(measure.mat$genus, as.character)
	measure.mat$genus[measure.mat$genus == ""] <- NA
	
	#### 
	#append regression catagories to each species
	####
	
	# reg.vec <- array(dim=nrow(measure.mat))
	# apply(measure.mat, 1, function(x) {
		# if (!is.na(x["family"]) & x["family"] %in% regMat$family) 
	
	#set whichever regression system to be the column reg
	if(regAuthor == "Janis") regMat$reg <- regMat$JanReg
	if(regAuthor == "Damuth") regMat$reg <- regMat$DamReg
	
	family.names <- uniqTax$family[match(x=rownames(measure.mat), table=uniqTax$taxon_name)]
	reg.vec <- regMat$reg[!is.na(regMat$family)][match(family.names, regMat$family[!is.na(regMat$family)])]  # this is the regression "labels" for the species from measure.mat in the correct order, based on family name
	
	genus.names <- uniqTax$genus[match(x=rownames(measure.mat), table=uniqTax$taxon_name)]
	reg.vec[is.na(reg.vec)] <- regMat$reg[!is.na(regMat$genus)][match(genus.names, regMat$genus[!is.na(regMat$genus)])][is.na(reg.vec)]   #this is the regression "labels" for the species from measure.mat in the correct order, based on genus name; it appears that having regMat$genus[!is.na(regMat$genus)] will cause the index to improperly assign regressions
	
	measure.mat$reg.vec <- reg.vec

	#check for taxa that are not recieving a regression
	#missingReg <- measure.mat[is.na(measure.mat $reg.vec),]
	#missingReg <- missingReg[!grepl("sp.",missingReg$species),]
	#missingReg <- missingReg[!grepl("indet.",missingReg$species),]
	#missingReg <- missingReg[!grepl("cf.",missingReg$species),]
	#write.csv(missingReg, '/Users/evandoughty/Dropbox/ungulate_RA/RCode/JonCode/2017_2_27_missingReg.csv')
	
	# nrow(measure.mat) #803 species remain after final removal 
	
	# rownames(measure.mat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(measure.mat))
	
	return(measure.mat)
}

############################################################################################################################################

approxBodyMass <- function(measure.mat, fill.missing=TRUE) {
	#Approximate body mass
	print("Building body mass estimates...")
	measure.mat[!is.na(measure.mat$family) & measure.mat$family=="Entelodontidae", c("P2_L", "P3_L", "p2_w", "m2_w", "m3_w")] <- NA	### entelodont tooth widths were generating >5 ton body masses, so not used for body mass estimation, here.
	measure.mat[,"bodyMass"] <- getBodyMassVectorFromMeasureMatAllMeasures(measure.mat = measure.mat, linked.files=FALSE)
	if (fill.missing) measure.mat <- fillMissingBodyMasses(measure.mat)	# this fills taxa missing their body mass with the average body mass of its cogeners
	# sort(unique(measure.mat[!sapply(measure.mat, function(x) is.character(x) | is.finite(x) | is.na(x))])) <- NA
	return(measure.mat)
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
