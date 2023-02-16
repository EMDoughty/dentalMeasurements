 #body mass estimation functions 

appendStDevToReg <- function(this.reg) {
	this.reg <- cbind(this.reg, stdev=log10(100 + this.reg["see"]) - 2)
	names(this.reg)[length(this.reg)] <- "stdev" 
	this.reg
}

getLikelihoodOneBodyMassFromMeasurementAndReg <- function(bm.estimate, this.species, this.reg) {
	sum(-dnorm(bm.estimate, mean=as.vector((log10(this.species) * this.reg[this.reg$m %in% colnames(this.species),]$slope) + this.reg[this.reg$m %in% colnames(this.species),]$intercept, mode="numeric"), sd=this.reg[this.reg$m %in% colnames(this.species),]$stdev, log=TRUE))
}

getMLbodyMassForOneSpecies <- function(this.species, this.reg, best.only=FALSE) {
	if (best.only) this.reg <- this.reg[this.reg$j %in% c("FLML","FLMA","SLML","SLMA","SUML","SUMA","TLMA","LMRL") ,]		#Janis's "best" variables Table 13.3, p.281
	this.species <- this.species[sapply(this.species, is.finite) & names(this.species) %in% this.reg$m]
	if (ncol(this.species)>0) {
		this.optim <- list(convergence=-1)
		while (this.optim$convergence !=0) {
			this.optim <- optim(par=runif(1, 0, 3), fn=getLikelihoodOneBodyMassFromMeasurementAndReg, this.species=data.matrix(this.species), this.reg=this.reg, method="L-BFGS-B", lower=-2, upper=5)
			# if (this.optim$convergence != 0) warning(paste("**** body mass optimization for", rownames(this.species), "did not converge", this.optim$convergence))
		}
		this.optim$par
	} else return (NA)
}

getMLBodyMasses <- function(this.m, reg.list, best.only=FALSE) {
	bm.vec <- array(NA, dim=c(nrow(this.m), 1), dimnames=list(rownames(this.m), "bodyMass"))
	for (i in seq_len(nrow(this.m))) {
		if (!is.na(this.m$reg[i]) & this.m$reg[i] %in% names(reg.list)) {
			# bm.vec[i] <- getMLbodyMassForOneSpecies_compiled(this.species=this.m[i,], this.reg=reg.list[[this.m$reg[i]]])
			bm.vec[i] <- getMLbodyMassForOneSpecies(this.species=this.m[i,], this.reg=reg.list[[this.m$reg[i]]])
		} else bm.vec[i] <- NA
	} 
	bm.vec
} 

getBodyMassVectorFromMeasureMatAllMeasures <- function(measure.mat, linked.files=FALSE) {
	#######################################################################################################################################
	##### read Janis/Damut 1990 regression parameters from file, and append standard deviations
	#######################################################################################################################################
	if (linked.files) {
		# reg.list <- list(ruminantia=read.csv("https://dl.dropbox.com/s/dcd0bs1x5v9e7lh/regRuminantia.csv"), perissodactyla=read.csv("https://dl.dropbox.com/s/04k387q7yh4wp9u/regPerissodactyla.csv"), ungulate=read.csv("https://dl.dropbox.com/s/310ayur1s1dc8sl/regAllUngulates.csv"))
	} else {
		reg.list <- list(ruminantia=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regRuminantia.csv"), 
						perissodactyla=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regPerissodactyla.csv"), 
						ungulate=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regAllUngulates.csv"),
						DamuthUngulate=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regAllArchaicUngulates.csv"),
						DamuthAllSelenodonts=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regArchaicAllSelenodonts.csv"),
						DamuthNonSelenodonts=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regArchaicNonselenodonts.csv"),
						DamuthSelenodontBrowsers=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regArchaicSelenodontBrowsers.csv"),
						DamuthSelenodontNonBrowsers=read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regArchaicSelenodontNonBrowsers.csv"))			
	}
	reg.list <- lapply(reg.list, appendStDevToReg)
	
	measure.mat$reg[measure.mat$reg==""] <- NA
	these.columns <- c(as.character(reg.list[["DamuthNonSelenodonts"]]$m)[as.character(reg.list[["DamuthNonSelenodonts"]]$m) %in% names(measure.mat)], "reg") # these.columns is the names of columns that are also in reg - published measuremnts
	these.taxa <- measure.mat$taxon[apply(X=data.matrix(measure.mat[,these.columns[seq_len(length(these.columns)-1)]]), MARGIN=1, FUN=function(x) any(sapply(x, is.finite)))]
	bm <- getMLBodyMasses(this.m=measure.mat[these.taxa, these.columns], reg.list=reg.list, best.only=FALSE)
	bm[match(measure.mat$taxon, rownames(bm))]	
}

fillMissingBodyMasses <- function(measure.mat) {
	# require(stringr)
	noMass <- measure.mat$taxon[!is.finite(measure.mat$bodyMass)]

		#### mean of congeners
		# noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) mean(measure.mat[grep(pattern=strsplit(x=x, split="_")[[1]][1], x=measure.mat$taxon),"bodyMass"], na.rm=TRUE)))
		# #### median of congeners
		noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) median(measure.mat[grep(pattern= strsplit(x=x, split="_")[[1]][1], x=measure.mat$taxon),"bodyMass"], na.rm=TRUE)))
		# #### randomly select a congener
		# noMass <- data.frame(species=noMass, bodyMass=sapply(X=noMass, FUN=function(x) sample(x=measure.mat[grep(pattern= strsplit(x=x, split="_")[[1]][1], x=measure.mat$taxon),"bodyMass"], size=1)))

	measure.mat$bodyMass[match(noMass$taxon, measure.mat$taxon)] <- noMass[,"bodyMass"]
	measure.mat
}

############################################################################################################################################

#Compile and format matrix of all measurements from multiple sources
appendRegressionCategories <- function(measure.mat, regMat, focal.tax) {
	
	regMat$family <- getCurrentTaxa(tax.vec = regMat$family)
	regMat$genus <- getCurrentTaxa(tax.vec = regMat$genus)
	regMat <- unique(regMat)
	
	uniqTax <- lapply(unlist(focal.tax), FUN=getTaxonomyForOneBaseTaxon)
	uniqTax <- makeMatrixFromList(uniqTax)
	uniqTax <- sapply(uniqTax, function(this.col) sapply(this.col, function(this.cell) if (grepl("_SPECIFIED", x=this.cell)) NA else this.cell))
	# uniqTax <- apply(as.matrix(uniqTax), 2, function(this.col) function(this.cell) ifgrepl("_SPECIFIED", x=this.col)] <- NA)
	uniqTax[uniqTax==""] <- NA
	uniqTax <- data.frame(uniqTax, stringsAsFactors=TRUE)
	uniqTax$taxon_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$taxon_name)
	uniqTax <- uniqTax[order(uniqTax$order, uniqTax$family, uniqTax$genus, uniqTax$taxon_name),]
	
	measure.mat$family <- uniqTax$family[match(x=measure.mat$taxon, table=uniqTax$taxon_name)]
	measure.mat$genus <- uniqTax$genus[match(x=measure.mat$taxon, table=uniqTax$taxon_name)]

	reg.vec <- regMat$DamReg[!is.na(regMat$family) & is.na(regMat$genus)][match(measure.mat$family, regMat$family[!is.na(regMat$family) & is.na(regMat$genus)])]  # this is the regression "labels" for the species from measure.mat in the correct order, based on family name
	
	special.genera <- sort(unique(measure.mat$genus[!is.na(measure.mat$genus) & measure.mat$genus %in% regMat$genus]))
	for (this.genus in special.genera) reg.vec[grepl(this.genus, measure.mat$genus)] <- regMat$DamReg[grepl(this.genus, regMat$genus)]
	
	measure.mat$reg.vec <- reg.vec

	#check for taxa that are not recieving a regression
	#missingReg <- measure.mat[is.na(measure.mat $reg.vec),]
	#missingReg <- missingReg[!grepl("sp.",missingReg$species),]
	#missingReg <- missingReg[!grepl("indet.",missingReg$species),]
	#missingReg <- missingReg[!grepl("cf.",missingReg$species),]
	#write.csv(missingReg, '/Users/evandoughty/Dropbox/ungulate_RA/RCode/JonCode/2017_2_27_missingReg.csv')
	
	# nrow(measure.mat) #803 species remain after final removal 
	
	# measure.mat$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = measure.mat$taxon)
	
	return(measure.mat)
}

############################################################################################################################################

approxBodyMass <- function(measure.mat, fill.missing=TRUE) {
	#Approximate body mass
	print("Building body mass estimates...")
	# measure.mat[!is.na(measure.mat$family) & measure.mat$family=="Entelodontidae", c("P2_L", "P3_L", "p2_w", "m2_w", "m3_w")] <- NA	### entelodont tooth widths were generating >5 ton body masses, so not used for body mass estimation, here.
	measure.mat[,"bodyMass"] <- getBodyMassVectorFromMeasureMatAllMeasures(measure.mat, linked.files=FALSE)
	if (fill.missing) measure.mat <- fillMissingBodyMasses(measure.mat)	# this fills taxa missing their body mass with the average body mass of its cogeners
	# sort(unique(measure.mat[!sapply(measure.mat, function(x) is.character(x) | is.finite(x) | is.na(x))])) <- NA
	return(measure.mat)
}

############################################################################################################################################

