# # setwd("C:/Users/Evan/Dropbox/ungulate_RA/RCode")

# #all ungulates
#load("~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata")

# #artio only
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/artio/handleyResult##------ Wed Nov 15 16:57:31 2017 ------##_artiodactyla.Rdata")

# #perisso only
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/perisso/handleyResult##------ Thu Nov 16 10:27:58 2017 ------##.Rdata")
# startTime <- Sys.time()
# load('~/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/repIntTaxa_SampleStandardized=TRUE##------ Fri Mar 29 21:20:21 2019 ------##.Rdata')

#sources for Jon Marcot's code and specimen measurements 
source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R") 
source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_makeRepIntOccs.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)

source("/Users/emdoughty/Dropbox/Code/R/dentalMeasurements/src/src_evanproposal.R")

require(stringr)
require(dplyr)
require(abind)

####################################################################################################################################
settings <- list()
settings$this.rank <- "species"
primary.workspace <- "~/Dropbox/Code/R/Files to Save outside git/Ch1_Analysis/Inputs/"

	options(timeout=300)
  # occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
	 occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
	  occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
	  occs <- occs[!occs$family %in% c("Allodelphinidae", "Allodesminae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
	  occs <- occs[!occs$genus %in% c("Enaliarctos", "Pteronarctos", "Kolponomos", "Pacificotaria", "Pinnarctidion", "Pteronarctos"), ]
	  occs <- occs[!occs$accepted_name %in% c("Archaeoceti", "Pinnipedia", "Imagotariinae"), ]
	  occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

	if (settings$this.rank %in% names(occs)) { occs <- occs[!(occs[, settings$this.rank] != "" & duplicated(occs[,c("collection_no", settings$this.rank)])),]
	} else if (settings$this.rank == "species") {
		occs <- occs[!(occs$accepted_rank=="species" & duplicated(occs[,c("collection_no", "accepted_name")])),]
	}

	 occs.filename <- paste0(primary.workspace,"occs_",timestamp(),".csv")
	 settings$occs.filename <- occs.filename
	 write.csv(occs, file = occs.filename)
	  
####################################################################################################################################

settings$focal.tax$clade <- c("Artiodactyla", "Perissodactyla", "Condylarthra", "Dinocerata", "Taeniodonta", "Pantodonta", "Tillodontia", "Arctocyonidae","Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae") #Including Proboscidea here causes functions that use focal.clade to query PBDB to break.  Doesn't return elephant-morphs just Hemiptera.
#settings$focal.tax$order <- c("Artiodactyla", "Perissodactyla", "Dinocerata", "Cimolesta") #Cimolesta is used as there are a few Taeniodont and tillodonts that lack a family designation
#settings$focal.tax$family <- c("Arctocyonidae", "Hyopsodontidae", "Phenacodontidae", "Periptychidae")
settings$bmBreaks_herb <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
	 
measure.mat <- getMeasureMatWithBodyMasses(settings)
measure.mat$SizeCat <- measure.mat$bodyMass
for(xx in seq(1, length(settings$bmBreaks_herb)-1, 1)){
  measure.mat$SizeCat[measure.mat$bodyMass > settings$bmBreaks_herb[xx] & measure.mat$bodyMass < settings$bmBreaks_herb[xx+1]] <- xx
} 

#append Proboscideans onto measure mat
probo.mat <- as.data.frame(matrix(nrow=length(unique(occs$accepted_name[occs$order %in% "Proboscidea" & occs$accepted_rank %in% settings$this.rank])), ncol=ncol(measure.mat))); colnames(probo.mat) <- colnames(measure.mat)
probo.mat$taxon <- unique(occs$accepted_name[occs$order %in% "Proboscidea" & occs$accepted_rank %in% settings$this.rank]); probo.mat <- probo.mat[!probo.mat$taxon %in% "",]
probo.mat$SizeCat <- 5
measure.mat <- rbind(measure.mat, probo.mat)
measure.mat <- measure.mat[order(measure.mat$taxon),]

if(settings$this.rank =="genus") measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat)

####################################################################################################################################
#### reduces matrix to just the focal order(s)
####################################################################################################################################
#need way to get around issue of species without higher taxonimic designations in occs despite being in a clade.  Could do work around where bigList is derived from getHigher taxon and then compared with occs
bigList <- unique(getTaxaInClade(clades = settings$focal.tax$clade, occs = occs, save.file = NULL))
bigList <- unique(bigList[,c("order","family", "genus", "accepted_name")])
proboList <- occs[occs$order %in% "Proboscidea",]; proboList <- unique(proboList[,c("order", "family", "genus","accepted_name")])
bigList <- rbind(bigList, proboList)

bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]

# bigList <- unheadique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% settings$focal.tax$order), c("order","family", "genus", "accepted_name")])
#bigList <- unique(occs[(occs$accepted_rank =="species" & (occs$order %in% settings$focal.tax$order | occs$family %in% settings$focal.tax$family)), c("order","family", "genus", "accepted_name")])
#bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]

shortFam <- sort(unique(bigList$family)) #[bigList$order %in% settings$focal.tax$clade]))
if (any(shortFam == "NO_FAMILY_SPECIFIED")) shortFam <- shortFam[-which(shortFam == "NO_FAMILY_SPECIFIED")]
if (any(shortFam == " ")) shortFam <- shortFam[-which(shortFam == " ")]

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
 all.taxa <- sort(unique(bigList$accepted_name))

measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name, ]

measure.mat.filename <- paste0(primary.workspace,"measure_mat_", settings$this.rank, "_", timestamp())
write.csv(measure.mat, file = paste0(measure.mat.filename,".csv"))
save(settings, bigList, shortFam, all.taxa, measure.mat, file = paste0(measure.mat.filename,".Rdata"))

####################################################################################################################################
settings$focal.tax$clade <- c("Carnivoramorpha", "Creodonta","Hyaenodonta", "Mesonychidae")
settings$pred.data.filename <- "~/Dropbox/Proposal/Proposal/predator_data_final.csv"
pred.data <- read.csv(settings$pred.data.filename)
colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl",	"bodyMass",	"Citation")

#append Mesonychidae from  Zhao?
meson.mat <- read.csv("/Users/emdoughty/Dropbox/Code/R/Files to Save outside git/Mesonychidae_BodySize.csv")
meson.mat <- meson.mat[meson.mat$family %in% "Mesonychidae", c("family", "accepted_name", "max_ma", "min_ma", "kg")]
meson.mat <- meson.mat[!meson.mat$accepted_name %in% "",]
meson.mat$m1L <-  meson.mat$rbl <- meson.mat$Citation <- NA
colnames(meson.mat) <- c("family", "taxon",	"max_ma","min_ma",	"bodyMass","m1L",	"rbl",	"Citation")
meson.mat <- meson.mat[, c("family", "taxon", "max_ma", "min_ma", "m1L", "rbl", "bodyMass", "Citation")]
pred.data <- rbind(pred.data, meson.mat)

pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)
rownames(pred.data) <- pred.data$taxon

pred.data[,c("bodyMass")] <- log10(pred.data[,c("bodyMass")])
pred.data <- pred.data[is.finite(pred.data$bodyMass),]

pred.data$genus <- unlist(lapply(strsplit(pred.data$taxon,"_"), function(x) x[1]))

if(settings$this.rank =="genus") pred.data <- makeOneGenusMatFromSpecimenMat(pred.data)

bigList <- unique(getTaxaInClade(clades = settings$focal.tax$clade, occs = occs, save.file = NULL))
bigList <- unique(bigList[,c("order","family", "genus", "accepted_name")])

bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]

shortFam <- sort(unique(bigList$family)) #[bigList$order %in% settings$focal.tax$clade]))
if (any(shortFam == "NO_FAMILY_SPECIFIED")) shortFam <- shortFam[-which(shortFam == "NO_FAMILY_SPECIFIED")]
if (any(shortFam == " ")) shortFam <- shortFam[-which(shortFam == " ")]

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
all.taxa <- sort(unique(bigList$accepted_name))

pred.data <- pred.data[pred.data$taxon %in% bigList$accepted_name, ]
###############

pred.diet <- read.csv("/Users/emdoughty/Dropbox/Proposal/Proposal/Diet Data PPP CJ.csv")
pred.diet$MASTER_LIST <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.diet$MASTER_LIST)

#pred.data.master  <- pred.data

pred.data <- pred.data[pred.data$taxon %in% pred.diet$MASTER_LIST,]
pred.data$Diet <- pred.diet[pred.diet$MASTER_LIST %in% pred.data$taxon, "PPP_diet_2"]

pred.data.hyper <- pred.data[pred.data$Diet %in% "hypercarnivore",]
pred.data.meso <- pred.data[pred.data$Diet %in% "mesocarnivore",]
pred.data.hypo <- pred.data[pred.data$Diet %in% "hypocarnivore",]
pred.data.hyper.meso <- pred.data[pred.data$Diet %in% c("hypercarnivore","mesocarnivore"),]
pred.data.hyper.hypo <- pred.data[pred.data$Diet %in% c("hypercarnivore","hypocarnivore"),]
pred.data.meso.hypo <- pred.data[pred.data$Diet %in% c("mesocarnivore","hypocarnivore"),]

pred.data.filename <- paste0(primary.workspace,"pred_data_",settings$this.rank,"_", timestamp())
write.csv(pred.data, file = paste0(pred.data.filename,".csv"))
save(settings, bigList, shortFam, all.taxa, 
     pred.data, pred.data.hyper, pred.data.meso, pred.data.hypo, pred.data.hyper.meso, pred.data.hyper.hypo, pred.data.meso.hypo,
     file = paste0(pred.data.filename,".Rdata"))

####################################################################################################################################

for(xx in c(3.5,4))
{
settings <- list()
settings$occs.filename <- occs.filename

settings$int_length <- 2
settings$min_age <- 1
settings$max_age <- 65

if(is.numeric(settings$int_length)) intervals <- makeIntervals(settings$min_age, settings$max_age, settings$int_length) else (intervals <- getBinsNALMA(settings))
intList <- listifyMatrixByRow(intervals)

settings$n.reps <- 10000
settings$do.subsample <- TRUE
settings$quota <- 0.4

settings$bootstrapSpecimens <- FALSE
settings$bootstrapSpecies <- FALSE
settings$bootstrapSpeciesWithinIntervals <- FALSE
settings$plotHist <- FALSE

settings$do.rangethrough <- FALSE

	if (settings$plotHist) {
		quartz("Guild Histograms")
		par(mfrow=c((nrow(intervals)), 3), mar=c(0,0,0.75,0), cex.axis=0.5, cex.main=0.75)
	}
  
	repIntOccs <- getRepIntOccs(settings = settings, intervals=intervals)

	repInt.folder <- paste0(primary.workspace,"repInt_numReps=", settings$n.reps,"_int_length=", settings$int_length,"_age=",settings$max_age,"to", settings$min_age,"_subsample=", settings$do.subsample,"/")
	if(!dir.exists(repInt.folder)) dir.create(repInt.folder)
	repIntOccs.filename <- paste0(repInt.folder, "repIntOccs", timestamp(), ".Rdata")
	save(settings, repIntOccs, file = repIntOccs.filename) 
	print("**** repIntOccs saved to file...")

	settings$repIntOccs.filename <- repIntOccs.filename
	
	repIntSp.start <- Sys.time()
	
	#####
	settings$this.rank <- "species"
	repIntTaxa <- getRepIntTaxaFromRepIntOccs(settings, repIntOccs, do.rangethrough=settings$do.rangethrough)
	save(settings, repIntTaxa, file = paste0(repInt.folder,"repIntTax_", settings$this.rank, "_do.rangethrough=",settings$do.rangethrough,"_", timestamp(),".Rdata"))

	repIntSp.end <- Sys.time()
	
	repIntGen.start <- Sys.time()
	
	settings$this.rank <- "genus"
	repIntTaxa <- getRepIntTaxaFromRepIntOccs(settings, repIntOccs, do.rangethrough=settings$do.rangethrough)
	save(settings, repIntTaxa, file = paste0(repInt.folder,"repIntTax_", settings$this.rank, "_do.rangethrough=",settings$do.rangethrough,"_", timestamp(),".Rdata"))
	
	repIntGen.end <- Sys.time()
	
	print("Completed getting taxa with intervals")
	

	repIntSp.end-repIntSp.start
	repIntGen.end-repIntGen.start
}


####################################################################################################################################
### Handley analysis of taxonomic distributions
print("Beginning median taxonomic Handley analysis...")

do.parallel <- TRUE
	if (do.parallel) require(parallel)
extra.intvs <- 0
do.disparity <- FALSE
do.heuristic <- FALSE

class.vec.tax <- sort(unique(unlist(repIntTaxa)))
class.vec.tax <- class.vec.tax[class.vec.tax %in% as.character(bigList$accepted_name)]
# this.class[this.class=="NO_FAMILY_SPECIFIED"] <- NA
class.vec.tax <- factor(sapply(class.vec.tax, function(this.taxon) shortFam[match(bigList$family[as.character(bigList$accepted_name) %in% this.taxon], shortFam)], simplify="array"), exclude="NO_FAMILY_SPECIFIED")

range.mat.median <- apply(taxRangeBox, c(1,2), function(x) as.logical(round(median(x, na.rm=TRUE))))

optList_tax_median <- doHandleyTest(master.class.list=class.vec.tax, master.range.mat=range.mat.median, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)

print("Beginning taxonomic Handley analysis for all n.reps...")
optList_tax_allReps <- list()
for (this.rep in seq_len(settings$n.reps)) {
	optList_tax_allReps[[this.rep]] <- doHandleyTest(master.class.list=class.vec.tax, master.range.mat=taxRangeBox[,,this.rep], do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)
	if (this.rep %% 100 == 0) cat("Taxonomic Handley Rep:", this.rep, "\n")
}

####################################################################################################################################
# if(Sys.info()["sysname"] == "Darwin"){
	save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, file=paste0("~/Desktop/EcologyResults/Taxon_handleyResult_SampleStandardized=", settings$do.subsample, timestamp(),".Rdata"))
	# save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=", settings$do.subsample, timestamp(),".Rdata"))
	# #load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
# } else if(Sys.info()["sysname"] == "Windows"){
	# save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=", settings$do.subsample, timestamp(),".Rdata"))
	# # load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
# }

####################################################################################################################################
### Handley analysis of body mass distributions
print("Beginning median body mass Handley analysis...")

# bmBreaks <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
bmBreaks <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
break.mat <- sapply(bmBreaks, function(x) x + runif(n=length(repIntOccs), min=-0.2, max=0.2))
class.vec.bm <- apply(break.mat, 1, function(this.bm.breaks) array(sapply(measure.mat$bodyMass, function(x) as.integer(max(which(x > this.bm.breaks)))), dimnames=list(rownames(measure.mat))))

# bmBreaks <- c(-Inf, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, Inf) #Badgely and Fox 2000
# bmBreaks <- hist(measure.mat$bodyMass, plot=FALSE)$breaks

# countCube <- sapply(repIntTaxa, function(this.rep) {
	# sapply(this.rep, function(this.intv, this.rep) {
		# hist(measure.mat$bodyMass[match(this.intv, measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts
		# }, this.rep=this.rep)
	# }, simplify = "array")
# dimnames(countCube)[[1]] <- 1:(length(bmBreaks)-1)
	
# countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 
optList_bm_median <- doHandleyTest(master.class.list=as.factor(apply(class.vec.bm, 1, median, na.rm=TRUE)), master.range.mat=range.mat.median, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)

# optList_bm_median <- doHandleyTest(this.counts=countBox[2,,], this.range.frame=appendBMToOneRanges(this.ranges=apply(taxRangeCube, c(1,2), median, na.rm=TRUE)), intervals=intervals, do.heuristic=do.heuristic, extra.intvs=extra.intvs)
	# quartz(width=3.3, height=9.8)
	# par(mfrow=c(ncol(countBox[2,,])/2,2), mar=c(0,3,0.5, 0.5), mfg=c(2,1))
	# for (i in seq(from=1, to=ncol(countBox[2,,]), by=1)) barplot(countBox[2,,][,i], width=c(0.68897, 0.68897, 0.778151, 0.522879, 1.18786),space=0, cex.axis=0.5, ylim=c(0,30))
	
	# quartz()
	# thisTab <- table(unlist(sapply (optList_bm_allReps, function(x) x[[(length(x) - 1)]]$optBreaks)))
	# thisTab <- array(thisTab[match(seq_len(nrow(intervals)), names(thisTab))], dimnames=list(rownames(intervals)))
	# barplot(rev(thisTab)/settings$n.reps, cex.names=0.5, ylim=c(0,1))
	# abline(h=c(0.95, 0.75, 0.5), lty=3, col="gray50")

class.vec.bm <- data.frame(apply(class.vec.bm, 2, function(x) as.factor(array(x[dimnames(taxRangeBox)[[1]]], dimnames=list(dimnames(taxRangeBox)[[1]])))), stringsAsFactors=TRUE)
print("Beginning body mass Handley analysis for all n.reps...")
optList_bm_allReps <- list()
for (this.rep in seq(from=1, to=100)) {
	# optList_bm_allReps[[this.rep]] <- doHandleyTest(this.counts=countCube[,,this.rep], this.range.frame=appendBMToOneRanges(this.ranges=taxRangeCube[,,this.rep]), intervals=intervals, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)
	optList_bm_allReps[[this.rep]] <- doHandleyTest(master.class.list=class.vec.bm[,this.rep], master.range.mat=taxRangeBox[,,this.rep], do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)
	if (this.rep %% 100 == 0) cat("Body Mass Handley Rep:",this.rep, "\n")
}

####################################################################################################################################
if(Sys.info()["sysname"] == "Darwin"){
	save(repIntTaxa, repIntOccs, optList_bm_median, optList_bm_allReps, file=paste0("~/Desktop/EcologyResults/BM_handleyResult_SampleStandardized=", settings$do.subsample, timestamp(),".Rdata"))
	# save(repIntTaxa, repIntOccs, optList_bm_median, optList_bm_allReps, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", settings$do.subsample, timestamp(),".Rdata"))
	#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
} else if(Sys.info()["sysname"] == "Windows"){
	save(repIntTaxa, repIntOccs, optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", settings$do.subsample, timestamp(),".Rdata"))
	# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
}

# save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/handleyResult", timestamp(),".Rdata"))
####################################################################################################################################

####################################################################################################################################
quartz(width=10)
	### Handley-block histogram series
	this.rep <- 1
	old.shift <- nrow(intervals)
	shift.ints <- optList_bm_median[[length(optList_bm_median) - 1]]$optBreaks
	# shift.ints <- rev(which(intervals$ageBase %in% optList_bm_allReps[[this.rep]][[length(optList_bm_allReps[[this.rep]]) - 1]]$optBreaks))
	par(mfrow=c(1,length(shift.ints)+1), mar=c(4,2,0.5,0.5), col.axis="gray50", col.lab="gray50", fg="gray50")
	hist.breaks <- c(min(measure.mat$bodyMass, na.rm=TRUE),bmBreaks[2:6],max(measure.mat$bodyMass, na.rm=TRUE))
	for (this.shift in c(shift.ints, 0)) {
		hist(measure.mat$bodyMass[unique(unlist(repIntTaxa[[this.rep]][seq(from=old.shift, to=(this.shift+1))]))], freq=FALSE, breaks=hist.breaks, col=rainbow(length(hist.breaks)), main="", xlab="log Body Mass", ylab="", ylim=c(0,1))
	}

	# Number of shifts per rep
	quartz()
	par(mfrow=c(2,1), mar=c(5, 4, 4, 1))
	hist(sapply(optList_tax_allReps, function(x) length(x)), breaks=seq(-0.5, 10, 1.0), col="orchid4", main="Number of taxonomic shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
	hist(sapply(optList_bm_allReps, function(x) length(x)), breaks=seq(-0.5, 11.5, 1.0), col="firebrick4", main="Number of body mass shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
	
	# quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
####################################################################################################################################
	### number of Replicates with a shift in that interval
	quartz()
	par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
	# optList_tax_allReps <- optList_tax_allReps_heuristic
	# optList_tax_allReps <- optList_tax_allReps_full
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_tax_allReps, function(x) x[[length(x)]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,5), xlim=rev(range(intervals)), xaxp =c(55,5,10), ylim=c(0,settings$n.reps), main="Number of Replicates with a Taxonomic Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="orchid4", border="orchid1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,10),  xlim=c(55, 0), ylim=c(0,settings$n.reps), add=TRUE, cex=0.5)

	# par(mfrow=c(2,1), mar=c(3.5, 3.5, 1, 1))
	# optList_bm_allReps <- optList_bm_allReps_heuristic
	# optList_bm_allReps <- optList_bm_allReps_full
	### number of Replicates with a shift in that interval
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_bm_allReps, function(x) unique(x[[length(x)]]$optBreaks)))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=rev(range(intervals)), ylim=c(0,settings$n.reps), main="Number of Replicates with a Body Mass Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="firebrick4", border="firebrick1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=c(55, 0), ylim=c(0,settings$n.reps), add=TRUE, cex=0.5)


####################################################################################################################################


####################################################################################################################################
######
# make three-panel figure
#####
	quartz(width=6.89)
		par(mfrow=c(4,1), mar=c(0,4,0.5,0.5), mgp=c(2, 1,0))
	### isotope panel
		if(Sys.info()["sysname"] == "Darwin"){
			source('~/Dropbox/code/R/common_src/isotopes.R', chdir = TRUE)
			} else if(Sys.info()["sysname"] == "Windows"){
				source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
			}
		
	#source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
	# source("~/Dropbox/code/R/common_src/isotopes.R")
		
	optList_topes <- doTopesRateAnalysis(intervals, do.parallel=do.parallel)
	plotTopesRateAnalysis(optList_topes, intervals, x.axis=FALSE) #
	box(lwd=1)
	# getAlroyStatistics(intervals)
			
	### taxonomy panel
		par(mar=c(0,4, 2.5,0.5))
		
		# taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% measure.mat$taxon[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
		#taxCubeG <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$genus) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
		taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
		dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
		
		# prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
		prop <- t(apply(taxCube, c(1,2), mean, na.rm=TRUE))
		colnames(prop)[colnames(prop)==""] <- "indeterminate"
		# dimnames(prop) <- list(rownames(intervals), shortFam)
		source("https://dl.dropbox.com/s/iy0tu983xesbig2/taxonomicEv.R")
		plotStackedRichness(this.box=prop, intervals=intervals, do.log=FALSE, overlay.labels=TRUE, numbers.only=TRUE, legend=FALSE, xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)))
		#med.n <- median(length(unique(unlist(sapply(repIntTaxa[[this.rep]], function(x) measure.mat$taxon[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
		# med.n <- median(sapply(repIntTaxa, function(x) length(unique(unlist(sapply(x, function(y) measure.mat$taxon[y]))))))
		# optList_tax <- doHandleyTest(this.counts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
		abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median)]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
		box(lwd=1)
	
	### body mass panel
		thisRanges <- getTaxonRangesFromOccs(occs=occs, this.rank=this.rank, random=FALSE)
		rownames(thisRanges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisRanges))
		measure.mat[,c("FO","LO")] <- thisRanges[match(measure.mat$taxon, rownames(thisRanges)),]

		par(mar=c(0,4,2.5,0.5))
		# quartz(width=12, height=6)
		plot(measure.mat$FO, measure.mat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(55,5,10), xlab="Time (Ma)", cex.axis=1.5, cex.lab=1.5)
		# plot(measure.mat$FO, measure.mat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(50,0,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, col="gray75", fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #alter xaxpto change x-axis values
		# rect(-10e6, -10e6, 10e6, 10e6, col="white")
		overlayCzTimescale(do.subepochs=TRUE)
		
		famColors <- rainbow(length(shortFam))
		colorList <- famColors[match(bigList$family[as.character(bigList$accepted_name) %in% measure.mat$taxon], shortFam)]
		colorList[is.na(colorList)] <- "gray25"
		
		orderColors <- array(NA, dim=nrow(measure.mat))
		# orderColors[bigList$order[match(measure.mat$taxon, bigList$accepted_name)]=="Perissodactyla"] <- "dodgerblue4"
		# orderColors[bigList$order[match(measure.mat$taxon, bigList$accepted_name)] =="Artiodactyla"] <- "deeppink4"

		for (i in seq_len(nrow(measure.mat))) {
			# lines(x=c(this["FO"], x["LO"]), y=c(x["bodyMass"], x["bodyMass"]), lwd=3, pch=21, col=famColors[match(bigList[match(measure.mat$taxon, bigList[,1]),2], shortFam)])
			# lines(x=c(measure.mat$FO[i], measure.mat$LO[i]), y=c(measure.mat$bodyMass[i], measure.mat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor(colorList[i], 0.75))
			# lines(x=c(thisRanges[match(measure.mat$taxon[i], rownames(thisRanges)),"FO"], thisRanges[match(measure.mat$taxon[i], rownames(thisRanges)),"LO"]), y=c(measure.mat$bodyMass[i], measure.mat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor("gray0", 0.75)) #
			if (is.finite(measure.mat$FO[i]) & is.finite(measure.mat$LO[i]) & measure.mat$FO[i] != measure.mat$LO[i]) lines(x=measure.mat[i,c("FO","LO")], y=c(measure.mat$bodyMass[i], measure.mat$bodyMass[i]), lwd=0.75, pch=21, col=alphaColor("gray0", 0.5)) #alphaColor(orderColors[i], 0.5)
		}
		points(measure.mat[complete.cases(measure.mat[ ,c("FO","LO")]) & measure.mat$FO == measure.mat$LO, c("FO","bodyMass")], pch=21, col=alphaColor("gray0", 0.5), cex=0.25) #this line is not generating the proper output for the final graph due to c("FO","bodyMass") causing a  "undefined columns selected" error
		
		# optList_bm <- doHandleyTest(this.counts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
		# optList_bm <- doHandleyTest(this.counts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
		abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median)]]$optBreaks,2], range(intervals))), lwd=1.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)]]$optBreaks) + 1)), pos=3, cex=0.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)]]$optBreaks,2])), "Ma"), adj=c(0,0),cex=0.5, col="firebrick4")
		
		quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
		polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
		polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor("darkorange4", 0.25), border="darkorange4")
		lines(rowMeans(intervals), quants[3,], col=alphaColor("goldenrod1", 0.5), lwd=5)
		lines(rowMeans(intervals), quants[3,], col=alphaColor("darkorange4", 1.0), lwd=3)
		points(rowMeans(intervals), quants[3,], col=alphaColor("darkorange1", 0.5), cex=0.5)
		box(lwd=1)

################################################################################################################################################################

		endTime <- Sys.time()
		RunTime2Ma <- endTime - startTime

		###Compile regime distribution histograms and net/change histograms using the range of taxon occurance
#		regimeHist(repIntTaxa = repIntTaxa, breaks = c(51,47,37,21, 5,2), optList=optList_bm_median, measure.mat = measure.mat, netFreq = TRUE, regimeFreq=FALSE,
#							 netPlotType = "pos/neg", plot.together = FALSE)
		
		###Compile regime distribution histograms and net/change histograms using the midpoint of taxon occurance
		##this was used to produce the regime and net change histograms for publication
		#get midpoint of each taxon 
		taxonranges <- getTaxonRangesFromOccs(occs = occs, random = FALSE)
		rownames(taxonranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(taxonranges))
		taxonranges <- cbind(taxonranges, (taxonranges[,"FO"]+taxonranges[,"LO"])/2); colnames(taxonranges)[3] <- "MP"
		occDatesMP <- taxonranges[,"MP"]
		intSpMP <- apply(intervals, 1, function(thisIntv) taxonranges[taxonranges[,"MP"] > thisIntv[1] & taxonranges[,"MP"] <= thisIntv[2],])
		countCubeMP <- sapply(intSpMP, function(x) hist(measure.mat$bodyMass[match(rownames(x), measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts, simplify = "array")
		#drop the intervals that include the Quaternary (<3 Ma)
		countCubeMP <- countCubeMP[,!as.double(str_remove(colnames(countCubeMP), " Ma")) < 3]
		optList_bm_medianMP <- doHandleyTest(countCubeMP, n=nrow(measure.mat), do.heuristic=do.heuristic, extra.intvs=extra.intvs)
		regimeHist_countBox(countBox = countCubeMP, breaks = c(51,47,37,21), optList=optList_bm_medianMP, 
												measure.mat = measure.mat, netFreq = TRUE, regimeFreq=FALSE,
												netPlotType = "pos/neg", plot.together = FALSE, plot.axes = TRUE,
												grayscale = FALSE)
		
		
		##Get frequency of combinations for breaks
		bm_allReps <- BreakComboFreq(optList_bm_allReps)   #cannot find function count error
		tax_allReps <- BreakComboFreq(optList_tax_allReps)
		bm_allReps[order(-bm_allReps$freq),]
		tax_allReps[order(-tax_allReps$freq),]
		
######################################################################################################
    #load final
    #RAW
    load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=FALSE##------ Sun Mar 31 15:24:34 2019 ------##.Rdata")
    load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=FALSE##------ Sun Mar 31 15:17:50 2019 ------##.Rdata")
    
    #Subsampled
    load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=TRUE##------ Fri Mar 29 22:06:27 2019 ------##.Rdata")
    load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=TRUE##------ Fri Mar 29 22:00:47 2019 ------##.Rdata")
    
    ##Get frequency of combinations for breaks
    bm_allReps <- BreakComboFreq(optList_bm_allReps)   #cannot find function count error
    tax_allReps <- BreakComboFreq(optList_tax_allReps)
    bm_allReps[order(-bm_allReps$freq),]
    tax_allReps[order(-tax_allReps$freq),]
    
    #taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$genus) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
    taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
    dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
    taxPlio_median <- apply(taxCube, c(1,2), median); taxPlio_median[,c("3.5 Ma", "4.5 Ma", "5.5 Ma", "6.5 Ma")]
    taxPlio_min <- apply(taxCube, c(1,2), min); taxPlio_min[,c("3.5 Ma", "4.5 Ma", "5.5 Ma", "6.5 Ma")]
    taxPlio_max <- apply(taxCube, c(1,2), max); taxPlio_max[,c("3.5 Ma", "4.5 Ma", "5.5 Ma", "6.5 Ma")]
    
    #1-Ma bins c(50,46,40,26,16,2)
    #2-Ma bins c(51,47,37,21, 5,2)
    regimeHist(repIntTaxa = repIntTaxa, breaks = c(51,47,37,21, 5,2), optList=optList_bm_median, measure.mat = measure.mat, netFreq = TRUE, regimeFreq=FALSE,
    					 netPlotType = "pos/neg", plot.together = FALSE)
    
    #get midpoint of each taxon 
    taxonranges <- getTaxonRangesFromOccs(occs = occs, random = FALSE)
    rownames(taxonranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(taxonranges))
    taxonranges <- cbind(taxonranges, (taxonranges[,"FO"]+taxonranges[,"LO"])/2); colnames(taxonranges)[3] <- "MP"
    occDatesMP <- taxonranges[,"MP"]
    intSpMP <- apply(intervals, 1, function(thisIntv) taxonranges[taxonranges[,"MP"] > thisIntv[1] & taxonranges[,"MP"] <= thisIntv[2],])
    #intSpMP <- lapply(intSpMP, function (y) gsub(pattern = "[[:space:]]", replacement = "_", x = y))
    
    #intSpMP <- sapply(intOccsMP, function(x) sort(unique(as.character(gsub(pattern = "[[:space:]]", 
    #																																			 replacement = "_", x = occs$accepted_name[occs$occurrence_no %in% x])))))
    countCubeMP <- sapply(intSpMP, function(x) hist(measure.mat$bodyMass[match(rownames(x), measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts, simplify = "array")
    #drop the intervals that include the Quaternary (<3 Ma)
    countCubeMP <- countCubeMP[,!as.double(str_remove(colnames(countCubeMP), " Ma")) < 3]
    optList_bm_medianMP <- doHandleyTest(countCubeMP, n=nrow(measure.mat), do.heuristic=do.heuristic, extra.intvs=extra.intvs)
    regimeHist_countBox(countBox = countCubeMP, breaks = c(51,47,37,21), optList=optList_bm_medianMP, 
    										measure.mat = measure.mat, netFreq = TRUE, regimeFreq=FALSE,
    										netPlotType = "pos/neg", plot.together = FALSE, plot.axes = TRUE,
    										grayscale = FALSE)
    
    #regimeHist_HistMedian(repIntTaxa = repIntTaxa, breaks = c(50,46,40,26,16,2), optList=optList_bm_median, measure.mat = measure.mat, netFreq = TRUE, regimeFreq=FALSE,
    #											netPlotType = "pos/neg", plot.together = FALSE)
   
    
    # write function to lot distribtion of body mass within each break interval
    #Jons code for histograms below: run tests to see how they change iteration to iteration#make into function 
    this.rep <- 1
    old.shift <- nrow(intervals)
    optList <- optList_bm_allReps
    shift.ints <- rev(which(intervals$ageBase %in% optList[[this.rep]][[length(optList[[this.rep]])]]$optBreaks))
    par(mfrow=c(1,length(shift.ints)+1), mar=c(4,2,0.5,0.5))
    for (this.shift in c(shift.ints, 0)) 
    {
    	hist(measure.mat$bodyMass[unique(unlist(repIntTaxa[[this.rep]][seq(from=old.shift, to=(this.shift+1))]))], freq=TRUE, breaks=c(min(measure.mat$bodyMass, na.rm=TRUE),bmBreaks[2:6],max(measure.mat$bodyMass, na.rm=TRUE)), col=rainbow(n=length(shift.ints)+1), main="", xlab="log Body Mass", ylab="", ylim=c(0,1))
    }
    
    
    
     
###
#Compiled Histogram of Proportional shifts in body mass
###

 #install.packages("ggplot2")
 #require(ggplot2)
 #install.packages("Rcmdr")
 #require(Rcmdr)
 #quartz(width=6.89)
	# par(mfrow=c(2,1), mar=c(4,4, 1,0.5), mgp=c(2, 1,0))
	#	plot(measure.mat$FO, measure.mat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="% Species Diversity", xaxp =c(55,5,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, fg="black", bg="black", col.axis="black", col.lab="black")
	#	rect(-10e6, -10e6, 10e6, 10e6, col="white")
	#	overlayCzTimescale(do.subepochs=TRUE)
		
		# #assign species bin from Janis 2000
		# #need a vector like countCube but made with raw data prior to distributions
		# #get list of taxa throughout all intervals by running single run of species interval and rangethrough
	#	countInt
				
		# #load("/Users/evandoughty/Dropbox/ungulate_RA/RCode/EcologyAnalysisResults/countInt.RData")
		# #countInt <- countInt_bin1myr
		# #countInt <- countInt_bin2myr		
			# y <- countInt
		# #countDist <- sapply(countInt, function(y) sapply(y, function(x) hist(this.column[x], breaks=breaks, plot=FALSE)$counts), simplify = "array")
		# countDist <- sapply(countInt, function(x) hist(this.column[x], breaks=breaks, plot=FALSE)$counts)
		# #dimnames(countDist)
		
		# r <- matrix()
		# r <- countDist
		# colnames(r) <- intervals$ageTop
		# a <- r[,ncol(r):1]
		# r <- a
				
		# #Need to Remove Intervals from Handley Bins and into actual intervals
		# cols <- palette(rainbow(6))
		# distBar <- barplot(r, ylim=c(0,100), xlab="Time(Ma)", ylab = "Number of taxa", col =cols, space = 0)

	# par(mar=c(6,4,2.5,0.5))	
		# summedDist <- colPercents(r)
		# sumBar <- rev(barplot(summedDist[-c(7:9),], col = cols, xlab="Time(Ma)", ylab = "% Species Diversity"))
		
	
# ####
# #NMDS
# ###
# install.packages("vegan")
# require(vegan)
# bmBox <- sapply(repIntTaxa, function(thisRep) t(sapply(thisRep, function(taxa) hist(measure.mat$bodyMass[taxa], freq=TRUE, xlim=c(min(bmBreaks), max(bmBreaks)), ylim=c(0,22), breaks=bmBreaks, col="firebrick1", ylab="", plot=settings$plotHist)$counts)), simplify="array") #, main=paste(intervals[intv, "ageTop"], "-", intervals[intv, "ageBase"], "Ma (n = ", length(intSp[[intv]]), ")")
# par(mfrow=c(2,1), mar=c(3,3,0.5,0.5), mgp=c(1.25, 0.5, 0))
# tax <- plotNMDS(thisBox=t(apply(taxCube, c(1,2), mean, na.rm=TRUE)), intervals, polygon.ints=optList_tax[[(length(optList_tax)-1)]]$optBreaks, scaler=5, title="taxonomy") #, filename="~/Desktop/NMDS_bodyMass.pdf"
# bm <- plotNMDS(thisBox=bmBox, intervals, polygon.ints=oneOpt[[(length(oneOpt)-1)]]$optBreaks, scaler=5, title="Body Mass") #, filename="~/Desktop/NMDS_bodyMass.pdf"
