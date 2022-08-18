# # setwd("C:/Users/Evan/Dropbox/ungulate_RA/RCode")

# #all ungulates
#load("~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata")

# #artio only
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/artio/handleyResult##------ Wed Nov 15 16:57:31 2017 ------##_artiodactyla.Rdata")

# #perisso only
# load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/perisso/handleyResult##------ Thu Nov 16 10:27:58 2017 ------##.Rdata")
 startTime <- Sys.time()
# load('~/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/repIntTaxa_SampleStandardized=TRUE##------ Fri Mar 29 21:20:21 2019 ------##.Rdata')

#sources for Jon Marcot's code and specimen measurements 
source("~/Dropbox/code/R/common_src/strat.R")
source("~/Dropbox/code/R/common_src/occFns.R")
source("~/Dropbox/code/R/common_src/sampling.R") 
source("~/Dropbox/code/R/common_src/utils_marcot.R")
source("~/Dropbox/code/R/common_src/CzTimescale.R") 

source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
source('~/Dropbox/code/R/dentalMeasurements/src/src_ecologyAnalysisFns.R', chdir = TRUE)

source('~/Dropbox/Code/Recreate Dental Plots Source 2021_6_2.R')

####################################################################################################################################

  occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
  #occs <- read.csv("/Users/emdoughty/Dropbox/Code/Occs_2021_9_5.csv")
  occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
  occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                   "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                   "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", 
                                   "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
  occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

################################################  
#Settings
run.taxon <- "herbivores" #"ungulates"
this.rank <- "genus" #"genus"
interval.type <- "bins" #"nalma"
do.parallel <- TRUE
if (do.parallel) require(parallel)
save.files <- TRUE

herbivore.size.cat <- c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, Inf) #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)

predator.size.cat  <- c(-Inf, 0, 0.845098, 1.322219, 2, Inf) #includes <1kg, adn 1-7 kg categories
#predator.size.cat  <- c(-Inf,0.845098,1.322219,2,Inf) #<7 kg as single category
pred.cat.name <- "allCateg" #"sub7kg"

#the primary location of where the output files will go.  Make sure this only includes the file destination as the rest of the filename is concatenated in its respective section.
save.pathname <- "~/Dropbox/Code/R/Results/"

#toggle so one can fire and forget without having to go back and forth hunting for each section individually
analysis.toggle <- c("bmHandley") #"taxHandley",

#if you want to load a repIntOccs or repIntTaxa from file put the pathname as this object.  otherwise keep as NUll to make a new repIntOccs and repIntTaxa using the settings below
repIntLoad <- NULL #"/Users/emdoughty/Dropbox/Code/R/Results/repIntMaster__this.rank=species_timebin=2Mabins_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Fri Mar  4 21:58:08 2022 ------##.Rdata"
  #NULL 
  #"/Users/emdoughty/Dropbox/Code/R/Results/repIntMaster__this.rank=species_timebin=nalma_SampleStandardized=TRUE_Reps=10000Jonathans_MBP.lan##------ Thu Mar  3 18:05:02 2022 ------##.Rdata"
  #NULL

################################################

######################################################################################################################################################
  #Ungualtes
##################
if(run.taxon == "ungulates")
{
  measure.mat <- getMeasureMatWithBodyMasses()
  
  ####################################################################################################################################
  #### reduces matrix to just the focal order(s)
  ####################################################################################################################################
  
  archaic.ung <- read.csv("/Users/emdoughty/Dropbox/Code/ArchaicUngulate_UploadFile_2021_4_29.csv")
  
  archaic.Mat.tot <- getMeasureMatCondylarths(data.raw = archaic.ung, occs = occs, 
                                              col.order = colnames(measure.mat), 
                                              all.bm = FALSE,
                                              regression = "ArchaicNonselenodonts")
  #test <- jon_archicBodyMassEstim()
  #test.bm <- 10^test$bodyMass
  #test.bm <- test.bm/1000
  #test.bm <- log10(test.bm)
  
  measure.mat <- rbind(measure.mat, archaic.Mat.tot)
  
  if (this.rank=="genus") measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat)
  
  focal.order <- c("Artiodactyla", "Perissodactyla")
  focal.family <- unique(occs[occs$order %in% focal.order,]$family)
  #search through those without order
  add.family <- c("Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
  focal.family <- c(as.character(focal.family), add.family)
  focal.family <- focal.family[!focal.family %in% ""]
  focal.family <- focal.family[order(focal.family)]
  
  bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
  bigList <- rbind(bigList,bigList.cond)
  
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  # bigList[order(bigList$family, bigList$accepted_name),]
  shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	
  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$family %in% shortFam], ]
}
#################

#################
#Predators
################
if(run.taxon == "carnivores")
{
  pred.data <- read.csv("~/Dropbox/Proposal/predator_data_final.csv")
  colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl", "bodyMass", 	"Citation") #"BM_all_carnivoran","BM_extant_reg")
  pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)
  rownames(pred.data) <- pred.data$taxon
  
  pred.data[,c("bodyMass")] <- log10(pred.data[,c("bodyMass")])
  pred.data <- pred.data[is.finite(pred.data$bodyMass),]
  
  if (this.rank=="genus") pred.data <- makeOneGenusMatFromSpecimenMat(pred.data)
  
  focal.order <- c("Carnivora", "Creodonta","Hyaenodonta")
  #occsPred <- occs[occs$order %in% focal.orderPred,]
  focal.family <- unique(occs[occs$order %in% focal.order,]$family)
  add.family <- c("Viverravidae")
  
  focal.family <- c(as.character(focal.family), add.family)
  focal.family <- focal.family[!focal.family %in% ""]
  focal.family <- focal.family[order(focal.family)]
  
  bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList.cond <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% add.family), c("order","family", "genus", "accepted_name")])
  bigList <- rbind(bigList,bigList.cond)
  
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  # bigList[order(bigList$family, bigList$accepted_name),]
  shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	
  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  measure.mat <- pred.data[pred.data$taxon %in% bigList$accepted_name[bigList$family %in% shortFam], ]
  measure.mat <- measure.mat[is.finite(measure.mat$bodyMass),]
}

###############################################################################################################################################

if(interval.type == "bins")
{
  int_length <- 2
  intervals <- makeIntervals(0, 64, int_length)
  intList <- listifyMatrixByRow(intervals)
  
  save.path.bins <- paste0(int_length,"Ma",interval.type)
}
#######################

if(interval.type == "nalma")
{
  nalma.mark <- read.csv("/Users/emdoughty/Dropbox/Proposal/NOW_intervals_edit.csv")
  nalma.mark <- nalma.mark[,1:3]
  #nalma.add  <- rbind(c("Aquilian", 84, 70),c("Lancian", 70, 66),c("Puercan", 66, 64.81)); colnames(nalma.add) <- colnames(nalma.mark)
  # nalma.mark <- rbind(nalma.mark, nalma.add)
  #  nalma.mark[,2] <- as.numeric(nalma.mark[,2])
  # nalma.mark[,3] <- as.numeric(nalma.mark[,3])
  # nalma.mark <- nalma.mark[order(as.numeric(nalma.mark$Max_age), decreasing = TRUE),]
  rownames(nalma.mark) <- nalma.mark$CHRON
  colnames(nalma.mark) <- c("NALMA_Subdivision", "ageBase","ageTop")
  nalma.mark <- nalma.mark[,-1]
  nalma.mark <- nalma.mark[,c(2,1)]
  intervals <- nalma.mark
  
  save.path.bins <- paste0(interval.type)
}

##############################################################################################################################################

if(is.null(repIntLoad))
{
  do.parallel <- TRUE
  	if (do.parallel) require(parallel)
  reps <- 10000
  do.subsample <- TRUE
  quota <- 0.4
  do.disparity <- FALSE
  bootstrapSpecimens <- FALSE
  bootstrapSpecies <- FALSE
  bootstrapSpeciesWithinIntervals <- FALSE
  plotHist <- FALSE
  do.heuristic <- TRUE
  	extra.intvs <- 0
  do.rangethrough <- TRUE
  save.files <- TRUE
  
  if (bootstrapSpecies) holderMat <- measure.mat
  
  if (plotHist) {
  	quartz("Guild Histograms")
  	par(mfrow=c((nrow(intervals)), 3), mar=c(0,0,0.75,0), cex.axis=0.5, cex.main=0.75)
  }
  
  repIntOccs <- list()
  
  ################################################################################################################################################
  #### get species within intervals
  ################################################################################################################################################
  
  for (rep in seq_len(reps)) {
  	cat("Beginning Rep", rep, "of", reps, "...\r")
  	##################################################We need to update this sbootstrap section
  	if (bootstrapSpecimens) {
  		measure.mat <- specimenMat[sample.int(nrow(specimenMat), size=nrow(specimenMat), replace=TRUE),]
  		measure.mat <- aggregate(measure.mat, by=list(taxon=specimenMat$taxon), mean, na.rm=TRUE)
  		# measure.mat <- measure.mat[,apply(!sapply(measure.mat, is.na), 2, any)]
  		rownames(measure.mat) <- measure.mat$taxon
  		measure.mat[sapply(measure.mat, is.nan)] <- NA
  		# measure.mat<-cbind(measure.mat, cbind(FO=vector(length=nrow(measure.mat), mode="numeric"), LO=vector(length=nrow(measure.mat), mode="numeric")))
  		measure.mat[,"reg"] <- as.character(famList$reg[match(measure.mat$taxon,famList$taxon)])
  		measure.mat[,"bodyMass"] <- makeBodyMasses(measure.mat, regList, best.only=TRUE)
  		# measure.mat[,"PC2"] <- pcVec[match(measure.mat$taxon, names(pcVec))]
  		# measure.mat[,"PC3"] <- pcaLo$x[match(measure.mat$taxon, rownames(pcaLo$x)),3]
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
  
  repIntTaxa <- getRepIntTaxaFromRepIntOccs(repIntOccs, this.rank=this.rank, do.rangethrough=do.rangethrough)
  
  print("Completed getting taxa with intervals")
  
  ###################################################################################################################################
  if(save.files)
  {
  	if(Sys.info()["sysname"] == "Darwin"){
  	  if(run.taxon == "carnivores") run.taxon <- paste0(run.taxon, "_",pred.cat.name)
  	  save(repIntTaxa, repIntOccs, intervals, reps, do.subsample, quota,do.disparity, 
  	       bootstrapSpecimens,bootstrapSpecies,bootstrapSpeciesWithinIntervals ,
  	       plotHist,do.heuristic,extra.intvs,do.rangethrough,
  	       file=paste0(save.pathname,"repIntMaster_",
  	                                           "_this.rank=", this.rank,
  	                                           "_timebin=", save.path.bins,
  	                                           "_SampleStandardized=", do.subsample, 
  	                                           "_Reps=", reps, gsub("-","_",Sys.info()["nodename"]),
  	                                           timestamp(),".Rdata"))
  		#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
  	} else if(Sys.info()["sysname"] == "Windows"){
  		save(repIntTaxa, repIntOccs, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/repIntTaxa_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
  		# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
  	}
  }
}

if(!is.null(repIntLoad)) load(repIntLoad)
####################################################################################################################################
### Handley analysis of taxonomic distributions
if("taxHandley" %in% analysis.toggle){
  print("Beginning median taxonomic Handley analysis...")
  
  # bigList <- bigList[bigList$order %in% focal.order,]
  # shortFam <- sort(unique(bigList$family))
  
  taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
  med.n <- median(sapply(repIntTaxa, function(x) length(unique(unlist(sapply(x, function(y) y))))))
  optList_tax_median <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	
  
  print("Beginning taxonomic Handley analysis for all reps...")
  optList_tax_allReps <- list()
  for (this.rep in seq_len(reps)) {
  	taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  	this.n <- length(unique(unlist(sapply(repIntTaxa [[this.rep]], function(x) x))))
  	optList_tax_allReps[[this.rep]] <- doHandleyTest(taxCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	
  	if(this.rep %% 100 == 0) cat("Taxonomic Handley Rep:", this.rep, "\n")
  }
  
  ####################################################################################################################################
  if(save.files)
  {
    if(Sys.info()["sysname"] == "Darwin"){
      if(run.taxon == "carnivores") run.taxon <- paste0(run.taxon, "_",pred.cat.name)
      save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, 
           file=paste0(save.pathname,"Taxon_handleyResult_", run.taxon,
                                               "_this.rank=", this.rank,
                                               "_timebin=", save.path.bins,
                                               "_SampleStandardized=", do.subsample, 
                                               "_Reps=", reps, gsub("-","_",Sys.info()["nodename"]),
                                               timestamp(),".Rdata"))
  
    	#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    } else if(Sys.info()["sysname"] == "Windows"){
    	save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, 
    	     file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_", run.taxon,
  		                 "_SampleStandardized=", do.subsample, 
  		                 "_Reps=", reps, timestamp(),".Rdata"))
    	                                                                       
    	# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    }
  }
}
####################################################################################################################################
### Handley analysis of body mass distributions
if("bmHandley" %in% analysis.toggle) {
  print("Beginning median body mass Handley analysis...")
  if(run.taxon == "ungulates") bmBreaks <- herbivore.size.cat
  if(run.taxon == "carnivores") bmBreaks <- predator.size.cat
  #bmBreaks <- c(-Inf, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, Inf) #Badgely and Fox 2000
  # bmBreaks <- hist(measure.mat$bodyMass, plot=FALSE)$breaks
  
  countCube <- sapply(repIntTaxa, function(this.rep) {
  	sapply(this.rep, function(this.intv, this.rep) {
  		hist(measure.mat$bodyMass[match(this.intv, measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts
  		}, this.rep=this.rep)
  	}, simplify = "array")
  	
  countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 
  
  optList_bm_median <- doHandleyTest(countBox[2,,], n=nrow(measure.mat), do.heuristic=do.heuristic, extra.intvs=extra.intvs)
  	# quartz(width=3.3, height=9.8)
  	# par(mfrow=c(ncol(countBox[2,,])/2,2), mar=c(0,3,0.5, 0.5), mfg=c(2,1))
  	# for (i in seq(from=1, to=ncol(countBox[2,,]), by=1)) barplot(countBox[2,,][,i], width=c(0.68897, 0.68897, 0.778151, 0.522879, 1.18786),space=0, cex.axis=0.5, ylim=c(0,30))
  	
  	# quartz()
  	# thisTab <- table(unlist(sapply (optList_bm_allReps, function(x) x[[(length(x) - 1)]]$optBreaks)))
  	# thisTab <- array(thisTab[match(seq_len(nrow(intervals)), names(thisTab))], dimnames=list(rownames(intervals)))
  	# barplot(rev(thisTab)/reps, cex.names=0.5, ylim=c(0,1))
  	# abline(h=c(0.95, 0.75, 0.5), lty=3, col="gray50")
  
  print("Beginning body mass Handley analysis for all reps...")
  optList_bm_allReps <- list()
  for (this.rep in seq_len(reps)) {
  	this.n <- length(unique(unlist(sapply(repIntTaxa [[this.rep]], function(x) x))))
  	optList_bm_allReps[[this.rep]] <- doHandleyTest(countCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)
  	if(this.rep %% 100 == 0) cat("Body Mass Handley Rep:",this.rep, "\n")
  }
  
  ####################################################################################################################################
  if(save.files)
  {
    if(Sys.info()["sysname"] == "Darwin"){
      if(run.taxon == "carnivores") run.taxon <- paste0(run.taxon, "_",pred.cat.name)
      save(repIntTaxa, repIntOccs, optList_bm_median, optList_bm_allReps,
           file=paste0(save.pathname,"BM_handleyResult_", run.taxon,
                                               "_this.rank=", this.rank,
                                               "_timebin=", save.path.bins,
                                               "_SampleStandardized=", do.subsample, 
                                               "_Reps=", reps, gsub("-","_",Sys.info()["nodename"]),
                                               timestamp(),".Rdata"))
      
    	#load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    	} else if(Sys.info()["sysname"] == "Windows"){
    		save(repIntTaxa, repIntOccs, optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
    		# load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    }
  }
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
	hist(sapply(optList_tax_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 10, 1.0), col="orchid4", main="Number of taxonomic shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
	hist(sapply(optList_bm_allReps, function(x) length(x) - 2), breaks=seq(-0.5, 11.5, 1.0), col="firebrick4", main="Number of body mass shifts in each rep", xlab="Number of Shifts", xlim=c(0,10))
	
  quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,"bodyMass"], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
####################################################################################################################################
	### number of Replicates with a shift in that interval
	quartz()
	par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
	# optList_tax_allReps <- optList_tax_allReps_heuristic
	# optList_tax_allReps <- optList_tax_allReps_full
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_tax_allReps, function(x) x[[length(x) - 1]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xlim=rev(range(intervals)), xaxp =c(65,5,10), ylim=c(0,reps), main="Number of Replicates with a Taxonomic Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="orchid4", border="orchid1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(65,5,10),  xlim=c(55, 0), ylim=c(0,reps), add=TRUE)

	# par(mfrow=c(2,1), mar=c(3.5, 3.5, 1, 1))
	# optList_bm_allReps <- optList_bm_allReps_heuristic
	# optList_bm_allReps <- optList_bm_allReps_full
	### number of Replicates with a shift in that interval
	breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList_bm_allReps, function(x) unique(x[[length(x) - 1]]$optBreaks)))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
	plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=rev(range(intervals)), ylim=c(0,reps), main="Number of Replicates with a Body Mass Distribution Shift", xlab="Time (Ma)")
	overlayCzTimescale(do.subepochs=TRUE)
	plot(breakHist, col="firebrick4", border="firebrick1", labels=TRUE, freq=TRUE, cex=0.3, xaxp =c(55,5,10), xlim=c(55, 0), ylim=c(0,reps), add=TRUE)


####################################################################################################################################


####################################################################################################################################
######
# make three-panel figure
#####
	quartz(width=6.89)
		par(mfrow=c(3,1), mar=c(0,4,0.5,0.5), mgp=c(2, 1,0))
	### isotope panel
		if(Sys.info()["sysname"] == "Darwin"){
			source('~/Dropbox/code/R/common_src/isotopes.R', chdir = TRUE)
			} else if(Sys.info()["sysname"] == "Windows"){
				source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
			}
		
	#source("C:/Users/Blaire/Dropbox/ungulate_RA/RCode/isotopes.R")
	# source("~/Dropbox/code/R/common_src/isotopes.R")
		
	optList_topes <- doTopesRateAnalysis(intervals)
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
		plotStackedRichness(this.box=prop, intervals=intervals, do.log=FALSE, overlay.labels=TRUE, numbers.only=TRUE, legend=TRUE , xlim=c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE)))
		#med.n <- median(length(unique(unlist(sapply(repIntTaxa[[this.rep]], function(x) measure.mat$taxon[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
		# med.n <- median(sapply(repIntTaxa, function(x) length(unique(unlist(sapply(x, function(y) measure.mat$taxon[y]))))))
		# optList_tax <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
		abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
		text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
		box(lwd=1)
	
	### body mass panel
		thisRanges <- getTaxonRangesFromOccs(occs=occs, random=FALSE)
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
		
		# optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
		# optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
		abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="firebrick4")
		text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0),cex=0.5, col="firebrick4")
		
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
		
if(FALSE)
{
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
# bmBox <- sapply(repIntTaxa, function(thisRep) t(sapply(thisRep, function(taxa) hist(measure.mat$bodyMass[taxa], freq=TRUE, xlim=c(min(bmBreaks), max(bmBreaks)), ylim=c(0,22), breaks=bmBreaks, col="firebrick1", ylab="", plot=plotHist)$counts)), simplify="array") #, main=paste(intervals[intv, "ageTop"], "-", intervals[intv, "ageBase"], "Ma (n = ", length(intSp[[intv]]), ")")
# par(mfrow=c(2,1), mar=c(3,3,0.5,0.5), mgp=c(1.25, 0.5, 0))
# tax <- plotNMDS(thisBox=t(apply(taxCube, c(1,2), mean, na.rm=TRUE)), intervals, polygon.ints=optList_tax[[(length(optList_tax)-1)]]$optBreaks, scaler=5, title="taxonomy") #, filename="~/Desktop/NMDS_bodyMass.pdf"
# bm <- plotNMDS(thisBox=bmBox, intervals, polygon.ints=oneOpt[[(length(oneOpt)-1)]]$optBreaks, scaler=5, title="Body Mass") #, filename="~/Desktop/NMDS_bodyMass.pdf"
}