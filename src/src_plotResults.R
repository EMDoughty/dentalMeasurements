plotTalkFigs_NAPC2019 <- function()
{
  #body mass and k/s set-up
  getEcoAnalSrc()
  occs <- getOccs()
  thisMat <- getMeasureMat()
  bigList <- getDataLists_Eco_OrigExt(object = "bigList", occs = occs, thisMat = thisMat)
  shortFam <- getDataLists_Eco_OrigExt(object = "shortFam", occs = occs, thisMat = thisMat)
  thisMat <- getDataLists_Eco_OrigExt(object = "thisMat", occs = occs, thisMat = thisMat)
  
  int_length <- 2
  intervals <- makeIntervals(1, 55, int_length)
  intList <- listifyMatrixByRow(intervals)
  
  #load repIntSp, taxon hadley, and body mass hadley results
  load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/RepIntSp_SampleStandardized=TRUE##------ Fri Mar 29 21:20:21 2019 ------##.Rdata")
  load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/Taxon_handleyResult_SampleStandardized=TRUE##------ Fri Mar 29 22:00:47 2019 ------##.Rdata")
  load("/Users/emdoughty/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/BM_handleyResult_SampleStandardized=TRUE##------ Fri Mar 29 22:06:27 2019 ------##.Rdata")
  
  
  load("~/Dropbox/ungulate_RA/EcologyResults/repIntTaxa_SampleStandardized=TRUE##------ Wed Jun 19 16:49:25 2019 ------##.Rdata")
  load("~/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=TRUE##------ Wed Jun 19 17:38:37 2019 ------##.Rdata")
  load("~/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=TRUE##------ Wed Jun 19 17:44:48 2019 ------##.Rdata")
  
  print("Setup Done")
  ############################################################################################################
  #3panel
  plot3panel(repIntSp=repIntOccs, bigList = bigList, shortFam = shortFam, intervals = intervals, optList_tax_median = optList_tax_median, optList_bm_median = optList_bm_median)
  
  print("3panel")
  #body mass sholder plot (no breaks added)
  plotSinglePanelPlot(panel="bodyMass", includeBreaks = FALSE,
                      repIntSp=repIntSp, bigList = bigList, shortFam = shortFam, intervals = intervals, optList_tax_median = optList_tax_median, optList_bm_median = optList_bm_median)
  
  print("1 panel")
  #Regime Distributions histgram
  RegimeNetdistribution_Midpoint(thisMat = thisMat, occs = occs, intervals = intervals)
  
  print("eco hists")
  
  #Origination Extinction
  ##K/S histograms
  source("~/Dropbox/ungulate_RA/RCode/Origination-Extinction/OrigExtinctScript_Draft_2018_1_9.R")
  source("~/Dropbox/ungulate_RA/RCode/Origination-Extinction/OrigExtinct_Main_2019_2_27.R")
  OrigExt_KSAnalysis(repIntSp = repIntSp,bigList = bigList, shortFam = shortFam, thisMat = thisMat)
  
  print("K/S")
  ############################################################################################################  
  #evolution analysis
  require(phytools)
  require(mvMORPH)
  require(stringr)
  require(parallel)
  
  evoSetup()
  
  #filenames.2Ma <- file.info(list.files("~/Dropbox/ungulate_RA/RCode/Results/TestDate/", pattern = "*.Rdata", full.names=TRUE))
  filenames.2Ma <- file.info(list.files("/Users/emdoughty/Google Drive/MarcLab/EvAnalysisResults20181023/", pattern = "*.Rdata", full.names=TRUE))
  #filenames.2Ma <- file.info(list.files("/Users/emdoughty/Google Drive/EvAnalysisResults20181017_dep/", pattern = "*.Rdata", full.names=TRUE))
  filenames.Trees <- file.info(list.files("~/Dropbox/ungulate_RA/RCode/NA_Ungulate_Trees/", pattern = "*.Rdata", full.names=TRUE))
  
  filenames.2Ma <- rownames(filenames.2Ma[with(filenames.2Ma, order(as.POSIXct(mtime))), ])
  filenames.Trees <- rownames(filenames.Trees[with(filenames.Trees, order(as.POSIXct(mtime))), ])
  
  Optbreaks <-vector()
  OptSigma <- list()
  optBM <- list()
  Analysis.trees <- list()
  Int.List <- list()
  
  for(ii in seq(1, length(filenames.2Ma),1)){ #rerun and save data file for sigma and separate for break.dates
    load(as.character(filenames.2Ma[[ii]]))
    #	load(as.character(filenames.Trees[[ii]]))
    #	Analysis.trees[[ii]] <- this.tree
    #	Int.List[[ii]] <- intervals
    #	load(as.character(file.list[[ii]][1,])) #to get this.tree
    #	resultsTemporalShifts
    
    bmBest <- getOptModels(opt.list = resultsTemporalShifts, model = "BM")
    optBM[[ii]] <-list(BM=bmBest[c("sigma","break.dates")], OU=NULL)
    
    #	Optbreaks <- append(Optbreaks, optBM[[ii]]$break.dates)
    
    print(ii)
  } #vector memory runs out at ~648 
  
  
  evoModelRates(evoResults = optBM, runOnRates = "BestRates")
  evoModelHists(optBM, intervals, runHistOn = "BestRates", 
                model = "BM", plotPercent = TRUE, ylim=c(0, 0.15), relativeFreq = TRUE)
  
  #evoModelHists(optBM, intervals, runHistOn = "BestRates", 
  #              model = "BM", plotPercent = FALSE)  #currently does not work
  
  plotRatesBM(this.rez=optBM, intervals = intervals, num.Trees = length(optBM), 
              PlotXmax = 56, PlotXmin = 4)
  
  RiseDeclineList <- lapply(optBM, getSigmaRateRiseDecline)
  plotBreaksRiseDecline(RiseDeclineList, intervals)
  plotBreaksRiseDecline(RiseDeclineList, intervals, plotPercent = TRUE)
  
  return()
}

plotTalkFigs_NAPC2019()