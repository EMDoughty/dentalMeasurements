require(ape)

sourceCall <- function(){
  
  # require(phytools)
  # require(mvMORPH)
  # require(stringr)
  # require(parallel)
  
  #sources for Jon Marcot's code and specimen measurements 
  # source("https://dl.dropbox.com/s/8jy9de5owxj72p7/strat.R")
  # source("https://dl.dropbox.com/s/253p4avcvb66795/occFns.R")
  # source("https://dl.dropbox.com/s/9gdafsqss2b586x/phy_dateTree.R")
  # source("https://dl.dropbox.com/s/9tdawj35qf502jj/amandaSrc.R")
  # source("https://dl.dropbox.com/s/rlof7juwr2q4y77/blasto_Birlenbach.R")
  # source("https://dl.dropbox.com/s/643op7ye4s49w8p/utils_marcot.R")
  # source("https://dl.dropbox.com/s/pxnbmroe2hgcgo3/kozak_src.R")
  if(Sys.info()["sysname"] == "Darwin"){
    ################MAC
    source('~/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call cource file for analysis functions
    source('~/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R') #call cource file for data functions
    source('~/Dropbox/ungulate_RA/RCode/EvAnalysesTreeSrc.R') #call cource file for tree functions
    source('~/Dropbox/ungulate_RA/RCode/EvAnalysesPlotSrc.R') #call cource file for tree functions
    
  } else if(Sys.info()["sysname"] == "Windows"){
    #################PC   #had to swap my name with Blaires in pathnames
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/2017_2_22_CopesRule_Source_Func_Clavel_ver1.R') #call source for use on Evan's PC
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesDataSrc.R') #call cource file for data functions
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesTreeSrc.R') #call cource file for tree functions
    source('C:/Users/Van Valkenburgh Lab/Dropbox/ungulate_RA/RCode/EvAnalysesPlotSrc.R') #call cource file for tree functions
    
  } else print("Mac or Windows operating systems are not detected")
}


getDat.vec <- function(tree.list, measureMat) {
  measureMatCrop <- matrix()
  dat.vec <- data.frame()
  if (reps == 1) { measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[1]], thisMat = measureMat)
  tree.list[[1]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[1]], thisMat = measureMat)
  } else { 
    for (ii in seq(1, length(tree.list))) { 
      measureMatCrop <- comp_TreeSpecMat_outMatrix(this.tree = tree.list[[ii]], thisMat = measureMat) 
      tree.list[[ii]] <- comp_TreeSpecMat_outTree(this.tree = tree.list[[ii]], thisMat = measureMat)
    }
  }
  dat.vec <- as.data.frame(data.matrix(frame = measureMatCrop[,"bodyMass"]))
  dimnames(dat.vec) <- list(rownames(measureMatCrop), "body.mass")
  names(dat.vec) <- "body.mass"
  
  print("Dat vector completed")
  return(dat.vec)
}

getOptModelValues <- function(fileList = NULL, ExtractAttributes = c("sigma","break.dates"))
{
  optBM <- list()
  for(ii in seq(1, length(filenames.2Ma),1)){ #rerun and save data file for sigma and separate for break.dates
    load(as.character(filenames.2Ma[[ii]]))
    #	load(as.character(filenames.Trees[[ii]]))

    bmBest <- getOptModels(opt.list = resultsTemporalShifts, model = "BM")
    optBM[[ii]] <-list(BM=bmBest[ExtractAttributes], OU=NULL)
    
    if(ii %% 10) print(ii)
  } 
  
  return(optBM)
}