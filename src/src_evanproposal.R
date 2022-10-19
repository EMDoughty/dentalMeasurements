get.speciesMeans <- function(data, 
                             measure.colnames, 
                             meta.colnames, 
                             species.Mat = NA, 
                             fam.name=NA,
                             genus.name = NA, 
                             reg.vec = NA, 
                             bodyMass = NA, 
                             col = NA, 
                             fam.symbol = 1,
                             append.Mat <- FALSE,
                             complete = TRUE)
{
 
  if(!is.na(species.Mat)) append.mat <- TRUE
  
  uniq.sp <- unique(data$taxon)
  
  species.Mat.new <- matrix(nrow = length(uniq.sp), ncol =length(c(measure.colnames)))
  
  #combine for species average
  for(xx in seq(1, length(uniq.sp),1))
  {
    #search main doc for specific species
    phen.species <- data[data$taxon %in% uniq.sp[xx],]
    
    #aggregate emasurements
    phen.reduced <- phen.species[,colnames(phen.species) %in% measure.colnames]
    
    phen.means <- colMeans(phen.reduced[,colnames(phen.reduced) %in% measure.colnames], na.rm = TRUE)
    
    species.Mat.new[xx,] <- phen.means
  }
  
  # attributes(measure.rem)
  #species.Mat.new$family <- fam.name
  #species.Mat.new$taxon <- rownames(species.Mat.nes)
  #species.Mat.new$genus <- genus.name
  #species.Mat.new$reg.vec <- reg.vec
  #species.Mat.new$bodyMass <- bodyMass
  #species.Mat.new$col <- col
  #species.Mat.new$symbol <- fam.symbol
  
  if(append.Mat == TRUE) rbind(species.Mat, species.Mat.new)
  if(append.Mat == FALSE) species.Mat <- species.Mat.new
  
  rownames(species.Mat) <- uniq.sp
  colnames(species.Mat) <- measure.colnames
  
  if(complete == TRUE) species.Mat <- as.data.frame(species.Mat[complete.cases(species.Mat),])
  
  return(species.Mat)
}


data.coverage <- function(collected.dat = archaic.merg,   #data.coverage(archaic.merg, bigList.check, focal.archaic, "family")
                          occs.list = bigList.check, 
                          taxon.vec = focal.archaic, 
                          column.taxon = "family",
                          this.rank)
{
  #remove genus only from occs
  # occs.list <- occs.list[!occs.list$accepted_name %in% occs.list$genus,]
  
  occs.list <- occs.list[occs.list$accepted_rank %in% this.rank,]
  total.cov <- length(unique(occs.list$accepted_name[occs.list[,column.taxon] %in% taxon.vec]))
  
  total.sampled <- length(unique(collected.dat$accepted_name[collected.dat[,column.taxon] %in% taxon.vec]))
  
  print(paste("Total Coverage: ",paste(total.sampled,total.cov, sep = "/")))
  
  # print("\n")
  
  for(xx in seq(1, length(taxon.vec),1))
  {
    
    family.sam <- nrow(collected.dat[collected.dat[,column.taxon] %in% taxon.vec[xx],])
    family.tot <- length(unique(occs.list$accepted_name[occs.list[,column.taxon] %in% taxon.vec[xx]]))
    print(paste(paste(taxon.vec[xx], ": ", sep=""),paste(family.sam,family.tot, sep = "/")))
    #  print("\n")
  }
  
  return()
}

getMeasureMatCondylarths<- function(data.raw, 
                                    occs, 
                                    col.order = NULL,
                                    all.bm = TRUE,
                                    regression = NULL)
{
  archaic.ung <- data.raw
  
  meta.colnames <- c("taxon","Order","Family", "Accepted.Genus","Accepted.Species")
  measure.colnames <- c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W",
                        "p2_l","p2_w","p3_l","p3_w","p4_l","p4_w","m1_l","m1_w","m2_l","m2_w","m3_l","m3_w")
  
  #remove genus only entries
  ##indet specimens will be fine due sp. designation
  archaic.ung <- archaic.ung[!archaic.ung$Verbatim.Species %in% "",]
  
  archaic.Mat <- getSingleSpeciesMatrix_Archaic(archaic.ung)
  
  rownames(archaic.Mat) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(archaic.Mat))
  
  #get the information of data intergrity and coverage
  #scrape PBDB for all "valid" archaic ungulate names in North America
  focal.archaic <- c("Condylarthra", "Arctocyonidae", "Chriacidae", "Hyopsodontidae","Periptychidae","Phenacodontidae") 
 
  bigList.check <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & (occs$family %in% focal.archaic) | occs$order %in% focal.archaic), 
                               c("order","family", "genus", "accepted_name")])
  bigList.check <- bigList.check[order(bigList.check$order, bigList.check$family, bigList.check$genus, bigList.check$accepted_name),]
  shortFam <- sort(unique(bigList.check$family[bigList.check$family %in% focal.archaic]))	
  
  
  archaic.Mat <- as.data.frame(archaic.Mat)
  archaic.Mat$accepted_name<- rownames(archaic.Mat)
  
  archaic.merg <- merge(archaic.Mat, bigList.check, by = "accepted_name", all.x = TRUE, sort = FALSE) #not sure what this line is for
  
  
  order.col.name <- c("accepted_name", "order", "family", "genus", 
                      "P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W",
                      "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w")
  
  archaic.merg <- archaic.merg[,order.col.name]
  
  #remove rows with no measurement data
  archaic.merg <- archaic.merg[rowSums(is.na(archaic.merg[,5:ncol(archaic.merg)])) != ncol(archaic.merg)-4,]
  
  #total vs per family coverage
  
 # data.coverage(archaic.merg,bigList.check,focal.archaic,"family")
  
  colnames(archaic.merg)[colnames(archaic.merg) %in% "accepted_name"] <- "taxon"
  order.col.name[order.col.name %in% "accepted_name"] <- "taxon"
  rownames(archaic.merg) <- archaic.merg$taxon
  
  archaic.merg <- getDasmuthMeasures(measure.mat = archaic.merg) # function gets area measurements so Dasmuth 1990 regressions work

  
  #if(all.bm) estimBmAllDasmuthRegress(measure.mat, fill.missing)
  #else { 
  archaic.merg$reg.vec <- regression
  #archaic.merg$bodyMass <- 1 #placeholder so other parts code works
  archaic.merg <- approxBodyMassDasmuth(measure.mat = archaic.merg, fill.missing = FALSE)
  #}
  
  #measure.culled <- merge(measure.culled, plot.det, by = "family", all = FALSE, no.dups = TRUE)
  #archaic.merg <- archaic.merg[,c(order.col.name, "reg.vec", "bodyMass")]
  
  archaic.merg <- archaic.merg[,col.order]
  rownames(archaic.merg) <- archaic.merg$taxon
  
  #log.colnames <- measure.colnames 
  #archaic.merg[,log.colnames] <- archaic.merg[,log.colnames]#/10 #get dental measures into centimeters
  
 # archaic.temp[,"bodyMass"] <- log10(archaic.temp[,"bodyMass"])
  
  return(archaic.merg)
}

get.plot.det.Settings <- function(data.mat, 
                                  archaic.fam = NULL,
                                  archaic.col = NULL,
                                  archaic.sym = NULL)
{
  #set color assignments for modern groups in measure.mat
  fam.uniq <- unique(data.mat$family)
  fam.uniq <- fam.uniq[order(fam.uniq)]
  
  col.fam <- palette(rainbow(length(unique(fam.uniq))))
  sym.fam <- rep(c(15,16,17,18,45,43,15,16),4)
  
  #append settings for archaic ungulates
  plot.det <- data.frame(family = c(as.character(fam.uniq),archaic.fam), 
                         col = c(col.fam, archaic.col),
                         symbol = c(sym.fam, archaic.sym))
  #if this ever goes official or need an overhaul 
  ###idea is to have an input matrtix of 1 to 2 rows that detail any changes that are need...
  ###granted you wouldbe makingbasically what this outputs....
  
  return(plot.det)
}

getDentalPCA <- function(data.mat, 
                         col.rem = NULL,
                         plot.pca = FALSE, 
                         this.x = 2,
                         this.y = 3,
                         scale.load = 1)
{
  measure.pca <- data.mat[,!colnames(data.mat) %in% col.rem]
  
  #clean dataset of incomplete values
  measure.pca <- measure.pca[complete.cases(measure.pca),]
  
  #plot PC2/PC3
  dental.pca <- prcomp(measure.pca[,3:(ncol(measure.pca)-5)])
  
  if(plot.pca) {
    #plot(NA, xlim = c(-4,6.5), ylim = c(-2,2.5), ylab = "PC3", xlab = "PC2")
    plot(NA, xlim = c(-3,3), ylim = c(-1.5,1.5), ylab = "PC3", xlab = "PC2")
    points(dental.pca$x[,this.x],dental.pca$x[,this.y], pch = as.numeric(measure.pca$symbol), cex = 0.5,
           col = as.character(measure.pca$col))
    text(dental.pca$x[,this.x],dental.pca$x[,this.y],labels = rownames(dental.pca$x), 
         col = as.character(measure.pca$col), pos = 4, cex = 0.1)
    
    require(dplyr)
    
    legend.info <- measure.pca %>% distinct(family, symbol,col)
    
    legend(x  = "bottomleft", legend = legend.info$family, pch = as.numeric(legend.info$symbol), 
           col = as.character(legend.info$col), cex = 0.5)
    
    for(xx in seq(1, nrow(dental.pca$rotation),1)) {
      lines(c(0,dental.pca$rotation[xx,this.x]*scale.load), c(0,dental.pca$rotation[xx,this.y]*scale.load), col = "red")
      text(x = dental.pca$rotation[xx,this.x]*scale.load, y = dental.pca$rotation[xx, this.y]*scale.load, 
           labels = rownames(dental.pca$rotation)[xx], col = "red")
    }
  }
  return(dental.pca)
}

plot.Dental.ratio <- function(plot.data, 
                              column.x, #make sure it the name of the column
                              column.y, #make sure it the name of the column
                              xlim = NULL,
                              ylim = NULL,
                              xlab=NULL, 
                              ylab=NULL, 
                              cex.plot = 1,
                              text.on = TRUE,
                              pos = 1,
                              adj = 1,
                              cex.text = 1,
                              legend.on = FALSE,
                              legend.x = "bottomleft",
                              cex.legend = 1)
{
  original.data <- plot.data[,c(column.x,column.y,"family","symbol","col")]
  plot.data <- original.data[complete.cases(original.data),]
  
 # if(plot.sym)) plot.sym <- rep(16,nrow(plot.data))
 # if(is.na(plot.col)) plot.col <- rep("black",nrow(plot.data))
  
  if(is.null(xlab)) xlab <- column.x
  if(is.null(ylab)) ylab <- column.y
  if(is.null(xlim)) xlim <- c(0,max(plot.data[,column.x]))
  if(is.null(ylim)) ylim <- c(0,max(plot.data[,column.y]))
  
  plot(NA, xlim = xlim, ylim = ylim, ylab = ylab , xlab = xlab)
  points(plot.data[,column.x],plot.data[,column.y], 
         pch = as.numeric(plot.data$symbol), cex = cex.plot,
         col = as.character(plot.data$col))
  if(text.on) text(x = plot.data[,column.x],
                   y = plot.data[,column.y],
                   labels = rownames(plot.data), 
                   col = as.character(plot.data$col), 
                   pos = pos, adj = adj, cex = cex.text)
  
  require(dplyr)
  legend.info <- original.data %>% distinct(family, symbol,col)
  
  if(legend.on)
  {
   legend(x  = legend.x, legend = legend.info[,1], pch = as.numeric(legend.info$symbol), 
          col = as.character(legend.info$col), cex = cex.legend)
  }
  
  return() 
}

getBodyMassVectorFromMeasureMatAllMeasuresDasmuth <- function(measure.mat, linked.files=FALSE) 
{
  #######################################################################################################################################
  ##### read Dasmuth 1990 regression parameters from file, and append standard deviations
  #######################################################################################################################################
  if (linked.files) {
    # regList <- list(ruminantia=read.csv("https://dl.dropbox.com/s/dcd0bs1x5v9e7lh/regRuminantia.csv"), perissodactyla=read.csv("https://dl.dropbox.com/s/04k387q7yh4wp9u/regPerissodactyla.csv"), ungulate=read.csv("https://dl.dropbox.com/s/310ayur1s1dc8sl/regAllUngulates.csv"))
  } else {
    regList <- list(ArchaicAllSelenodonts=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicAllSelenodonts.csv"), 
                    ArchaicNonselenodonts=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicNonselenodonts.csv"), 
                    ArchaicSelenodontNonBrowsers=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicSelenodontNonBrowsers.csv"), 
                    ArchaicSelenodontBrowsers=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regArchaicSelenodontBrowsers.csv"), 
                    DasmuthAllUngulates=read.csv("~/Dropbox/Code/R/DentalMeasurements/dat/regAllArchaicUngulates.csv"))
  }
  regList <- lapply(regList, appendStDevToReg)
  
  #######################################################################################################################################
  ##### get body mass for only those taxa that have all (i.e., are not missing any) of the published measurements (about _____ species)
  #######################################################################################################################################
  measure.mat$reg[measure.mat$reg==""] <- NA
  #theseColumns <- c(as.character(regList[[1]]$m)[-which(as.character(regList[[1]]$m) %in% c("loP", "loM"))], "reg") # theseColumns is the names of columns that are also in reg - published measuremnts
  theseColumns <- as.character(regList[[1]]$m)
  theseColumns <- theseColumns[which(theseColumns %in% colnames(measure.mat))] #remove measues not present in dataset or those not
  
 # bm <- getMLBodyMasses_compiled(measure.mat[, theseColumns], regList, best.only=FALSE) # for some reason the use of these columns subset causes the is.na() check to fail.  Check works when code is ran manually. 
  
  
  #######################################################################################################################################
  ##### get regression parameters for measurements not in published regression
  #######################################################################################################################################
 # other_m <- c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M3_L", "M3_W")
  # 	other_m <- colnames(measure.mat[,sapply(measure.mat, is.numeric)])[!colnames(measure.mat[,sapply(measure.mat, is.numeric)]) %in% theseColumns]
  
  #otherReg <- lapply(X=names(regList), FUN=function(this.group) {
  #  shortMat <- measure.mat[rownames(measure.mat) %in% rownames(bm) & measure.mat$reg==this.group, other_m]
  #  short.bm <- bm[rownames(bm) %in% rownames(shortMat),]
  #  apply(log10(shortMat), MARGIN=2, FUN=function(x, bm) {
  #    lm(bm ~ x) } , bm=short.bm ) 
  #} )
  
  #names(otherReg) <- names(regList)
  
  #otherList <- lapply(otherReg, function(thisReg) t(sapply(thisReg, function(x) c(coef(x), summary(x)$sigma))))
  
  #######################################################################################################################################
  ##### merge regression parameters for unpublished measurements (otherList) with those from published (regList)
  #######################################################################################################################################
  #for (i in seq_along(regList)) {
  #  colnames(otherList[[i]]) <- c("intercept", "slope", "stdev")
  #  otherList[[i]] <- data.frame(m=rownames(otherList[[i]]), otherList[[i]], stringsAsFactors=FALSE)
  #  regList[[i]] <- merge(regList[[i]], otherList[[i]], all=TRUE, sort=FALSE)
  #}
  
  #######################################################################################################################################
  ##### recalculate body masses of all taxa with all (merged) measurements
  #######################################################################################################################################
  bm <-  getMLBodyMasses_compiled(measure.mat = measure.mat, regList = regList, best.only=FALSE)

  #convert g to kg
  bm.g <- 10^bm
  bm.g <- bm.g/1000
  bm <- log10(bm.g)
  
  bm[match(measure.mat$taxon, rownames(bm))]	
}

approxBodyMassDasmuth <- function(measure.mat, fill.missing=TRUE) 
{
  #Approximate body mass
  print("Building body mass estimates...")
  measure.mat[!is.na(measure.mat$family) & measure.mat$family=="Entelodontidae", c("P2_L", "P3_L", "p2_w", "m2_w", "m3_w")] <- NA	### entelodont tooth widths were generating >5 ton body masses, so not used for body mass estimation, here.
  #compile dentla measures to get composite measures
  measure.mat[,"bodyMass"] <- getBodyMassVectorFromMeasureMatAllMeasuresDasmuth(measure.mat = measure.mat, linked.files=FALSE) # bodyMass is in grams
  if (fill.missing) measure.mat <- fillMissingBodyMasses(measure.mat)	# this fills taxa missing their body mass with the average body mass of its cogeners
  # sort(unique(measure.mat[!sapply(measure.mat, function(x) is.character(x) | is.finite(x) | is.na(x))])) <- NA
  return(measure.mat)
}

getDasmuthMeasures <- function(measure.mat)
{
  measure.mat$P3_A <- measure.mat$P3_L*measure.mat$P3_W
  measure.mat$P4_A <- measure.mat$P4_L*measure.mat$P4_W
  measure.mat$M1_A <- measure.mat$M1_L*measure.mat$M1_W
  measure.mat$M2_A <- measure.mat$M2_L*measure.mat$M2_W
  measure.mat$M3_A <- measure.mat$M3_L*measure.mat$M3_W
  measure.mat$p3_a <- measure.mat$p3_l*measure.mat$p3_w
  measure.mat$p4_a <- measure.mat$p4_l*measure.mat$p4_w
  measure.mat$m1_a <- measure.mat$m1_l*measure.mat$m1_w
  measure.mat$m2_a <- measure.mat$m2_l*measure.mat$m2_w
  measure.mat$m3_a <- measure.mat$m3_l*measure.mat$m3_w
  
  measure.mat$M1_3_A <- measure.mat$m1_3_a <- measure.mat$M1_2_A <- measure.mat$m1_2_a <- measure.mat$UpM <- measure.mat$loM <- measure.mat$P4_M3_A <- 
  measure.mat$p4_m3_a <- measure.mat$P4_M3_L <- measure.mat$p4_m3_l <- measure.mat$P4_M2_A <- measure.mat$p4_m2_a  <-
  measure.mat$P4_M2_L <- measure.mat$p4_m2_l <- measure.mat$P3_M3_A <- measure.mat$p3_m3_a <- NA
    
  for(ii in seq_len(nrow(measure.mat))) measure.mat[ii,] <- calcCompositeDasmuthMeasures(measure.mat = measure.mat[ii,])
  
  return(measure.mat)
}
  
calcCompositeDasmuthMeasures <- function(measure.mat)
{
  if(!is.na(measure.mat$M1_A) & !is.na(measure.mat$M2_A) & !is.na(measure.mat$M3_A)) measure.mat$M1_3_A <- measure.mat$M1_A + measure.mat$M2_A + measure.mat$M3_A
  if(!is.na(measure.mat$m1_a) & !is.na(measure.mat$m2_a) & !is.na(measure.mat$m3_a)) measure.mat$m1_3_a <- measure.mat$m1_a + measure.mat$m2_a + measure.mat$m3_a
  if(!is.na(measure.mat$M1_A) & !is.na(measure.mat$M2_A)) measure.mat$M1_2_A <- measure.mat$M1_A + measure.mat$M2_A
  if(!is.na(measure.mat$m1_a) & !is.na(measure.mat$m2_a)) measure.mat$m1_2_a <- measure.mat$m1_a + measure.mat$m2_a
  if(!is.na(measure.mat$M1_L) & !is.na(measure.mat$M2_L) & !is.na(measure.mat$M3_L)) measure.mat$UpM <- measure.mat$M1_L + measure.mat$M2_L + measure.mat$M3_L
  if(!is.na(measure.mat$m1_l) & !is.na(measure.mat$m2_l) & !is.na(measure.mat$m3_l)) measure.mat$loM <- measure.mat$m1_l + measure.mat$m2_l + measure.mat$m3_l
  if(!is.na(measure.mat$P4_A) & !is.na(measure.mat$M1_A) & !is.na(measure.mat$M2_A) & !is.na(measure.mat$M3_A)) measure.mat$P4_M3_A <- measure.mat$P4_A + measure.mat$M1_A + measure.mat$M2_A + measure.mat$M3_A
  if(!is.na(measure.mat$p4_a) & !is.na(measure.mat$m1_a) & !is.na(measure.mat$m2_a) & !is.na(measure.mat$m3_a)) measure.mat$p4_m3_a <- measure.mat$p4_a + measure.mat$m1_a + measure.mat$m2_a + measure.mat$m3_a
  if(!is.na(measure.mat$P4_L) & !is.na(measure.mat$M1_L) & !is.na(measure.mat$M2_L) & !is.na(measure.mat$M3_L)) measure.mat$P4_M3_L <- measure.mat$P4_L + measure.mat$M1_L + measure.mat$M2_L + measure.mat$M3_L
  if(!is.na(measure.mat$p4_l) & !is.na(measure.mat$m1_l) & !is.na(measure.mat$m2_l) & !is.na(measure.mat$m3_l)) measure.mat$p4_m3_l <- measure.mat$p4_l + measure.mat$m1_l + measure.mat$m2_l + measure.mat$m3_l
  if(!is.na(measure.mat$P4_A) & !is.na(measure.mat$M1_A) & !is.na(measure.mat$M2_A)) measure.mat$P4_M2_A <- measure.mat$P4_A + measure.mat$M1_A + measure.mat$M2_A
  if(!is.na(measure.mat$P4_L) & !is.na(measure.mat$M1_L) & !is.na(measure.mat$M2_L)) measure.mat$P4_M2_L <- measure.mat$P4_L + measure.mat$M1_L + measure.mat$M2_L
  if(!is.na(measure.mat$p4_a) & !is.na(measure.mat$m1_a) & !is.na(measure.mat$m2_a)) measure.mat$p4_m2_a  <- measure.mat$p4_a + measure.mat$m1_a + measure.mat$m2_a
  if(!is.na(measure.mat$p4_l) & !is.na(measure.mat$m1_l) & !is.na(measure.mat$m2_l)) measure.mat$p4_m2_l <- measure.mat$p4_l + measure.mat$m1_l + measure.mat$m2_l
  if(!is.na(measure.mat$P3_A) & !is.na(measure.mat$P4_A) & !is.na(measure.mat$M1_A) & !is.na(measure.mat$M2_A) & !is.na(measure.mat$M3_A)) measure.mat$P3_M3_A <- measure.mat$P3_A + measure.mat$P4_A + measure.mat$M1_A + measure.mat$M2_A + measure.mat$M3_A
  if(!is.na(measure.mat$p3_a) & !is.na(measure.mat$p4_a) & !is.na(measure.mat$m1_a) & !is.na(measure.mat$m2_a) & !is.na(measure.mat$m3_a)) measure.mat$p3_m3_a <- measure.mat$p3_a + measure.mat$p4_a + measure.mat$m1_a + measure.mat$m2_a + measure.mat$m3_a
  
  return(measure.mat)
}

getMeasureMatWithBodyMassesDasmuth <- function() {
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

estimBmAllDasmuthRegress <- function(measure.mat, fill.missing = TRUE) #function is not functional
{
  DasRegress <- c("ArchaicAllSelenodonts", "ArchaicNonselenodonts", "ArchaicSelenodontNonBrowsers", "ArchaicSelenodontBrowsers","DasmuthAllUngulates")
  
  measure.mat$reg.vec <- DasRegress[1]
  
  colname.bm <- paste("bm", DasRegress, sep=".")
  
  measure.mat <- approxBodyMassDasmuth(measure.mat = measure.mat, filling.missing = fill.missing)
  
  for(ii in seq(2,length(DasRegress)))
  {
    measure.mat$reg.vec <- DasRegress[ii]
    
    temp.mat <- approxBodyMassDasmuth(measure.mat = measure.mat, fill.missing = fill.missing)
    measure.mat$DasRegress[ii]
    
  }
  
}

#check that no duplicate species names in rows


getDateTaxa <- function(measure.mat, occs, this.rank = "species")
{
  thisRanges <- getTaxonRangesFromOccs(occs = occs, rank = this.rank, random=FALSE)
  rownames(thisRanges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisRanges))
  measure.mat[,c("FO","LO")] <- thisRanges[match(measure.mat$taxon, rownames(thisRanges)),]
  return(measure.mat)
}

shoulderPlot <- function(measure.mat, plot.y, intervals, occs, bigList, shortFam, repIntTaxa = NULL, quants = NULL,
                         optList_bm_median = NULL, specOcc.col = "gray0", specOcc.alpha = 0.5,
                         plot.breaks = FALSE, manual.breaks = NULL, break.text = TRUE, break.col = "firebrick4",
                         this.rank = "species",
                         main = NULL, ylab = NULL, ylim = NULL, yaxp = NULL, xlab = NULL, xlim = NULL, xaxp = c(55,5,10), cex.axis = 1.5, cex.lab = 1.5, 
                         do.subepochs=TRUE, overlay.labels = FALSE, overlay.color=TRUE, thisAlpha.intervals=0.33, thisAlpha.text = 0.33, borderCol="white", invertTime=FALSE, scale.cex=0.75, scale.headers = 0.95, text.offset = 0.025,
                         do.quants = FALSE, col.axis = "black", col.lab = "black", poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
{

  #for later versions make plot.x so user can define x axis but for now keep as $FO
  
  #check that measure.mat has FO and LO
  if(!"FO" %in% colnames(measure.mat) || !"LO" %in% colnames(measure.mat)) measure.mat <- getDateTaxa(measure.mat = measure.mat, occs = occs, this.rank = this.rank)
  
  #par(mar=c(4,4,2.5,0.5))
  # quartz(width=12, height=6)
  if(is.null(xlim)) xlim <- c(max(intervals), min(intervals))
  
  plot(measure.mat$FO, measure.mat[,plot.y], type="n", xlim= xlim, ylim=ylim, yaxp = yaxp, ylab = ylab, xaxp = xaxp, xlab = xlab, main = main, cex.axis = cex.axis, cex.lab = cex.lab, col.axis = col.axis, col.lab = col.lab)
  # plot(measure.mat$FO, measure.mat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(50,0,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, col="gray75", fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #alter xaxpto change x-axis values
  # rect(-10e6, -10e6, 10e6, 10e6, col="white")
  overlayCzTimescale(do.subepochs= do.subepochs, color = overlay.color, thisAlpha.text = thisAlpha.text, thisAlpha.intervals = thisAlpha.intervals, borderCol = borderCol, invertTime = invertTime, scale.cex = scale.cex, scale.headers = scale.headers, text.offset = text.offset)
  
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
    if (is.finite(measure.mat$FO[i]) & is.finite(measure.mat$LO[i]) & measure.mat$FO[i] != measure.mat$LO[i]) lines(x=measure.mat[i,c("FO","LO")], y=c(measure.mat[,plot.y][i], measure.mat[,plot.y][i]), lwd=0.75, pch=21, col=alphaColor(specOcc.col, specOcc.alpha)) #alphaColor(orderColors[i], 0.5)
  }
  points(measure.mat[complete.cases(measure.mat[ ,c("FO","LO")]) & measure.mat$FO == measure.mat$LO, c("FO", plot.y)], pch=21, col=alphaColor(specOcc.col, specOcc.alpha), cex=0.25) #this line is not generating the proper output for the final graph due to c("FO","bodyMass") causing a  "undefined columns selected" error
  
  if(plot.breaks)
  {
    if(is.null(manual.breaks))
    {
      # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
      # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
      abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col = break.col)
      if(break.text)
      {
        text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col=break.col)
        text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0),cex=0.5, col=break.col)
      }
    }
    
    if(!is.null(manual.breaks))
    {
      abline(v= manual.breaks, lwd=1.5, col = break.col)
      if(break.text)
      {
        text(x= manual.breaks - 0.35, y=par()$usr[3], labels= manual.breaks, pos=3, cex=0.5, col = break.col)
        text(x= manual.breaks - 0.35, y=par()$usr[3], labels= paste(manual.breaks, "Ma"), adj=c(0,0),cex=0.5, col = break.col)
      }
    }
  }
  
  if(do.quants)
  {
    #make function to run Handley method and run the quants
    if(is.null(quants)) quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,plot.y], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
    polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor(poly.col, 0.25), border=poly.col)
    polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor(poly.col, 0.25), border=poly.col)
    lines(rowMeans(intervals), quants[3,], col=alphaColor(median.col[1], 0.5), lwd=5)
    lines(rowMeans(intervals), quants[3,], col=alphaColor(median.col[2], 1.0), lwd=3)
    points(rowMeans(intervals), quants[3,], col=alphaColor(median.col[3], 0.5), cex=0.5)
    box(lwd=1)
  }
}



#work in progress
getRepIntData <- function(measure.mat, 
                          col.nam = NULL,
                          occs,
                          intervals = NULL,
                          int_length = 2,
                          int_min = 1,
                          int_max = 56,
                          do.parallel = FALSE, 
                          reps = 10,
                          do.subsample = FALSE, 
                          sub.type = "global", #c("global", "interval)
                          quota = NULL,
                          intOccs2 = NULL,
                          do.disparity = FALSE,
                          bootstrapSpecimens = FALSE, #require specimenMat which may be from now depricated versions this code
                          bootstrapSpecies = FALSE,
                          bootstrapSpeciesWithinIntervals = FALSE,
                          do.heuristic = TRUE,
                          extra.intvs = 0,
                          this.rank = "species",
                          do.rangethrough = TRUE,
                          save.reps = FALSE,
                          plotHist = FALSE)
{
  if(is.null(intervals))
  {
    intervals = makeIntervals(int_min, int_max, int_length)
    intList <- listifyMatrixByRow(intervals)
  }
  
  if (do.parallel) require(parallel)
  
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
      measure.mat[,col.nam] <- makeBodyMasses(measure.mat, regList, best.only=TRUE)
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
      
      if(sub.type %in% "global") intOccs <- subsampleReps.Global(intOccs = intOccs, nTaxa = 0,quota = quota)
      
      if(sub.type %in% "interval")
      {
        if(!length(intOccs) == length(intOccs2))
        {
          cat("Objects are of different lengths. Subsample not completed \n")
        } else {
           intOccs <- subsampleReps.WithinInterval(intOccs = intOccs, intOccs2 = intOccs2)
        }
      }
     # nOccs <- sapply(intOccs, length)
      # nTaxa <- sapply(intSp, length)			### if you want to set the quota no lower than the maximum number of SIB taxa; intSp is required for this to work, so has to be done above
      #nTaxa <- 0									### set to zero to simply set the quota to the minimum number of occurrences
      #quota <- max(c(max(nTaxa), min(nOccs)))		### quota is either the maximum number of observed taxa, or the minimum number of occurrences
      #cat("Subsampling quota set to", quota, "occurrences")
      
      #intOccs <- lapply(X=intOccs, FUN=sample, size=quota)
    }
    
    repIntOccs[[rep]] <- intOccs 
  }
  
  repIntTaxa <- getRepIntTaxaFromRepIntOccs(repIntOccs, this.rank=this.rank, do.rangethrough=do.rangethrough)
  
  print("Completed getting taxa with intervals")
  
  if(save.reps)
  {
    if(Sys.info()["sysname"] == "Darwin"){
      save(repIntTaxa, repIntOccs, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/repIntTaxa_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      #load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    } else if(Sys.info()["sysname"] == "Windows"){
      save(repIntTaxa, repIntOccs, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/repIntTaxa_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      # load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    }
  }
  
  repIntAll <- list(RepsTaxa = repIntTaxa, RepOccs = repIntOccs)
  
  return(repIntAll)
}


taxHadley <- function(repIntTaxa, 
                      occs,
                      shortFam,
                      bigList,
                      intervals,
                      reps = 10,
                      extra.intvs = 0, 
                      do.parallel=FALSE, 
                      do.heuristic = FALSE,
                      do.save = FALSE,
                      run.update = 100) #runs needed to list how many runs completed to date
{
  
####################################################################################################################################
  ### Handley analysis of taxonomic distributions
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
    if(this.rep %% run.update == 0) cat("Taxonomic Handley Rep:", this.rep, "\n")
  }
  
  ####################################################################################################################################
  if(do.save)
  {
    if(Sys.info()["sysname"] == "Darwin"){
      save(optList_tax_median, optList_tax_allReps, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      #load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    } else if(Sys.info()["sysname"] == "Windows"){
      save(optList_tax_median, optList_tax_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/Taxon_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      # load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    }
  }
  
  repIntTaxaAll <- list(OptListTaxMedian = optList_tax_median, OptListTaxAllReps = optList_tax_allReps, taxCube = taxCube)
 
  return(repIntTaxaAll)
}

traitHadley <- function(repIntTaxa,
                        trait.col,
                        measure.mat,
                        occs,
                        shortFam,
                        bigList, 
                        intervals, 
                        reps = 10,
                        bmBreaks = c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf), #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
                        extra.intvs = 0, 
                        do.parallel=FALSE, 
                        do.heuristic = FALSE,
                        do.save = FALSE,
                        run.update = 100)
{
  ####################################################################################################################################
  ### Handley analysis of body mass distributions
  print("Beginning median body mass Handley analysis...")
  
  countCube <- sapply(repIntTaxa, function(this.rep) {
    sapply(this.rep, function(this.intv, this.rep) {
      hist(measure.mat[,trait.col][match(this.intv, measure.mat$taxon)], breaks=bmBreaks, plot=FALSE)$counts
    }, this.rep=this.rep)
  }, simplify = "array")
  
  countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 
  
  optList_bm_median <- doHandleyTest(countBox[2,,], n=nrow(measure.mat), do.heuristic=do.heuristic, extra.intvs=extra.intvs)
  
  print("Beginning body mass Handley analysis for all reps...")
  optList_bm_allReps <- list()
  for (this.rep in seq_len(reps)) {
    this.n <- length(unique(unlist(sapply(repIntTaxa [[this.rep]], function(x) x))))
    optList_bm_allReps[[this.rep]] <- doHandleyTest(countCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)
    if(this.rep %% run.update == 0) cat("Body Mass Handley Rep:",this.rep, "\n")
  }
  
  if(do.save)
  {
    ####################################################################################################################################
    if(Sys.info()["sysname"] == "Darwin"){
      save(optList_bm_median, optList_bm_allReps, file=paste0("~/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      #load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    } else if(Sys.info()["sysname"] == "Windows"){
      save(optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/BM_handleyResult_SampleStandardized=", do.subsample, timestamp(),".Rdata"))
      # load('~/Dropbox/ungulate_RA/EcologyResults/allUngulates/handleyResult##------ Thu Nov  9 02:12:20 2017 ------##_allUngulates.Rdata')
    }
    
    # save(repIntTaxa, repIntOccs, optList_tax_median, optList_tax_allReps, optList_bm_median, optList_bm_allReps, file=paste0("C:/Users/Blaire/Dropbox/ungulate_RA/EcologyResults/handleyResult", timestamp(),".Rdata"))
    ####################################################################################################################################
  }
  
  repIntTraitAll <- list(OptListTraitMedian = optList_bm_median, OptListTraitAllReps = optList_bm_allReps, countCube = countCube, countBox = countBox)
  
  return(repIntTraitAll)
}



hadleyHisto <- function(measure.mat, bmBreaks, repIntTaxa, optList_bm_median, intervals)
{
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
  return()
}

plotNumberShiftPerRep <- function(hist1.dat = NULL, hist1.breaks = seq(-0.5, 10, 1.0), hist1.col = "orchid4", hist1.main = "Number of taxonomic shifts in each rep", hist1.xlab = "Number of Shifts", hist1.xlim = c(0,10),
                                  plothist2 = FALSE, hist2.dat = NULL, hist2.breaks = seq(-0.5, 11.5, 1.0), hist2.col = "firebrick4", hist2.main = "Number of body mass shifts in each rep", hist2.xlab = "Number of Shifts", hist2.xlim = c(0,10))
{
  # Number of shifts per rep
  hist(sapply(hist1.dat, function(x) length(x) - 2), breaks = hist1.breaks, col=hist1.col, main=hist1.main , xlab = hist1.xlab, xlim = hist1.xlim)
  if(plothist2) hist(sapply(hist2.dat, function(x) length(x) - 2), breaks = hist2.breaks, col=hist2.col, main=hist2.main , xlab = hist2.xlab, xlim = hist2.xlim)
} 
  
getTraitQuants <- function(measure.mat, 
                           traitCol = "bodyMass", 
                           repIntTaxa)
{
  quants <- apply(sapply(repIntTaxa, function(y) sapply(y, function(x) quantile(measure.mat[x,traitCol], probs=c(0, 0.25, 0.5, 0.75, 1.0), na.rm=TRUE)), simplify = "array"), c(1,2), median, na.rm=TRUE)
  return(quants)
}  

plotHistShiftinInterval <- function(optList,
                                    intervals,
                                    reps = 10,
                                    xaxp = c(65,5,10),
                                    main = "Number of Replicates with a Taxonomic Distribution Shift",
                                    xlab = "Time (Ma)",
                                    xlim = c(55,0),
                                    cex = 0.3,
                                    col = "orchid4",
                                    border = "orchid1",
                                    do.subepochs = TRUE,
                                    labels = TRUE,
                                    freq=TRUE)
{
####################################################################################################################################
  ### number of Replicates with a shift in that interval
  # optList_tax_allReps <- optList_tax_allReps_heuristic
  # optList_tax_allReps <- optList_tax_allReps_full
  breakHist <- hist(rowMeans(intervals)[unlist(sapply(optList, function(x) x[[length(x) - 1]]$optBreaks))], breaks=sort(unique(unlist(intervals))), plot=FALSE)
  plot(breakHist, col=NA, border=NA, labels=FALSE, freq=TRUE, cex=cex, 
       xlim=rev(range(intervals)), xaxp =xaxp, ylim=c(0,reps), main=main, xlab=xlab)
  overlayCzTimescale(do.subepochs=do.subepochs)
  plot(breakHist, col=col, border=border, freq=freq, cex=cex, 
       xaxp = xaxp,  xlim=xlim, ylim=c(0,reps), add=TRUE)
 if(labels == TRUE) text(x = breakHist$mids, breakHist$counts, labels = breakHist$counts, cex = cex, pos = 3)
}

plotZachosCurve <- function(intervals)
{
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
}

plotTaxonThroughTime <- function(repIntTaxa, bigList, shortFam, taxCube = NULL, 
                                 intervals,optList_tax_median, plot.breaks = FALSE, 
                                 manual.breaks = NULL, xlim = NULL,
                                 legend = TRUE)
{
  # taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% measure.mat$taxon[x]], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  #taxCubeG <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$genus) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  if(is.null(taxCube))
  {
   taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
  }
   dimnames(taxCube) <- list(shortFam, rownames(intervals), NULL)
  # prop <- t(apply(taxCube, c(1,2), median, na.rm=TRUE))
  prop <- t(apply(taxCube, c(1,2), mean, na.rm=TRUE))
  colnames(prop)[colnames(prop)==""] <- "indeterminate"
  # dimnames(prop) <- list(rownames(intervals), shortFam)
  source("~/Dropbox/Code/R/common_src/taxonomicEv.R")
 
  if(is.na(xlim)) xlim <- c(max(intervals, na.rm=TRUE),min(intervals, na.rm=TRUE))
  plotStackedRichness(this.box=prop, intervals=intervals, do.log=FALSE, 
                      overlay.labels=TRUE, numbers.only=TRUE, 
                      legend=legend ,xlim=xlim)
  #med.n <- median(length(unique(unlist(sapply(repIntTaxa[[this.rep]], function(x) measure.mat$taxon[x]))))) #what is this.rep set to during this function?  variable is used in for loop in Handley
  # med.n <- median(sapply(repIntTaxa, function(x) length(unique(unlist(sapply(x, function(y) measure.mat$taxon[y]))))))
  # optList_tax <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
  if(plot.breaks)
  {
    if(is.null(manual.breaks))
    {
      abline(v=sort(c(intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2], range(intervals))), lwd=1.5, col="darkorchid4")
      text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels=rev(seq_len(length(optList_tax_median[[length(optList_tax_median)-1]]$optBreaks) + 1)), pos=3, cex=0.5, col="darkorchid4")
      text(x= sort((c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3], labels= paste(sort(c(max(intervals), intervals[optList_tax_median[[length(optList_tax_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0), cex=0.5, col="darkorchid4")
    }
  
    if(!is.null(manual.breaks))
    {
      abline(v= manual.breaks, lwd=1.5, col="firebrick4")
      text(x= manual.breaks - 0.35, y=par()$usr[3], labels= manual.breaks, pos=3, cex=0.5, col="firebrick4")
      text(x= manual.breaks - 0.35, y=par()$usr[3], labels= paste(manual.breaks, "Ma"), adj=c(0,0),cex=0.5, col="firebrick4")
    }
  }
  box(lwd=1)
}

plot3Panel <- function(intervals,repIntTaxa, bigList, shortFam, measure.mat, plot.y, occs,
                       taxCube = NULL,  optList_tax_median, optList_bm_median = NULL,
                       plot.breaks = FALSE, this.rank = "species",
                       ylab = NULL, xlab = NULL, xaxp = c(65,5,10), cex.axis = 1.5, cex.lab = 1.5, 
                       do.subepochs = TRUE,  do.quants = FALSE)
{
  quartz(width=6.89)
  par(mfrow=c(4,1), mar=c(0,4,0.5,0.5), mgp=c(2, 1,0))
  plotZachosCurve(intervals = intervals)
  plotTaxonThroughTime(repIntTaxa = repIntTaxa, bigList = bigList, shortFam = shortFam, 
                       taxCube = taxCube, intervals = intervals,optList_tax_median = optList_tax_median, 
                       plot.breaks = plot.breaks)
  shoulderPlot(measure.mat = measure.mat, plot.y = plot.y, intervals = intervals, occs = occs, 
               bigList = bigList, shortFam = shortFam, repIntTaxa, quants = quants,
               optList_bm_median = optList_bm_median, plot.breaks = plot.breaks, this.rank = this.rank,
               ylab = ylab, xlab = xlab, xaxp = xaxp, cex.axis = cex.axis, cex.lab = cex.lab, 
               do.subepochs = do.subepochs,  do.quants = do.quants)
}

firstDifferencesBinned <- function(data.input) #need to test
{
  fd <- vector()
  fd.name <- vector()
  for(ii in seq(2, length(data.input)-1,1)) #intervals and quants will be in different orders?
  {
    fd[ii-1] <- data.input[ii] - data.input[ii-1]
    fd.name[ii-1] <- paste(names(data.input[ii-1]),names(data.input[ii]),sep="to")
  }
  
  names(fd) <- fd.name
  
  return(fd)
}

plotRegimeDistribution <- function(occs) #not functional
{
  ###Compile regime distribution histograms and net/change histograms using the range of taxon occurance
  #		regimeHist(repIntTaxa = repIntTaxa, breaks = c(51,47,37,21, 5,2), optList=optList_bm_median, measure.mat = measure.mat, netFreq = TRUE, regimeFreq=FALSE,
  #							 netPlotType = "pos/neg", plot.together = FALSE)
  
  ###Compile regime distribution histograms and net/change histograms using the midpoint of taxon occurance
  ##this was used to produce the regime and net change histograms for publication
  #get midpoint of each taxon 
  require(stringr)
  taxonranges <- getTaxonRangesFromOccs(occs = occs, random = FALSE)
  rownames(taxonranges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(taxonranges))
  taxonranges <- cbind(taxonranges, (taxonranges[,"FO"]+taxonranges[,"LO"])/2); colnames(taxonranges)[3] <- "MP"
  occDatesMP <- taxonranges[,"MP"]
  intSpMP <- apply(intervals, 1, function(thisIntv) taxonranges[taxonranges[,"MP"] > thisIntv[1] & taxonranges[,"MP"] <= thisIntv[2],])
  countCubeMP <- sapply(intSpMP, function(x) hist(measure.mat$bodyMass[match(rownames(x), measure.mat$taxon)], breaks=c(-Inf, 0, 0.845098, 1.322219, 2, Inf), plot=FALSE)$counts, simplify = "array")
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
}


#regimeHist is meant to parse through the repIntTaxa list object to print the distribution of traits
#within the regimes demarcated by breaks or designated by optList (if breaks are not provided).
#netFreq and regime Freq sets whether to use proportions or actually counts for the net gain/loss and regime 
#distribution histograms, respectively.
regimeHist_v2 <- function(repIntTaxa = NULL, 
                          trait.Col = "bodyMass", 
                          breaks = NULL,
                          trait.breaks = c(-2, 0, 0.845098, 1.322219, 2, 4), 
                          optList, 
                          thisMat, 
                          netFreq = TRUE, 
                          regimeFreq = FALSE, 
                          netPlotType = "absolute", 
                          plot.together = FALSE,
                          plot.quartz = TRUE,
                          plot.histo = c("all", "withinregime", "btwnregime")) 
{
  #function to generate histograms of species within each regime
  #get list of intervals that comprise a regime
  repIntTaxa_regimes <- repIntTaxa
 
  regimeBM.list <- list()
  regimeSp.List <- list()
  
  for (ii in seq(1, length(repIntTaxa_regimes), 1)) {
    names(repIntTaxa_regimes[[ii]]) <- str_remove(names(repIntTaxa_regimes[[ii]]), " Ma")
  }
  
  regimeSp <- unique(unlist(repIntTaxa_regimes[[1]][which(as.double(names(repIntTaxa_regimes[[ii]])) > breaks[1])]))
  regimeSp <- regimeSp[regimeSp %in% thisMat$taxon]
  regimeBM <- thisMat[thisMat$taxon %in% regimeSp, trait.Col]
  
  regimeSp.List[[1]] <- regimeSp
  regimeBM.list[[1]] <- regimeBM
  
  rm(regimeSp, regimeBM)
  
  for (mm in seq(2, length(breaks), 1)) {
    #get regimes for remaining sections
    
    #seq(min(which(as.double(names(repIntTaxa_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntTaxa_regimes[[ii]])) > breaks[mm])),
    #max(which(as.double(names(repIntTaxa_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntTaxa_regimes[[ii]])) > breaks[mm])),1))
    regimeSp <- unique(unlist(repIntTaxa_regimes[[1]][which(as.double(names(repIntTaxa_regimes[[ii]])) < breaks[mm - 1] & as.double(names(repIntTaxa_regimes[[ii]])) > 
                                                              breaks[mm])]))
    regimeSp <- regimeSp[regimeSp %in% thisMat$taxon]
    regimeBM <- thisMat[thisMat$taxon %in% regimeSp, trait.Col]
    
    regimeSp.List[[mm]] <- regimeSp
    regimeBM.list[[mm]] <- regimeBM
    
    rm(regimeSp, regimeBM)
  }
  
  regimeSp <- unique(unlist(repIntTaxa_regimes[[1]][which(as.double(names(repIntTaxa_regimes[[ii]])) < breaks[length(breaks)])]))
  regimeSp <- regimeSp[regimeSp %in% thisMat$taxon]
  regimeBM <- thisMat[thisMat$taxon %in% regimeSp, trait.Col]
  
  regimeSp.List[[length(breaks) + 1]] <- regimeSp
  regimeBM.list[[length(breaks) + 1]] <- regimeBM
  
  rm(regimeSp, regimeBM)
  
  regimeNames <- vector()
  regimeNames[1] <- paste(">", max(breaks, sep = ""))
  for (ii in seq(2, length(regimeBM.list) - 1, 1)) {
    regimeNames[ii] <- paste(paste(breaks[ii - 1], " to ", sep = ""), breaks[ii], sep = "")
  }
  regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list) - 1], ">", sep = "")
  names(regimeBM.list) <- regimeNames
  bmBreaks <- trait.breaks
  breakCol <- rainbow(length(bmBreaks))
  
  plot.histo.withinregime<- FALSE
  if(plot.histo == "all" || plot.histo == "withinregime") {plot.histo.withinregime <- TRUE}
  
    hist.list <- list()
    if(plot.quartz) quartz()
    par(mfrow = c(2, length(regimeBM.list)/2 + 1))
    for (jj in seq(1, length(regimeBM.list), 1)) 
    {
      hist.list[[jj]] <- hist(regimeBM.list[[jj]], main = names(regimeBM.list[jj]), 
                              breaks = bmBreaks, las = 1, col = breakCol, freq = regimeFreq, 
                              xlab = "log Body Mass (kg)", ylim = c(0, 1), plot = plot.histo.withinregime)
   
    }
  
  
  #these are here to set up a list to contain how regimes change between breaks
  regimeNetChange <- hist.list
  regimeNetChange[[length(regimeNetChange)]] <- NULL
  for(ii in seq(1,length(regimeNetChange),1))
  {
    regimeNetChange[[ii]]$xname <- breaks[ii]
  }

  
  for (tt in seq(2, length(hist.list), 1)) {
    regimeNetChange[[tt - 1]]$counts <- hist.list[[tt]]$counts - hist.list[[tt - 1]]$counts
    tt <- tt + 1
  }
  
  if(plot.histo == "all" || plot.histo == "btwnregime")
  {
    if (netPlotType == "pos/neg") {
      #net change between regimes
      if(plot.quartz) quartz()
      par(mfrow = c(2, length(hist.list)/2))
      for (hh in seq(1, length(regimeNetChange), 1)) {
        plot(regimeNetChange[[hh]], main = regimeNetChange[[hh]]$xname, freq = netFreq, las = 1, col = breakCol, xlab = "log Body Mass (kg)") #####for some reason is identical to regular hist
      }
    }
    if (netPlotType == "absolute") {
      #absolute net change between regimes
      quartz()
      par(mfrow = c(2, length(hist.list)/2))
      for (dd in seq(1, length(regimeNetChange), 1)) regimeNetChange[[dd]]$counts <- abs(regimeNetChange[[dd]]$counts)
      for (hh in seq(1, length(regimeNetChange), 1)) {
        plot(regimeNetChange[[hh]], main = names(regimeNetChange[hh]), freq = netFreq, las = 1, col = breakCol, xlab = "log Body Mass (kg)") #####for some reason is identical to regular hist
      }
    }
  }
  ####Get a hsitogram for each entry of repIntTaxa to avoid missing unique species in each regime
  ####and then take median across all of the histograms.	
}

subsampleReps.Global <- function(intOccs, 
                          nTaxa = 0, ### set to zero to simply set the quota to the minimum number of occurrences
                          quota = NULL)
{
  nOccs <- sapply(intOccs, length)
  # nTaxa <- sapply(intSp, length)			### if you want to set the quota no lower than the maximum number of SIB taxa; intSp is required for this to work, so has to be done above
  
  if(is.null(quota)) quota <- max(c(max(nTaxa), min(nOccs, na.rm = TRUE), na.rm = TRUE)) ### quota is either the maximum number of observed taxa, or the minimum number of occurrences
  cat("Subsampling quota set to", quota, "occurrences\n")
  
  intOccs <- lapply(X=intOccs, FUN=sample, size=quota)
  
  return(intOccs)
}

subsampleReps.WithinInterval <- function(intOccs, 
                                         intOccs2)
{
  nOccs <- sapply(intOccs, length)
  #nOccs2 <- sapply(intOccs2, length)
 
 
  for(ii in seq(1, length(intOccs), 1))
  {
    quota <- max(c(0, intOccs2[[ii]], na.rm = TRUE)) ### quota is either the maximum number of observed taxa, or the minimum number of occurrences
    cat("Subsampling quota set to", quota, "occurrences\n")
    
    intOccs[[ii]] <- sample(intOccs[[ii]],size = quota)
  }
  
  return(intOccs)
}

#how wnat do it up to me...doesnt...conceptually sample ungualte at carnivoran level, more prey than preds, 
#artificallly crush curve to flat line
#way now most agnostic
#sample all mammals then cull taxons we want
#do our subsample for ungulates and carnivores trea equally

#don't have different quota for each interval
##when per interval not standardized per interval, not standard variation over time, i dont know what, think on what standardize,, what keep evem

#interesting sample standardized=carnivores more sensitive, all occurences 30k mammals 60% ungulates then 20% canrivorans , 20% rodents, then %> everything else
#standardize all mammal set by ungualtes most likely

#makr doingn=need standardizing=lotta variation sampling carnivores

#estim smapling prob over time=pyrates go to....=estim some sampling/preservation bias vs all mammals=see how change over time

#subsampling shift thing be a neogene story

getNALMAIntervals <- function(startDate = 86, endDate = 0)
{
  #Janis 1998 volume designations
  
  NALMA <- rev(c("Aquilian","Judithian", "Lancian", "Puercan", "Torrejonian", "Tiffanian", "Clarkforkian", "Wasatchian", "Bridgerian", 
                 "Uintan", "Duchesnean", "Chadronian", "Orellan","Whitneyan", "Arikareean", 
                 "Hemingfordian", "Barstovian", "Clarendonian", "Hemphillian", "Blancan", 
                 "Irvingtonian", "Rancholabrean"))
  
  interval.start <- rev(c(86, 84, 70, 66.043, 63.3, 61.7, 56.8, 55.8, 50.3, 46.2, 40.4, 37.2, 33.9, 33.3, 30.8,
                      20.43, 15.97, 13.6, 10.3, 4.9, 1.8, 0.3))
  interval.end <- rev(c(84, 70, 66.043, 63.3, 61.7, 56.8, 55.8, 50.3, 46.2, 40.4, 37.2, 33.9,33.3, 30.8, 20.43, 15.97,
                    13.6, 10.3, 4.9, 1.8, 0.3, 0.012))
  nalma.intervals <- as.data.frame(cbind(agetop = interval.end, ageBase = interval.start))
  
  rownames(nalma.intervals) <- NALMA
  
  #Remove entries not within bounds of dates set
  if(is.null(startDate)) 
    
  nalma.intervals <- nalma.intervals[nalma.intervals$ageBase <= startDate & nalma.intervals$agetop >= endDate,]
  
  #need way to handle dates falling within intervals
  ####rbind these onto one another
  #nalma.intervals[nalma.intervals$agetop <= startDate & nalma.intervals$ageBase >= endDate,]
  
  return(nalma.intervals)
}

getBigList <- function(focal.order = NULL, focal.family.manual = NULL)
{
  
  focal.family <- unique(occs[occs$order %in% focal.order,]$family)
  focal.family <- c(as.character(focal.family),focal.family.manual) #"Arctocyonidae", "Hyopsodontidae","Periptychidae","Phenacodontidae",
  #"Hapalodectidae","Mesonychidae","Prodinoceratidae","Uintatheriidae","Triisodontidae"
  focal.family <- focal.family[!focal.family %in% ""]
  focal.family <- focal.family[order(focal.family)]
  
  bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% focal.family & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
  bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
  # bigList[order(bigList$family, bigList$accepted_name),]

  
  bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)
  
  #bigList <- bigList[bigList$family %in% shortFam & bigList$order %in% focal.order,]
  
  return(bigList)
}

getDiversity1stDiff <- function(data.mat,intervals)
{
  div.diff <- matrix(nrow = nrow(data.mat), ncol = ncol(data.mat)-1)
  for(xx in seq(1, nrow(data.mat),1))
  {
    div.diff[xx,] <- diff(data.mat[xx,], lag = 1)
  }
  
  NALMA_rownames <- vector()
  for(xx in seq(1, nrow(intervals),1))
  {
    NALMA_rownames[xx] <- paste(rownames(intervals)[xx], 
                                rownames(intervals)[xx + 1], sep="-")
  }
  colnames(div.diff) <- NALMA_rownames[1:length(NALMA_rownames)-1]
  rownames(div.diff) <- rownames(data.mat)
  
  return(div.diff)
}
###############################################

jon_archicBodyMassEstim <- function()
{

source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)

dat.condy <- read.csv(file="~/Dropbox/code/R/dentalMeasurements/dat/ArchaicUngulate_UploadFile_2021_4_29.csv", strip.white=TRUE)
this.reg <- read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regArchaicNonselenodonts.csv", strip.white=TRUE)
this.reg <- cbind(this.reg, stdev=(log10(100 + this.reg["see"]) - 2))
colnames(this.reg)[ncol(this.reg)] <- "stdev"
# head(dat.condy)
# head(measure.mat)
# sapply(dat.condy, class)

colnames(dat.condy)[colnames(dat.condy)=="Family"] <- "family"
colnames(dat.condy)[colnames(dat.condy)=="Accepted.Genus"] <- "genus"

dat.condy <- dat.condy[!dat.condy$Accepted.Species %in% c("", "sp.") ,]				#removes indeterminate species
dat.condy$species <- apply(dat.condy[,c("genus", "Accepted.Species")], 1, function(x) paste(x[1], x[2], sep="_"))
dat.condy <- dat.condy[apply(dat.condy[, names(dat.condy) %in% this.reg$m], 1, function(x) any(!is.na(x))),] #retain only rows with data

# head(dat.condy[, names(dat.condy) %in% names(measure.mat)])

dat.condy$P4_M2_L <- dat.condy$P4_L + dat.condy$M1_L + dat.condy$M2_L
dat.condy$P4_M3_L <- dat.condy$P4_L + dat.condy$M1_L + dat.condy$M2_L + dat.condy$M3_L
dat.condy$UpM <- dat.condy$M1_L + dat.condy$M2_L + dat.condy$M3_L
dat.condy$p4_m2_l <- dat.condy$p4_l + dat.condy$m1_l + dat.condy$m2_l
dat.condy$p4_m3_l <- dat.condy$p4_l + dat.condy$m1_l + dat.condy$m2_l + dat.condy$m3_l
dat.condy$loM <- dat.condy$m1_l + dat.condy$m2_l + dat.condy$m3_l

dat.condy$p3_a <- dat.condy$p3_l * dat.condy$p3_w
dat.condy$p4_a <- dat.condy$p4_l * dat.condy$p4_w
dat.condy$m1_a <- dat.condy$m1_l * dat.condy$m1_w
dat.condy$m2_a <- dat.condy$m2_l * dat.condy$m2_w
dat.condy$m3_a <- dat.condy$m3_l * dat.condy$m3_w

dat.condy$P3_A <- dat.condy$P3_L * dat.condy$P3_W
dat.condy$P4_A <- dat.condy$P4_L * dat.condy$P4_W
dat.condy$M1_A <- dat.condy$M1_L * dat.condy$M1_W
dat.condy$M2_A <- dat.condy$M2_L * dat.condy$M2_W
dat.condy$M3_A <- dat.condy$M3_L * dat.condy$M3_W

dat.condy$M1_2_A <- dat.condy$M1_A + dat.condy$M2_A
dat.condy$M1_3_A <- dat.condy$M1_A + dat.condy$M2_A + dat.condy$M3_A
dat.condy$P4_M2_A <- dat.condy$P4_A + dat.condy$M1_A + dat.condy$M2_A
dat.condy$P4_M3_A <- dat.condy$P4_A + dat.condy$M1_A + dat.condy$M2_A + dat.condy$M3_A
dat.condy$P3_M3_A <- dat.condy$P3_A + dat.condy$P4_A + dat.condy$M1_A + dat.condy$M2_A + dat.condy$M3_A

dat.condy$m1_2_a <- dat.condy$m1_a + dat.condy$m2_a
dat.condy$m1_3_a <- dat.condy$m1_a + dat.condy$m2_a + dat.condy$m3_a
dat.condy$p4_m2_a <- dat.condy$p4_a + dat.condy$m1_a + dat.condy$m2_a
dat.condy$p4_m3_a <- dat.condy$p4_a + dat.condy$m1_a + dat.condy$m2_a + dat.condy$m3_a
dat.condy$p3_m3_a <- dat.condy$p3_a + dat.condy$p4_a + dat.condy$m1_a + dat.condy$m2_a + dat.condy$m3_a

dat.condy <- aggregate(dat.condy[,names(dat.condy) %in% this.reg$m], by = list(taxon = dat.condy$species), median, na.rm = TRUE)
rownames(dat.condy) <- dat.condy$taxon

dat.condy <- dat.condy[, names(dat.condy) %in% this.reg$m]
# sapply(dat.condy, function(x) any(!is.na(x)))

bm.vec <- array(NA, dim=c(nrow(dat.condy), 1), dimnames=list(rownames(dat.condy), "bodyMass"))
for (i in seq_len(nrow(dat.condy))) {
  bm.vec[i] <- getMLbodyMassForOneSpecies(this=dat.condy[i,], thisReg=this.reg)
} 
10^(bm.vec-3)

dat.condy$bodyMass <- bm.vec

return(dat.condy)
}

jon_archicBodyMassEstim_edited <- function()
{
  
  source('~/Dropbox/code/R/dentalMeasurements/src/src_bodyMassEstimation.R', chdir = TRUE)
  
  dat.condy <- read.csv(file="~/Dropbox/code/R/dentalMeasurements/dat/ArchaicUngulate_UploadFile_2021_4_29.csv", strip.white=TRUE)
  this.reg <- read.csv("~/Dropbox/code/R/dentalMeasurements/dat/regArchaicNonselenodonts.csv", strip.white=TRUE)
  this.reg <- cbind(this.reg, stdev=(log10(100 + this.reg["see"]) - 2))
  colnames(this.reg)[ncol(this.reg)] <- "stdev"
  # head(dat.condy)
  # head(measure.mat)
  # sapply(dat.condy, class)
  
  colnames(dat.condy)[colnames(dat.condy)=="Family"] <- "family"
  colnames(dat.condy)[colnames(dat.condy)=="Accepted.Genus"] <- "genus"
  
  dat.condy <- dat.condy[!dat.condy$Accepted.Species %in% c("", "sp.") ,]				#removes indeterminate species
  dat.condy$species <- apply(dat.condy[,c("genus", "Accepted.Species")], 1, function(x) paste(x[1], x[2], sep="_"))
  dat.condy <- dat.condy[apply(dat.condy[, names(dat.condy) %in% this.reg$m], 1, function(x) any(!is.na(x))),] #retain only rows with data
  
  # head(dat.condy[, names(dat.condy) %in% names(measure.mat)])
  
  dat.condy$P4_M2_L <- dat.condy$P4_L + dat.condy$M1_L + dat.condy$M2_L
  dat.condy$P4_M3_L <- dat.condy$P4_L + dat.condy$M1_L + dat.condy$M2_L + dat.condy$M3_L
  dat.condy$UpM <- dat.condy$M1_L + dat.condy$M2_L + dat.condy$M3_L
  dat.condy$p4_m2_l <- dat.condy$p4_l + dat.condy$m1_l + dat.condy$m2_l
  dat.condy$p4_m3_l <- dat.condy$p4_l + dat.condy$m1_l + dat.condy$m2_l + dat.condy$m3_l
  dat.condy$loM <- dat.condy$m1_l + dat.condy$m2_l + dat.condy$m3_l
  
  dat.condy$p3_a <- dat.condy$p3_l * dat.condy$p3_w
  dat.condy$p4_a <- dat.condy$p4_l * dat.condy$p4_w
  dat.condy$m1_a <- dat.condy$m1_l * dat.condy$m1_w
  dat.condy$m2_a <- dat.condy$m2_l * dat.condy$m2_w
  dat.condy$m3_a <- dat.condy$m3_l * dat.condy$m3_w
  
  dat.condy$P3_A <- dat.condy$P3_L * dat.condy$P3_W
  dat.condy$P4_A <- dat.condy$P4_L * dat.condy$P4_W
  dat.condy$M1_A <- dat.condy$M1_L * dat.condy$M1_W
  dat.condy$M2_A <- dat.condy$M2_L * dat.condy$M2_W
  dat.condy$M3_A <- dat.condy$M3_L * dat.condy$M3_W
  
  measure.mat$M1_3_A <- measure.mat$m1_3_a <- measure.mat$M1_2_A <- measure.mat$m1_2_a <- measure.mat$UpM <- measure.mat$loM <- measure.mat$P4_M3_A <- 
    measure.mat$p4_m3_a <- measure.mat$P4_M3_L <- measure.mat$p4_m3_l <- measure.mat$P4_M2_A <- measure.mat$p4_m2_a  <-
    measure.mat$P4_M2_L <- measure.mat$p4_m2_l <- measure.mat$P3_M3_A <- measure.mat$p3_m3_a <- NA
  
  dat.condy$M1_2_A <- dat.condy$M1_A + dat.condy$M2_A
  dat.condy$M1_3_A <- dat.condy$M1_A + dat.condy$M2_A + dat.condy$M3_A
  dat.condy$P4_M2_A <- dat.condy$P4_A + dat.condy$M1_A + dat.condy$M2_A
  dat.condy$P4_M3_A <- dat.condy$P4_A + dat.condy$M1_A + dat.condy$M2_A + dat.condy$M3_A
  dat.condy$P3_M3_A <- dat.condy$P3_A + dat.condy$P4_A + dat.condy$M1_A + dat.condy$M2_A + dat.condy$M3_A
  
  dat.condy$m1_2_a <- dat.condy$m1_a + dat.condy$m2_a
  dat.condy$m1_3_a <- dat.condy$m1_a + dat.condy$m2_a + dat.condy$m3_a
  dat.condy$p4_m2_a <- dat.condy$p4_a + dat.condy$m1_a + dat.condy$m2_a
  dat.condy$p4_m3_a <- dat.condy$p4_a + dat.condy$m1_a + dat.condy$m2_a + dat.condy$m3_a
  dat.condy$p3_m3_a <- dat.condy$p3_a + dat.condy$p4_a + dat.condy$m1_a + dat.condy$m2_a + dat.condy$m3_a
  
  dat.condy <- aggregate(dat.condy[,names(dat.condy) %in% this.reg$m], by = list(taxon = dat.condy$species), median, na.rm = TRUE)
  rownames(dat.condy) <- dat.condy$taxon
  
  dat.condy <- dat.condy[, names(dat.condy) %in% this.reg$m]
  # sapply(dat.condy, function(x) any(!is.na(x)))
  
  bm.vec <- array(NA, dim=c(nrow(dat.condy), 1), dimnames=list(rownames(dat.condy), "bodyMass"))
  for (i in seq_len(nrow(dat.condy))) {
    bm.vec[i] <- getMLbodyMassForOneSpecies(this=dat.condy[i,], thisReg=this.reg)
  } 
  10^(bm.vec-3)
  
  dat.condy$bodyMass <- bm.vec
  
  return(dat.condy)
}

getCurrentHigherTaxonomy <- function(archaic.ung, save.file=NULL) { #this function is meant to update the .csv files (or dataframe) of dental measures and update the taxonomy so that the accepted genus and species and higher taxonomy is correct (as much as possible).  Made for use with my Conylarth dental measurement sheet.
  require(stringr)
  #remove question marks and quotes form verbatim names temporarily
  
  Verbatim.Genus <- gsub("\"", "", archaic.ung$Verbatim.Genus)
  Verbatim.Species<- gsub("\"", "", archaic.ung$Verbatim.Species)
  
  Verbatim.Genus <- gsub("[?]", "", Verbatim.Genus)
  Verbatim.Species<- gsub("[?]", "", Verbatim.Species)
  
  Verbatim.Genus <- gsub("cf. ", "", Verbatim.Genus)
  Verbatim.Species<- gsub("cf. ", "", Verbatim.Species)
  
  Verbatim.Genus <- gsub("aff. ", "", Verbatim.Genus)
  Verbatim.Species<- gsub("aff. ", "", Verbatim.Species)
  
  archaic.ung$Verbatim.taxon <- paste(Verbatim.Genus, Verbatim.Species, sep=" ")
 # archaic.ung$Verbatim.taxon <- paste(archaic.ung$Verbatim.Genus, archaic.ung$Verbatim.Species, sep=" ")
  
  Accepted.Name <- getCurrentTaxa(tax.vec = archaic.ung$Verbatim.taxon) #this does not work on most Genus sp. designations or entries that are in quotes.  Is hit or miss on entries with question marks
  
  Accepted.Genus <- str_split_fixed(string = Accepted.Name, " ", n = Inf)[,1]
  Accepted.Species <- paste0(str_split_fixed(string = Accepted.Name, " ", n = Inf)[,2],str_split_fixed(string = Accepted.Name, " ", n = Inf)[,3])
  
  archaic.ung$Accepted.Name <- Accepted.Name
  archaic.ung$Accepted.Name <- gsub(pattern = "[[:space:]]", replacement = "_", x = archaic.ung$Accepted.Name)
  
  archaic.ung$Accepted.Genus <- Accepted.Genus
  archaic.ung$Accepted.Species <- Accepted.Species
  
  uniqTax <- lapply(c("Mammalia"), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  uniqTax <- rbind(uniqTax[[1]])
  uniqTax$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$accepted_name)
  
  archaic.ung$Family <- uniqTax$family[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)]
  archaic.ung$Family[is.na(archaic.ung$Family)] <- uniqTax$family[match(x=archaic.ung$Accepted.Genus[is.na(archaic.ung$Family)], table=uniqTax$taxon_name)]

  archaic.ung$Order <- uniqTax$order[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)]
  archaic.ung$Order[is.na(archaic.ung$Order)] <- uniqTax$order[match(x=archaic.ung$Accepted.Genus[is.na(archaic.ung$Order)], table=uniqTax$taxon_name)]
  
  
  if(!is.null(save.file)) write.csv(archaic.ung, file = save.file)
     
  return(archaic.ung)
}

getTaxaInClade <- function(clades, occs, save.file=NULL) {
  uniqTax <- lapply(c(clades), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  if(length(uniqTax) == 1) { uniqTax <- rbind(uniqTax[[1]])
    }
  if(length(uniqTax) == 2) { uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]])
  }
  if(length(uniqTax) > 2) { 
    taxaMat <- rbind(uniqTax[[1]], uniqTax[[2]])
    for(xx in seq(3, length(uniqTax),1))
    {
      taxaMat <- rbind(taxaMat, uniqTax[[xx]])
    }
    uniqTax <- taxaMat
  }
  
  uniqTax$accepted_species <- str_split_fixed(string = uniqTax$accepted_name, " ", n = Inf)[,2]
  
  uniqTax$verbatim_genus <- str_split_fixed(string = uniqTax$taxon_name, " ", n = Inf)[,1]
  uniqTax$verbatim_species <- str_split_fixed(string = uniqTax$taxon_name, " ", n = Inf)[,2]
  
  uniqTax <- uniqTax[,c("phylum", "class", "order", "family", "genus", "accepted_species", "verbatim_genus", "verbatim_species","accepted_name", "taxon_name")]
  
  if(!is.null(save.file)) write.csv(uniqTax, file = save.file)
  
  return(uniqTax)
}

getTaxonomyForOneBaseTaxon_AcceptedName <- function(this.taxon) {
  this.URL <- URLencode(paste0(server, "taxa/list.csv?base_name=", this.taxon, "&taxon_status=all&show=class"))
  this.names <- getStuffWithoutWarnings(this.URL)
  this.names[,c("phylum", "class", "order", "family", "genus", "accepted_name","taxon_name")]
}

getTargetTaxa<- function(measure.mat, uniqTax, occs, save.file = NULL){ # get a list of taxa that indicates #occs it has and whether its in North America if it does
  
  uniqTax$NoOccs <- 0 
  
  uniqTax$InNorAmer <- FALSE
  
  uniqTax$DataCollected <- 0
  
  measure.mat$accepted_name <- paste(measure.mat$Accepted.Genus, measure.mat$Accepted.Species, sep = " ")
  
  dental.col <-  c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W",
  "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w")
  
  for(xx in uniqTax$accepted_name)
  {
    uniqTax$NoOccs[uniqTax$accepted_name %in% xx] <- nrow(occs[occs$accepted_name %in% xx,])
    

    if(!uniqTax$NoOccs[uniqTax$accepted_name %in% xx] <= 0) 
    {
      if(any(occs$cc[occs$accepted_name %in% xx] %in% c("US", "CA", "MX"))) uniqTax$InNorAmer[uniqTax$accepted_name %in% xx] <- TRUE
    }
    
    #remove blank rows so they are not skwing the target list.  This will remove entries that are present due to images only or compilation measures of the tooth row.  Likely make sthe Accepted name portion of line below redundant
    measure.mat <- measure.mat[!apply(measure.mat[,dental.col] == "", 1, all),]
    #get number of entries present in measure.mat
    uniqTax$DataCollected[uniqTax$accepted_name %in% xx] <- nrow(measure.mat[measure.mat$accepted_name %in% xx & !measure.mat$Catalog.Number %in% "Accepted Names PBDB",])
  }
 
  #likely need to use the higher taxonomy columns to get actual counts for occurences in families and genera.  the for loop above only grabs those entered into the database at the family or generic level.

  if(!is.null(save.file)) write.csv(uniqTax, file = save.file)
}

data.coverage_v2 <- function(clades, clade.level = "family", data.mat, occs, measure.colnames, save.file = NULL)
{
  output.list <- list(AllTaxa = NA, OccsOnly = NA, MissingTaxa = NA) #, NoOccs = NA)
  
  #need this function to 
  ##1) read in a given taxon and then collate the genera and species within
  uniqTax <- lapply(c(clades), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  if(length(uniqTax) == 1) { uniqTax <- rbind(uniqTax[[1]])
  }
  if(length(uniqTax) == 2) { uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]])
  }
  if(length(uniqTax) > 2) { 
    taxaMat <- rbind(uniqTax[[1]], uniqTax[[2]])
    for(xx in seq(3, length(uniqTax),1))
    {
      taxaMat <- rbind(taxaMat, uniqTax[[xx]])
    }
    uniqTax <- taxaMat
  }
  
  tax.mat <- uniqTax[uniqTax$taxon_name %in% clades, c("order","family")]
  
  tax.mat$Percent.Species <- tax.mat$No.Species <- tax.mat$Sampled.Species <- tax.mat$Percent.Genera <- tax.mat$No.Genera <- tax.mat$Sampled.Genera  <-  NA
  
  for(xx in tax.mat$family)
  {
    uniqTax.temp <- uniqTax[uniqTax$family %in% xx,]
    tax.mat$No.Genera[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$genus[uniqTax.temp$family %in% xx])) - 1 #-1 is used to remove the entry for blank genera (usually from higher level taxa in database)
                                                                
    tax.mat$No.Species[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$genus != "" & uniqTax.temp$accepted_name != uniqTax.temp$genus]))
  }
  
  ##2) compare that list versus a list for the taxa that I have measurements for individual teeth (excludes image only or specimens with tooth row measurements but not individual teeth)
  #### cleanup the data.mat to remove extra columns
  data.mat <- data.mat[, colnames(data.mat)[1:66]]
  data.mat <- data.mat[!data.mat$Catalog.Number %in% "Accepted Names PBDB",]
  
  data.mat$accepted_name <- paste(data.mat$Accepted.Genus, data.mat$Accepted.Species, sep = " ")
  
  ####remove rows that lack measurements
  data.mat <- data.mat[!apply(is.na(data.mat[,measure.colnames]) | data.mat[,measure.colnames] == "", 1, all),]
 
  sample.genera <- unique(c(data.mat$Accepted.Genus, 
                          data.mat$Accepted.Genus[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus]))
  
  sample.species <- unique(data.mat$accepted_name[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus])
 
  for(xx in tax.mat$family)
  {
    uniqTax.temp <- uniqTax[uniqTax$family %in% xx,]
    tax.mat$Sampled.Genera[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$genus[uniqTax.temp$genus %in% sample.genera]))-1
    
    tax.mat$Sampled.Species[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$accepted_name %in% sample.species]))
  }

  tax.mat$Percent.Genera <- (tax.mat$Sampled.Genera/tax.mat$No.Genera)*100
  tax.mat$Percent.Species <- (tax.mat$Sampled.Species/tax.mat$No.Species)*100
  
  output.list$AllTaxa <- tax.mat
  
  ###########################################################################################################################################################################
  ##3) Generate a table for species with occs in database
  tax.mat$Percent.Species <- tax.mat$No.Species <- tax.mat$Sampled.Species <- tax.mat$Percent.Genera <- tax.mat$No.Genera <- tax.mat$Sampled.Genera  <-  NA
  
  for(xx in tax.mat$family)
  {
    uniqTax.temp <- uniqTax[uniqTax$family %in% xx & (uniqTax$accepted_name %in% occs$accepted_name | uniqTax$accepted_name %in% occs$genus),]
    tax.mat$No.Genera[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$genus[uniqTax.temp$family %in% xx])) #-1 is used to remove the entry for blank genera (usually from higher level taxa in database)
    
    tax.mat$No.Species[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$genus != "" & uniqTax.temp$accepted_name != uniqTax.temp$genus]))
  }
   
  sample.genera <- unique(data.mat$Accepted.Genus[data.mat$accepted_name %in% occs$accepted_name])
  
  sample.species <- unique(data.mat$accepted_name[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus])
  sample.species <- sample.species[sample.species %in% occs$accepted_name]
  
  for(xx in tax.mat$family)
  {
    uniqTax.temp <- uniqTax[uniqTax$family %in% xx & (uniqTax$accepted_name %in% occs$accepted_name | uniqTax$accepted_name %in% occs$genus),]
    tax.mat$Sampled.Genera[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$genus[uniqTax.temp$genus %in% sample.genera]))
    
    tax.mat$Sampled.Species[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$accepted_name %in% sample.species]))
  }
  
  tax.mat$Percent.Genera <- (tax.mat$Sampled.Genera/tax.mat$No.Genera)*100
  tax.mat$Percent.Species <- (tax.mat$Sampled.Species/tax.mat$No.Species)*100
  
  output.list$OccsOnly <- tax.mat
  
  #get taxa that are missing
  
  UnSampled.Taxa<- list()
  
  for(xx in tax.mat$family)
  {
    uniqTax.temp <- uniqTax[uniqTax$family %in% xx & (uniqTax$accepted_name %in% occs$accepted_name) & !(uniqTax$accepted_name %in% data.mat$accepted_name | uniqTax$accepted_name %in% data.mat$Accepted.Genus),]
    UnSampled.Genera <- uniqTax.temp[!uniqTax.temp$genus %in% unique(sample.genera),]
    
    UnSampled.Species <- uniqTax.temp[!uniqTax.temp$accepted_name %in% unique(sample.species),]
    
    UnSampled.Taxa[[xx]] <- unique(rbind(UnSampled.Genera, UnSampled.Species))
  }
  
  if(length(UnSampled.Taxa) == 1) { UnSampled.Taxa <- rbind(UnSampled.Taxa[[1]])
  }
  if(length(UnSampled.Taxa) == 2) { UnSampled.Taxa <- rbind(UnSampled.Taxa[[1]], UnSampled.Taxa[[2]])
  }
  if(length(UnSampled.Taxa) > 2) { 
    taxaMat <- rbind(UnSampled.Taxa[[1]], UnSampled.Taxa[[2]])
    for(xx in seq(3, length(UnSampled.Taxa),1))
    {
      taxaMat <- rbind(taxaMat, UnSampled.Taxa[[xx]])
    }
    UnSampled.Taxa <- taxaMat
  }
  
  if(!is.null(save.file)) write.csv(UnSampled.Taxa, file = save.file)
  
  output.list$MissingTaxa <- UnSampled.Taxa
  
  ###########################################################################################################################################################################
  ##4) Generate a table for species without occs in database
 # tax.mat$Percent.Species <- tax.mat$No.Species <- tax.mat$Sampled.Species <- tax.mat$Percent.Genera <- tax.mat$No.Genera <- tax.mat$Sampled.Genera  <-  NA
  
#  for(xx in tax.mat$family)
#  {
#    prev.uniqTax.temp <- uniqTax[uniqTax$family %in% xx & (uniqTax$accepted_name %in% occs$accepted_name | uniqTax$accepted_name %in% occs$genus),]
#    uniqTax.temp <- uniqTax[uniqTax$family %in% xx & !(uniqTax$accepted_name %in% occs$accepted_name | uniqTax$accepted_name %in% occs$genus) & !(uniqTax$genus %in% prev.uniqTax.temp$genus),]
    
#    tax.mat$No.Genera[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$genus[uniqTax.temp$family %in% xx])) #-1 is used to remove the entry for blank genera (usually from higher level taxa in database)
    
#    tax.mat$No.Species[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$genus != "" & uniqTax.temp$accepted_name != uniqTax.temp$genus]))
#  }
  
#  sample.genera <- unique(data.mat$Accepted.Genus[!data.mat$accepted_name %in% occs$accepted_name])
  
#  sample.species <- unique(data.mat$accepted_name[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus])
#  sample.species <- sample.species[!sample.species %in% occs$accepted_name]
  
#  for(xx in tax.mat$family)
#  {
#    uniqTax.temp <- uniqTax[uniqTax$family %in% xx,]
#    tax.mat$Sampled.Genera[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$genus[uniqTax.temp$genus %in% sample.genera]))
    
#    tax.mat$Sampled.Species[tax.mat$family %in% xx] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$accepted_name %in% sample.species]))
 # }
  
#  tax.mat$Percent.Genera <- (tax.mat$Sampled.Genera/tax.mat$No.Genera)*100
#  tax.mat$Percent.Species <- (tax.mat$Sampled.Species/tax.mat$No.Species)*100
  
#  output.list$NoOccs <- tax.mat
  
  return(output.list)
}

getSingleSpeciesMatrix_Archaic <- function(specimen.mat = NULL) {
  #compile and label dental measurments for specimens
 
  specimen.mat$species <- paste(specimen.mat$Accepted.Genus, specimen.mat$Accepted.Species, sep= " ")
  
  #need to drop duplicate entries from the dataset
  specimen.mat <- specimen.mat[!specimen.mat$Include %in% FALSE,]
  
  specimen.mat[sapply(specimen.mat, is.nan)] <- NA
  specimen.mat$species <- getCurrentTaxa(tax.vec = specimen.mat$species)
  specimen.mat$species <- gsub(pattern = "[[:space:]]", replacement = "_", x = specimen.mat$species)
  
  upLabels <- c("P2_L", "P2_TrigW", "P2_TalW", "P2_W", 
                "P3_L", "P3_TrigW", "P3_TalW", "P3_W", 
                "P4_L","P4_TrigW", "P4_TalW", "P4_W", 
                "M1_L", "M1_TrigW", "M1_TalW","M1_W", 
                "M2_L", "M2_TrigW", "M2_TalW", "M2_W", 
                "M3_L", "M3_TrigW", "M3_TalW", "M3_W") #\P2_L\",\"P2_W\","
  loLabels <- casefold(upLabels)
  
  #remove symbols from measurement section (?, *, etc.)
  for(xx in c(upLabels, loLabels))
  {
    specimen.mat[,xx] <- as.numeric(gsub("\\*", "",  as.character(specimen.mat[,xx])))
    specimen.mat[,xx] <- as.numeric(gsub("*", "",  as.character(specimen.mat[,xx])))
    specimen.mat[,xx] <- as.numeric(str_remove(as.character(specimen.mat[,xx]),"\\*"))
    specimen.mat[,xx] <- as.numeric(gsub("\\?", "",  as.character(specimen.mat[,xx])))
  }
  
  ############################
  #### note that specimens are aggregated by their medians, so as to minimize the effect of outlier measurements
  ############################
  
  # measure.mat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(taxon = specimen.mat$species), mean, na.rm = TRUE)
  measure.mat <- aggregate(specimen.mat[, c(upLabels, loLabels)], by = list(taxon = specimen.mat$species), median, na.rm = TRUE)
  
  #measure.mat[, sapply(measure.mat, is.numeric)] <- measure.mat[, sapply(measure.mat, is.numeric)]#/10 #converts mm measurements to cm for compatibility with Janis regressions
  #measure.mat <- transform(measure.mat, p4_a = p4_l * p4_w, m1_a = m1_l * m1_w, m2_a = m2_l * m2_w, m3_a = m3_l * m3_w, M2_A = M2_L * M2_W)
  measure.mat[sapply(measure.mat, is.nan)] <- NA
  
  rownames(measure.mat) <- measure.mat$taxon
  
  return(measure.mat)
}

sensitivity.numberReps <- function(countCube_herb = NULL, countCube_pred = NULL, number.reps = 1000)
{
  
  #look at how # reps alters
  ##1) the number of species in each category per bin
  ##2) the variance of occs in bin when doing subsampling
  ##3) how the correlation coef. of the median assemblage changes with reps
  
  
  countBox <- list()
  prop <- list()
  
  rep.test <- seq(number.reps,dim(countCube_herb)[3], number.reps)
  
  for(xx in seq(1, length(rep.test), 1))
  {
    countBox_herb <- apply(countCube_herb[,,c(1,rep.test[xx])], c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 
    countBox_pred <- apply(countCube_pred[,,c(1,rep.test[xx])], c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE) 
    
    singleBox <- list(herb=countBox_herb, pred=countBox_pred)
    
    countBox[[xx]] <- singleBox
    
    ##############################################################
    prop_herb <- t(apply(countCube_herb[,,c(1,rep.test[xx])], c(1,2), median, na.rm=TRUE))
    colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"
    
    prop_pred <- t(apply(countCube_pred[,,c(1,rep.test[xx])], c(1,2), median, na.rm=TRUE))
    colnames(prop_pred)[colnames(prop_pred)==""] <- "indeterminate"
    
    singleProp <- list(herb=prop_herb, pred=prop_pred)
    
    prop[[xx]] <- singleProp
  }
  
  ### raw
  par(mfrow = c(1,2))
  ungulates <- prop_herb
  predators <- prop_pred
  
  corr.results.both <- cor(predators, ungulates, method = "spearman")
 
 #1st Dif
  pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred), intervals = intervals)
  UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb), intervals = intervals)
  
  ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
  predators <- t(as.matrix(pred.prey.DivFirstDiff))

  corr.results.both <- cor(predators, ungulates, method = "spearman")
  
  
  
   
  return()
}

sensitivity.temporalBinning <- function() {
  
  return()
}

sensitivity.bodyMassBinning <- function() {
  
  return()
}

estiamteMissingDiversity_Alroy <- function() {
  
  return()
}



