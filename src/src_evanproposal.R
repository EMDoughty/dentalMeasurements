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
                             append.Mat = FALSE,
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

getDiversity1stDiff <- function(data.mat, output.rownames = NULL)
{
  div.diff <- matrix(nrow = nrow(data.mat), ncol = ncol(data.mat)-1)
  for(xx in seq(1, nrow(data.mat),1))
  {
    div.diff[xx,] <- diff(data.mat[xx,], lag = 1)
  }
  
  dat.rownames <- vector()
  for(xx in seq(1, length(output.rownames),1))
  {
   dat.rownames[xx] <- paste(output.rownames[xx], 
                             output.rownames[xx + 1], sep="-")
  }
  colnames(div.diff) <- dat.rownames[1:length(dat.rownames)-1]
  rownames(div.diff) <- rownames(data.mat)
  
  return(div.diff)
}
###############################################
getCurrentHigherTaxonomy <- function(archaic.ung, save.file=NULL) { #this function is meant to update the .csv files (or dataframe) of dental measures and update the taxonomy so that the accepted genus and species and higher taxonomy is correct (as much as possible).  Made for use with my Conylarth dental measurement sheet.
  require(stringr)
  #remove question marks and quotes form verbatim names temporarily
  
  Verbatim.Genus <- gsub("\"", "", archaic.ung$Verbatim.Genus)
  Verbatim.Species <- gsub("\"", "", archaic.ung$Verbatim.Species)
  
  Verbatim.Genus <- gsub("[?]", "", Verbatim.Genus)
  Verbatim.Species <- gsub("[?]", "", Verbatim.Species)
  
  Verbatim.Genus <- gsub("cf. ", "", Verbatim.Genus)
  Verbatim.Species <- gsub("cf. ", "", Verbatim.Species)
  
  Verbatim.Genus <- gsub("aff. ", "", Verbatim.Genus)
  Verbatim.Species <- gsub("aff. ", "", Verbatim.Species)
  
  archaic.ung$Verbatim.taxon <- NA
  
  archaic.ung$Verbatim.taxon[!Verbatim.Species %in% ""] <- paste(Verbatim.Genus[!Verbatim.Species %in% ""], Verbatim.Species[!Verbatim.Species %in% ""], sep=" ") #7/13/2023 genera only entries have added space = conflict with PBDB
  archaic.ung$Verbatim.taxon[Verbatim.Species %in% "sp." | Verbatim.Species %in% ""] <- Verbatim.Genus[Verbatim.Species %in% "sp." | Verbatim.Species %in% ""] #this sets Genus only desingations (i.e. Genus sp.) serpeately to avoid a space after he genus nae
  
  archaic.ung$Accepted.Name <- getCurrentTaxa(tax.vec = archaic.ung$Verbatim.taxon) #this does not work on most Genus sp. designations or entries that are in quotes.  Is hit or miss on entries with question marks
  archaic.ung$Accepted.Name <- gsub(pattern = "[[:space:]]", replacement = "_", x = archaic.ung$Accepted.Name)
  
  uniqTax <- lapply(c("Mammalia"), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  uniqTax <- rbind(uniqTax[[1]])
  uniqTax$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$accepted_name)
  
  archaic.ung$Accepted.Genus <- uniqTax$genus[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)] 
  
  #need to add a check for dealing with sub genera and subspecies
  archaic.ung$Accepted.Species <- unlist(lapply(strsplit(uniqTax$accepted_name[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)], "_"), function(x) x[2])) 
  archaic.ung$Accepted.Species[is.na(archaic.ung$Accepted.Species)] <- ""
  
  uniqTax$accepted_name[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)]
  
  archaic.ung$Family <- uniqTax$family[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)]

  archaic.ung$Order <- uniqTax$order[match(x=archaic.ung$Accepted.Name, table=uniqTax$accepted_name)]

  if(!is.null(save.file)) write.csv(archaic.ung, file = save.file)
     
  return(archaic.ung)
}

getTaxaInClade <- function(clades, occs, save.file=NULL) {
  uniqTax <- lapply(c(clades), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  uniq.len <- length(uniqTax)
 if(uniq.len == 1) { uniqTax <- rbind(uniqTax[[1]])
    }
  if(uniq.len == 2) { uniqTax <- rbind(uniqTax[[1]], uniqTax[[2]])
  }
  if(uniq.len > 2) 
  { 
    taxaMat <- rbind(uniqTax[[1]], uniqTax[[2]])
    for(xx in seq(3, length(uniqTax),1))
    {
      taxaMat <- rbind(taxaMat, uniqTax[[xx]])
    }
    uniqTax <- taxaMat
  }
  
 # uniqTax$accepted_species <- str_split_fixed(string = uniqTax[,"accepted_name"], " ", n = Inf)[,2]
  uniqTax$accepted_species <- str_split_fixed(string = uniqTax$accepted_name, " ", n = Inf)[,2]
  
  uniqTax$verbatim_genus <- str_split_fixed(string = uniqTax$taxon_name, " ", n = Inf)[,1]
  uniqTax$verbatim_species <- str_split_fixed(string = uniqTax$taxon_name, " ", n = Inf)[,2]
  
  uniqTax$accepted_name <- gsub(" ","_", uniqTax$accepted_name)
  
  uniqTax <- uniqTax[,c("phylum", "class", "order", "family", "genus", "accepted_species", "verbatim_genus", "verbatim_species","accepted_name", "taxon_name")]
  
  if(!is.null(occs)) uniqTax <- unique(uniqTax[uniqTax$accepted_name %in% occs$accepted_name[occs$accepted_rank =="species"],]) #This should remove non mammal Proboscidea designation in the PBDB query since occs is for all mammals
  
  if(!is.null(save.file)) write.csv(uniqTax, file = save.file)
  
  return(uniqTax)
}

getTaxonomyForOneBaseTaxon_AcceptedName <- function(this.taxon) {
  this.URL <- URLencode(paste0(server, "taxa/list.csv?base_name=", this.taxon, "&taxon_status=all&show=class"))
  this.names <- getStuffWithoutWarnings(this.URL)
  this.names[,c("phylum", "class", "order", "family", "genus", "accepted_name","taxon_name")]
}

getTargetTaxa<- function(measure.mat, uniqTax, occs, uniqOnly = FALSE, species.only = FALSE, save.file = NULL){ # get a list of taxa that indicates #occs it has and whether its in North America if it does
  
  uniqTax$NoOccs <- 0 
  
  uniqTax$InNorAmer <- FALSE
  
  uniqTax$DataCollected <- 0
  
  measure.mat$accepted_name <- paste(measure.mat$Accepted.Genus, measure.mat$Accepted.Species, sep = " ")
  
  dental.col <-  c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W",
  "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w")
  
  if(uniqOnly == TRUE) uniqTax <- unique(uniqTax[,c(1:6,9)])
  if(species.only == TRUE) uniqTax <- uniqTax[uniqTax$accepted_name %in% occs$accepted_name[occs$accepted_rank %in% "species"],]
  
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


data.coverage_v3 <- function(clades, clade.level = "family", clade.ranked = TRUE, data.mat, occs, measure.colnames, save.file = NULL)
{
  output.list <- list(AllTaxa = NA, OccsOnly = NA, MissingTaxa = NA) #, NoOccs = NA)
  
  #data.mat <- data.mat[, colnames(data.mat)[1:66]]
  data.mat <- data.mat[!data.mat$Catalog.Number %in% "Accepted Names PBDB",]
  data.mat$accepted_name <- paste(data.mat$Accepted.Genus, data.mat$Accepted.Species, sep = " ")
  ####remove rows that lack measurements
  data.mat <- data.mat[!apply(is.na(data.mat[,measure.colnames]) | data.mat[,measure.colnames] == "", 1, all),]
  
  #need this function to 
  ##1) read in a given taxon and then collate the genera and species within
  uniqTax <- lapply(c(clades), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  
  tax.mat <- matrix(nrow = 1, ncol = 6)
  colnames(tax.mat) <- c("Percent.Species", "No.Species", "Sampled.Species", "Percent.Genera", "No.Genera", "Sampled.Genera")
  
  for(xx in seq(1, length(uniqTax),1))
  {
    
    fam.list <- unique(uniqTax[[xx]]$family)[!unique(uniqTax[[xx]]$family) %in% ""]
    
    tax.mat <- matrix(nrow = length(fam.list), ncol = 6)
    rownames(tax.mat) <- c(fam.list)
    colnames(tax.mat) <- c("Percent.Species", "No.Species", "Sampled.Species", "Percent.Genera", "No.Genera", "Sampled.Genera")
    
    for(yy in fam.list)
    {
      uniqTax.temp <- uniqTax[[xx]][uniqTax[[xx]]$family %in% yy,]
      uniqTax.temp <- uniqTax.temp[!uniqTax.temp$genus %in% "",]
      
      tax.mat[rownames(tax.mat) %in% yy, c("No.Genera")] <- length(unique(uniqTax.temp$genus[uniqTax.temp$family %in% yy]))
      tax.mat[rownames(tax.mat) %in% yy, c("No.Species")] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$genus != "" & uniqTax.temp$accepted_name != uniqTax.temp$genus]))
      
      sample.genera <- unique(c(data.mat$Accepted.Genus, 
                                data.mat$Accepted.Genus[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus]))
      sample.genera <- sample.genera[sample.genera %in% unique(uniqTax.temp$genus[uniqTax.temp$family %in% yy])]
      tax.mat[rownames(tax.mat) %in% yy, c("Sampled.Genera")] <- length(sample.genera)
      
      
      sample.species <- unique(data.mat$accepted_name[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus])
      sample.species <- sample.species[sample.species%in% unique(uniqTax.temp$accepted_name[uniqTax.temp$family %in% yy])]
      tax.mat[rownames(tax.mat) %in% yy, c("Sampled.Species")] <-  length(sample.species)
    }
    
    tax.mat[, c("Percent.Genera")] <- (tax.mat[,"Sampled.Genera"]/tax.mat[,"No.Genera"])*100
    tax.mat[, c("Percent.Species")] <- (tax.mat[,"Sampled.Species"]/tax.mat[,"No.Species"])*100
   
    if(xx == 1) { output.list$AllTaxa <- list(tax.mat)
    } else { output.list$AllTaxa <- append(output.list$AllTaxa, list(tax.mat)) }
    
  }
  
 names(output.list$AllTaxa) <- clades

  ###########################################################################################################################################################################
  ##3) Generate a table for species with occs in database
  for(xx in seq(1, length(uniqTax),1))
  {
    fam.list <- unique(uniqTax[[xx]]$family)[!unique(uniqTax[[xx]]$family) %in% ""]
    
    tax.mat <- matrix(nrow = length(fam.list), ncol = 6)
    rownames(tax.mat) <- c(fam.list)
    colnames(tax.mat) <- c("Percent.Species", "No.Species", "Sampled.Species", "Percent.Genera", "No.Genera", "Sampled.Genera")
    
    for(yy in fam.list)
    {
      uniqTax.temp <- uniqTax[[xx]][uniqTax[[xx]]$family %in% yy & (uniqTax[[xx]]$accepted_name %in% occs$accepted_name | uniqTax[[xx]]$accepted_name %in% occs$genus),]
      uniqTax.temp <- uniqTax.temp[!uniqTax.temp$genus %in% "",]
      
      tax.mat[rownames(tax.mat) %in% yy, c("No.Genera")] <- length(unique(uniqTax.temp$genus[uniqTax.temp$family %in% yy]))
      tax.mat[rownames(tax.mat) %in% yy, c("No.Species")] <- length(unique(uniqTax.temp$accepted_name[uniqTax.temp$genus != "" & uniqTax.temp$accepted_name != uniqTax.temp$genus]))
      
      sample.genera <- unique(c(data.mat$Accepted.Genus, 
                                data.mat$Accepted.Genus[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus]))
      sample.genera <- sample.genera[sample.genera %in% unique(uniqTax.temp$genus[uniqTax.temp$family %in% yy])]
      tax.mat[rownames(tax.mat) %in% yy, c("Sampled.Genera")] <- length(sample.genera)
      
      sample.species <- unique(data.mat$accepted_name[!gsub(pattern = "[[:space:]]", replacement = "", x = data.mat$accepted_name) %in% data.mat$Accepted.Genus])
      sample.species <- sample.species[sample.species%in% unique(uniqTax.temp$accepted_name[uniqTax.temp$family %in% yy])]
      tax.mat[rownames(tax.mat) %in% yy, c("Sampled.Species")] <-  length(sample.species)

    }
    tax.mat[, c("Percent.Genera")] <- (tax.mat[,"Sampled.Genera"]/tax.mat[,"No.Genera"])*100
    tax.mat[, c("Percent.Species")] <- (tax.mat[,"Sampled.Species"]/tax.mat[,"No.Species"])*100
    
    if(xx == 1) { output.list$OccsOnly <- list(tax.mat)
    } else { output.list$OccsOnly <- append(output.list$OccsOnly, list(tax.mat)) }
  }

 names(output.list$OccsOnly) <- clades
  
  #get taxa that are missing
  
  UnSampled.Taxa<- list()
  
  for(xx in seq(1, length(uniqTax),1))
  {
    fam.list <- unique(uniqTax[[xx]]$family)[!unique(uniqTax[[xx]]$family) %in% ""]
    
    uniqTax.temp <- uniqTax[[xx]][uniqTax[[xx]]$family %in% fam.list & (uniqTax[[xx]]$accepted_name %in% occs$accepted_name | uniqTax[[xx]]$accepted_name %in% occs$genus),]
    uniqTax.temp <- uniqTax.temp[!(uniqTax.temp$accepted_name %in% data.mat$accepted_name),]
    uniqTax.temp <- uniqTax.temp[!uniqTax.temp$genus %in% "",]
      
    UnSampled.Genera <- uniqTax.temp[!uniqTax.temp$genus %in% unique(sample.genera),]
      
    UnSampled.Species <- uniqTax.temp[!uniqTax.temp$accepted_name %in% unique(sample.species),]
      
    UnSampled.Taxa[[xx]] <- unique(rbind(UnSampled.Genera, UnSampled.Species))
  }
  
  output.list$MissingTaxa <- UnSampled.Taxa
  names(output.list$MissingTaxa) <- clades
  
  if(!is.null(save.file)) write.csv(UnSampled.Taxa, file = save.file)
  
  return(output.list)
}

pred_herb_correl <- function(prop_herb = NULL, prop_pred = NULL, output.filename = NULL)
{
  png(output.filename, width = 7, height = 3, units = "in", res = 200)
  par(mfrow = c(1,2))
  
  all.val <- cbind(prop_pred, prop_herb)
  
  corr.results.both <- cor(prop_pred, prop_herb, method = "spearman")
  cor.p <- cor.mtest(all.val)$p
  
  cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(prop_herb)]
  cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(prop_pred),]
  
  corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                     cl.align.text = 'l', 
                     #addCoef.col = 'black',
                     addgrid.col = 'black',
                     method = "color",
                     na.label = "-",
                     sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                     pch.cex = 1,
                     tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
  
  #First Differences Spearman Correlation
  pred.prey.DivFirstDiff <- getDiversity1stDiff(data.mat = t(prop_pred), output.rownames = rownames(prop_pred))
  
  UngualteBMGroupDiv_FirstDiff <- getDiversity1stDiff(data.mat = t(prop_herb), output.rownames = rownames(prop_herb))
  
  ungulates <- t(as.matrix(UngualteBMGroupDiv_FirstDiff))
  predators <- t(as.matrix(pred.prey.DivFirstDiff))
  
  all.val <- cbind(predators, ungulates)
  
  corr.results.both <- cor(predators, ungulates, method = "spearman")
  cor.p <- cor.mtest(all.val)$p
  
  cor.p.culled <- cor.p[,colnames(cor.p) %in% colnames(ungulates)]
  cor.p.culled <- cor.p.culled[rownames(cor.p.culled) %in% colnames(predators),]
  
  corrplot::corrplot(corr.results.both,  mar = c(0,0,0,0), p.mat = cor.p.culled,
                     cl.align.text = 'l', 
                     #addCoef.col = 'black',
                     addgrid.col = 'black',
                     method = "color",
                     na.label = "-",
                     sig.level = c(0.05, 0.01, 0.001), insig = 'label_sig', 
                     pch.cex = 1,
                     tl.pos = c("lt"), tl.offset = 1, tl.cex = 1) -> p1
  
  dev.off()
}

plot.medianRichnessInBin <- function(countCube = NULL, rep.test = 1, 
                                     xlim = c(66,0), ylim = c(0,50), xaxt = "s", xlab = "Time (Ma)", ylab = "Number of Taxa", main = "Median Number of Taxa Per Time Bin",
                                     axis.label = TRUE)
{
  plot(0,0, col = "red", type = "n", xaxt = xaxt, xlim = c(66,0), ylim = ylim,
       xlab = xlab, ylab = ylab, main = main)
  axis(side = 1, at = seq(xlim[2], xlim[1], 5), labels = axis.label)
  
  for(xx in seq(1, rep.test, 1))
  {
    prop <- t(apply(countCube[,,c(1,rep.test[xx])], c(1,2), median, na.rm=TRUE))
    colnames(prop)[colnames(prop)==""] <- "indeterminate"
    
    lines(rowMeans(intervals), prop[,1], col = alphaColor("red", 0.05))
    lines(rowMeans(intervals), prop[,2], col = alphaColor("orange", 0.05))
    lines(rowMeans(intervals), prop[,3], col = alphaColor("green", 0.05))
    lines(rowMeans(intervals), prop[,4], col = alphaColor("blue", 0.05))
    lines(rowMeans(intervals), prop[,5], col = alphaColor("purple", 0.05))
  }
  #median across all replicates
  lines(rowMeans(intervals), prop[,1], col = alphaColor("red",1))
  lines(rowMeans(intervals), prop[,2], col = alphaColor("orange", 1))
  lines(rowMeans(intervals), prop[,3], col = alphaColor("green", 1))
  lines(rowMeans(intervals), prop[,4], col = alphaColor("blue", 1))
  lines(rowMeans(intervals), prop[,5], col = alphaColor("purple", 1))
}

estiamteMissingDiversity_Alroy <- function() {
  getThreeTimerRichnessFromOneOccList()
  getThreeTimerRatesFromCounts()
  
  return()
}

getBinsNALMA <- function(settings)
{
  nalma.mark <- read.csv("/Users/emdoughty/Dropbox/Code/R/dentalMeasurements/dat/NOW_intervals_edit.csv")
  nalma.mark <- nalma.mark[,1:3]; rownames(nalma.mark) <- nalma.mark$CHRON; colnames(nalma.mark) <- c("NALMA_Subdivision", "ageBase","ageTop")
  nalma.mark <- nalma.mark[,-1]
  nalma.mark <- nalma.mark[,c(2,1)]
  intervals <- nalma.mark[nrow(nalma.mark):1,]
  
  intervals <- intervals[intervals$ageTop >= settings$min_age & intervals$ageTop <= settings$max_age | intervals$ageBase >= settings$min_age & intervals$ageBase <= settings$max_age,]
  
  return(intervals)
}

getCountCube <- function(repIntTaxa, measure.mat, target.column = "bodyMass", bmBreaks, sizecateg, intervals)
{
  countCube <- sapply(repIntTaxa, function(this.rep) {
    sapply(this.rep, function(this.intv, this.rep) {
      hist(measure.mat[,target.column][match(this.intv, measure.mat$taxon)], 
           breaks= bmBreaks, plot=FALSE)$counts
    }, this.rep=this.rep)
  }, simplify = "array")
  
  #countCube <- countCube[,,1]
  dimnames(countCube) <- list(sizecateg, rownames(intervals), NULL)
  
  return(countCube)
}

sensitivity.bodymass.kmeans <- function(mom.data = NULL, k.max = 10)
{
  kmeans_results <- list()
  kmeans_breaks <- list()
  
  wss <- vector(length = k.max-1)
  sil <- vector(length = k.max-1)
  cal.hara <- vector(length = k.max-1)
  davies <- vector(length = k.max-1)
  
  for(ii in seq(2,k.max,1))
  {
    set.seed(1)
    kmeans_results[[ii]] <- kmeans(mom.data$LogMass..kg., centers = ii, iter.max = 10000, nstart = 50)
    
    wss[ii-1] <- kmeans_results[[ii]]$tot.withinss
    sil[ii-1] <- mean(silhouette(kmeans_results[[ii]]$cluster, dist(mom.data$LogMass..kg.))[,3])
    cal.hara[ii-1] <-  calinhara(mom.data$LogMass..kg., kmeans_results[[ii]]$cluster)
    davies[ii-1] <- index.DB(mom.data$LogMass..kg., kmeans_results[[ii]]$cluster, centrotypes="centroids")$DB
    
    kmean_Onebreaks <- matrix(ncol = 3, nrow = ii)
    colnames(kmean_Onebreaks) <- c("k", "min.kg", "max.kg")
    
    for(jj in seq(1, max(unique(kmeans_results[[ii]]$cluster)),1)) 
    {
      min.categ <- min(mom.data[kmeans_results[[ii]]$cluster == jj,"LogMass..kg."])
      max.categ <- max(mom.data[kmeans_results[[ii]]$cluster == jj,"LogMass..kg."])
      
      kmean_Onebreaks[jj,] <- c(jj, min.categ, max.categ)
    }
    kmeans_breaks[[ii]] <- kmean_Onebreaks[order(kmean_Onebreaks[,c("min.kg")]),] 
  }
  
  kmeans_out <- list(results = kmeans_results, breaks = kmeans_breaks, 
                     wss = wss, sil = sil, cal.hara = cal.hara)
    
  return(kmeans_out)
}

overlayCzTimescale_evan <- function(do.subepochs=FALSE, color=TRUE, 
                                    thisAlpha.intervals=0.33, thisAlpha.text = 0.33, 
                                    borderCol="white", invertTime=FALSE, 
                                    scale.cex=0.75, scale.headers = 0.95, text.offset = 0.025, 
                                    top = NULL, bottom = NULL, epoch.label.vert.move = 0, subepoch.label.vert.move = 0) {
  textCol <- rgb(0,0,0,thisAlpha.text)
  # textShadowCol<-"gray50"
  old.cex<-par()$cex
  par(cex=old.cex * scale.cex)
  epochs=data.frame(name=c("Camb",
                           "eO",
                           "mO",
                           "lO",
                           "eS",
                           "mS",
                           "lS",
                           "eD",
                           "mD",
                           "lD",
                           "Miss",
                           "Penn",
                           "eP",
                           "mP",
                           "lP",
                           "eT", 
                           "mT",
                           "lT",
                           "eJ",
                           "mJ",
                           "lJ",
                           "eK",
                           "lK",
                           "Paleo",
                           "Eo",
                           "Oligo",
                           "Mio",
                           "Plio",
                           "Pl.",
                           "Recent"),
                    ageBase =c(542,
                               486,
                               472,
                               461,
                               444,
                               428,
                               423,
                               416,
                               398,
                               385,
                               359,
                               318,
                               299,
                               271,
                               260,
                               
                               251.0,
                               245.0,
                               235.0,
                               201.6,
                               176.0,
                               161.0,
                               145.5,
                               99.6,
                               65.5,
                               55.8,
                               33.9,
                               23.03,
                               5.33,
                               2.58,
                               0))
  if (color) { epochs<-data.frame(epochs, rgb =c(
    rgb(0.533, 0.671, 0.494, thisAlpha.intervals),
    rgb(0.0, 0.686, 0.565, thisAlpha.intervals),
    rgb(0.118, 0.737, 0.624, thisAlpha.intervals),
    rgb(0.459, 0.796, 0.690, thisAlpha.intervals),
    rgb(0.565, 0.835, 0.788, thisAlpha.intervals),
    rgb(0.675, 0.871, 0.831, thisAlpha.intervals),
    rgb(0.725, 0.894, 0.867, thisAlpha.intervals),
    rgb(0.906, 0.698, 0.471, thisAlpha.intervals),
    rgb(0.953, 0.788, 0.557, thisAlpha.intervals),
    rgb(0.953, 0.875, 0.710, thisAlpha.intervals),
    rgb(0.427, 0.624, 0.533, thisAlpha.intervals),
    rgb(0.584, 0.769, 0.780, thisAlpha.intervals),
    rgb(0.937, 0.463, 0.404, thisAlpha.intervals),
    rgb(0.988, 0.549, 0.478, thisAlpha.intervals),
    rgb(0.996, 0.702, 0.647, thisAlpha.intervals),
    
    rgb(0.643, 0.365, 0.627, thisAlpha.intervals),
    rgb(0.718, 0.510, 0.71, thisAlpha.intervals),
    rgb(0.745, 0.616, 0.776, thisAlpha.intervals),
    rgb(0.0, 0.718, 0.906, thisAlpha.intervals),
    rgb(0.392, 0.816, 0.918, thisAlpha.intervals),
    rgb(0.647, 0.882, 0.973, thisAlpha.intervals),
    rgb(0.5803922, 0.7960784, 0.4745098, thisAlpha.intervals),
    rgb(0.7803922, 0.8784314, 0.6156863, thisAlpha.intervals),
    rgb(0.9803922, 0.6980392, 0.4862745, thisAlpha.intervals),
    rgb(0.9843137, 0.7372549, 0.5294118, thisAlpha.intervals),
    rgb(0.9960784, 0.8588235, 0.6745098, thisAlpha.intervals),
    rgb(1, 0.945098, 0, thisAlpha.intervals),
    rgb(1, 0.9764706, 0.6823529, thisAlpha.intervals),
    rgb(1, 0.9411765, 0.7490196, thisAlpha.intervals),
    rgb(1, 0.9529412, 0.9333333, thisAlpha.intervals)), stringsAsFactors = FALSE) 
  } else { epochs<-data.frame(epochs, rgb =c(
    rgb(0.57, 0.57, 0.57, thisAlpha.intervals),
    rgb(0.51, 0.51, 0.51, thisAlpha.intervals),
    rgb(0.59, 0.59, 0.59, thisAlpha.intervals),
    rgb(0.68, 0.68, 0.68, thisAlpha.intervals),
    rgb(0.79, 0.79, 0.79, thisAlpha.intervals),
    rgb(0.79, 0.79, 0.79, thisAlpha.intervals),
    rgb(0.79, 0.79, 0.79, thisAlpha.intervals),
    rgb(0.69, 0.69, 0.69, thisAlpha.intervals),
    rgb(0.77, 0.77, 0.77, thisAlpha.intervals),
    rgb(0.85, 0.85, 0.85, thisAlpha.intervals),
    rgb(0.51, 0.51, 0.51, thisAlpha.intervals),
    rgb(0.68, 0.68, 0.68, thisAlpha.intervals),
    rgb(0.53, 0.53, 0.53, thisAlpha.intervals),
    rgb(0.61, 0.61, 0.61, thisAlpha.intervals),
    rgb(0.73, 0.73, 0.73, thisAlpha.intervals),
    
    rgb(0.39, 0.39, 0.39, thisAlpha.intervals),
    rgb(0.51, 0.51, 0.51, thisAlpha.intervals),
    rgb(0.59, 0.59, 0.59, thisAlpha.intervals),
    rgb(0.58, 0.58, 0.58, thisAlpha.intervals),
    rgb(0.71, 0.71, 0.71, thisAlpha.intervals),
    rgb(0.81, 0.81, 0.81, thisAlpha.intervals),
    rgb(0.69, 0.69, 0.69, thisAlpha.intervals),
    rgb(0.81, 0.81, 0.81, thisAlpha.intervals),
    rgb(0.71, 0.71, 0.71, thisAlpha.intervals),
    rgb(0.75, 0.75, 0.75, thisAlpha.intervals),
    rgb(0.85, 0.85, 0.85, thisAlpha.intervals),
    rgb(0.93, 0.93, 0.93, thisAlpha.intervals),
    rgb(0.96, 0.96, 0.96, thisAlpha.intervals),
    rgb(0.93, 0.93, 0.93, thisAlpha.intervals),
    rgb(1, 1, 1, thisAlpha.intervals)), stringsAsFactors = FALSE) }
  # } else { epochs <- data.frame(epochs, rgb = rep("000000", times=nrow(epochs)), stringsAsFactors = FALSE) }
  
  if (do.subepochs) subepochs=data.frame(modifier=c("e", "m", "l", "e", "m", "l", "e", "l", "e", "m", "l", "PP"), ageBase=c(65.5, 61.7, 58.7, 55.8, 48.6, 37.2, 33.9, 28.4, 23.03, 15.97, 11.61, 5.33), stringsAsFactors = FALSE)
  
  if (invertTime) {
    epochs$ageBase<- par()$usr[2]/1.137059-epochs$ageBase
    if (do.subepochs) subepochs$ageBase<- par()$usr[2]/1.137059-subepochs$ageBase
  }
  
  if(is.null(top)) top=par()$usr[4]
  if(is.null(bottom)) bottom=par()$usr[3]
  
  #lays down colors
  for (i in 1:(nrow(epochs)-1)) {
    polygon(c(epochs[i,2], epochs[i,2], epochs[(i+1),2], epochs[(i+1),2]), c(bottom, top, top, bottom), col=epochs$rgb[i], border=borderCol, lty=0)
  }
  
  #lays down Era boundaries
  abline(v=251,lwd=1.5, col=borderCol)
  abline(v=65.5,lwd=1.5, col=borderCol)
  
  if (do.subepochs) {
    top=(scale.headers*(top - bottom)) + bottom
    for (i in 1:(nrow(subepochs)-1)) {
      polygon(cbind(c(subepochs[i,2], subepochs[i,2], subepochs[(i+1),2], subepochs[(i+1),2]), c(bottom, top, top, bottom)), lwd=0.25, border=borderCol)
      text (x=mean(c(as.numeric(subepochs[i,2]), as.numeric(subepochs[(i+1),2]))), y=(((scale.headers + text.offset)*(top - bottom)) + bottom + subepoch.label.vert.move), label=subepochs[i,1], col=textCol)
    }
  }
  
  #lays down borders and text
  if(is.null(top)) top=par()$usr[4]
  for (i in 1:(nrow(epochs)-1)) {
    polygon(cbind(c(epochs[i,2], epochs[i,2], epochs[(i+1),2], epochs[(i+1),2]), c(bottom, top, top, bottom)), col=NA, border=borderCol, lwd=0.5)
    # text (mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), (0.975*top), label=epochs[i,1], col=textShadowCol, adj=c(0.475,0.525), font=2) #cex=par()$cex+0.05
    # text (mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), (0.975*top), label=epochs[i,1], col=textShadowCol, font=2) #cex=par()$cex+0.05
    text (x=mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), y=(((scale.headers + text.offset)*(top - bottom)) + bottom + epoch.label.vert.move), label=epochs[i,1], col=textCol, font=2)
  }
  abline(h=((scale.headers*(top - bottom)) + bottom), lwd=1.0, col=borderCol)
  
  par(cex=old.cex)
}

sensitivity.MedianRichnessInBin <- function(countCube_herb = NULL, countCube_pred = NULL, intervals = NULL, number.reps = 1000, 
                                            ylim = c(0,60),
                                            output.filename = NULL, output.filepath = NULL)
{
  if(!is.null(output.filename)) { png(output.filename, width = 10, height = 6, units = "in", res = 200)
  } else {png(paste0(output.filepath, "medianRichnessInBin_", number.reps, "reps_", timestamp(),".png"), width = 11, height = 6, units = "in", res = 200)}

  rep.test <- seq(number.reps,dim(countCube_herb)[3], number.reps)
  
  par(mfrow = c(2,1), mar = c(0.5,1,1,1), oma = c(3,2,0,0))
    
  plot.medianRichnessInBin(countCube = countCube_herb, rep.test = rep.test,
                           xlim = c(66,0), ylim = ylim, xaxt = "n", xlab = "", ylab = "Number of Taxa", main = "Median Number of Taxa Per Time Bin",
                           axis.label = FALSE)
  overlayCzTimescale_evan(do.subepochs = TRUE, color = TRUE, borderCol = "black")
  plot.medianRichnessInBin(countCube = countCube_pred, rep.test = rep.test,
                           xlim = c(66,0), ylim = ylim/2, xlab = "Time (Ma)", ylab = "Number of Taxa", main = "")
  overlayCzTimescale_evan(do.subepochs = TRUE, color = TRUE, borderCol = "black")
  dev.off()
    
  return()
}


sensitivity.var.CorrelCoef.Single <- function(countCube_herb = NULL, countCube_pred = NULL, intervals = NULL, number.reps = 1000, 
                                              var.CorrelCoef.Single = c(1,1),
                                              ylim = NULL,
                                              output.filename = NULL, output.filepath = NULL)
{
  if(!is.null(output.filename)) { png(output.filename, width = 10, height = 6, units = "in", res = 200)
  } else {png(paste0(output.filepath,"var.CorrelCoef_single_", number.reps, "reps_", timestamp(),".png"), width = 11, height = 6, units = "in", res = 200)}
  
  rep.test <- seq(number.reps,dim(countCube_herb)[3], number.reps)
  
    corr.results.both <- vector()
    var.corr <- vector()
    
    for(xx in seq(1, rep.test, 1))
    {
      ##############################################################
      prop_herb <- t(apply(countCube_herb[,,seq(1, xx,1)], c(1,2), median, na.rm=TRUE))
      colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"
      
      prop_pred <- t(apply(countCube_pred[,,seq(1, xx,1)], c(1,2), median, na.rm=TRUE))
      colnames(prop_pred)[colnames(prop_pred)==""] <- "indeterminate"
      
      ungulates <- prop_herb
      predators <- prop_pred
      
      all.val <- cbind(predators, ungulates)
      
      corr.results.both[xx] <- cor(predators, ungulates, method = "spearman")[var.CorrelCoef.Single[1],var.CorrelCoef.Single[2]] #for antelope vs wolf sized critters
      
      var.corr[xx] <- var(corr.results.both)
    }
    plot(seq(1, dim(countCube_herb)[3],1), var.corr, col = alphaColor("gray75",0.5),
         xlab = "Cumulative Number of Replicates", ylab = "Variance of Correlation Coefficient of Median Assemblage")
    lines(seq(1, dim(countCube_herb)[3],1), var.corr, col = "black")
  
  dev.off()
  return()
}

sensitivity.var.CorrelCoef.All <- function(countCube_herb = NULL, countCube_pred = NULL, intervals = NULL, number.reps = 1000, 
                                           ylim = NULL, FD.correl = FALSE,
                                           output.filename = NULL, output.filepath = NULL)
{
  if(!is.null(output.filename)) { png(output.filename, width = 10, height = 6, units = "in", res = 200)
  } else {png(paste0(output.filepath, "var.CorrelCoef.All_", number.reps, "reps_", timestamp(),".png"), width = 11, height = 6, units = "in", res = 200)}
  
  rep.test <- seq(number.reps,dim(countCube_herb)[3], number.reps)
  
  corr.results.both <- array(numeric(0),dim=c(nrow(countCube_herb),nrow(countCube_pred),dim(countCube_herb)[3]))
  var.corr <- vector()
  for(zz in seq(1, rep.test, 1))
  {
     ##############################################################
    prop_herb <- t(apply(countCube_herb[,,seq(1, zz, 1)], c(1,2), median, na.rm=TRUE))
    colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"
      
    prop_pred <- t(apply(countCube_pred[,,seq(1, zz, 1)], c(1,2), median, na.rm=TRUE))
    colnames(prop_pred)[colnames(prop_pred)==""] <- "indeterminate"
      
    if(FD.correl) {
      predators <- getDiversity1stDiff(data.mat = t(prop_pred), output.rownames = rownames(prop_pred))
      ungulates <- getDiversity1stDiff(data.mat = t(prop_herb), output.rownames = rownames(prop_herb))
      
    } else{
      ungulates <- prop_herb
      predators <- prop_pred
    }
      
    all.val <- cbind(predators, ungulates)
    
    corr.results.both[,,zz] <- cor(predators, ungulates, method = "spearman") 
  }  
  dimnames(corr.results.both) <- dimnames(cor(predators, ungulates, method = "spearman"))
  
  if(length(ylim) == 2) ylim <- array(c(rep(ylim[1], times = nrow(countCube_pred)), rep(ylim[2],times = nrow(countCube_pred))), dim = c(nrow(countCube_pred), 2, nrow(countCube_herb))) # = c(2, nrow(countCube_pred), nrow(countCube_herb)))
    
  par(mfrow=c(nrow(countCube_herb),nrow(countCube_pred)), mar = c(1,4.5,0,0), oma = c(4,4,3,3))
  for(xx in seq(1, nrow(countCube_pred),1))
  {
    for(yy in seq(1, nrow(countCube_herb),1))
    {
        
      for(zz in seq(1, rep.test, 1))
      {
        var.corr[zz] <- var(corr.results.both[xx,yy,c(seq(1,zz,1))])
      }
        
       if(!is.null(ylim)){
         plot(0, 0, type = "n", xlim = c(0, dim(countCube_herb)[3]), ylim = ylim[yy,,xx],
              main = "", xlab = "", ylab = "", axes= FALSE, yaxt = "n", xaxt = "n")
          
        axis(2, labels = FALSE)
        axis(1, labels = FALSE)
        if(yy == 1) axis(2, labels = TRUE)
        if(xx == max(nrow(countCube_pred))) axis(1, labels = TRUE)
          
        points(seq(1, dim(countCube_herb)[3],1), var.corr, col = alphaColor("gray75",0.5))
        lines(seq(1, dim(countCube_herb)[3],1), var.corr, col = "black")
          
        if(yy == 1 & xx == ceiling(nrow(countCube_pred)/2)){ mtext("Variance of Correlation Coefficient of Median Assemblage", side = 2, line = 3, cex = 1)}
        if(yy == ceiling(nrow(countCube_herb)/2) & xx == max(nrow(countCube_pred))){ mtext("Number of Cumulative Replicates", side = 1, line = 2, cex = 1)}
        if(xx == 1) { mtext(rownames(countCube_herb)[yy], side = 3, line = 1, cex = 0.5)}
        if(yy == 1) { mtext(rownames(countCube_pred)[xx], side = 2, line = 2, cex = 0.5)}
          
      } else {
        plot(0, 0, type = "n", xlim = c(0, dim(countCube_herb)[3]), ylim = c(0, max(var.corr, na.rm = T)),
             main = "", xlab = "", ylab = "", axes= FALSE, yaxt = "n", xaxt = "n")
          
        options(scipen = -2, digits = 3)
        axis(2, labels = TRUE, lwd = 1, las = 1, line = 0.5)
        options(scipen = 0, digits = 7)
          
        axis(1, labels = FALSE)
        if(xx == max(nrow(countCube_pred))) axis(1, labels = TRUE)
          
        points(seq(1, dim(countCube_herb)[3],1), var.corr, col = alphaColor("gray75",0.5))
        lines(seq(1, dim(countCube_herb)[3],1), var.corr, col = "black")
          
        if(yy == 1 & xx == ceiling(nrow(countCube_pred)/2)){ mtext("Variance of Correlation Coefficient of Median Assemblage", side = 2, line = 6, cex = 1)}
        if(yy == ceiling(nrow(countCube_herb)/2) & xx == max(nrow(countCube_pred))){ mtext("Number of Cumulative Replicates", side = 1, line = 2.5, cex = 1)}
        if(xx == 1) { mtext(rownames(countCube_herb)[yy], side = 3, line = 1, cex = 0.75)}
        if(yy == 1) { mtext(rownames(countCube_pred)[xx], side = 2, line = 4.5, cex = 0.75)}
      }
        
      print(paste0(yy," ",xx," ",zz))
    }
  }
  
  dev.off()
  return()
}

sensitivity.CorrelCoef.All <- function(countCube_herb = NULL, countCube_pred = NULL, intervals = NULL, number.reps = 1000, 
                                       ylim = NULL, FD.correl = FALSE,
                                       output.filename = NULL, output.filepath = NULL)
{
  if(!is.null(output.filename)) { png(output.filename, width = 10, height = 6, units = "in", res = 200)
  } else {png(paste0(output.filepath, "CorrelCoef.All_", number.reps, "reps_", timestamp(),".png"), width = 11, height = 6, units = "in", res = 200)}
  
  rep.test <- seq(number.reps,dim(countCube_herb)[3], number.reps)
  
  corr.results.both <- array(numeric(0),dim=c(nrow(countCube_herb),nrow(countCube_pred),dim(countCube_herb)[3]))
  corr.coef <- vector()
  for(zz in seq(1, rep.test, 1))
  {
    ##############################################################
    prop_herb <- t(apply(countCube_herb[,,seq(1, zz, 1)], c(1,2), median, na.rm=TRUE))
    colnames(prop_herb)[colnames(prop_herb)==""] <- "indeterminate"
    
    prop_pred <- t(apply(countCube_pred[,,seq(1, zz, 1)], c(1,2), median, na.rm=TRUE))
    colnames(prop_pred)[colnames(prop_pred)==""] <- "indeterminate"
    
    if(FD.correl) {
      predators <- getDiversity1stDiff(data.mat = t(prop_pred), output.rownames = rownames(prop_pred))
      ungulates <- getDiversity1stDiff(data.mat = t(prop_herb), output.rownames = rownames(prop_herb))
      
    } else{
      ungulates <- prop_herb
      predators <- prop_pred
    }
    all.val <- cbind(predators, ungulates)
    
    corr.results.both[,,zz] <- cor(predators, ungulates, method = "spearman")
    
  }
  dimnames(corr.results.both) <- dimnames(cor(predators, ungulates, method = "spearman"))
  
  if(length(ylim) == 2) ylim <- array(c(rep(ylim[1], times = nrow(countCube_pred)), rep(ylim[2],times = nrow(countCube_pred))), dim = c(nrow(countCube_pred), 2, nrow(countCube_herb))) # = c(2, nrow(countCube_pred), nrow(countCube_herb)))
  
  par(mfrow=c(nrow(countCube_herb),nrow(countCube_pred)), mar = c(1,4.5,0,0), oma = c(4,4,3,3))
  for(xx in seq(1, nrow(countCube_pred),1))
  {
    for(yy in seq(1, nrow(countCube_herb),1))
    {
      
      if(!is.null(ylim)){
        plot(0, 0, type = "n", xlim = c(0, dim(countCube_herb)[3]), ylim = ylim[yy,,xx],
             main = "", xlab = "", ylab = "", axes= FALSE, yaxt = "n", xaxt = "n")
        
        axis(2, labels = FALSE)
        axis(1, labels = FALSE)
        if(yy == 1) axis(2, labels = TRUE)
        if(xx == max(nrow(countCube_pred))) axis(1, labels = TRUE)
        
        points(seq(1, dim(countCube_herb)[3],1), apply(corr.results.both,c(3), function(x) x[xx,yy]), col = alphaColor("gray75",0.5))
        lines(seq(1, dim(countCube_herb)[3],1), apply(corr.results.both,c(3), function(x) x[xx,yy]), col = "black")
        
        if(yy == 1 & xx == ceiling(nrow(countCube_pred)/2)){ mtext("Correlation Coefficient of Median Assemblage", side = 2, line = 6, cex = 1)}
        if(yy == ceiling(nrow(countCube_herb)/2) & xx == max(nrow(countCube_pred))){ mtext("Number of Cumulative Replicates", side = 1, line = 2.5, cex = 1)}
        if(xx == 1) { mtext(rownames(countCube_herb)[yy], side = 3, line = 1, cex = 0.75)}
        if(yy == 1) { mtext(rownames(countCube_pred)[xx], side = 2, line = 4.5, cex = 0.75)}
        
      } else {
        plot(0, 0, type = "n", xlim = c(0, dim(countCube_herb)[3]), ylim = c(min(apply(corr.results.both,c(3), function(x) x[xx,yy]), na.rm = T), max(apply(corr.results.both,c(3), function(x) x[xx,yy]), na.rm = T)),
             main = "", xlab = "", ylab = "", axes= FALSE, yaxt = "n", xaxt = "n")
        
        options(scipen = -2, digits = 3)
        axis(2, labels = TRUE, lwd = 1, las = 1, line = 0.5)
        options(scipen = 0, digits = 7)
        
        axis(1, labels = FALSE)
        if(xx == max(nrow(countCube_pred))) axis(1, labels = TRUE)
        
        points(seq(1, dim(countCube_herb)[3],1), apply(corr.results.both,c(3), function(x) x[xx,yy]), col = alphaColor("gray75",0.5))
        lines(seq(1, dim(countCube_herb)[3],1), apply(corr.results.both,c(3), function(x) x[xx,yy]), col = "black")
        
        if(yy == 1 & xx == ceiling(nrow(countCube_pred)/2)){ mtext("Correlation Coefficient of Median Assemblage", side = 2, line = 6, cex = 1)}
        if(yy == ceiling(nrow(countCube_herb)/2) & xx == max(nrow(countCube_pred))){ mtext("Number of Cumulative Replicates", side = 1, line = 2.5, cex = 1)}
        if(xx == 1) { mtext(rownames(countCube_herb)[yy], side = 3, line = 1, cex = 0.75)}
        if(yy == 1) { mtext(rownames(countCube_pred)[xx], side = 2, line = 4.5, cex = 0.75)}
      }
      
      print(paste0(xx," ",yy," ",zz))
    }
  }
  
  dev.off()
  return()
}

plotSpeciesInGenera <- function(repIntTaxa, measure.mat, bmBreaks, sizecateg, intervals, do.reps=100)
{
  ##############This is a work in progress
  
  #check repIntTaxa to determine how many species are to a genera within a given bin.
  ##need to relegate this to only sampling congenerics so that speciose genera with longer duration (canis) don't overestimate by taking all at once
  ##get genus name, compare with PBDB for species filtered for time bin is in, if valid count them
  
  #remove non target taxa
  
  #compare interval (species level) with PBDB occs

  rownames(measure.mat) <- measure.mat$taxon
  measure.mat.genus <- makeOneGenusMatFromSpecimenMat(measure.mat)
  
  for(xx in seq(1, length(settings$bmBreaks_herb)-1, 1)){
    measure.mat.genus$SizeCat[measure.mat.genus$bodyMass > bmBreaks[xx] & measure.mat.genus$bodyMass < bmBreaks[xx+1]] <- xx
  } 
  
  measure.mat.genus$SizeCat[rownames(measure.mat.genus) %in% unique(occs$genus[occs$order %in% "Proboscidea"])] <- 5
  
  gen.count.inrep.list <- list()
  gen.count.alrep.list <- list()
  
  max.rep <- do.reps
  for(zz in seq(1, max.rep,1))
  {
    for(yy in seq(1,length(repIntTaxa[[zz]]),1))
    {
      genus.vec <- unique(occs$genus[occs$accepted_name %in% repIntTaxa[[zz]][[yy]] & occs$accepted_name %in% measure.mat$taxon]) # get list of genera
      genus.vec <- genus.vec[order(genus.vec)]
      if(length(genus.vec)>0)
      {
        genus.count.mat <- data.frame(genus = genus.vec[order(genus.vec)], sp.per.gen.in.bin = NA, size.cat = measure.mat.genus$SizeCat[match(genus.vec[order(genus.vec)],measure.mat.genus$taxon)]) #setup matrix to store values
      } else {
        genus.count.mat <- data.frame(genus = NA, sp.per.gen.in.bin = NA, size.cat = NA)
      }
      
      for(xx in seq(1, nrow(genus.count.mat))) 
      {
        #add code to seperate into size cats too
        sp.in.gen <- unique(occs$accepted_name[occs$genus %in% genus.count.mat[xx,1] & occs$accepted_rank %in% "species"]) # get unique species in a genus
        genus.count.mat[xx,2] <- length(sp.in.gen[sp.in.gen %in% repIntTaxa[[zz]][[yy]]]) #check that species are in the bin (may be missing as species duration not in bin, not sampled, etc.)
      }
      
      print(genus.count.mat)
      
      gen.count.list[[yy]] <- data.frame(mean=mean(genus.count.mat$sp.per.gen.in.bin), 
                                         median=median(genus.count.mat$sp.per.gen.in.bin),
                                         min=min(genus.count.mat$sp.per.gen.in.bin),
                                         max=max(genus.count.mat$sp.per.gen.in.bin))
    }
    names(gen.count.list) <- names(repIntTaxa[[1]])
    gen.count.alrep.list[[zz]] <- gen.count.list
    
    print(paste0(zz," out of ", max.rep))
  }
  
  gen.count.alrep.list[[1]]

#need way to reliably condense into median
 genusCube.mean <- rowMeans(sapply(gen.count.alrep.list, function (this.rep) {
                                   sapply(this.rep, function(this.intv, this.rep) {
                                          as.numeric(this.intv$mean)
                                   }, this.rep=this.rep)
                                 },simplify = "array"))
  
 genusCube.median <- sapply(gen.count.alrep.list, function (this.rep) {
                            sapply(this.rep, function(this.intv, this.rep) {
                                   as.numeric(this.intv$median)
                            }, this.rep=this.rep)
                           },simplify = "array")
 genusCube.median <- apply(genusCube.median, c(1), median, na.rm=FALSE)
 
 genusCube.min <- sapply(gen.count.alrep.list, function (this.rep) {
                         sapply(this.rep, function(this.intv, this.rep) {
                                as.numeric(this.intv$min)
                         }, this.rep=this.rep)
                        },simplify = "array")
 genusCube.min <- apply(genusCube.min, c(1), min, na.rm=FALSE)
 
 genusCube.max <- sapply(gen.count.alrep.list, function (this.rep) {
                         sapply(this.rep, function(this.intv, this.rep) {
                                as.numeric(this.intv$max)
                         }, this.rep=this.rep)
                        },simplify = "array")
 genusCube.max <- apply(genusCube.max, c(1), max, na.rm=FALSE)

 plot(rowMeans(intervals),genusCube.median, xlim = c(66,0), ylim = c(0,15), type = "n",
      xlab = "Time (Ma)", ylab = "Species Per Genera") 
 lines(rowMeans(intervals), genusCube.median, col ="black", lwd = 2)
 lines(rowMeans(intervals), genusCube.mean, col ="red", lwd = 2)
 lines(rowMeans(intervals), genusCube.max, col = "darkgreen", lwd = 2)
 legend(x = 20, y = 15, legend = c("median", "mean", "max"), fill = c("black","red", "darkgreen"))
  
 return()
}
