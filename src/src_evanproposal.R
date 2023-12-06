data.coverage <- function(clades, clade.level = "family", clade.ranked = TRUE, data.mat, occs, measure.colnames, save.file = NULL)
{
  output.list <- list(AllTaxa = NA, OccsOnly = NA, MissingTaxa = NA) #, NoOccs = NA)
  
  #data.mat <- data.mat[, colnames(data.mat)[1:66]]
  #data.mat <- data.mat[!data.mat$Catalog.Number %in% "Accepted Names PBDB",]
  #data.mat$accepted_name <- paste(data.mat$Accepted.Genus, data.mat$Accepted.Species, sep = " ")
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
    
    uniqTax[[xx]]$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax[[xx]]$accepted_name)
    
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

getCorrelationPlot <- function(prop1, prop2, correlation.type = "spearman",mar = c(0,0,0,0), cl.pos = NULL, tl.pos = c("lt"), sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', number.digits = 3)
{
  corr.results.both <- cor.p <- matrix(nrow = ncol(prop1), ncol= ncol(prop2)) 
  dimnames(corr.results.both) <- dimnames(cor.p) <- list(colnames(prop1), colnames(prop2))
  
  if(correlation.type == "pearson") finite.obs <- 3
  if(correlation.type == "spearman") finite.obs <- 2
  
  for (xx in seq_len(ncol(prop1))) {
    for (yy in seq_len(ncol(prop2))) {
      if(sum(complete.cases(prop1[,xx], prop2[,yy]), na.rm = TRUE) < finite.obs) { corr.results.both[xx,yy]  <- NA } else {
        
        corr.results.temp <- cor.test(prop1[,xx], prop2[,yy], method = correlation.type)
        
        corr.results.both[xx,yy] <- corr.results.temp$estimate
        cor.p[xx,yy] <- corr.results.temp$p.value
      }
    }
  }
  
  corrplot::corrplot(corr.results.both,  mar = mar, p.mat = cor.p,
                     cl.align.text = 'l', cl.pos = cl.pos, cl.cex = 1,
                     #addCoef.col = 'black',
                     addgrid.col = 'black',
                     method = "color",
                     na.label = "-",
                     sig.level = sig.level, insig = insig, number.digits = number.digits,
                     pch.cex = 1, number.cex = 1,
                     tl.pos = tl.pos, tl.offset = 1, tl.cex = 1) -> p1
  
  return(list(correl.coef = corr.results.both, p_value = cor.p))
}

getCountCube <- function(repIntTaxa, measure.mat, target.column = "bodyMass", bmBreaks, sizecateg, intervals)
{
  countCube <- sapply(repIntTaxa, function(this.rep) {
    sapply(this.rep, function(this.intv, this.rep) {
      hist(measure.mat[,target.column][match(this.intv, measure.mat$taxon)], 
           breaks= bmBreaks, plot=FALSE)$counts
    }, this.rep=this.rep)
  }, simplify = "array")
  
  dimnames(countCube) <- list(sizecateg, rownames(intervals), NULL)
  
  return(countCube)
}

getCurrentHigherTaxonomy <- function(archaic.ung, focal.taxa = "Mammalia", save.file=NULL) { #this function is meant to update the .csv files (or dataframe) of dental measures and update the taxonomy so that the accepted genus and species and higher taxonomy is correct (as much as possible).  Made for use with my Conylarth dental measurement sheet.
  require(stringr)
  #remove question marks and quotes form verbatim names
  Verbatim.Genus <- gsub("\"", "", archaic.ung$Verbatim.Genus)
  Verbatim.Species <- gsub("\"", "", archaic.ung$Verbatim.Species)
  
  Verbatim.Genus <- gsub("[?]", "", Verbatim.Genus)
  Verbatim.Species <- gsub("[?]", "", Verbatim.Species)
  
  Verbatim.Genus <- gsub("cf[[:punct:]] ", "", Verbatim.Genus)
  Verbatim.Species <- gsub("cf[[:punct:]] ", "", Verbatim.Species)
  
  Verbatim.Genus <- gsub("aff[[:punct:]] ", "", Verbatim.Genus)
  Verbatim.Species <- gsub("aff[[:punct:]] ", "", Verbatim.Species)
  
  archaic.ung$verbatim.name <- NA
  
  archaic.ung$verbatim.name[!Verbatim.Species %in% ""] <- paste(Verbatim.Genus[!Verbatim.Species %in% ""], Verbatim.Species[!Verbatim.Species %in% ""], sep=" ") 
  archaic.ung$verbatim.name[Verbatim.Species %in% "sp." | Verbatim.Species %in% ""] <- Verbatim.Genus[Verbatim.Species %in% "sp." | Verbatim.Species %in% ""] #this sets Genus only designations (i.e. Genus sp.) separately to avoid a space after he genus name.  Also prevents sp. from conflicting with PBDB.
  
  archaic.ung$accepted_name <- getCurrentTaxa(tax.vec = archaic.ung$verbatim.name)
  archaic.ung$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = archaic.ung$accepted_name)
  
  if(length(focal.taxa) == 1) 
  {
    uniqTax <- lapply(c(focal.taxa), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
    uniqTax <- uniqTax[[1]]
  } else {
    uniqTax <- lapply(unlist(focal.taxa), FUN=getTaxonomyForOneBaseTaxon_AcceptedName) 
    uniqTax <- makeMatrixFromList(uniqTax)
  }
  uniqTax$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$accepted_name)
  
  #need to add a check for dealing with sub genera and subspecies
  #get the genus and species from $accepted_name so that it relates to occurrences.  This would otherwise break if genus was taken from uniqTax as the generic designation sometimes differs from the genus within the $accepted_name binomial (e.g., Tricentes subtrigonus will be under Chriacus)
  archaic.ung$Accepted.Genus <- unlist(lapply(strsplit(uniqTax$accepted_name[match(x=archaic.ung$accepted_name, table=uniqTax$accepted_name)], "_"), function(x) x[1])) 
  archaic.ung$Accepted.Species <- unlist(lapply(strsplit(uniqTax$accepted_name[match(x=archaic.ung$accepted_name, table=uniqTax$accepted_name)], "_"), function(x) x[2])) 
  archaic.ung$Accepted.Species[is.na(archaic.ung$Accepted.Species)] <- ""
  
  uniqTax$accepted_name[match(x=archaic.ung$accepted_name, table=uniqTax$accepted_name)]
  
  #get higher taxonomy
  archaic.ung$Genus <- uniqTax$genus[match(x=archaic.ung$accepted_name, table=uniqTax$accepted_name)] #some generic designations differ between $accepted_name and the genus that is returned in uniqTax (e.g., Tricentes subtrigonus will be under Chriacus)
  archaic.ung$Family <- uniqTax$family[match(x=archaic.ung$accepted_name, table=uniqTax$accepted_name)]
  archaic.ung$Order <- uniqTax$order[match(x=archaic.ung$accepted_name, table=uniqTax$accepted_name)]
  
  #names(archaic.ung)[colnames(archaic.ung) %in% "accepted_name"] <- "accepted_name"
  names(archaic.ung)[colnames(archaic.ung) %in% "verbatim.name"] <- "identified.name"
  
  if(!is.null(save.file)) write.csv(archaic.ung, file = save.file)
  
  return(archaic.ung)
}

getDiversity1stDiff <- function(data.mat, output.rownames = NULL)
{
  div.diff <- matrix(nrow = nrow(data.mat)-1, ncol = ncol(data.mat))
  for(xx in seq(1, ncol(data.mat),1))
  {
    #div.diff[xx,] <- diff(data.mat[xx,], lag = 1) 
    div.diff[,xx] <- rev(diff(rev(data.mat[,xx])))
  }
  
  dat.rownames <- vector()
  for(xx in seq(2, length(output.rownames),1))
  {
    dat.rownames[xx-1] <- paste(output.rownames[xx-1], 
                                output.rownames[xx], sep="-")
  }
  rownames(div.diff) <- dat.rownames
  colnames(div.diff) <- colnames(data.mat)
  
  return(div.diff)
}

getDateTaxa <- function(measure.mat, occs, this.rank = "species")
{
  thisRanges <- getTaxonRangesFromOccs(occs = occs, rank = this.rank, random=FALSE)

  rownames(thisRanges) <- gsub(pattern = "[[:space:]]", replacement = "_", x = rownames(thisRanges))
  measure.mat[,c("FO","LO")] <- thisRanges[match(measure.mat$taxon, rownames(thisRanges)),]
  return(measure.mat)
}

getNALMAIntervals <- function(startDate = 1, endDate = 86)
{
  if (startDate>endDate) { 
    holder<-startDate
    startDate<-endDate
    endDate<-holder
  }
  #Janis 1998 volume designations
  NALMA <- rev(c("Aquilian","Judithian", "Lancian", "Puercan", "Torrejonian", "Tiffanian", "Clarkforkian", "Wasatchian", "Bridgerian", 
                 "Uintan", "Duchesnean", "Chadronian", "Orellan","Whitneyan", "Arikareean", 
                 "Hemingfordian", "Barstovian", "Clarendonian", "Hemphillian", "Blancan", 
                 "Irvingtonian", "Rancholabrean"))
  
  interval.start <- rev(c(86, 84, 70, 66.043, 63.3, 61.7, 56.8, 55.8, 50.3, 46.2, 40.4, 37.2, 33.9, 33.3, 30.8,
                          20.43, 15.97, 13.6, 10.3, 4.9, 1.8, 0.3))
  interval.end <- rev(c(84, 70, 66.043, 63.3, 61.7, 56.8, 55.8, 50.3, 46.2, 40.4, 37.2, 33.9,33.3, 30.8, 20.43, 15.97,
                        13.6, 10.3, 4.9, 1.8, 0.3, 0.012))
  nalma.intervals <- as.data.frame(cbind(ageTop = interval.end, ageBase = interval.start))
  
  rownames(nalma.intervals) <- NALMA
  
  #Remove entries not within bounds of dates set
  #if(is.null(startDate)) 
  
  nalma.intervals <- nalma.intervals[nalma.intervals$ageBase >= startDate & nalma.intervals$ageTop <= endDate,]
  
  #need way to handle dates falling within intervals
  ####rbind these onto one another
  #nalma.intervals[nalma.intervals$agetop <= startDate & nalma.intervals$ageBase >= endDate,]
  
  return(nalma.intervals)
}

getSizeCategString <- function(bmBreaks, sig.fig = 2)
{
  sizeCateg <- vector()
  for(xx in seq(1, length(bmBreaks)-1,1))
  {
    if(is.infinite(bmBreaks[xx]) & !is.infinite(bmBreaks[xx+1])) { sizeCateg[xx] <- paste0("<",round(10^bmBreaks[xx+1],sig.fig)," kg") } else if(!is.infinite(bmBreaks[xx]) & is.infinite(bmBreaks[xx+1])) { sizeCateg[xx] <- paste0(">",round(10^bmBreaks[xx],sig.fig), " kg")} else {sizeCateg[xx] <- paste0(round(10^bmBreaks[xx],sig.fig),"-",round(10^bmBreaks[xx+1],sig.fig), " kg")}
  }
  return(sizeCateg)
}

getTargetTaxa<- function(measure.mat, uniqTax, occs, measure.columns = c("P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W",
                                                                         "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w",
                                                                         "M1.3_L", "m1.3_L", "P4.M3_L", "p4.m3_l", "P3.M3_L", "p3.m3_l", "PCRU", "PCRL"), 
                         accepted_nameOnly = FALSE, species.only = FALSE, save.file = NULL)
{ # get a list of taxa that indicates #occs it has and whether its in North America if it does
  
  uniqTax$NoOccs <- 0 
  
  uniqTax$InNorAmer <- FALSE
  
  uniqTax$DataCollected <- 0
  
  if(!"accepted_name" %in% colnames(measure.mat)) {print("Function requires an accepted_name column");return()}
  
  dental.col <- measure.columns
  
  if(accepted_nameOnly) uniqTax <- unique(uniqTax[,c(1:5)])
  if(species.only) uniqTax <- uniqTax[uniqTax$accepted_name %in% occs$accepted_name[occs$accepted_rank %in% "species"],]
  
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

getTaxaInClade <- function(clades, occs, save.file=NULL) 
{
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
  
  uniqTax$accepted_name <- gsub(" ","_", uniqTax$accepted_name)
  
  uniqTax <- uniqTax[,c("phylum", "class", "order", "family", "accepted_name", "taxon_name")]
  
  uniqTax <- uniqTax[order(uniqTax$phylum, uniqTax$class, uniqTax$order, uniqTax$family, uniqTax$accepted_name),]
  
  if(!is.null(occs)) 
  {
    occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)
    uniqTax <- unique(uniqTax[uniqTax$accepted_name %in% occs$accepted_name[occs$accepted_rank %in% c("species", "genus")],]) #This should remove non mammal Proboscidea designation in the PBDB query since occs is for all mammals
  }
  if(!is.null(save.file)) write.csv(uniqTax, file = save.file)
  
  return(uniqTax)
}

getTaxonomyForOneBaseTaxon_AcceptedName <- function(this.taxon) 
{
  this.URL <- URLencode(paste0(server, "taxa/list.csv?base_name=", this.taxon, "&taxon_status=all&show=class"))
  this.names <- getStuffWithoutWarnings(this.URL)
  this.names[,c("phylum", "class", "order", "family", "genus", "accepted_name","taxon_name")]
}

correl.scatter.plot <- function(prop1, prop2, point.labels,
                                mar = c(2,2,1,1), oma = c(2,4,4,2),
                                axis.lims = FALSE, xlim = NULL, xlim.nudge = c(0,0), ylim = NULL, ylim.nudge=c(0,0),
                                col = alphaColor(rev(COL2("RdYlBu", n = length(prop1[,yy]))), alpha = 1),
                                plot.single.correl = NULL)

{
  if(!is.null(plot.single.correl)){
    par.mfrow.rows <- 1
    par.mfrow.cols <- 1
    
    for.xx <- plot.single.correl[1]
    for.yy <- plot.single.correl[2]
  } else {
    par.mfrow.rows <- ncol(prop2) #y axis within plot but is rows on array
    par.mfrow.cols <- ncol(prop1) #x axis within plot but is cols on array
    
    for.xx <- seq_len(ncol(prop2))
    for.yy <- seq_len(ncol(prop1))
  }
  
  par(mfrow = c(par.mfrow.rows,par.mfrow.cols), mar = mar, oma = oma)
  for(xx in for.xx)
  {
    for(yy in for.yy)
    {
      if(!any(complete.cases(cbind(prop1[,yy],prop2[,xx])))) 
      { 
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      } else {
        if(axis.lims & is.null(xlim)) {
          xlim.this.plot <- c(round_any(min(prop1[complete.cases(cbind(prop1[,yy],prop2[,xx])),yy], na.rm = TRUE),5, f= floor) + xlim.nudge[1],
                              round_any(max(prop1[complete.cases(cbind(prop1[,yy],prop2[,xx])),yy], na.rm = TRUE),5, f= ceiling)+xlim.nudge[2]) 
        } else if(axis.lims & !is.null(xlim)) {
          xlim.this.plot <- xlim
        } else {xlim.this.plot <- NULL}
        
        if(axis.lims & is.null(ylim)) {
          ylim.this.plot <- c(round_any(min(prop2[complete.cases(cbind(prop1[,yy],prop2[,xx])),xx], na.rm = TRUE),5, f = floor) + ylim.nudge[1],
                              round_any(max(prop2[complete.cases(cbind(prop1[,yy],prop2[,xx])),xx], na.rm = TRUE),5, f = ceiling) + ylim.nudge[2])
        } else if(axis.lims & !is.null(ylim)) {
          ylim.this.plot <- ylim
        } else {ylim.this.plot <- NULL}
       
        plot(prop1[,yy],
             prop2[,xx],
             xlim=xlim.this.plot, ylim=ylim.this.plot,
             ylab = "Predator Richness", xlab = "Herbivores Richness",
             col = col,
             pch = 16)
        points(prop1[,yy], prop2[,xx], pch = 21, cex = 1)
        if(!is.null(point.labels)) text(prop1[,yy], prop2[,xx], point.labels, cex = 0.75, adj = -0.25)
        
        prop.linear <- lm(prop2[,xx] ~ prop1[,yy])
        if(is.finite(prop.linear$coefficients[1]) & is.finite(prop.linear$coefficients[2])) abline(lm(prop2[,xx] ~ prop1[,yy]))
        
        if(xx == 3 & yy == 1) mtext("Predator Richness", side = 2, line = 4)
        if(xx == 1 & yy == 3) mtext("Herbivore Richness", side = 3, line = 2.5)
        if(yy == 1) mtext(colnames(prop2)[xx], side = 2, line =2)
        if(xx == 1) mtext(colnames(prop1)[yy], side = 3, line =0.5)
      }
    }
  }
}

#function to remove zeros
removeZeros <- function(dat.sub, dat.master, remove.zeros = c("all", "false.zero","leading","trailing"))
{
  if(any(remove.zeros %in% "false.zero"))
  {
    for(xx in seq_len(length(dat.sub)))
    {
      #remove false zeros
      if(dat.sub[xx] == 0 & dat.master[xx] > 0) dat.sub[xx] <- NA
    }
  }
  if(any(remove.zeros %in% "leading"))
  {
    run.while <- TRUE
    yy <- length(dat.sub)
    while(run.while)
    {
      if(!is.na(dat.master[yy]) & dat.master[yy] == 0)
      {
        dat.sub[yy] <- dat.master[yy] <- NA
        yy <- yy - 1
        if(yy == 0) run.while <- FALSE
      } else {
        run.while <- FALSE
      }
    }
  }
  if(any(remove.zeros %in% "trailing"))
  {
    run.while <- TRUE
    yy <- 1
    while(run.while)
    {
      if(!is.na(dat.master[yy]) & dat.master[yy] == 0)
      {
        dat.sub[yy] <- dat.master[yy] <- NA
        yy <- yy + 1
        if(yy == 0) run.while <- FALSE
      } else {
        run.while <- FALSE
      }
    }
  }
  if(any(remove.zeros %in% "all"))
  {
    dat.sub[dat.sub == 0] <- NA
    dat.master[dat.master == 0] <- NA
  } 
  
  return(cbind(dat.sub=dat.sub, dat.master=dat.master))
}

transpose.array <- function(countCube)
{
  countCube.flip <- array(NA, dim = c(ncol(countCube),nrow(countCube),dim(countCube)[3]))
  for(zz in seq_len(dim(countCube)[3]))
  {
    for(xx in seq_len(nrow(countCube[,,zz])))
    {
      for(yy in seq_len(ncol(countCube[,,zz])))
      {
        countCube.flip[yy,xx,zz] <- countCube[xx,yy,zz]
      }
    }
  }
  dimnames(countCube.flip) <- list(colnames(countCube), rownames(countCube), NULL)
  countCube.flip
}

make_mock_median_prop <- function(nrow = 30, ncol = 5)
{
  mock.dat <- matrix(nrow = nrow, ncol = ncol)
  mock.vec <- c(rep(x = 10, times = nrow/2), rep(x = 0, times = nrow/2))
  mock.dat[,1] <- mock.dat[,2] <- mock.vec
  mock.dat[,3] <- c(rep(x = 0, times = nrow/3), rep(x = 10, times = nrow/3), rep(x = 0, times = nrow/3))
  mock.dat[,4] <- mock.dat[,5] <- rev(mock.vec)
   
  return()
}

taxHandley <- function(repIntTaxa, 
                       measure.mat,
                      occs,
                      #shortFam, #should be list of families in analysis.  Need to determine if those taxa without family will be aggregated together or aggregated by higher level taxonomy (e.g. order or other clade designation)
                      bigList, #should be list of genus or species in analysis
                      intervals,
                      whichHandley = c("median", "allReps"),
                      extra.intvs = 0, 
                      do.parallel=FALSE, 
                      this.cores = NULL,
                      do.heuristic = FALSE,
                      do.save = FALSE,
                      run.update = 100, #runs needed to list how many runs completed to date
                      file.path = NULL,
                      filename = NULL)
{
  
  ####################################################################################################################################
  ### Handley analysis of taxonomic distributions
  if(any(whichHandley %in% "median"))
  {
    print("Beginning median taxonomic Handley analysis...")
    
    # bigList <- bigList[bigList$order %in% focal.order,]
    # shortFam <- sort(unique(bigList$family))
    
#    taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
    taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], unique(bigList$family)), nbins=length(unique(bigList$family))), simplify="array"), simplify="array")
    
#   taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(bin=as.factor(bigList$family[as.character(bigList$accepted_name) %in% x]), nbins=length(bigList$family)), simplify="array"), simplify="array")
    
    
    dimnames(taxCube) <- list(unique(bigList$family), rownames(intervals), NULL)
    #med.n <- median(sapply(repIntTaxa, function(x) length(unique(unlist(sapply(x, function(y) y))))))
    med.n <- nrow(measure.mat)
    optList_tax_median <- doHandleyTest(thisCounts=apply(taxCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel, this.cores = this.cores)	
    
    if(do.save)
    {
      save(optList_tax_median, taxCube, bigList, file=paste0(file.path, filename,"_optList_tax_median.Rdata"))
    }
    beep(3)
  
  } else if(any(whichHandley %in% "allReps")) {
    print("Beginning taxonomic Handley analysis for all reps...")
    optList_tax_allReps <- list()
    for (this.rep in seq_len(reps)) {
     # taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array") #duplicate code, can be removed if taxcube is made above
      this.n <- length(unique(unlist(sapply(repIntTaxa[[this.rep]], function(x) x))))
      optList_tax_allReps[[this.rep]] <- doHandleyTest(taxCube[,,this.rep], n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel, this.cores = this.cores)	
      if(this.rep %% run.update == 0) cat("Taxonomic Handley Rep:", this.rep, "\n")
    }
    
    if(do.save)
    {
      save(optList_tax_allReps, taxCube, bigList, file=paste0(file.path, filename,"_optList_bm_allReps.Rdata"))
    }
    beep(5)
  } else {
    print("Unable to perform. Please select whether the Hadley analysis will be performed over the median species richness (median) and/or all pseudoreplicates (allReps)")
  }
  
#  repIntTaxaAll <- list(OptListTaxMedian = optList_tax_median, OptListTaxAllReps = optList_tax_allReps, taxCube = taxCube)
 # return(repIntTaxaAll)
  return()
}

traitHandley <- function(countCube,
                         measure.mat,
                         whichHandley = c("median", "allReps"),
                         extra.intvs = 0, 
                         do.parallel=FALSE,
                         this.cores = NULL,
                         do.heuristic = FALSE,
                         do.save = FALSE,
                         run.update = 100,
                         file.path = NULL,
                         filename = NULL)
{
  ####################################################################################################################################
  ### Handley analysis of body mass distributions
  if(any(whichHandley %in% "median"))
  {
    print("Beginning median body mass Handley analysis...")
    
    countBox <- apply(countCube, c(1,2), quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)
    
    start.test <- Sys.time()
    optList_bm_median <- doHandleyTest(thisCounts = t(countBox[2,,]), n=nrow(measure.mat), do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel = do.parallel, this.cores = this.cores, 
                                       objectSize.filepath = paste0(file.path, "objectFileSize.txt"))
    end.test <- Sys.time()
    print(end.test - start.test)
    beep(3)
    ##Median Herbivores (subsampled; optimal is 7 breaks, so needs to run through 8 breaks)
    #runs in about 38.04854 minutes with 2 cores
    #runs in about 25.51563 minutes with 3 cores
    #runs in about 21.86931 minutes with 4 cores
    #runs in about 22.05333 minutes with 5 cores
    #runs in about 24.65575 minutes with 6 cores 
    
    if(do.save)
    {
      save(optList_bm_median, file=paste0(file.path, filename,"_optList_bm_median.Rdata"))
    }
    beep(3)
    
  } else if(any(whichHandley %in% "allReps")) {
    start.test <- Sys.time()
    print("Beginning body mass Handley analysis for all reps...")
    optList_bm_allReps <- list()
    for (this.rep in seq_len(dim(countCube)[3])) {
      this.n <- length(unique(unlist(sapply(repIntTaxa[[this.rep]], function(x) x))))
      optList_bm_allReps[[this.rep]] <- doHandleyTest(thisCounts = t(countCube[,,this.rep]), n=this.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel, this.cores = this.cores)
      if(this.rep %% run.update == 0) cat("Body Mass Handley Rep:",this.rep, "\n")
    }
    end.test <- Sys.time()
    print(end.test - start.test)
    
    ##Herbivores (___ subsampled reps; optimal is 7 breaks, so needs to run through 8 breaks)
    #runs in about ___ minutes with 2 cores
    #runs in about ___ minutes with 3 cores ~3 GB at start of 8 rate,
    #runs in about ___ minutes with 4 cores 
    #runs in about ___ minutes with 5 cores
    #runs in about ___ minutes with 6 cores ~1.5 GB per core at 8 rates, ~2.5 per core at 9 rates before killing program
    
    if(do.save)
    {
      save(optList_bm_allReps, file=paste0(file.path, filename,"_optList_bm_allReps.Rdata"))
    }
    beep(5)
  } else {
    print("Unable to perform. Please select whether the Hadley analysis will be performed over the median species richness (median) and/or all pseudoreplicates (allReps)")
    }
  
  #repIntTraitAll <- list(OptListTraitMedian = optList_bm_median, OptListTraitAllReps = optList_bm_allReps, countBox = countBox)
  #return(repIntTraitAll)
  return()
}

shoulderPlot <- function(measure.mat, plot.y, intervals, occs, bigList, repIntTaxa = NULL, optList_bm_median = NULL, quants = NULL, this.rank = settings$this.rank,
                         main = NULL, ylab = NULL, ylim = NULL, yaxp = NULL, yaxt = "s", xlab = NULL, xlim = NULL, xaxp = c(55,5,10), xaxt = "s", cex.axis = 1.5, cex.lab = 1.5, 
                         col.axis = "black", col.lab = "black",
                         specOcc.col = "gray0", specOcc.alpha = 0.5,
                         plot.breaks = FALSE, manual.breaks = NULL, break.text = TRUE, break.col = "firebrick4", breaks.lwd = 1.5, breaks.cex = 0.5, nudge.break.txt = 0,
                         do.subepochs=TRUE, overlay.labels = FALSE, overlay.color=TRUE, thisAlpha.intervals=0.33, thisAlpha.text = 0.33, borderCol="white", invertTime=FALSE, scale.cex=0.75, scale.headers = 0.95, text.offset = 0.025,
                         do.quants = FALSE, poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"))
{
  
  #for later versions make plot.x so user can define x axis but for now keep as $FO
  
  #check that measure.mat has FO and LO
  if(!"FO" %in% colnames(measure.mat) || !"LO" %in% colnames(measure.mat)) measure.mat <- getDateTaxa(measure.mat = measure.mat, occs = occs, this.rank = this.rank)
  
  #par(mar=c(4,4,2.5,0.5))
  # quartz(width=12, height=6)
  if(is.null(xlim)) xlim <- c(max(intervals), min(intervals))
  
  plot(measure.mat$FO, measure.mat[,plot.y], type="n", xlim= xlim, xaxp = xaxp, xaxt = xaxt, xlab = xlab, ylim=ylim, yaxp = yaxp, yaxt = yaxt, ylab = ylab, main = main, cex.axis = cex.axis, cex.lab = cex.lab, col.axis = col.axis, col.lab = col.lab)
  # plot(measure.mat$FO, measure.mat$bodyMass, xlim=c(max(intervals), min(intervals)), type="n", ylab="log-Body Mass (kg)", xaxp =c(50,0,5), xlab="Time (Ma)", cex.axis=1, cex.lab=1, col="gray75", fg="gray75", bg="gray75", col.axis="gray75", col.lab="gray75") #alter xaxpto change x-axis values
  # rect(-10e6, -10e6, 10e6, 10e6, col="white")
  overlayCzTimescale(do.subepochs= do.subepochs, color = overlay.color, thisAlpha.text = thisAlpha.text, thisAlpha.intervals = thisAlpha.intervals, borderCol = borderCol, invertTime = invertTime, scale.cex = scale.cex, scale.headers = scale.headers, text.offset = text.offset)
  
# famColors <- rainbow(length(bigList$family))
#  colorList <- famColors[match(bigList$family[as.character(bigList$accepted_name) %in% measure.mat$taxon], bigList$family)]
# colorList[is.na(colorList)] <- "gray25"
  
#  orderColors <- array(NA, dim=nrow(measure.mat))
  # orderColors[bigList$order[match(measure.mat$taxon, bigList$accepted_name)]=="Perissodactyla"] <- "dodgerblue4"
  # orderColors[bigList$order[match(measure.mat$taxon, bigList$accepted_name)] =="Artiodactyla"] <- "deeppink4"
  
  for (i in seq_len(nrow(measure.mat))) {
    # lines(x=c(this["FO"], x["LO"]), y=c(x["bodyMass"], x["bodyMass"]), lwd=3, pch=21, col=famColors[match(bigList[match(measure.mat$taxon, bigList[,1]),2], shortFam)])
    # lines(x=c(measure.mat$FO[i], measure.mat$LO[i]), y=c(measure.mat$bodyMass[i], measure.mat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor(colorList[i], 0.75))
    # lines(x=c(thisRanges[match(measure.mat$taxon[i], rownames(thisRanges)),"FO"], thisRanges[match(measure.mat$taxon[i], rownames(thisRanges)),"LO"]), y=c(measure.mat$bodyMass[i], measure.mat$bodyMass[i]), lwd=0.5, pch=21, col=alphaColor("gray0", 0.75)) #
    if (is.finite(measure.mat$FO[i]) & is.finite(measure.mat$LO[i]) & measure.mat$FO[i] != measure.mat$LO[i]) lines(x=measure.mat[i,c("FO","LO")], y=c(measure.mat[,plot.y][i], measure.mat[,plot.y][i]), lwd=0.75, pch=21, col=alphaColor(specOcc.col, specOcc.alpha)) #alphaColor(orderColors[i], 0.5)
  }
  points(measure.mat[complete.cases(measure.mat[ ,c("FO","LO")]) & measure.mat$FO == measure.mat$LO, c("FO", plot.y)], pch=21, col=alphaColor(specOcc.col, specOcc.alpha), cex=0.25) #this line is not generating the proper output for the final graph due to c("FO","bodyMass") causing a  "undefined columns selected" error
  
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
  
  if(plot.breaks)
  {
    if(is.null(manual.breaks))
    {
      # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), n=med.n, sig=0.01, do.heuristic=do.heuristic, extra.intvs=extra.intvs, do.parallel=do.parallel)	# based on means
      # optList_bm <- doHandleyTest(thisCounts=apply(countCube, c(1,2), median, na.rm=TRUE), sig=0.01, do.heuristic=TRUE, do.parallel=do.parallel)	# based on median
      abline(v=sort(c(intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2], range(intervals))), lwd= breaks.lwd, col = break.col) #ranges used to get start and stop of intervals since this method will find them as distinct changes (e.g, nothing to something and vice versa)
      if(break.text)
      {
        # text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3] + nudge.break.txt, labels=rev(seq_len(length(optList_bm_median[[length(optList_bm_median)-1]]$optBreaks) + 1)), pos=3, cex= breaks.cex, col=break.col)
        text(x= sort((c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2]) - 0.35)), y=par()$usr[3] + nudge.break.txt, labels= paste(sort(c(max(intervals), intervals[optList_bm_median[[length(optList_bm_median)-1]]$optBreaks,2])), "Ma"), adj=c(0,0),cex= breaks.cex, col=break.col)
      }
    }
    
    if(!is.null(manual.breaks))
    {
      abline(v= manual.breaks, lwd= breaks.lwd, col = break.col)
      if(break.text)
      {
        #  text(x= manual.breaks - 0.35, y=par()$usr[3] + nudge.break.txt, labels= manual.breaks, pos=3, cex= breaks.cex, col = break.col)
        text(x= manual.breaks - 0.35, y=par()$usr[3] + nudge.break.txt, labels= paste(manual.breaks, "Ma"), adj=c(0,0), cex= breaks.cex, col = break.col)
      }
    }
  }
}


plotRegimeBoxplot <- function(optList, prop, intervals, settings,
                              ylab = "Species Richness", ylim = NULL,
                              plot.legend = TRUE, legend.lwd = 5, legend.cex = 1, legend.seg.len = 1)
{
  require(stringr)
 
  regime.list <- getRegimeList(optList = optList, intervals = intervals)
  
  if(length(ylim) < 2 & ylim[1] %in% "rounded") ylim <- c(0, round(max(unlist(lapply(regime.list, function(x) max(prop[x,]))))/5)*5)
  
  par(mfrow = c(1,length(regime.list)),
      oma = c(3,3,0,0), mar = c(2,1,2,2))
  for(xx in seq_len(length(regime.list)))
  {
    boxplot(prop[regime.list[[xx]],], 
            ylab = ylab, ylim = ylim, xaxt = "n", 
            col = rainbow(ncol(prop)))
    if(xx == 1) 
    {
      mtext(paste0(">", intervals[optList[[length(optList)-1]]$optBreaks[xx],2], " Ma"), line = 0.75)
    } else if(xx == max(seq_len(length(regime.list))))
    {
      mtext(paste0("<", intervals[optList[[length(optList)-1]]$optBreaks[xx-1],2], " Ma"), line = 0.75)
    } else {
      mtext(paste0(intervals[optList[[length(optList)-1]]$optBreaks[xx-1],2],
                   " to ",
                   intervals[optList[[length(optList)-1]]$optBreaks[xx],2],
                   " Ma"), line = 0.75)
    }
    if(xx == 1 & !is.null(ylab)) 
    {
      mtext(paste0(" Median ", str_to_title(settings$this.rank)," Richness"), line = 2.5, side = 2)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom', y = 1, legend = colnames(prop), title = "Body Mass (kg)", col = rainbow(ncol(prop)), lwd = legend.lwd, xpd = TRUE, horiz = TRUE, cex = legend.cex, seg.len = legend.seg.len, bty = 'n')         
} 


plotRegimeBarplot <- function(optList, prop, intervals, settings,
                              ylab = "Species Richness", ylim = NULL,
                              plot.legend = TRUE, legend.lwd = 5, legend.cex = 1, legend.seg.len = 1)
{
  require(stringr)
  
  regime.list <- getRegimeList(optList = optList, intervals = intervals)
    
  if(length(ylim) < 2 & ylim[1] %in% "rounded") ylim <- c(0, round(max(unlist(lapply(regime.list, function(x) max(prop[x,]))))/5)*5)
  
  par(mfrow = c(1,length(regime.list)),
        oma = c(3,3,0,0), mar = c(2,1,2,2))
  for(xx in seq_len(length(regime.list)))
  {
    barplot(height =  apply(prop[regime.list[[xx]],], c(2), median, na.rm = TRUE),
            ylab = ylab, ylim = ylim, xaxt = "n", 
            col = rainbow(ncol(prop)))
    if(xx == 1) 
    {
      mtext(paste0(">", intervals[optList[[length(optList)-1]]$optBreaks[xx],2]," Ma"), line = 0.75)
    } else if(xx == max(seq_len(length(regime.list))))
    {
      mtext(paste0("<", intervals[optList[[length(optList)-1]]$optBreaks[xx-1],2], " Ma"), line = 0.75)
    } else {
    mtext(paste0(intervals[optList[[length(optList)-1]]$optBreaks[xx-1],2],
                 " to ",
                 intervals[optList[[length(optList)-1]]$optBreaks[xx],2],
                 " Ma"), line = 0.75)
    }
    if(xx == 1 & !is.null(ylab)) 
    {
      mtext(paste0(" Median ", str_to_title(settings$this.rank)," Richness"), line = 2.5, side = 2)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom', y = 1, legend = colnames(prop), title = "Body Mass (kg)", col = rainbow(ncol(prop)), lwd = legend.lwd, xpd = TRUE, horiz = TRUE, cex = legend.cex, seg.len = legend.seg.len, bty = 'n')         
} 
  

netChangeBarplot <- function(optList, prop, intervals, settings, 
                             ylab = "Species Richness", ylim = NULL,
                             plot.legend = TRUE, legend.lwd = 5, legend.cex = 1, legend.seg.len = 1)
{
  require(stringr)
  
  regime.list <- getRegimeList(optList = optList, intervals = intervals)
  
  if(length(ylim) < 2 & ylim[1] %in% "rounded")
  {
    temp.test <- list()                                                                                                         
    for(xx in seq_len(length(regime.list)-1))
    {
      temp.test[[xx]] <- apply(prop[regime.list[[xx+1]],], c(2), median, na.rm = TRUE) - apply(prop[regime.list[[xx]],], c(2), median, na.rm = TRUE)
    }
    ylim <- c(round(min(unlist(temp.test))/5)*5, round(max(unlist(temp.test)/5))*5)
  }
  
  par(mfrow = c(1,length(optList[[length(optList)-1]]$optBreaks)),
      oma = c(2,3,0,0), mar = c(2,1,2,2))
  for(xx in seq_len(length(regime.list)-1))
  {
    barplot(height =  apply(prop[regime.list[[xx+1]],], c(2), median, na.rm = TRUE) - apply(prop[regime.list[[xx]],], c(2), median, na.rm = TRUE),
            ylab = ylab, ylim = ylim, xaxt = "n", 
            col = rainbow(ncol(prop)))
    mtext(paste0(intervals[optList[[length(optList)-1]]$optBreaks[xx],2], " Ma"), line = 0.75)
    if(xx == 1 & is.null(ylab)) 
    {
      mtext(paste0("\u394", " Median ", str_to_title(settings$this.rank)," Richness"), line = 2.5, side = 2)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom', y = 1, legend = colnames(prop), title = "Body Mass (kg)", col = rainbow(ncol(prop)), lwd = legend.lwd, xpd = TRUE, horiz = TRUE, cex = legend.cex, seg.len = legend.seg.len, bty = 'n')         
}

getRegimeList <- function(optList, intervals)
{
  regime.list <- list()
  for(xx in seq_len(length(optList[[length(optList)-1]]$optBreaks)))
  {
    if(xx == 1) 
    {
      regime.list[[1]] <- seq(optList[[length(optList)-1]]$optBreaks[xx]+1, nrow(intervals),1)
    } else {
      regime.list[[xx]]<- seq(optList[[length(optList)-1]]$optBreaks[xx]+1, #its +1 so that the regime ends at the break e.g.(break at 41 Ma will have last bin be 41-43Ma with 42 Ma label)
                              optList[[length(optList)-1]]$optBreaks[xx-1],
                              1)
    }
    if(xx == length(optList[[length(optList)-1]]$optBreaks)) 
    {
      regime.list[[xx+1]]<- seq(1, optList[[length(optList)-1]]$optBreaks[xx],1)
    }
  }
  return(regime.list)
}

comparePBDB_NOW_collectionDurations <- function(occs)
{
  pbdb_NOW_map <- read.csv("~/Dropbox/Code/R/dentalMeasurements/dat/pbdb_NOW_map.csv")
  NOW_loc_dates <- read.csv("~/Dropbox/Code/R/dentalMeasurements/dat/NOW_loc_dates.csv")
  
  #remove PBDB rows that are empty
  pbdb_NOW_map <- pbdb_NOW_map[rowSums(is.na(pbdb_NOW_map)) != ncol(pbdb_NOW_map),]
  
  #remove rows without age values
  NOW_loc_dates <- NOW_loc_dates[!is.na(NOW_loc_dates$MAX_AGE) | !is.na(NOW_loc_dates$MIN_AGE),]

  #some of NOW_loc_dates are in the wrong date column
  #temp.max <- NOW_loc_dates[NOW_loc_dates$MAX_AGE < NOW_loc_dates$MIN_AGE & is.finite(NOW_loc_dates$MAX_AGE) & is.finite(NOW_loc_dates$MIN_AGE), "MAX_AGE"] 
  #temp.min <- NOW_loc_dates[NOW_loc_dates$MAX_AGE < NOW_loc_dates$MIN_AGE & is.finite(NOW_loc_dates$MAX_AGE) & is.finite(NOW_loc_dates$MIN_AGE), "MIN_AGE"]
  
  #NOW_loc_dates[NOW_loc_dates$MAX_AGE < NOW_loc_dates$MIN_AGE & is.finite(NOW_loc_dates$MAX_AGE) & is.finite(NOW_loc_dates$MIN_AGE), "MAX_AGE"]  <- temp.min
  #NOW_loc_dates[NOW_loc_dates$MAX_AGE < NOW_loc_dates$MIN_AGE & is.finite(NOW_loc_dates$MAX_AGE) & is.finite(NOW_loc_dates$MIN_AGE), "MIN_AGE"] <- temp.max
  
  NOW_loc_dates$MAX_AGE[NOW_loc_dates$MAX_AGE==64.81] <- 63.81	### fixes erroneous date of Pu3/To1 boundary in NOW databasae
  NOW_loc_dates$MIN_AGE[NOW_loc_dates$MIN_AGE==64.81] <- 63.81	### fixes erroneous date of Pu3/To1 boundary in NOW databasae

  NOW_loc_dates$diffNOW <- NOW_loc_dates$MAX_AGE - NOW_loc_dates$MIN_AGE

  compareDates <- cbind(pbdb_NOW_map[,c("collection_no", "NOW_loc","early_interval", "late_interval", "max_ma", "min_ma", "diff")],
                        NOW_loc_dates[match(pbdb_NOW_map$NOW_loc, NOW_loc_dates$LOC_SYNONYMS),c(c("BFA_MAX", "BFA_MIN", "MAX_AGE", "MIN_AGE","diffNOW", "CHRON"))])
  
  #empty rows are present after this step but are not getting removed
  #compareDates <- compareDates[rowSums(is.na(compareDates)) != ncol(compareDates),]
  #compareDates <- compareDates[is.na(compareDates$MAX_AGE) | is.na(compareDates$MIN_AGE),] 

  par(mfrow = c(2,2))
  hist(compareDates$diff, breaks = c(seq(0,25,1)), main = "PBDB Collections and Dates", ylim = c(0, 3500), xlim = c(0,25), xlab = "Collection Duration (Ma)")
  hist(compareDates$diffNOW, breaks = c(seq(0,25,1)), main = "PBDB Collections and NOW Dates", ylim = c(0,3500), xlim = c(0,25), xlab = "Collection Duration (Ma)")
  
  # how do occurrences change
  #par(mfrow = c(2,2))
  hist(occs$max_ma - occs$min_ma, breaks = c(seq(0,30,1)), main = "PBDB Occs and Dates", ylim = c(0, 15000), xlim = c(0,25), xlab = "Occurrence Duration (Ma)")
  
  occs$NOW_loc <- getNOWLocalityCodesFromPBDBCollectionNo(occs$collection_no)
  occs <- data.frame(occs, getNOWDatesMatFromLocCodes(occs$NOW_loc))
  occs$max_ma[is.finite(occs$MAX_AGE)] <- occs$MAX_AGE[is.finite(occs$MAX_AGE)]
  occs$min_ma[is.finite(occs$MIN_AGE)] <- occs$MIN_AGE[is.finite(occs$MIN_AGE)]
  
  hist(occs$max_ma - occs$min_ma, breaks = c(seq(0,30,1)), main = "PBDB Occs and NOW Dates", ylim = c(0, 25000), xlim = c(0,25), xlab = "Occurrence Duration (Ma)")
  
  return()
}

writeObjectSizeToFile <- function(nrates, object, object_name, file.path)
{
  #check if file has been created if not begin it
  write(paste0(quote(nrates), ":", nrates," ", object_name,":", format(object.size(object), units = "auto", standard = "SI")), file = file.path, append = TRUE)
  
  return()
}
  
  
  
  
  