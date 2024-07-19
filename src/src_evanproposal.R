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

getCorrelationPlot <- function(prop1, prop2, correlation.type = "spearman",mar = c(0,0,0,0), cl.pos = NULL, tl.pos = c("lt"), sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', number.digits = 3, p_val_correct_type = FALSE, p_val_correct_ntest = NULL)
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
  
  if(p_val_correct_type == "bonferroni") 
  {
    if(is.null(p_val_correct_ntest)) cor.p <- cor.p*(ncol(prop1)*ncol(prop2)) else cor.p <- cor.p*(p_val_correct_ntest)
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
                                round.x = 1, round.y = 1,
                                xlab = "Predator Richness", ylab = "Herbivore Richness",
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
  
  par(mfrow = c(par.mfrow.rows,par.mfrow.cols), mar = mar, oma = oma, xpd = FALSE)
  for(xx in for.xx)
  {
    for(yy in for.yy)
    {
      if(!any(complete.cases(cbind(prop1[,yy],prop2[,xx])))) 
      { 
        plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
      } else {
        if(axis.lims & is.null(xlim)) {
          xlim.this.plot <- c(round_any(min(prop1[complete.cases(cbind(prop1[,yy],prop2[,xx])),yy], na.rm = TRUE),round.x, f= floor) + xlim.nudge[1],
                              round_any(max(prop1[complete.cases(cbind(prop1[,yy],prop2[,xx])),yy], na.rm = TRUE),round.x, f= ceiling)+xlim.nudge[2]) 
        } else if(axis.lims & !is.null(xlim)) {
          xlim.this.plot <- xlim + xlim.nudge
        } else {xlim.this.plot <- NULL}
        
        if(axis.lims & is.null(ylim)) {
          ylim.this.plot <- c(round_any(min(prop2[complete.cases(cbind(prop1[,yy],prop2[,xx])),xx], na.rm = TRUE),round.y, f = floor) + ylim.nudge[1],
                              round_any(max(prop2[complete.cases(cbind(prop1[,yy],prop2[,xx])),xx], na.rm = TRUE),round.y, f = ceiling) + ylim.nudge[2])
        } else if(axis.lims & !is.null(ylim)) {
          ylim.this.plot <- ylim + ylim.nudge
        } else {ylim.this.plot <- NULL}
        
        plot(prop1[,yy],
             prop2[,xx],
             xlim=xlim.this.plot, ylim=ylim.this.plot,
             ylab = ylab, xlab = xlab,
             col = col,
             pch = 16)
        points(prop1[,yy], prop2[,xx], pch = 21, cex = 1)
        if(!is.null(point.labels)) text(prop1[,yy], prop2[,xx], point.labels, cex = 0.75, adj = -0.25)
        
        prop.linear <- lm(prop2[,xx] ~ prop1[,yy])
        if(is.finite(prop.linear$coefficients[1]) & is.finite(prop.linear$coefficients[2])) abline(lm(prop2[,xx] ~ prop1[,yy]))
        
        if(xx == 3 & yy == 1) mtext(ylab, side = 2, line = 4)
        if(xx == 1 & yy == 3) mtext(xlab, side = 3, line = 2.5)
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
        #countCube.flip[,,zz] <- as.data.frame(t(countCube[,,zz]))
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
    
#   taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], shortFam), nbins=length(shortFam)), simplify="array"), simplify="array")
    taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(match(bigList$family[as.character(bigList$accepted_name) %in% x], unique(bigList$family)), nbins=length(unique(bigList$family))), simplify="array"), simplify="array")
    
#   taxCube <- sapply(repIntTaxa, function(y) sapply(y, function(x) tabulate(bin=as.factor(bigList$family[as.character(bigList$accepted_name) %in% x]), nbins=length(bigList$family)), simplify="array"), simplify="array")
    
    if(dim(taxCube)[2] != nrow(intervals)) taxCube <- taxCube[,dimnames(taxCube)[[2]] %in% rownames(intervals),]
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
                         main = NULL, ylab = NULL, ylim = NULL, yaxp = NULL, yaxs = "i", yaxt = "s", xlab = NULL, xlim = NULL, xaxp = c(55,5,10), xaxs = "i", xaxt = "s", cex.axis = 1.5, cex.lab = 1.5,
                         col.axis = "black", col.lab = "black",
                         specOcc.col = "gray0", specOcc.alpha = 0.5,
                         plot.breaks = FALSE, manual.breaks = NULL, break.text = TRUE, break.col = "firebrick4", breaks.lwd = 1.5, breaks.cex = 0.5, nudge.break.txt = 0,
                         do.subepochs=TRUE, overlay.labels = FALSE, overlay.color=TRUE, thisAlpha.intervals=0.33, thisAlpha.text = 0.33, borderCol="black", invertTime=FALSE, scale.cex=0.75, scale.headers = 0.95, text.offset = 0.025,
                         do.quants = FALSE, poly.col = "darkorange4", median.col = c("goldenrod1", "darkorange4", "darkorange1"), poly.alpha = 0.25, median.alpha = 0.5)
{
  
  #for later versions make plot.x so user can define x axis but for now keep as $FO
  
  #check that measure.mat has FO and LO
  if(!"FO" %in% colnames(measure.mat) || !"LO" %in% colnames(measure.mat)) measure.mat <- getDateTaxa(measure.mat = measure.mat, occs = occs, this.rank = this.rank)
  
  #par(mar=c(4,4,2.5,0.5))
  # quartz(width=12, height=6)
  if(is.null(xlim)) xlim <- c(max(intervals), min(intervals))
  
  plot(measure.mat$FO, measure.mat[,plot.y], type="n", xlim= xlim, xaxp = xaxp, xaxs = xaxs, xaxt = xaxt, xlab = xlab, ylim=ylim, yaxp = yaxp, yaxs = yaxs, yaxt = yaxt, ylab = ylab, main = main, cex.axis = cex.axis, cex.lab = cex.lab, col.axis = col.axis, col.lab = col.lab)
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
    polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[1,], rev(quants[5,])), col=alphaColor(poly.col, poly.alpha), border=poly.col)
    polygon(c(rowMeans(intervals), rev(rowMeans(intervals))), c(quants[2,], rev(quants[4,])), col=alphaColor(poly.col, poly.alpha), border=poly.col)
    lines(rowMeans(intervals), quants[3,], col=alphaColor(median.col[1], median.alpha), lwd=5)
    lines(rowMeans(intervals), quants[3,], col=alphaColor(median.col[2], 1.0), lwd=3)
    points(rowMeans(intervals), quants[3,], col=alphaColor(median.col[3], median.alpha), cex=0.5)
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
                              use_par = TRUE, regime.lab.cex = 1,
                              ylab = "Species Richness", ylab.line = 2.5, ylab.cex = 1,
                              ylim = NULL, yaxt = "n", 
                              axis.at = NULL, cex.axis = 1,
                              text.col = "black",
                              axis.col.ticks="black", axis.col.axis = "black",
                              plot.legend = TRUE, legend.lwd = 5, legend.cex = 1, legend.seg.len = 1,
                              legend.text.col = "black")
{
  require(stringr)
 
  regime.list <- getRegimeList(optList = optList, intervals = intervals)
  
  if(length(ylim) < 2 & ylim[1] %in% "rounded") ylim <- c(0, round(max(unlist(lapply(regime.list, function(x) max(prop[x,]))))/5)*5)
  
  if(!is.null(use_par))
  {
    par(mfrow = c(1,length(regime.list)),
      oma = c(3,3,0,0), mar = c(2,1,2,2))
  }
  for(xx in seq_len(length(regime.list)))
  {
      boxplot(prop[regime.list[[xx]],], 
              ylab = ylab, ylim = ylim, xaxt = "n", yaxt = yaxt, cex.axis = cex.axis,
              col = rainbow(ncol(prop)))
     if(yaxt == "n" & xx == 1) axis(side = 2, at = axis.at, labels = TRUE, col.ticks= axis.col.ticks, col.axis = axis.col.axis, cex.axis = cex.axis)
     if(yaxt == "n" & xx != 1) axis(side = 2, at = axis.at, labels = FALSE, col.ticks= axis.col.ticks, col.axis = axis.col.axis, cex.axis = cex.axis)
    
    if(xx == 1) 
    {
      mtext(paste0(">", intervals[optList[[length(optList)-1]]$optBreaks[xx],2], " Ma"), line = 0.75, cex = regime.lab.cex, col = text.col)
    } else if(xx == max(seq_len(length(regime.list))))
    {
      mtext(paste0("<", intervals[optList[[length(optList)-1]]$optBreaks[xx-1],2], " Ma"), line = 0.75, cex = regime.lab.cex, col = text.col)
    } else {
      mtext(paste0(intervals[optList[[length(optList)-1]]$optBreaks[xx-1],2],
                   " to ",
                   intervals[optList[[length(optList)-1]]$optBreaks[xx],2],
                   " Ma"), line = 0.75, cex = regime.lab.cex, col = text.col)
    }
    if(xx == 1 & !is.null(ylab) & yaxt == "n") 
    {
      mtext(text = ylab, line = ylab.line, cex = ylab.cex, side = 2, col = text.col)
    }
  }
  if(plot.legend)
  {
    par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottom', y = 1, legend = colnames(prop), title = "Body Mass (kg)", col = rainbow(ncol(prop)), lwd = legend.lwd, xpd = TRUE, horiz = TRUE, cex = legend.cex, seg.len = legend.seg.len, bty = 'n', text.col = legend.text.col)         
  }
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
            cex.axis = cex.axis, 
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
  if(plot.legend)
  {
    par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottom', y = 1, legend = colnames(prop), title = "Body Mass (kg)", col = rainbow(ncol(prop)), lwd = legend.lwd, xpd = TRUE, horiz = TRUE, cex = legend.cex, seg.len = legend.seg.len, bty = 'n')         
  }
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
  if(plot.legend)
  {
    par(fig = c(0, 1, 0, 1), oma = c(1, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend('bottom', y = 1, legend = colnames(prop), title = "Body Mass (kg)", col = rainbow(ncol(prop)), lwd = legend.lwd, xpd = TRUE, horiz = TRUE, cex = legend.cex, seg.len = legend.seg.len, bty = 'n')         
  }
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

  NOW_loc_dates$MAX_AGE[NOW_loc_dates$MAX_AGE==64.1] <- 65.83	### makes any collection from Pu3 from combined Pu2/Pu3. as of 1/12/2024 the bound is set to 65.84 
  NOW_loc_dates$MIN_AGE[NOW_loc_dates$MIN_AGE==64.1] <- 64.81	### makes any collection from Pu2 from combined Pu2/Pu3. as of 1/12/2024 the bound is set to 64.75 

  NOW_loc_dates$diffNOW <- NOW_loc_dates$MAX_AGE - NOW_loc_dates$MIN_AGE

  compareDates <- cbind(pbdb_NOW_map[,c("collection_no", "NOW_loc","early_interval", "late_interval", "max_ma", "min_ma", "diff")],
                        NOW_loc_dates[match(pbdb_NOW_map$NOW_loc, NOW_loc_dates$LOC_SYNONYMS),c(c("BFA_MAX", "BFA_MIN", "MAX_AGE", "MIN_AGE","diffNOW", "CHRON"))])
  
  #empty rows are present after this step but are not getting removed
  #compareDates <- compareDates[rowSums(is.na(compareDates)) != ncol(compareDates),]
  #compareDates <- compareDates[is.na(compareDates$MAX_AGE) | is.na(compareDates$MIN_AGE),] 

  png(paste0(save.folder, "PBDB versus NOW durations histograms ", settings$this.rank,".png"),
      width = 8, height = 4, units = "in", res = 200)
  par(mfrow = c(1,2), mar = c(2,2,1,2), oma = c(1,1,1,1))
#  hist(compareDates$diff, breaks = c(seq(0,25,1)), main = "PBDB Collections and Dates", ylim = c(0, 3500), xlim = c(0,25), xlab = "Collection Duration (Ma)")
  hist(compareDates$diff, breaks = c(seq(0,25,1)), main = "", ylim = c(0, 3500), xlim = c(0,25), xlab = "Collection Duration (Ma)")
  
  mtext("Duration (Myr)", side = 1, line = 2, cex = 1)
  mtext("Frequency", side = 2, line = 2, cex = 1)
  mtext(text= "A", side = 3, line = 0.2, cex = 2, adj = 0)
#  hist(compareDates$diffNOW, breaks = c(seq(0,25,1)), main = "PBDB Collections and NOW Dates", ylim = c(0,3500), xlim = c(0,25), xlab = "Collection Duration (Ma)")
  hist(compareDates$diffNOW, breaks = c(seq(0,25,1)), main = "", ylim = c(0,3500), xlim = c(0,25), xlab = "Collection Duration (Ma)")
  mtext("Duration (Myr)", side = 1, line = 2, cex = 1)
  mtext(text= "B", side = 3, line = 0.2, cex = 2, adj = 0)
  
  dev.off()
  
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
  
generatePyrateAnalysisCode <- function(filepathname, python.version = 3,
                                       iter.num = 10000000, sample.freq = 1000,
                                       main.algorithm = "RJMCMC", #""
                                       preserv.model = "NHPP", epoch.boundary.file, TPPgamma.dist = c(1.5,1.5),
                                       Gamma.Model = FALSE, save.per.lin.rates = FALSE,
                                       check.rep = NULL)
  
{
  pythonscript <- paste0("python", python.version," PyRate.py ", filepathname)
  
  if(main.algorithm == "BDMCMC") {main.al.out <- " -A 2"} else {main.al.out <- ""}
  pythonscript <- paste0(pythonscript, main.al.out)
  
  pythonscript <- paste0(pythonscript, " -n ", format(iter.num, scientific = F), " -s ", format(sample.freq, scientific = F))
  
  if(preserv.model == "HPP") {pres.mod.out <- " -mHPP"} else if(preserv.model == "TPP") {pres.mod.out <- paste0(" -qShift ",epoch.boundary.file, " -pP ", TPPgamma.dist[1], " ", TPPgamma.dist[2])} else {pres.mod.out <- ""}
  pythonscript <- paste0(pythonscript, pres.mod.out)
  
  if(Gamma.Model) pythonscript <- paste0(pythonscript, " -mG")
  
  if(save.per.lin.rates) pythonscript <- paste0(pythonscript, " -log_sp_q_rates")
  
  if(!is.null(check.rep)) pythonscript <- paste0(pythonscript, " -j ", check.rep)
  
  return(pythonscript)
}

generatePyratePlotCode <- function(filepathname, python.version = 3,
                                   iter.num = 10000000, sample.freq = 1000,
                                   main.algorithm = "RJMCMC", #""
                                   preserv.model = "NHPP", epoch.boundary.file, TPPgamma.dist = c(1.5,1.5),
                                   Gamma.Model = FALSE, save.per.lin.rates = FALSE,
                                   check.rep = NULL)
  
{
  pythonscript <- paste0("python", python.version," PyRate.py ", filepathname)
  
  #work in progress
  
  return(pythonscript)
}
  

test_func <- function(num_plots = 4, plot_id_start = 1, ncol = 11, plot_alignment = c("left"), widen_multiplier = 1)
{
  layout_seq <- seq(plot_id_start, plot_id_start + (num_plots-1),1)
  if(widen_multiplier > 1) 
  {
    temp_seq <- integer()
    for(xx in layout_seq)
    {
      temp_seq <- c(temp_seq, rep(xx, widen_multiplier))
    }
    layout_seq <- unlist(temp_seq)
  }
  if(length(layout_seq) > ncol) {print(paste("There are more plots than columns.  Please increase the number of columns"))
  } else if(length(layout_seq) < ncol){
      seq_diff <- ncol - length(layout_seq)
      if(plot_alignment == "left")   {layout_seq_out <- c(layout_seq, rep(0, seq_diff))}
      if(plot_alignment == "right")  {layout_seq_out <- c(rep(0, seq_diff), layout_seq)}
      if(plot_alignment == "middle") {
        if((seq_diff %% 2) == 0) {layout_seq_out <- c(rep(0, seq_diff/2), layout_seq, rep(0, seq_diff/2))}
        if((seq_diff %% 2) != 0) {
          seq_front <- floor(seq_diff/2)
          seq_back  <- ceiling(seq_diff/2)
          if(seq_diff != seq_front + seq_back) {print(paste("The function is having trouble aligning the plots.  Please use left or right"))}
          layout_seq_out <- c(rep(0, floor(seq_diff/2)), layout_seq, rep(0, ceiling(seq_diff/2)))
        }
      } 
    } else {layout_seq_out <- layout_seq}
  
  return(layout_seq_out)
}  


plotCzEnviroEvents <- function(pch.y.pos = 1.075, pch = 25, pch.cex = 2.5, text.cex = 1, lwd = 1, pch.bg = NA)
{
  #climatic event dates taken form Zachos etal 2008 and Westerhold etal 2020
  #Vegetative transitions taken from Stromberg etal 2011, Stromberg and McInerney 2011, and Townsend et al. 2010. 
  par(xpd = TRUE)
  #PETM
  points(56, pch.y.pos*par()$usr[4], pch = pch, cex = pch.cex, col = "red", bg = pch.bg)
  text(x = 56, y = pch.y.pos*par()$usr[4], labels = "1", cex = text.cex, col = "red")
  #EECO
  points(c(53,49), rep(pch.y.pos*par()$usr[4],2), pch = pch, cex = pch.cex, col = "red", bg = pch.bg)
  text(c(53,49), rep(pch.y.pos*par()$usr[4],2), labels = c("2","2"), cex = text.cex, col = "red")
  lines( x =c(52.5,49.5), y = rep(pch.y.pos*par()$usr[4],2), col = "red", lwd = lwd)
  #Earliest Inferred Open Habitat (Townsend et al. 2010)
  points(42, pch.y.pos*par()$usr[4], pch = pch, cex = pch.cex, bg=pch.bg)
  text(x = 42, y = pch.y.pos*par()$usr[4], labels = "3", cex = text.cex)
  #Middle Eocene Climatic Optimum
  points(40.25, pch.y.pos*par()$usr[4], pch = pch, cex = pch.cex, col = "red", bg = pch.bg) #average of 40.5-40.1
  text(x = 40.25, y = pch.y.pos*par()$usr[4], labels = "4", cex = text.cex, col = "red")
  #E/O Boundary
  points(33.9, pch.y.pos*par()$usr[4], pch = pch, cex = pch.cex, col = "blue", bg = pch.bg)
  text(x = 33.9, y = pch.y.pos*par()$usr[4], labels = "5", cex = text.cex, col = "blue")
  #Onset of Widespread Grass-Dominated Habitats (Stromberg etal 2011)
  points(26, pch.y.pos*par()$usr[4], pch = pch, cex = pch.cex, bg = pch.bg)
  text(x = 26, y = pch.y.pos*par()$usr[4], labels = "6", cex = text.cex)
  #lines( x =c(26,6), y = rep(pch.y.pos*par()$usr[4],2), lwd = lwd)
 # line.alpha <- rev(seq( 1-(13*0.05),1,0.05))
  #line.alpha <- alphaColor("black", line.alpha)
  #segments(x0 = rev(seq(19.5, 25.5, 0.5)), y0 = pch.y.pos*par()$usr[4],
   #        x1 = rev(seq(19, 25, 0.5)), y1 = pch.y.pos*par()$usr[4],
    #       col = line.alpha)
  arrows(x0 = 25.5, y0 = pch.y.pos*par()$usr[4],
         x1 = 19, y1 = pch.y.pos*par()$usr[4],
         length = 0.1, angle = 25)
  #MMCO
  points(c(17,14), rep(pch.y.pos*par()$usr[4],2), pch = pch, cex = pch.cex, col = "red", bg = pch.bg)
  text(c(17,14), rep(pch.y.pos*par()$usr[4],2), labels = c("7","7"), cex = text.cex, col = "red")
  lines( x =c(16.5,14.5), y = rep(pch.y.pos*par()$usr[4],2), col = "red", lwd = lwd)
  #C3/C4 transition (Stromberg and McInerney 2011)
  points(c(8,5.5), rep(pch.y.pos*par()$usr[4],2), pch = pch, cex = pch.cex, bg = pch.bg)
  text(c(8,5.5), rep(pch.y.pos*par()$usr[4],2), labels = c("8","8"), cex = text.cex)
  lines( x =c(7.5,6), y = rep(pch.y.pos*par()$usr[4],2), lwd = lwd)
  par(xpd = FALSE)
}

flatten_occs <- function(strat_uncert, intervals)
{
  PBDB_flat <- as.data.frame(matrix(nrow = nrow(strat_uncert), ncol = nrow(intervals)))
  rownames(PBDB_flat) <- rownames(strat_uncert)
  colnames(PBDB_flat) <- rownames(intervals)
  
  for(yy in rownames(PBDB_flat))
  {
    #need to prevent taxa that exist outside our binning interval from causing errors
    if(!strat_uncert[rownames(strat_uncert) %in% yy,"FO"] <= min(intervals) | strat_uncert[rownames(strat_uncert) %in% yy,"LO"] >= max(intervals))
    {
      #need to handle with an FO before our first interval, an LO after our last interval, and normal taxa with range contained by our bins
      if(strat_uncert[rownames(strat_uncert) %in% yy,"FO"] >= max(intervals))
      { 
        taxon_intervals <- rownames(intervals)[intervals$ageBop >= strat_uncert[rownames(strat_uncert) %in% yy,"LO"]]
      } else if (strat_uncert[rownames(strat_uncert) %in% yy,"LO"] <= min(intervals)) { taxon_intervals <- rownames(intervals)[intervals$ageBop >= strat_uncert[rownames(strat_uncert) %in% yy,"FO"]]
      } else {taxon_intervals <- rownames(intervals)[intervals$ageBase >= strat_uncert[rownames(strat_uncert) %in% yy,"LO"] & intervals$ageTop <= strat_uncert[rownames(strat_uncert) %in% yy,"FO"]]
      }
      
      PBDB_flat[rownames(PBDB_flat) %in% yy, colnames(PBDB_flat) %in% taxon_intervals] <- yy
    }
  }
  
  #Make each column a list so that I can use our functions to get the distribution of body mass/families
  PBDB_flat_list <- list()
  
  for(xx in seq_len(ncol(PBDB_flat)))
  {
    PBDB_flat_list[[xx]] <- unique(PBDB_flat[,xx])[!is.na(unique(PBDB_flat[,xx]))]
  }
  names(PBDB_flat_list) <- colnames(PBDB_flat)
  
  return(PBDB_flat_list)
}

getBigList <- function(focal.tax = NULL, 
                       #rank.vec = c("subspecies","species"), 
                       OccsCols = c("class","order","family","genus", "accepted_name"), 
                       seperateNoFamilyTaxa = FALSE)
{
  uniqTax <- lapply(unlist(focal.tax), FUN=getTaxonomyForOneBaseTaxon_AcceptedName)
  uniqTax <- makeMatrixFromList(uniqTax)
  uniqTax <- uniqTax[!uniqTax$class %in% "Insecta",] #remove entries for Insecta due to Proboscidea being a synonym for Hemiptera
  uniqTax <- uniqTax[!uniqTax$family %in% "",] #remove entries with empty family assignments.  generally entries reported at the family level or higher
  uniqTax$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = uniqTax$accepted_name)
  uniqTax <- unique(uniqTax[,OccsCols])
 
  #uniqTax <- uniqTax[uniqTax$accepted_name %in% occs$accepted_name[occs$accepted_rank %in% rank.vec],] #remove taxa without occurrences
  
  bigList_herb <- uniqTax[order(uniqTax$family),] #4 occurrences for Pleistocene Suidae are in occs.  Likely just modern specimens or mistaken Tayassuids.
  
  if(seperateNoFamilyTaxa)
  {
    no_fam_Taxa <- bigList_herb$accepted_name[bigList_herb$family %in% "NO_FAMILY_SPECIFIED" & !bigList_herb$order %in% "NO_ORDER_SPECIFIED"]
    bigList_herb$family[bigList_herb$accepted_name %in% no_fam_Taxa] <- paste(bigList_herb$order[bigList_herb$accepted_name %in% no_fam_Taxa], bigList_herb$family[bigList_herb$accepted_name %in% no_fam_Taxa], sep = "_")
    
    #some runs may not contain taxa that lack a order designation (e.g., run only Perissodactyls) so need to check to avoid errors
    if(length(bigList_herb$accepted_name[bigList_herb$family %in% "NO_FAMILY_SPECIFIED" & bigList_herb$order %in% "NO_ORDER_SPECIFIED"])>0)
    {
      no_ord_no_fam_Taxa <- bigList_herb$accepted_name[bigList_herb$family %in% "NO_FAMILY_SPECIFIED" & bigList_herb$order %in% "NO_ORDER_SPECIFIED"]
      bigList_herb$family[bigList_herb$accepted_name %in% no_ord_no_fam_Taxa] <- paste(bigList_herb$class[bigList_herb$accepted_name %in% no_ord_no_fam_Taxa], bigList_herb$family[bigList_herb$accepted_name %in% no_ord_no_fam_Taxa], sep = "_")
    }
  }
  return(bigList_herb)
}

getIntMeasure.mat <- function(measure.mat, 
                              #settings, 
                              breaks, 
                              uniqIntTaxa)
{
  if(settings$this.rank != "species"){
    measure.mat_Int <- measure.mat[measure.mat[,settings$this.rank] %in% uniqIntTaxa,] # need to make sure that there is an argument to restrict collection to those in this interval and not all of them
    
    intSizeCat <- as.data.frame(matrix(data=NA, nrow = length(unique(measure.mat_Int[,settings$this.rank])), ncol = 2, 
                                       dimnames = list(unique(measure.mat_Int[,settings$this.rank]), c("bodyMass","SizeCat"))))
  } else {measure.mat_Int <- measure.mat[measure.mat[,"taxon"] %in% uniqIntTaxa,]
  intSizeCat <- as.data.frame(matrix(data=NA, nrow = length(unique(measure.mat_Int[,"accepted_name"])), ncol = 2, 
                                     dimnames = list(unique(measure.mat_Int[,"accepted_name"]), c("bodyMass","SizeCat"))))
  }
  
  for(mm in unique(measure.mat_Int[,settings$this.rank]))
  {
    intSizeCat[mm,"bodyMass"] <- mean(measure.mat_Int$bodyMass[measure.mat_Int[,settings$this.rank] %in% mm])
  }
  
  ##bin into size categories
  for(nn in seq(1, length(settings$bmBreaks_herb)-1, 1)){
    intSizeCat$SizeCat[intSizeCat$bodyMass >= settings$bmBreaks_herb[nn] & intSizeCat$bodyMass < settings$bmBreaks_herb[nn+1]] <- nn
  } 
  
  #set Proboscideans to sizeCat 5 manually since they lack body mass measures
  for(zz in unique(measure.mat_Int[,settings$this.rank]))
  {
    if(zz %in% unique(occs$genus[occs$order %in% "Proboscidea"])) {intSizeCat[zz,"SizeCat"] <- 5 
    } else {intSizeCat[zz,"bodyMass"] <- mean(measure.mat_Int$bodyMass[measure.mat_Int[,settings$this.rank] %in% zz])}
  }
  
  return(intSizeCat)
}


testFactorial4OccupancyProb <- function()
{
  n <- 100
  y <- 50
  p <- 0.58
  
  factSeries_n <- seq_len(n)
  factSeries_y <- seq_len(y)
  factSeries_diffny <- seq_len(n-y)
  
  probFunc <- prod(factSeries_n[!factSeries_n %in% factSeries_diffny])/prod(factSeries_y)
  probFunc*(p^y)*((1-p)^(n-y))
  fact(170)
}

checkifGenusNamesValid <- function() #for checking if a genus is valid.  Not all genera have their own entry in accepted_names.  this checks the genus in occs to see if valid species are present
{
  genus.vec <- c("Navahoceros", "Ottoceros", "Konobelodon", "Yumaceras", "Heteropliohippus", "Colbertchoerus", "Procranioceras", "Acritohippus",
                  "Drepanomeryx", "Mediochoerus", "Stirtonhyus", "Tedfordhyus", "Lucashyus", "Mediochoerus", "Nexuotapirus", "Paramerychyus", 
                  "Problastomeryx", "Pseudoblastomeryx", "Stuckyhyus", "Wrightohyus", "Desmatochoerus", "Floridaceras", "Nexuotapirus", 
                  "Paramerychyus", "Pseudoblastomeryx")
  
  genus.vec[!genus.vec %in% occs$genus] #check to see if all 
  
  sp.per.gen <- vector()
  for(xx in seq_len(length(genus.vec)))
  {
    sp.per.gen[xx] <- length(unique(occs$accepted_name[occs$genus %in% genus.vec[xx]]))
  }
  
  cbind(genus.vec, sp.per.gen)
  
}

getOccupancyOneInterval_evan <- function(this.intv, settings, intervals)
{
  ####################################################################
  #setup a data.frame to be filled in with presence/absence values
  if(!nrow(this.intv) == 0)
  {
    uniqSiteVec <- unique(this.intv[,settings$occSiteName])
    if(settings$occSiteName != "collection_no") 
    {
      presenceMat <- as.data.frame(matrix(data=NA, nrow = length(uniqSiteVec), ncol = max(table(this.intv[,settings$occSiteName])), dimnames = list(sort(uniqSiteVec),c(paste0("Sample_",seq_len(max(table(this.intv[,settings$occSiteName]))),sep="")))))
    } else {
      presenceMat <- as.data.frame(matrix(data=NA, nrow = length(uniqSiteVec), ncol = 1, dimnames = list(sort(uniqSiteVec),c("Sample_1"))))
    }
    
    #make data.frame to hold the occupancy values for individual taxa within interval
    if(settings$this.rank %in% "species")
    {
      uniqIntTaxa <- sort(unique(occs$accepted_name[occs$collection_no %in% this.intv$collection_no]))
    } else {
      uniqIntTaxa <- sort(unique(occs[occs$collection_no %in% this.intv$collection_no, settings$this.rank]))
    } 
    
    ####################################################################
    #restrict occupancy to just herbivore taxa to cut down on run time
    uniqIntTaxa <- uniqIntTaxa[uniqIntTaxa %in% measure.mat[,settings$this.rank]]
    occupancyMatTaxa <- as.data.frame(matrix(data = NA, nrow = length(uniqIntTaxa), ncol = 3, 
                                             dimnames = list(uniqIntTaxa,  c("NaiveOcc","NumSitesDetected", "NumCollsDetected"))))
    
    #get measure.mat for this interval and append bodymass and sizeCat to occupancyMatTaxa.  This will require the aggregation of con generics and reassignment of Proboscideans.
    intSizeCat <- getIntMeasure.mat(measure.mat = measure.mat, settings = settings, breaks = settings$bmBreaks_herb, uniqIntTaxa = uniqIntTaxa)
    occupancyMatTaxa <- cbind(occupancyMatTaxa, intSizeCat)
    
    ####################################################################
    #make data.frame to hold the occupancy values for body mass categories within interval
    size.categ <- c("Chevrotain-size (<5 kg)", "Javelina-size (5-25 kg)", "Antelope-size (25-150 kg)", "Horse-size (150-500 kg)", "Rhino-size (>500 kg)")
    occupancyMatBM <- as.data.frame(matrix(data = NA, nrow = length(size.categ), ncol = 3, 
                                           dimnames = list(size.categ, c("NaiveOcc","NumSitesDetected", "NumCollsDetected"))))
    
    ####################################################################
    #make data.frame to hold the occupancy values for body mass categories within interval
    occupancyMatAllHerb <- as.data.frame(matrix(data = NA, nrow = 1, ncol = 3, 
                                                dimnames = list("AllHerb", c("NaiveOcc","NumSitesDetected", "NumCollsDetected"))))
    
    ####################################################################
    #Get occupancy for individual taxa
    ####################################################################
    for(zz in seq_len(nrow(occupancyMatTaxa)))
    {
      presenceMat.temp <- presenceMat
      for(pp in seq_len(nrow(presenceMat.temp)))
      {
        #get collection within site
        collAtSite <- this.intv[this.intv[,settings$occSiteName] %in% rownames(presenceMat.temp)[pp], "collection_no"]
        
        for(mm in seq_len(length(collAtSite)))
        {
          if(rownames(occupancyMatTaxa)[zz] %in% occs[occs$collection_no %in% collAtSite[mm], settings$this.rank])
          {
            presenceMat.temp[pp,mm] <- 1
          } else {presenceMat.temp[pp,mm] <- 0}
        }
      }
      
      ###################################################
      #Naive Occupancy
      occupancyMatTaxa[zz, "NaiveOcc"] <- sum(apply(presenceMat.temp, 1 , max, na.rm = TRUE))/nrow(presenceMat.temp)
      
      occupancyMatTaxa$NaiveOcc[occupancyMatTaxa$NaiveOcc %in% 0] <- NA
      
      occupancyMatTaxa[zz,]$NumSitesDetected <- sum(apply(presenceMat.temp, 1 , max, na.rm = TRUE))
      occupancyMatTaxa[zz,]$NumCollsDetected <- sum(presenceMat.temp, na.rm = TRUE)
      occupancyMatTaxa
    }
    ####################################################################
    #Get occupancy for body mass categories
    ####################################################################
    for(zz in seq_len(nrow(occupancyMatBM)))
    {
      presenceMat.temp <- presenceMat
      for(pp in seq_len(nrow(presenceMat.temp)))
      {
        #get collection within site
        collAtSite <- this.intv[this.intv[,settings$occSiteName] %in% rownames(presenceMat.temp)[pp], "collection_no"]
        
        for(mm in seq_len(length(collAtSite)))
        {
          if(zz %in% measure.mat$SizeCat[measure.mat[,settings$this.rank] %in% occs[occs$collection_no %in% collAtSite[mm], settings$this.rank]])
          {
            presenceMat.temp[pp,mm] <- 1
          } else {presenceMat.temp[pp,mm] <- 0}
        }
      }
      
      ###################################################
      #Naive Occupancy
      occupancyMatBM[zz, "NaiveOcc"] <- sum(apply(presenceMat.temp, 1 , max, na.rm = TRUE))/nrow(presenceMat.temp)
      
      occupancyMatBM$NaiveOcc[occupancyMatBM$NaiveOcc %in% 0] <- NA
      
      occupancyMatBM[zz,]$NumSitesDetected <- sum(apply(presenceMat.temp, 1 , max, na.rm = TRUE))
      occupancyMatBM[zz,]$NumCollsDetected <- sum(presenceMat.temp, na.rm = TRUE)
    }
    ####################################################################
    #Get occupancy for herbviores as a whole
    ####################################################################
    presenceMat.temp <- presenceMat
    for(pp in seq_len(nrow(presenceMat.temp)))
    {
      #get collection within site
      collAtSite <- this.intv[this.intv[,settings$occSiteName] %in% rownames(presenceMat.temp)[pp], "collection_no"]
      
      for(mm in seq_len(length(collAtSite)))
      {
        if(any(rownames(occupancyMatTaxa) %in% occs[occs$collection_no %in% collAtSite[mm], settings$this.rank]))
        {
          presenceMat.temp[pp,mm] <- 1
        } else {presenceMat.temp[pp,mm] <- 0}
      }
    }
    ###################################################
    #Naive Occupancy
    occupancyMatAllHerb[1,"NaiveOcc"] <- sum(apply(presenceMat.temp, 1 , max, na.rm = TRUE))/nrow(presenceMat.temp)
    
    occupancyMatAllHerb$NaiveOcc[occupancyMatAllHerb$NaiveOcc %in% 0] <- NA
    
    occupancyMatAllHerb[1,]$NumSitesDetected <- sum(apply(presenceMat.temp, 1 , max, na.rm = TRUE))
    occupancyMatAllHerb[1,]$NumCollsDetected <- sum(presenceMat.temp, na.rm = TRUE)
    
    ###################################################
    #output
    
    this.intv.OccProb <- list(OccupancyTaxa = occupancyMatTaxa, OccupancyBM = occupancyMatBM, OccupancyAllHerb = occupancyMatAllHerb) #, Sites=repIntColl[[xx]])
  } else {this.intv.OccProb <- list(OccupancyTaxa = NULL, OccupancyBM = NULL, OccupancyAllHerb = NULL)} #, Sites=NULL)}
  this.intv.OccProb
}

getOccupancyOneRep_evan <- function(this.rep, settings, intervals)
{
  this.rep.OccProb <- lapply(this.rep, function(this.intv) getOccupancyOneInterval_evan(this.intv = this.intv, settings = settings, intervals = intervals))
  this.rep.OccProb
}

getOccList2Mat <- function(Occs.list, col.label)
{
  outMat <- matrix(data = NA, nrow = length(Occs.list), ncol = c(nrow(intervals)))
  colnames(outMat) <- rownames(intervals)
  
  this.rep.vec <- vector()
  for(xx in seq_len(length(Occs.list)))
  {
    for(yy in seq_len(length(Occs.list[[xx]])))
    {
      this.rep.vec[yy] <- Occs.list[[xx]][[yy]][col.label]
    }
    outMat[xx,] <- unlist(this.rep.vec)
  }
  outMat
}
