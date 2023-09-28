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

getCorrelationPlot <- function(prop1, prop2, mar = c(0,0,0,0), cl.pos = NULL, tl.pos = c("lt"))
{
  corr.results.both <- cor.p <- matrix(nrow = ncol(prop1), ncol= ncol(prop2)) 
  dimnames(corr.results.both) <- dimnames(cor.p) <- list(colnames(prop1), colnames(prop2))
  
  for (xx in seq_len(ncol(prop1))) {
    for (yy in seq_len(ncol(prop2))) {
      if(sum(complete.cases(prop1[,xx], prop2[,yy]), na.rm = TRUE) < 2) { corr.results.both[xx,yy]  <- NA } else {
        
        corr.results.temp <- cor.test(prop1[,xx], prop2[,yy], method = "spearman")
        
        corr.results.both[xx,yy] <- corr.results.temp$estimate
        cor.p[xx,yy] <- corr.results.temp$p.value
      }
    }
  }
  
  corrplot::corrplot(corr.results.both,  mar = mar, p.mat = cor.p,
                     cl.align.text = 'l', cl.pos = cl.pos,
                     #addCoef.col = 'black',
                     addgrid.col = 'black',
                     method = "color",
                     na.label = "-",
                     sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig', 
                     pch.cex = 1,
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


