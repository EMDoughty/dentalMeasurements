#replicate jon's dental PCA plots

#read in/make dental file in measuremat
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

#occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&lngmin=-125.98&lngmax=-93.40&latmin=27&latmax=55.7&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                 "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                 "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", 
                                 "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

####################################################################################################################################

#dental measures are in centimeters while body mass is logged kg
measure.mat <- getMeasureMatWithBodyMasses()
this.rank <- "species"
if (this.rank=="genus") measure.mat <- makeOneGenusMatFromSpecimenMat(measure.mat)

####################################################################################################################################
#### reduces matrix to just the focal order(s)
####################################################################################################################################

# focal.order <- "Artiodactyla"
# focal.order <- "Perissodactyla"
focal.order <- c("Artiodactyla", "Perissodactyla")
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
focal.family <- c(as.character(focal.family),"Arctocyonidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.family <- focal.family[!focal.family %in% ""]
focal.family <- focal.family[order(focal.family)]

bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% focal.family), c("order","family", "genus", "accepted_name")])
bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFam <- sort(unique(bigList$family[bigList$family %in% focal.family]))	

bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)

# matrix(measure.mat$taxon[!measure.mat$taxon %in% bigList$accepted_name[bigList$order %in% focal.order]], ncol=1)
# matrix(bigList$accepted_name[bigList$order %in% focal.order][! bigList$accepted_name[bigList$order %in% focal.order] %in% measure.mat$taxon], ncol=1)
measure.mat <- measure.mat[measure.mat$taxon %in% bigList$accepted_name[bigList$order %in% focal.order], ]

#remove added columns...cus why I set it to not do that but still gets through....
#fix.columns<- c("taxon", "P2_L", "P2_W", "P3_L", "P3_W", "P4_L", "P4_W", "M1_L", "M1_W", "M2_L", "M2_W", "M3_L", "M3_W",
# "p2_l", "p2_w", "p3_l", "p3_w", "p4_l", "p4_w", "m1_l", "m1_w", "m2_l", "m2_w", "m3_l", "m3_w",
# "p4_a", "m1_a", "m2_a", "m3_a", "M2_A",
# "family", "genus", "reg.vec", "bodyMass", "col.fam", "sym.fam")
##
#measure.mat <- measure.mat[,colnames(measure.mat) %in% fix.columns]

#remove Megacerops curtus due to wierd placement on pca
measure.culled <- measure.mat[!rownames(measure.mat) %in% "Megacerops_curtus",]


archaic.ung <- read.csv("/Users/emdoughty/Dropbox/Code/ArchaicUngulate_UploadFile_2021_4_29.csv")

archaic.Mat.tot <- getMeasureMatCondylarths(data.raw = archaic.ung, occs = occs, 
                                            col.order = colnames(measure.culled), 
                                            all.bm = FALSE,
                                            regression = "ArchaicNonselenodonts")

write.csv(archaic.Mat.tot,"/Users/emdoughty/Dropbox/Code/R/Carnivore Paleobio Seminar/ArchaicUngulateBodyMass.csv")

measure.culled <- rbind(measure.culled, archaic.Mat.tot)

#using measure.mat so modern families are set. 
#im directly adding in settings for the archics. likely not best way probably but hey it works until it don't
plot.det <- get.plot.det.Settings(measure.mat,
                                  archaic.fam = c("Arctocyonidae", "Hyopsodontidae","Periptychidae","Phenacodontidae"),
                                  archaic.col = c("darkgreen", "magenta","darkblue", "darkgrey"),
                                  archaic.sym = c(42,42,42,42)) 

measure.culled <- merge(measure.culled, plot.det, by = "family", all = FALSE, no.dups = TRUE)

rownames(measure.culled) <- measure.culled$taxon

#write.csv(measure.culled, file = "/Users/emdoughty/Dropbox/Code/CheckMeasureCulled.csv")

#select dental measures
#("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W",
#  "p2_l","p2_w","p3_l","p3_w","p4_l","p4_w","m1_l","m1_w","m2_l","m2_w","m3_l","m3_w")

col.rem <- c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W",
             "p2_l","p2_w","p4_a","m1_a","m2_a","m3_a","M2_A")

getDentalPCA(measure.culled, 
             col.rem = col.rem,
             plot.pca = TRUE, 
             this.x = 2,
             this.y = 3,
             scale.load = 3)

#begin setting up code to make plots I discussed with mark on 6/7/2021

#make both surface area (L*W) and relative shape ration (L/W)
#rains plots of both

#need area for p3
measure.culled$p2_a <- measure.culled$p2_l*measure.culled$p2_w
measure.culled$p3_a <- measure.culled$p3_l*measure.culled$p3_w

#get molarization ratios L/W
measure.culled$P2_R <- measure.culled$P2_L/measure.culled$P2_W
measure.culled$P3_R <- measure.culled$P3_L/measure.culled$P3_W
measure.culled$P4_R <- measure.culled$P4_L/measure.culled$P4_W
measure.culled$M1_R <- measure.culled$M1_L/measure.culled$M1_W
measure.culled$M2_R <- measure.culled$M2_L/measure.culled$M2_W
measure.culled$M3_R <- measure.culled$M3_L/measure.culled$M3_W

measure.culled$p2_r <- measure.culled$p2_l/measure.culled$p2_w
measure.culled$p3_r <- measure.culled$p3_l/measure.culled$p3_w
measure.culled$p4_r <- measure.culled$p4_l/measure.culled$p4_w
measure.culled$m1_r <- measure.culled$m1_l/measure.culled$m1_w
measure.culled$m2_r <- measure.culled$m1_l/measure.culled$m2_w
measure.culled$m3_r <- measure.culled$m3_l/measure.culled$m3_w

quartz(height = 22, width= 8)
par(mfrow = c(5,2))
plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "p3_r", column.y = "p3_a", 
                  text.on = FALSE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = TRUE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "p3_r", column.y = "p3_a", 
                  text.on = TRUE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = FALSE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "p4_r", column.y = "p4_a", 
                  text.on = FALSE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = TRUE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "p4_r", column.y = "p4_a", 
                  text.on = TRUE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = FALSE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "m1_r", column.y = "m1_a", 
                  text.on = FALSE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = TRUE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "m1_r", column.y = "m1_a", 
                  text.on = TRUE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = FALSE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "m2_r", column.y = "m2_a", 
                  text.on = FALSE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = TRUE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "m2_r", column.y = "m2_a", 
                  text.on = TRUE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = FALSE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "m3_r", column.y = "m3_a", 
                  text.on = FALSE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = TRUE, legend.x = "bottomleft", cex.legend = 0.45)

plot.Dental.ratio(plot.data = measure.culled, 
                  column.x = "m3_r", column.y = "m3_a", 
                  text.on = TRUE,
                  cex.plot = 0.25, pos = 4, adj = 0.5, cex.text = 0.10,
                  legend.on = FALSE, legend.x = "bottomleft", cex.legend = 0.45)


#normalize between dentitions

attributes(measure.culled)$names

log.colnames <- c("P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W",
                  "p2_l","p2_w","p3_l","p3_w","p4_l","p4_w","m1_l","m1_w","m2_l","m2_w","m3_l","m3_w", 
                  "p4_a", "m1_a", "m2_a", "m3_a", "M2_A")

measure.Dental.log <- measure.culled

measure.Dental.log[,log.colnames] <- log10(measure.Dental.log[,log.colnames])

{
lowers.col.name <- c("accepted_name", "order", "family", "genus", 
                   # "p2_l", "p2_w", 
                    "p3_l", "p3_w", 
                    "p4_l", "p4_w", 
                    "m1_l", "m1_w",
                    "m2_l", "m2_w", 
                   "m3_l", "m3_w")

archaic.lowers <- archaic.merg[,colnames(archaic.merg) %in% lowers.col.name]
archaic.lowers <- archaic.lowers[complete.cases(archaic.lowers),]
data.coverage(archaic.lowers,bigList.check,focal.order,"family")

#temporal coverage
##make/add to body mass shoulder plot


#plots to illistrate molarization

#plot p4 l vs width....should see these equal out over time if molarization is present
plot(archaic.lowers$p4_l, archaic.lowers$p4_w, xlab = "p4_l", ylab = "p4_w")
text(archaic.lowers$p4_l, archaic.lowers$p4_w, label = archaic.lowers$accepted_name,
     cex = 0.35, adj = 1.15)

par(mfrow = c(ceiling((ncol(archaic.lowers)-4)/4),2))
for(xx in seq(1, ncol(archaic.lowers)-4,2))
  length(lowers.col.name)
{
  plot(archaic.lowers[,4+xx], archaic.lowers[,4+xx+1], 
       xlab = colnames(archaic.lowers)[4+xx], ylab = colnames(archaic.lowers)[4+xx+1])
  text(archaic.lowers[,4+xx], archaic.lowers[,4+xx+1], label = archaic.lowers$accepted_name,
       cex = 0.35, adj = 1.15)
}

#also do width/length as that can maybe better illistrate molarization
plot(archaic.lowers$p4_w/archaic.lowers$p4_l,archaic.lowers$m1_w/archaic.lowers$m1_l, 
     xlab = "p4_w/p4_l", ylab = "m1_w/m1_l")
text(archaic.lowers$p4_w/archaic.lowers$p4_l,archaic.lowers$m1_w/archaic.lowers$m1_l, 
     label = archaic.lowers$accepted_name,
     cex = 0.35, adj = 1.15)
 

archaic.area <- cbind(archaic.lowers[,1:4],p4_A = archaic.p4_A, m1_A = archaic.m1_A)
}





#incorporate in(to PCA of existing ungulate data

archaic.PCA <- (archaic.red-dental.pca$center) * dental.pca$rotation

points(archaic.PCA[,this.x], archaic.PCA[,this.y], 
       pch = archaic.merg$symbol[archaic.merg$accepted_name %in% rownames(archaic.PCA)], 
       col = archaic.merg$col[archaic.merg$accepted_name %in% rownames(archaic.PCA)])
text(archaic.PCA[,this.x], archaic.PCA[,this.y],labels = rownames(archaic.Mat), 
     col = archaic.merg$col[archaic.merg$accepted_name %in% rownames(archaic.PCA)], 
     cex= 0.25, adj = -0.5)

#try version where dataset are all together before pca
#append


#run pca





plot(measure.culled$p4_l, measure.culled$p4_w, xlab = "p4_l", ylab = "p4_w")
text(measure.culled$p4_l, measure.culled$p4_w, label = measure.culled$taxon,
     cex = 0.35, adj = 1.15)

plot.colnames <- c(#"P2_L","P2_W","P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W",
                   #"p2_l","p2_w",
                   #"p3_l","p3_w",
                   "p4_l","p4_w",
                   "m1_l","m1_w")#,
                   #"m2_l","m2_w",
                   #"m3_l","m3_w",
                    #M2_A,
                    #p3_a, 
                    #p4_a, 
                    #m1_a, 
                    #m2_a, 
                    #m3_a
    

plot.data <- measure.culled[, c("family", "taxon", plot.colnames)]
plot.data <- plot.data[complete.cases(plot.data),]

point.col <- merge(plot.data, plot.det, by = "family", all = FALSE, no.dups = TRUE)

par(mfrow = c(ceiling((length(plot.colnames))/4),2))
for(xx in seq(1, ncol(plot.data)-2,2))
  #length(lowers.col.name)-4
{
  plot(plot.data[,2+xx], plot.data[,2+xx+1], 
       xlab = colnames(plot.data)[2+xx], ylab = colnames(plot.data)[2+xx+1])
  text(plot.data[,2+xx], plot.data[,2+xx+1], label = plot.data$taxon,
       cex = 0.15, adj = 1.15)
}


#DON"T USE LOG FOR THIS PLOT
dev.off()
#also do width/length as that can maybe better illistrate molarization
plot(plot.data$p4_w/plot.data$p4_l,plot.data$m1_w/plot.data$m1_l, 
     xlab = "p4_w/p4_l", ylab = "m1_w/m1_l", pch = point.col$symbol)
point.col <- merge(plot.data, plot.det, by = "family", all = TRUE, no.dups = TRUE)
points(plot.data$p4_w/plot.data$p4_l,plot.data$m1_w/plot.data$m1_l, 
      col = as.character(point.col$col), pch = point.col$symbol)
text(plot.data$p4_w/plot.data$p4_l, plot.data$m1_w/plot.data$m1_l, 
     label = plot.data$taxon,
     cex = 0.15, adj = 1.15)

legend.info <- point.col %>% distinct(family, symbol,col)
legend(x  = "topleft", legend = legend.info$family, 
       pch = as.numeric(legend.info$symbol), 
       col = as.character(legend.info$col), cex = 0.5)

#estimate body size for acrhaics
#Dasmuth 1990 chapter indicates that Paleogene taxa have larger molar relative to body size=overestimate if use living species equations
#mentions molar thinner over time=shift in diet form omnivore/browser to grazer
##may be best to use nonselenodont dentition rather than all ungualte to avoid/limit the biases/error

##Dasmuth alos indicates on page 242 that grazers have thinner teeth than browsers when comparing body masses
####these grazers tend to be larger so reletively speaking thinner teeth despite the modification for chewing area

#combine archaics with Artio/Perisso dataset and set to all ungulate for now

#stepwise multi reggression over 54 dental measures shows each one capture soemthing of body mass not captured by body length



#####################################################################################################################
#Conduct code for Pred/Prey analysis
#####################################################################################################################

#plot and analyse archic and modern ungulate bm and compare to predator 
#may be best to compile the dataset for both archaic and modern ungulates and run through Handley method.  Do same for predators for rbl and body mass.

#1) update the bm through time shoulder plots using the archaics

#2) handley method to discern shifts in regime for bm and rbl compare when and where shifts occur.

#3) mark mentioned using first-difference tests?

#4)......

#5) profit.....

measure.culled <- getDateTaxa(measure.mat = measure.culled, occs=occs, this.rank = "species")

shoulderPlot(measure.mat = measure.culled, plot.y = "bodyMass", occs = occs, this.rank = "species", 
             bigList = bigList, shortFam = shortFam, repIntTaxa = repIntAll$RepsTaxa, quants = bm_quants,
             intervals = makeIntervals(1, 66, 2), optList_bm_median = NULL, plot.breaks = FALSE,
             ylab = "log bodymass (kg)", xlab = "Time (Ma)", xaxp = c(65,5,10), cex.axis = 1.5, cex.lab = 1.5, 
             do.subepochs = TRUE,  do.quants = TRUE)

#do hadley anaylsis on extended ungulate dataset
repIntAll <- getRepIntData(measure.mat = measure.culled,
                        col.nam = "bodyMass",
                        int_length = 2, 
                        int_min = 1,
                        int_max = 66,
                        do.parallel = FALSE, 
                        reps = 10,
                        do.subsample = FALSE, 
                        quota = 0.4,
                        do.disparity = FALSE,
                        bootstrapSpecimens = FALSE, #require specimenMat which may be from now depricated versions this code
                        bootstrapSpecies = FALSE,
                        bootstrapSpeciesWithinIntervals = FALSE,
                        do.heuristic = TRUE,
                        extra.intvs = 0,
                        do.rangethrough = TRUE,
                        save.reps = FALSE,
                        plotHist = FALSE)

int_length <- 2
intervals <- makeIntervals(1, 66, int_length)

optListTaxAll <-taxHadley(repIntTaxa = repIntAll$RepsTaxa,
                      occs = occs,
                      shortFam = shortFam,
                      bigList = bigList,
                      intervals = intervals,
                      extra.intvs = 0, 
                      do.parallel=FALSE, 
                      do.heuristic = FALSE,
                      do.save = FALSE)

save(measure.culled, repIntAll, optListTaxAll, file = "~/Dropbox/Code/R/Carnivore Paleobio Seminar/UngulateRepIntBackups_2021_7_12.RData")
load(file = "~/Dropbox/Code/R/Carnivore Paleobio Seminar/UngulateRepIntBackups_2021_7_12.RData")

optListBMAll <- traitHadley(repIntTaxa = repIntAll$RepsTaxa,
                            trait.col = "bodyMass",
                            measure.mat = measure.culled,
                            occs = occs,
                            shortFam = shortFam,
                            bigList = bigList, 
                            intervals = intervals, 
                            bmBreaks = c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf), #Janis 2000  max(measure.mat$bodyMass, na.rm=TRUE)
                            extra.intvs = 0, 
                            do.parallel=TRUE, 
                            do.heuristic = FALSE,
                            do.save = FALSE)


#get histos of where breaks occur (may be best to add that into the Hadley functions)

hadleyHisto(measure.mat = measure.culled, 
            bmBreaks =c(-Inf, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, Inf),
            repIntTaxa = repIntAll$RepsTaxa,
            optList_bm_median = optListBMAll$OptListTraitMedian,
            intervals = intervals)

plotNumberShiftPerRep(optList_tax_allReps = optListTaxAll$OptListTaxAllReps, 
                      optList_bm_allReps = optListBMAll$OptListTraitAllReps)

bm_quants <- getTraitQuants(measure.mat = measure.culled, traitCol = "bodyMass", repIntTaxa = repIntAll$RepsTaxa)

plotHistShiftinInterval(optList_tax_allReps = optListTaxAll$OptListTaxAllReps,
                        optList_bm_allReps = optListBMAll$OptListTraitAllReps,
                        intervals = intervals,
                        reps = 10)

####################################################################################
#plot to compare pred/prey runs
####################################################################################
#readin data

##Ungulates
#load()

#read in older standardized ungulate runs
load("~/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/BM_handleyResult_SampleStandardized=TRUE##------ Thu Jun 20 12:14:30 2019 ------##.Rdata")
load("~/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=TRUE/Taxon_handleyResult_SampleStandardized=TRUE##------ Thu Jun 20 12:12:50 2019 ------##.Rdata")

#read in older non standardized ungulate runs
load("~/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=FALSE/BM_handleyResult_SampleStandardized=FALSE##------ Sun Mar 31 15:24:34 2019 ------##.Rdata")
load("~/Dropbox/ungulate_RA/EcologyResults/Intervals=2Ma_Reps=1000_Subsampled=FALSE/Taxon_handleyResult_SampleStandardized=FALSE##------ Sun Mar 31 15:17:50 2019 ------##.Rdata")

##Carnivores
load("~/Dropbox/Code/CarnivoreHadleyOutput_2021_7_14.Rdata")
load("~/Dropbox/Code/Carnivore100RepHadleyOutput_2021_7_15.Rdata")

load("~/Dropbox/Code/AllCarniv100RepHadleyOutput_2021_7_16.Rdata")

pred.data <- read.csv("~/Dropbox/Code/R/Carnivore Paleobio Seminar/Predator_data_all.csv")

#get column names to be comparable to ungulate data so not break anything
colnames(pred.data) <- c("family", "taxon",	"max_ma","min_ma",	"m1L",	"rbl",	"BM_all_carnivoran",	"BM_extant_reg")

pred.data$taxon <- gsub(pattern = "[[:space:]]", replacement = "_", x = pred.data$taxon)

rownames(pred.data) <- pred.data$taxon

#purge species wihtout dates or body mass
pred.data <- pred.data[complete.cases(pred.data),]

#log bodymasses
pred.data[,c("BM_all_carnivoran",	"BM_extant_reg")] <- log10(pred.data[,c("BM_all_carnivoran",	"BM_extant_reg")])

occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&max_ma=66&min_ma=0&timerule=overlap&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", 
                                 "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", 
                                 "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", 
                                 "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
occs <- occs[!occs$family %in% c("Desmatophocidae","Panotariidae","Phocidae","Odobenidae"),]
occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

focal.order <- c("Carnivora", "Creodonta")
occs <- occs[occs$order %in% focal.order,]
focal.family <- unique(occs[occs$order %in% focal.order,]$family)
#focal.familyPred <- c(as.character(focal.family),"Arctocyonidae", "Hyopsodontidae","Periptychidae","Phenacodontidae")
focal.family <- focal.family[!focal.family %in% ""]
focal.family <- focal.family[order(focal.family)]

bigListPred <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$family %in% focal.family), c("order","family", "genus", "accepted_name")])
bigListPred <- bigListPred[order(bigListPred$order, bigListPred$family, bigListPred$genus, bigListPred$accepted_name),]
# bigList[order(bigList$family, bigList$accepted_name),]
shortFamPred <- sort(unique(bigListPred$family[bigListPred$family %in% focal.familyPred]))	

bigListPred$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigListPred$accepted_name)
pred.data <- getDateTaxa(measure.mat = pred.data, occs=occs, this.rank = "species")

########################
# Number breaks

##Ungulates
quartz()
par(mfrow=c(2,1), mar=c(5,4,4,1))
plotNumberShiftPerRep(hist1.dat = optList_tax_allReps, hist1.breaks = seq(-0.5, 10, 1.0), hist1.col = "orchid4", hist1.main = "Number of taxonomic shifts in each rep (ungulate)", hist1.xlab = "Number of Shifts", hist1.xlim = c(0,10),
                      plothist2 = TRUE, hist2.dat = optList_bm_allReps, hist2.breaks = seq(-0.5, 11.5, 1.0), hist2.col = "firebrick4", hist2.main = "Number of body mass shifts in each rep (ungulate)", hist2.xlab = "Number of Shifts", hist2.xlim = c(0,10))
#plotNumberShiftPerRep(hist1.dat = optList_bm_allReps, hist1.breaks = seq(-0.5, 11.5, 1.0), hist1.col = "firebrick4", hist1.main = "Number of body mass shifts in each rep (extent_reg)", hist1.xlab = "Number of Shifts", hist1.xlim = c(0,10),
#                      plothist2 = FALSE)

##Carnivores
quartz()
par(mfrow=c(3,1), mar=c(5,4,4,1))
plotNumberShiftPerRep(hist1.dat = optListTaxAllPred$OptListTaxAllReps, hist1.breaks = seq(-0.5, 10, 1.0), hist1.col = "orchid4", hist1.main = "Number of taxonomic shifts in each rep", hist1.xlab = "Number of Shifts", hist1.xlim = c(0,10),
                      plothist2 = TRUE, hist2.dat = optListBMAllCarniv$OptListTraitAllReps, hist2.breaks = seq(-0.5, 11.5, 1.0), hist2.col = "firebrick4", hist2.main = "Number of body mass shifts in each rep (all_Carnivoran)", hist2.xlab = "Number of Shifts", hist2.xlim = c(0,10))
plotNumberShiftPerRep(hist1.dat = optListBMAllextantreg$OptListTraitAllReps, hist1.breaks = seq(-0.5, 11.5, 1.0), hist1.col = "firebrick4", hist1.main = "Number of body mass shifts in each rep (extent_reg)", hist1.xlab = "Number of Shifts", hist1.xlim = c(0,10),
                      plothist2 = FALSE)

plotNumberShiftPerRep(hist1.dat = optListBMAllCarniv$OptListTraitAllReps, hist1.breaks = seq(-0.5, 11.5, 1.0), hist1.col = "firebrick4", hist1.main = "Number of body mass shifts in each rep (extent_reg)", hist1.xlab = "Number of Shifts", hist1.xlim = c(0,10),
                      plothist2 = FALSE)

#######################
#location of breaks between runs

##Ungulates
quartz()
par(mfrow=c(2,1), mar=c(4, 4, 1, 1))
plotHistShiftinInterval(optList = optList_tax_allReps,
                        intervals = intervals,
                        reps = 1000,
                        xaxp = c(70,0,14),
                        main = "Number of Replicates with a Taxonomic Distribution Shift (ungulate)",
                        xlab = "Time (Ma)",
                        xlim = c(65,0),
                        cex = 0.3,
                        col = "orchid4",
                        border = "orchid1",
                        do.subepochs = TRUE,
                        labels = TRUE,
                        freq=TRUE)
plotHistShiftinInterval(optList = optList_bm_allReps,
                        intervals = intervals,
                        reps = 1000,
                        xaxp = c(70,0,14),
                        main = "Number of Replicates with a Body Mass Distribution Shift (ungulate)",
                        xlab = "Time (Ma)",
                        xlim = c(65,0),
                        cex = 0.3,
                        col = "firebrick4",
                        border = "firebrick1",
                        do.subepochs = TRUE,
                        labels = TRUE,
                        freq=TRUE)

##Carnivores
quartz()
par(mfrow=c(3,1), mar=c(4, 4, 1, 1))
plotHistShiftinInterval(optList = optListTaxAll$OptListTaxAllReps,
                        intervals = intervals,
                        reps = 10,
                        xaxp = c(70,0,14),
                        main = "Number of Replicates with a Taxonomic Distribution Shift",
                        xlab = "Time (Ma)",
                        xlim = c(65,0),
                        cex = 0.3,
                        col = "orchid4",
                        border = "orchid1",
                        do.subepochs = TRUE,
                        labels = TRUE,
                        freq=TRUE)
plotHistShiftinInterval(optList = optListBMAllCarniv$OptListTraitAllReps,
                        intervals = intervals,
                        reps = 100,
                        xaxp = c(70,0,14),
                        main = "Number of Replicates with a Body Mass Distribution Shift (allCarniv)",
                        xlab = "Time (Ma)",
                        xlim = c(65,0),
                        cex = 0.3,
                        col = "firebrick4",
                        border = "firebrick1",
                        do.subepochs = TRUE,
                        labels = TRUE,
                        freq=TRUE)
plotHistShiftinInterval(optList = optListBMAllextantreg$OptListTraitAllReps,
                        intervals = intervals,
                        reps = 10,
                        xaxp = c(70,0,14),
                        main = "Number of Replicates with a Body Mass Distribution Shift (extant reg)",
                        xlab = "Time (Ma)",
                        xlim = c(65,0),
                        cex = 0.3,
                        col = "firebrick4",
                        border = "firebrick1",
                        do.subepochs = TRUE,
                        labels = TRUE,
                        freq=TRUE)

##############################
#shoulder plots with breaks

##Old Ungulate Data
shoulderPlot(measure.mat = measure.culled, plot.y = "bodyMass", occs = occs, this.rank = "species", 
             bigList = bigList, shortFam = shortFam, repIntTaxa = repIntAll$RepsTaxa, quants = NULL,
             intervals = makeIntervals(1, 66, 2), optList_bm_median = optList_bm_median, plot.breaks = TRUE,
             ylab = "log bodymass (kg) (ungulate)", xlab = "Time (Ma)", xaxp = c(70,0,14), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = FALSE)

##Carnivores
par(mar=c(0,4, 2.5,0.5))
quartz(width=6.89)
par(mfrow=c(3,1), mar=c(3.9,4,0.5,0.5))#, mgp=c(2, 1,0))
plotTaxonThroughTime(repIntTaxa = PredrepIntAll$RepsTaxa, bigList = bigListPred, shortFam = shortFamPred, 
                                 taxCube = optListTaxAllPred$taxCube, intervals = intervals,
                                 optList_tax_median = optListTaxAllPred$OptListTaxMedian, plot.breaks = TRUE)
shoulderPlot(measure.mat = pred.data, plot.y = "BM_all_carnivoran", intervals = intervals, occs = occs, 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = PredrepIntAll$RepsTaxa, quants = bm_quantsPred_allcarniv,
             optList_bm_median = optListBMAllCarniv$OptListTraitMedian, plot.breaks = TRUE, this.rank = "species",
             ylab = "log Bodymass (kg) (all_carnivoran)", xlab = "Time (Ma)", xaxp = c(75,0,15), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = FALSE)
shoulderPlot(measure.mat = pred.data, plot.y = "BM_extant_reg", intervals = intervals, occs = occs, 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = PredrepIntAll$RepsTaxa, quants = bm_quantsPred_extantreg,
             optList_bm_median = optListBMAllextantreg$OptListTraitMedian, plot.breaks = TRUE, this.rank = "species",
             ylab = "log Bodymass (kg) (extant_reg)", xlab = "Time (Ma)", xaxp = c(75,0,15), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = FALSE)

##Ungulate vs Carnivores
quartz(width=6.89)
par(mfrow=c(3,1), mar=c(3.9,4,0.5,0.5))#, mgp=c(2, 1,0))
shoulderPlot(measure.mat = measure.culled, plot.y = "bodyMass", occs = occs, this.rank = "species", 
             bigList = bigList, shortFam = shortFam, repIntTaxa = repIntAll$RepsTaxa, quants = NULL,
             intervals = makeIntervals(1, 66, 2), optList_bm_median = optList_bm_median, plot.breaks = TRUE,
             ylab = "log bodymass (kg) (ungulate)", xlab = "Time (Ma)", xaxp = c(70,0,14), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = TRUE)
shoulderPlot(measure.mat = pred.data, plot.y = "BM_all_carnivoran", intervals = intervals, occs = occs, 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = PredrepIntAll$RepsTaxa, quants = bm_quantsPred_allcarniv,
             optList_bm_median = optListBMAllCarniv$OptListTraitMedian, plot.breaks = TRUE, this.rank = "species",
             ylab = "log Bodymass (kg) (all_carnivoran)", xlab = "Time (Ma)", xaxp = c(75,0,15), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = TRUE)
shoulderPlot(measure.mat = pred.data, plot.y = "BM_extant_reg", intervals = intervals, occs = occs, 
             bigList = bigListPred, shortFam = shortFamPred, repIntTaxa = PredrepIntAll$RepsTaxa, quants = bm_quantsPred_extantreg,
             optList_bm_median = optListBMAllextantreg$OptListTraitMedian, plot.breaks = TRUE, this.rank = "species",
             ylab = "log Bodymass (kg) (extant_reg)", xlab = "Time (Ma)", xaxp = c(75,0,15), cex.axis = 1, cex.lab = 1, 
             do.subepochs = TRUE,  do.quants = FALSE)

####################################
#distribution within regimes

##old ungulate data

##carnivores
regimeHist_v2(repIntTaxa = PredrepIntAll$RepsTaxa, 
              trait.Col = "BM_all_carnivoran",
              breaks = c(41,27,21,11), 
              trait.breaks= ,
              optList=optListBMAllextantreg$OptListTraitMedian, 
              thisMat = pred.data, 
              netFreq = TRUE, 
              regimeFreq=FALSE,
              netPlotType = "pos/neg", 
              plot.together = FALSE,
              plot.quartz = TRUE)
regimeHist_v2(repIntTaxa = PredrepIntAll$RepsTaxa, 
              trait.Col = "BM_extant_reg",
              breaks = c(39,33,23,17,7), 
              trait.breaks= ,
              optList=optListBMAllextantreg$OptListTraitMedian, 
              thisMat = pred.data, 
              netFreq = TRUE, 
              regimeFreq=FALSE,
              netPlotType = "pos/neg", 
              plot.together = FALSE,
              plot.quartz = TRUE)

###################
fd.allcarnivMean <- firstDifferencesBinned(data.input = bm_quantsPred_allcarniv[3,])
fd.extantregMean <- firstDifferencesBinned(data.input = bm_quantsPred_extantreg[3,])





