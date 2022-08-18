source('~/Dropbox/code/R/dentalMeasurements/src/src_dentalDataFns.R', chdir = TRUE)
	
####################################################################################################################################
#### read occurrence data
####################################################################################################################################

  occs <- read.csv("http://paleobiodb.org/data1.2/occs/list.csv?base_name=Mammalia&continent=NOA&show=full&limit=all", stringsAsFactors=TRUE, strip.white=TRUE)
  occs <- occs[!occs$order %in% c("Cetacea", "Desmostylia", "Sirenia"), ]
  occs <- occs[!occs$family %in% c("Allodelphinidae", "Balaenidae", "Balaenopteridae", "Delphinidae", "Desmatophocidae", "Desmostylidae", "Didelphidae","Dugongidae","Enaliarctidae", "Eschrichtiidae","Iniidae", "Kentriodontidae", "Kogiidae", "Odobenidae", "Otariidae", "Paleoparadoxiidae", "Phocidae", "Physeteridae", "Platanistidae", "Pontoporiidae", "Protocetidae", "Squalodontidae", "Ziphiidae"), ]
  occs$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = occs$accepted_name)	#replace spaces with underscores

####################################################################################################################################
#### read measurement data
####################################################################################################################################

	measure.mat <- getSingleSpeciesMatrix()
	measure.mat <- measure.mat[-which(measure.mat$taxon=="Merycoidodon_(Merycoidodon)_presidioensis"),]		### this is a big outlier on PC2 x PC3 plot

#	measure.mat <- measure.mat[-which(measure.mat$taxon=="Cephalophus_silvicultor"),]						### this is a big outlier on PC2 x PC3 plot

	matrix(sort(unique(measure.mat$taxon[!measure.mat$taxon %in% occs$accepted_name])), ncol=1)
	matrix(sort(unique(occs$accepted_name[!occs$accepted_name %in% measure.mat$taxon & occs$order %in% focal.order & occs$accepted_rank=="species"])), ncol=1)

#PCA by specimen
	# pcaUp <- prcomp(log(specimenMat[complete.cases(specimenMat[,upLabels]),upLabels]))
	# pcaLo <- prcomp(log(specimenMat[complete.cases(specimenMat[,loLabels]),loLabels]))
	# pcaAll <- prcomp(log(specimenMat[complete.cases(specimenMat[,c(upLabels, loLabels)]),c(upLabels, loLabels)]))
	
#PCA by species
	upLabels<-c("P3_L","P3_W","P4_L","P4_W","M1_L","M1_W","M2_L","M2_W","M3_L","M3_W") #"P2_L","P2_W",
	loLabels <- casefold(upLabels)
	pcaUp <- prcomp(log(measure.mat[complete.cases(measure.mat[,upLabels]),upLabels]))
	pcaLo <- prcomp(log(measure.mat[complete.cases(measure.mat[,loLabels]),loLabels]))
	pcaAll <- prcomp(log(measure.mat[complete.cases(measure.mat[,c(upLabels, loLabels)]),c(upLabels, loLabels)]))

	# # # # # quartz()
	pc_h<- 2
	pc_v<- 3
	# par(mfrow=c(3,1), mar=c(4, 4, 1, 1))
	# famList <- read.csv("~/Dropbox/code/common_dat/taxonomy.csv")
	

####################################################################################################################################
#### set up colors and symbols for plots
####################################################################################################################################

  #by family
	# focal.order <- "Artiodactyla"
	# focal.order <- "Perissodactyla"
	# focal.order <- c("Artiodactyla", "Perissodactyla")
	# bigList <- unique(occs[((occs$accepted_rank =="species" | occs$accepted_rank =="genus") & occs$order %in% focal.order), c("order","family", "genus", "accepted_name")])
	# bigList <- bigList[order(bigList$order, bigList$family, bigList$genus, bigList$accepted_name),]
	# bigList[order(bigList$family, bigList$accepted_name),]

	m <- getTaxonomyForTaxa(tax.vec=sort(unique(gsub(pattern = "_", replacement = " ", x = measure.mat$taxon) )))
	m$species <- gsub(pattern = " ", replacement = "_", x = m$species)

	shortFam <- sort(unique(m$family))	
	
	# bigList$accepted_name <- gsub(pattern = "[[:space:]]", replacement = "_", x = bigList$accepted_name)

   #by parent
	# bigList <- cbind(as.character(famList$taxon[famList$taxon %in% unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]), as.character(famList$occurrences.parent_name[famList$taxon%in%unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))]))
	# bigList[bigList[,2] =="",2] <- as.character(famList$occurrences.family_name[famList$taxon %in% unique(c(rownames(pcaUp$x),rownames(pcaLo$x)))])[bigList[,2] ==""]

	famColors <- array(rainbow(length(shortFam)), dimnames=list(shortFam))
	symbolVec <- array(data=c(15:17, 21:22,24), dim=length(shortFam), dimnames=list(shortFam))

####################################################################################################################################
#### plot PCA using both uppers and lowers
####################################################################################################################################
  quartz()
	plot(pcaAll$x[,pc_h], pcaAll$x[,pc_v], xlim=c(min(pcaAll$x[,pc_h]),1.5*max(pcaAll$x[,pc_h])), xlab=paste("PC",pc_h," (",round(100*(pcaAll$sdev[pc_h]^2/sum(pcaAll$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaAll$sdev[pc_v]^2/sum(pcaAll$sdev^2)), digits=1),"%)", sep=""), type="n", main="Upper and Lower P3-M3")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(pcaAll$x[,pc_h], pcaAll$x[,pc_v], labels=rownames(pcaAll$x), cex=0.5, col= famColors[match(m$family[match(rownames(pcaAll$x), m$species)])

	text(pcaAll$x[,pc_h], pcaAll$x[,pc_v], labels=rownames(pcaAll$x), pos=4, cex=0.3, 
	     col=famColors[as.character(m$family[match(rownames(pcaAll$x), m$species)], shortFam)])
	points(pcaAll$x[,pc_h], pcaAll$x[,pc_v], cex=1.0, 
	       pch=symbolVec[match(m$family[match(rownames(pcaAll$x), m$species)])], 
	                           col=famColors[match(m$family[match(rownames(pcaAll$x), m$species)])])
	legend("bottomright", legend=shortFam[shortFam%in%bigList[bigList$accepted_name %in% rownames(pcaAll$x),2]], 
	       pch=symbolVec[shortFam%in%bigList[bigList$genus%in%rownames(pcaAll$x),2]], 
	       col=famColors[shortFam%in%bigList[bigList$genus%in%rownames(pcaAll$x),2]], 
	       box.col="gray50", bg="white", cex=0.55)
	
# # # # # pcaAll$rotation

####################################################################################################################################
#### plot PCA using only uppers
####################################################################################################################################

	plot(pcaUp$x[,pc_h], pcaUp$x[,pc_v], xlim=c(min(pcaUp$x[,pc_h]),1.5*max(pcaUp$x[,pc_h])), xlab=paste("PC",pc_h, " (",round(100*(pcaUp$sdev[pc_h]^2/sum(pcaUp$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaUp$sdev[pc_v]^2/sum(pcaUp$sdev^2)), digits=1),"%)", sep=""), type="n", main="Upper P3-M3", bg="gray90")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(pcaUp$x[,pc_h], pcaUp$x[,pc_v], labels=rownames(pcaUp$x), cex=0.33, col=famColors[match(bigList[match(rownames(pcaUp$x), bigList$genus),2], shortFam)])
	# text(pcaUp$x[,pc_h], pcaUp$x[,pc_v], labels=gsub("[[:space:]]", "\n", rownames(pcaUp$x)), cex=0.3, col=famColors[match(bigList[match(sapply(strsplit(rownames(pcaUp$x),"[[:punct:]]"), function(x) x[1]), bigList$genus),2], shortFam)])
	text(pcaUp$x[,pc_h], pcaUp$x[,pc_v], labels=rownames(pcaUp$x), pos=4, cex=0.3, col=famColors[match(bigList[match(sapply(strsplit(rownames(pcaUp$x),"[[:punct:]]"), function(x) x[1]), bigList$genus),2], shortFam)])
	points(pcaUp$x[,pc_h], pcaUp$x[,pc_v], cex=0.5, pch=symbolVec[match(bigList[match(rownames(pcaUp$x), bigList$accepted_name),2], shortFam)], col=famColors[match(bigList[match(rownames(pcaUp$x), bigList$accepted_name),2], shortFam)])
	legend("bottomright", legend=shortFam[shortFam %in% bigList[bigList$accepted_name %in% rownames(pcaUp$x),2]], pch=symbolVec[shortFam %in% bigList[bigList$accepted_name %in% rownames(pcaUp$x),2]], col=famColors[shortFam%in%bigList[bigList$accepted_name %in% rownames(pcaUp$x),2]], box.col="gray50", bg="white", cex=0.55)

 # # # pcaUp$rotation

####################################################################################################################################
#### plot PCA using only lowers
####################################################################################################################################

	plot(pcaLo$x[,pc_h], pcaLo$x[,pc_v], xlim=c(min(pcaLo$x[,pc_h]),1.5*max(pcaLo$x[,pc_h])), xlab=paste("PC",pc_h," (",round(100*(pcaLo$sdev[pc_h]^2/sum(pcaLo$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC",pc_v," (",round(100*(pcaLo$sdev[pc_v]^2/sum(pcaLo$sdev^2)), digits=1),"%)", sep=""), type="n", main="Lower P3-M3")
	# polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
	lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
	lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
	# text(pcaLo$x[,pc_h], pcaLo$x[,pc_v], labels=gsub("[[:space:]]", "\n", rownames(pcaLo$x)), cex=0.3, col=famColors[match(bigList[match(rownames(pcaLo$x), bigList$genus),2], shortFam)])
	text(pcaLo$x[,pc_h], pcaLo$x[,pc_v], labels=rownames(pcaLo$x), pos=4, cex=0.2, col=famColors[as.character(m$family[match(sapply(rownames(pcaLo$x), function(x) strsplit(x, split="_")[[1]][1]), m$genus)])])
	points(pcaLo$x[,pc_h], pcaLo$x[,pc_v], cex=0.5, pch=symbolVec[as.character(m$family[match(sapply(rownames(pcaLo$x), function(x) strsplit(x, split="_")[[1]][1]), m$genus)])], col=famColors[as.character(m$family[match(sapply(rownames(pcaLo$x), function(x) strsplit(x, split="_")[[1]][1]), m$genus)])])
	legend("bottomright", legend=shortFam, pch=symbolVec[shortFam], col=famColors[shortFam], box.col="gray50", bg="white", cex=0.55)

# # pcaLo$rotation

####################################################################################################################################
#### check correlations between PC's from uppers to PCs from lowers
####################################################################################################################################

pca.mat <- merge(data.frame(species=rownames(pcaLo$x), pcaLo$x), data.frame(species=rownames(pcaUp$x), pcaUp$x), all=TRUE, by="species", suffixes=c(".lo", ".up"))
pca.mat <- merge(pca.mat, data.frame(species=rownames(pcaAll$x), pcaAll$x), all=TRUE, by="species")

plot(pca.mat$PC2, pca.mat$PC2.lo, xlim=range(c(pca.mat$PC2, pca.mat$PC2.lo), na.rm=TRUE), ylim=range(c(pca.mat$PC2, pca.mat$PC2.lo), na.rm=TRUE))
abline(a=0, b=1, lty=3, col="gray50")
plot(pca.mat$PC2, pca.mat$PC2.up, xlim=range(c(pca.mat$PC2, pca.mat$PC2.up), na.rm=TRUE), ylim=range(c(pca.mat$PC2, pca.mat$PC2.up), na.rm=TRUE))
abline(a=0, b=1, lty=3, col="gray50")
plot(pca.mat$PC2.lo, pca.mat$PC2.up, xlim=range(c(pca.mat$PC2.lo, pca.mat$PC2.up), na.rm=TRUE), ylim=range(c(pca.mat$PC2.lo, pca.mat$PC2.up), na.rm=TRUE))
abline(a=0, b=1, lty=3, col="gray50")
cor.test(pca.mat$PC2.lo, pca.mat$PC2.up, method="pearson")

plot(pca.mat$PC3.lo, -pca.mat$PC3.up, xlim=range(c(pca.mat$PC3.lo, -pca.mat$PC3.up), na.rm=TRUE), ylim=range(c(pca.mat$PC3.lo, -pca.mat$PC3.up), na.rm=TRUE))
abline(a=0, b=1, lty=3, col="gray50")
cor.test(pca.mat$PC3.lo, pca.mat$PC3.up, method="pearson")

plot(pcaLo$rotation[, pc_v], -pcaUp$rotation[, pc_v], xlim=range(c(pcaLo$rotation[, pc_v], -pcaUp$rotation[, pc_v]), na.rm=TRUE), ylim=range(c(pcaLo$rotation[, pc_v], -pcaUp$rotation[, pc_v]), na.rm=TRUE))
abline(a=0, b=1, lty=3, col="gray50")
cor.test(pcaLo$rotation[, pc_v], pcaUp$rotation[, pc_v], method="pearson")

####################################################################################################################################
#### check correlation between PC3 and hypsodonty (horse data)
####################################################################################################################################

dat.hyps <- read.csv("~/Desktop/horse_isotopes_2006_all.csv")
names(dat.hyps)[names(dat.hyps)=="species"] <- "sp_name"
dat.hyps$species <- paste(dat.hyps$Genus, dat.hyps$sp_name, sep="_")

dat <- merge(data.frame(species=rownames(pcaLo$x), pcaLo$x), dat.hyps, all.x=FALSE, by="species")
plot(dat$HI, dat$PC3)
abline(a=0, b=1, lty=3, col="gray50")
abline(lm(dat$PC3 ~ dat$HI))
cor.test(dat$HI, dat$PC3, method="spearman")

####################################################################################################################################
#### check correlation between PC3 and hypsodonty (mostly artiodactyl data from Mendoza and Palmquist 2008)
####################################################################################################################################

dat.hyps <- read.csv("~/Dropbox/code/common_dat/Mendoza&Palmqvist_2008.csv")
dat.hyps$binomial <- gsub(pattern = " ", replacement = "_", x = dat.hyps$binomial)
names(dat.hyps)[names(dat.hyps)=="species"] <- "sp_name"
names(dat.hyps)[names(dat.hyps)=="binomial"] <- "species"

dat <- merge(data.frame(species=rownames(pcaLo$x), pcaLo$x), dat.hyps, all.x=FALSE, by="species")
plot(dat$HI, dat$PC3, type="n")
text(dat$HI, dat$PC3, labels=dat$species, cex=0.5)
abline(a=0, b=1, lty=3, col="gray50")
abline(lm(dat$PC3 ~ dat$HI))
cor.test(dat$HI, dat$PC3, method="spearman")

####################################################################################################################################
#### [WARNING: VERY INCOMPLETE]
#### Begin analyses to determine which measurements most strongly correlated with hypsodont 
####################################################################################################################################


dat.l <- log(measure.mat[complete.cases(measure.mat[,loLabels]),loLabels])
dat.l <- t(apply(dat.l, 1, function(x) x-sum(x)))
dat <- merge(data.frame(species=rownames(dat.l), dat.l), dat.hyps, all.x=FALSE, by="species")
fit <- lm(formula=dat$HI ~ dat$p3_l + dat$p3_w + dat$p4_l + dat$p4_w + dat$m1_l + dat$m1_w + dat$m2_l + dat$m2_w + dat$m3_l + dat$m3_w)
fit.anova <- anova(fit)
fit.anova <- aov(formula=HI ~ (p3_l + p3_w + p4_l + p4_w + m1_l + m1_w + m2_l + m2_w + m3_l + m3_w), data=dat)

