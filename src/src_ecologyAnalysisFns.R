 # require(gss)
makeRangeThroughOneRep <- function(this.rep) {
	tax.vec <- sort(unique(unlist(this.rep)))
	tax.ranges <- cbind(tax.vec, t(sapply(tax.vec, function(taxon) rev(range(which(sapply(this.rep, function(y, taxon) taxon %in% y, taxon)))))))
	for (i in seq_len(nrow(tax.ranges))) {
		for (i in tax.ranges[i, 2]:tax.ranges[i, 3]) this.rep[[i]] <- sort(unique(c(this.rep[[i]], tax.ranges[i, 1])))
	}
	this.rep
}

getMPWD<-function(datMat) {
	if (nrow(datMat)<2) return(NA)
	# mean(unlist(lapply(1:(nrow(datMat)-1), vGetEuclidianDist, datMat)), na.rm=TRUE)
	# system.time(mean(mapply(oged, rep(1:(nrow(datMat)-1), (nrow(datMat)-1):1), unlist(sapply(1:(nrow(datMat)-1), function(x) (x+1):nrow(datMat))), MoreArgs=list(datMat)), na.rm=TRUE))
	mean(apply(cbind(rep(1:(nrow(datMat)-1), (nrow(datMat)-1):1), unlist(sapply(1:(nrow(datMat)-1), function(x) (x+1):nrow(datMat)))), 1, function(y) cged(datMat[y[1],], datMat[y[2],])), na.rm=TRUE)
}

getMPWD_hard<-function(intSp, thisMat) {
	mean(mapply(cged_h, rownames(thisMat)[intSp[rep(1:(length(intSp)-1), (length(intSp)-1):1)]], rownames(thisMat)[intSp[unlist(sapply(1:(length(intSp)-1), function(x) (x+1):length(intSp)))]], MoreArgs=list(datMat=thisMat)), na.rm=TRUE)
	# mean(mcmapply(cged_h, rownames(thisMat)[intSp[rep(1:(length(intSp)-1), (length(intSp)-1):1)]], rownames(thisMat)[intSp[unlist(sapply(1:(length(intSp)-1), function(x) (x+1):length(intSp)))]], MoreArgs=list(datMat=thisMat), mc.cores=detectCores()-2), na.rm=TRUE)
}

plotUnivStats<-function(statbox, intervals, dispMat=NULL, new.window=TRUE, thisLab=NULL) {
	if (new.window) quartz("Univariate Statistics")
	bigStat<-apply(statbox, c(1,2), median, na.rm=TRUE)
	thisColors<-rainbow(n=dim(bigStat)[2]-1)
	if (new.window) par(mfrow=c(1,dim(bigStat)[2]))
	# titleList<-c("Mean", "Median", "Kurtosis", "Skewness", "Range", "MPWD")
	for (stat in seq_len(dim(bigStat[,!colnames(bigStat)%in%c("range_min", "range_max")])[2])) {
		if (stat>1) thisLab<-""
		plot(rowMeans(intervals), bigStat[,stat], xlim=c(max(rowMeans(intervals)), min(rowMeans(intervals))), ylim=c(min(bigStat[,stat], na.rm=TRUE), max(bigStat[,stat], na.rm=TRUE)), type="l", col=thisColors[stat], main=colnames(bigStat)[stat], xlab="Time (Ma)", ylab=thisLab)
		overlayCzTimescale()
		polygon(x=c(rowMeans(intervals)[1:nrow(intervals)], rowMeans(intervals)[nrow(intervals):1]), y=c(apply(statbox[,stat,], 1, quantile, probs=0.975, na.rm=TRUE), apply(statbox[nrow(statbox[,stat,]):1,stat,], 1, quantile, probs=0.025, na.rm=TRUE)), border=NA, col=rgb(t(col2rgb(thisColors[stat])), alpha=64, maxColorValue=255))
		lines(rowMeans(intervals), bigStat[,stat], col=thisColors[stat], lwd=2)
	}
	plot(rowMeans(intervals), bigStat[,"range_max"]-bigStat[,"range_min"], xlim=c(max(intervals), min(intervals)), type="n", main="Range", xlab="Time (Ma)", ylab=thisLab)
	overlayCzTimescale()
	lines(rowMeans(intervals), bigStat[,"range_max"]-bigStat[,"range_min"], col=thisColors[stat+1], lwd=2)

	plot(rowMeans(intervals), apply(dispMat, 2, median, na.rm=TRUE), xlim=c(max(rowMeans(intervals)), min(rowMeans(intervals))), ylim=c(min(dispMat, na.rm=TRUE), max(apply(dispMat, 2, median, na.rm=TRUE), na.rm=TRUE)), type="l", lwd=2, main="disparity (MPWD)", xlab="Time (Ma)", ylab="Disparity (MPWD)")
	overlayCzTimescale()
	polygon(x=c(rowMeans(intervals), rev(rowMeans(intervals))), y=c(apply(dispMat, 2, quantile, probs=0.975, na.rm=TRUE), rev(apply(dispMat, 2, quantile, probs=0.025, na.rm=TRUE))), border=NA, col=rgb(t(col2rgb("black")), alpha=64, maxColorValue=255))
	lines(rowMeans(intervals), apply(dispMat, 2, median, na.rm=TRUE), lwd=2)
}

getIntervalColors<-function(intervals) {
	intColors<-vector(length=nrow(intervals))
	for (i in seq_len(nrow(intervals))) {
		if (rowMeans(intervals)[i]>48.6) { intColors[i]<-"tomato4"
		} else if (rowMeans(intervals)[i]<48.6 & rowMeans(intervals)[i]>37.2) { intColors[i] <- "firebrick2"
		} else if (rowMeans(intervals)[i]<37.2 & rowMeans(intervals)[i]>33.9) { intColors[i] <- "orange"
		} else if (rowMeans(intervals)[i]<33.9 & rowMeans(intervals)[i]>28.4) { intColors[i] <- "darkolivegreen4"	
		} else if (rowMeans(intervals)[i]<28.4 & rowMeans(intervals)[i]>23.03) { intColors[i] <- "seagreen3"
		} else if (rowMeans(intervals)[i]<23.03 & rowMeans(intervals)[i]>20.43) { intColors[i] <- "dodgerblue2"
		} else if (rowMeans(intervals)[i]<20.43 & rowMeans(intervals)[i]>15.97) { intColors[i] <- "slateblue3"
		} else intColors[i] <- "orchid4"
	}
	intColors
}

plotCluster<-function(histMat, intervals) {
	# plot(bigStat[,3], apply(dispMat, 2, median, na.rm=TRUE)[2:length(intervals)])
	histMat<-histMat[apply(histMat, 1, function(x) { if (all(x==0)) return(FALSE) else return(TRUE) }),]
	histDist<-vegdist(histMat, method="bray", diag=FALSE, upper=FALSE, na.rm=FALSE)
	raw_cluster<-hclust(histDist, method="ward") #same as previous
	quartz("Cluster Analysis")
	plot(raw_cluster, frame.plot=TRUE, cex=0.5, col=getIntervalColors(intervals))
}

plotNMDS<-function(thisBox, intervals, scaler=1, polygon.ints=NULL, title=NULL, filename=NULL) {
	if (length(dim(thisBox)) > 2) thisBox <- apply(thisBox, c(1,2), mean, na.rm=TRUE)
	# thisBox<-apply(thisBox, c(1,2), median, na.rm=TRUE)
	x <- metaMDS(thisBox, distance = "bray", k = 2, trymax = 100, autotransform =TRUE,
	        noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
	        plot = FALSE, old.wa = FALSE, zerodist="add")

	intColors <- getIntervalColors(intervals)
	# ordiplot(x, display = "sites", choices = c(1, 2), type="n", xlim=c(min(x$species, na.rm=TRUE), max(x$species, na.rm=TRUE)), ylim=c(min(x$species, na.rm=TRUE), max(x$species, na.rm=TRUE)))
	# if (!is.null(title)) thisTitle <- paste("NMDS (", title, ")", sep="") else thisTitle<-"NMDS"
	# if (is.null(filename)) { quartz(title=thisTitle, width=scaler*(max(x$points[,1])-min(x$points[,1])),height=scaler*(max(x$points[,2])-min(x$points[,2])))
	# } else pdf(file=filename, width=scaler*(max(x$points[,1])-min(x$points[,1])),height=scaler*(max(x$points[,2])-min(x$points[,2])))
	ordiplot(x, display = "sites", choices = c(1, 2), type="n")
	# arrows(0,0,x$species[,1]*0.5, x$species[,2]*0.5, lwd=1.5, length=0.1, col="gray75")
	# text(x$species[,1]*0.5, x$species[,2]*0.5, labels = round(as.numeric(rownames(x$species)), digits=2), cex=0.5)

	if (!is.null(polygon.ints)) {
		thisBreaks <- sort(unique(c(0,polygon.ints, nrow(intervals))))
		for (i in 2:length(thisBreaks)) {
			shortBreaks <- x$points[thisBreaks[i]:(thisBreaks[(i-1)]+1),]
			polygon(shortBreaks[chull(shortBreaks),], border=alphaColor("black", 0.5), col=alphaColor("gray75", 0.5))
		}
	}
	
	points(x, display = "sites", choices = c(1, 2), col=alphaColor(intColors,0.5), pch=19)
	
	if (!is.null(polygon.ints)) {
		for (i in 2:length(thisBreaks)) {
			shortBreaks <- x$points[thisBreaks[i]:(thisBreaks[(i-1)]+1),]
			points(shortBreaks[chull(shortBreaks),], bg=NA, col=intColors[match(rownames(shortBreaks)[chull(shortBreaks)], rownames(intervals))])
		}
	}

	# shift <- 0.01*(max(x$points[,1])-min(x$points[,1]))
	# text(x$points[,1], x$points[,2], dimnames(thisBox)[[1]], pos=1, col=intColors, cex=0.5)
	legend("topleft", legend=c("late Miocene", "middle Miocene", "early Miocene", "late Oligocene", "early Oligocene", "late Eocene", "middle Eocene", "early Eocene"), col=unique(intColors), pt.bg=alphaColor(unique(intColors),0.5), pch=21, cex=1)
	if (!is.null(filename)) dev.off()
	x$points
	return(x)
}

NMDSDist<-function(thisBox) {
	y<-vector()
	for (i in seq_len(dim(thisBox)[3])) {
		x <- metaMDS(thisBox[,,i], distance = "bray", k = 2, trymax = 100, autotransform =TRUE,
	        noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
	        plot = FALSE, old.wa = FALSE, zerodist="add")
		y<-cbind(y, pairwiseNMDSDistSubsequent(x))	
	}
	cbind(intervals, t(apply(y, 1, quantile, probs=c(0.025, 0.5, 0.975), na.rm=TRUE)))
}

pairwiseNMDSDistSubsequent<-function(x) {
	y<-vector()
	for (i in 2:nrow(x$points)) y<-c(y, sqrt(((x$points[i-1,1]-x$points[i,1])^2)+((x$points[i-1,2]-x$points[i,2])^2)))
	c(y,NA)
}

pairwiseKSTestsSubsequent<-function(thisBox) {
	bigHist <- apply(thisBox, c(1,2), median, na.rm=TRUE)
	y <- list()
	for (i in 2:nrow(bigHist)) y[[(length(y)+1)]] <- ks.test(bigHist[i-1,], bigHist[i,])
	yy <- sapply(y, function(x) { c(statistic=x$statistic, p.value=x$p.value) })
	cbind(intervals, rbind(t(yy), c(NA, NA)))
}


ksMatrix <- function(thisBox) {
	bigHist<-apply(thisBox, c(1,2), mean, na.rm=TRUE)
	y <- matrix(NA, nrow=nrow(bigHist), ncol=nrow(bigHist))
	for (i in seq_len((nrow(bigHist)-1))) {
		for (j in (i+1):nrow(bigHist)) {
			y[i,j]<-ks.test(bigHist[i,],bigHist[j,])$p.value
		}
	}
	# yy<-as.matrix(sapply(y, function(x) { return(x$p.value) }), ncol=1)
	# cbind(intervals, rbind(yy,NA))
	rownames(y)<-rownames(bigHist)
	colnames(y)<-rownames(bigHist)
	y
} 

# src <- '
	# Rcpp::NumericVector n_(n);
	# Rcpp::NumericVector t_(t);
	
	# int nRateClasses = n_[0];
	# int totalClasses = t_[0];

	# arma::rowvec rateClass(arma::zeros(totalClasses));
	# arma::mat klondike(arma::zeros(2,totalClasses));

	# return Rcpp::wrap(klondike);
	
# //	classList = getRateClassListMachine(1, 1, nRateClasses, cartoon, klondike);

# '

# getRateClassList_cpp <- cxxfunction(signature(n="numeric", t="numeric"), src, plugin = "RcppArmadillo", includes="#include \"/Users/jmarcot/Dropbox/code/C++/rateClass/rateClass.h\"\n#include \"/Users/jmarcot/Dropbox/code/C++/rateClass/rateClass.cpp\"")

findOptimalPartitions<-function(thisCol, dframe, intervals, breaks, minTimes=5, minPartDiff=0.5, enforce.normal=FALSE) {
	require(nloptr)
	require(RcppArmadillo)
	require(inline)
	src<-'
		Rcpp::NumericVector partitions(p);
		Rcpp::NumericVector dat(d);
		Rcpp::NumericVector FO(fo);
		Rcpp::NumericVector LO(lo);
		Rcpp::NumericMatrix intervals(v);
		Rcpp::NumericVector br(b);
		arma::vec breaks(br.begin(), br.size(), false);
		arma::vec thisSample;
		arma::vec counts;
		double lnl=0.0;
		int i, j;
		
/*
		for (i=0;i<partitions.length()-1;i++) {
			if (partitions[i]>=intervals(intervals.nrow()-1,1) || partitions[i]<=intervals(0,0)) return Rcpp::wrap(1e10);	// discard if partitions are outside window of observation (intervals)
			for (j=i+1;j<partitions.length();j++) if (abs(partitions[i]-partitions[j])<0.5) return Rcpp::wrap(2e10);		// discard if any partitions are within a minimum time interval of one another
		}
*/
		partitions.push_back(intervals(intervals.nrow()-1,1));		// add the maximum age of the window of observation to the partition list
	
		double top = intervals(0,0) - 1e-5;			// set the "top" of the intitial partition interval to just less than youngest age of the window of observation - the 1e-05 is to allow taxa whose age == the min age of the window
		for (i=0;i<partitions.length();i++) {
			thisSample.reset();
			for (j=0;j<FO.length();j++) if (FO[j]>top && LO[j]<partitions[i]) {
				thisSample.resize(thisSample.n_elem+1);
				thisSample[thisSample.n_elem-1]=dat[j];
			}
			if (thisSample.n_elem<2) return Rcpp::wrap(3e10);							// discard if there are fewer than two taxa in the partition
	
			counts=arma::conv_to<arma::vec>::from(arma::histc(thisSample, breaks));		// histogram counts converted to a double vector (vec)
			counts=counts.elem(arma::find(counts));										// discards histogram categories that have 0 members
			lnl += -sum(counts%(log(counts)-log(thisSample.n_elem)));
	
			top=partitions[i];															// top of next partition interval is the base of this one
		}
		return Rcpp::wrap(lnl);
	'
	lkPartitionModel_cpp <- cxxfunction(signature(p="numeric", d="numeric", fo="numeric", lo="numeric", v="numeric", b="numeric", minPartDiff="numeric"), src, plugin = "RcppArmadillo")
	# lkPartitionModel_cpp(p=partitions, d=thisMat[(is.finite(thisMat[,"bodyMass"])&is.finite(thisMat[,"FO"])),"bodyMass"], fo=thisMat[(is.finite(thisMat[,"bodyMass"])&is.finite(thisMat[,"FO"])),"FO"], lo=thisMat[(is.finite(thisMat[,"bodyMass"])&is.finite(thisMat[,"LO"])),"LO"], v=data.matrix(intervals), b=bmBreaks, minPartDiff=0.5)
	enforce.normal=FALSE
	dframe <- thisMat
	thisCol <- which(colnames(dframe)=="bodyMass")
	# thisCol <- which(colnames(dframe)=="PC3")
	dframe <- dframe[is.finite(dframe[,thisCol]) & is.finite(dframe$FO) & is.finite(dframe$LO),]
	breaks <- bmBreaks
	# breaks<-pc3Breaks
	minTimes <- 5
	minPartDiff <- 0.5
	bigList <- list(objective=lkPartitionModel_cpp(p=vector(), d=data.matrix(dframe[,thisCol]), fo=data.matrix(dframe[,"FO"]), lo=data.matrix(dframe[,"LO"]), v=data.matrix(intervals), b=breaks, minPartDiff=minPartDiff))
	for (i in 1:4) {
		repList<-list()
		eval_g0 <- function( x, minPartDiff, d, fo, lo, v, b) {
			if (length(x)<2) return(0.0)
			out<-vector()
			for (part in 2:(length(x))) {
				out<-c(out, (x[part-1]+minPartDiff)-x[part])
			}
			return(out)
		}
		while(length(repList) < minTimes) { 
			partitions<-sort(runif(i, min=min(intervals)+(2*minPartDiff), max=max(intervals)-(2*minPartDiff)))
			# if (length(bigList)>1) partitions<-sort(c(sample(bigList[[i]], 1)[[1]]$par, runif(1, min=min(intervals), max=max(intervals)))) else partitions<-sort(runif(i, min=min(intervals), max=max(intervals)))
			if (enforce.normal) { thisTry<-optim(partitions, lkpm_n, thisCol=thisCol, dframe=dframe, intervals=intervals, method="L-BFGS-B",lower=rep((min(intervals)+0.5), i), upper=rep((max(intervals)-0.5), i), control=list(maxit=2500, trace=0))#, )) , control=list(), , 
			# } else thisTry<-optim(partitions, lkPartitionModel_cpp, d=data.matrix(dframe[,thisCol]), fo=data.matrix(dframe[,"FO"]), lo=data.matrix(dframe[,"LO"]), v=data.matrix(intervals), b=breaks, method="L-BFGS-B",lower=rep((min(intervals)+minPartDiff), i), upper=rep((max(intervals)-minPartDiff), i), control=list(maxit=2500, trace=0, parscale=rep(31.25, length(partitions)), fnscale=850.0))
			# } else thisTry <- optim(partitions, lkPartitionModel_cpp, d=data.matrix(dframe[,thisCol]), fo=data.matrix(dframe[,"FO"]), lo=data.matrix(dframe[,"LO"]), v=data.matrix(intervals), b=breaks, method="SANN", control=list(maxit=10000, trace=0, parscale=rep(5, length(partitions)), fnscale=850.0))
			# } else thisTry<-nlminb(start=partitions, objective=lkPartitionModel_cpp, gradient=NULL, hessian=NULL, d=data.matrix(dframe[,thisCol]), fo=data.matrix(dframe[,"FO"]), lo=data.matrix(dframe[,"LO"]), v=data.matrix(intervals), b=breaks, control=list(trace=0, eval.max=2000, iter.max=1500), lower=rep(min(intervals)+minPartDiff, length(partitions)), upper=rep(max(intervals)-minPartDiff, length(partitions)))
			} else thisTry<-nloptr(x0=partitions, eval_f=lkPartitionModel_cpp, lb=rep(min(intervals)+minPartDiff, length(partitions)), ub=rep(max(intervals)-minPartDiff, length(partitions)), eval_g_ineq=eval_g0, opts=list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel"=1e-12, "maxeval"=1e06), d=data.matrix(dframe[,thisCol]), fo=data.matrix(dframe[,"FO"]), lo=data.matrix(dframe[,"LO"]), v=data.matrix(intervals), b=breaks, minPartDiff=minPartDiff)
			# if (length(repList)==0 || thisTry$value<repList[[1]]$value) { repList<-list(thisTry)
			# } else if (thisTry$value==repList[[1]]$value) repList[[length(repList)+1]]<-thisTry
			# cat("***This logLikelihood =", thisTry$value, "- found", length(repList), "of", minTimes, i, "-partition models (", repList[[length(repList)]]$value, "[", paste(repList[[length(repList)]]$par), "])\r", sep=" ")
			if (length(repList)==0 || thisTry$objective<repList[[1]]$objective) { repList<-list(thisTry)
			} else if (thisTry$objective==repList[[1]]$objective) repList[[length(repList)+1]]<-thisTry
			cat("***This logLikelihood =", thisTry$objective, "- found", length(repList), "of", minTimes, i, "-partition models (", repList[[length(repList)]]$objective, "[", paste(repList[[length(repList)]]$solution), "])\r", sep=" ")
		}
		cat("\n")
		# lk<-sapply(repList, function(x) x$objective)
		# bigList[[length(bigList)+1]]<-repList[which(lk==min(lk))]
		bigList[[length(bigList)+1]]<-repList
	}
	# while(currentAIC>bestAIC) {
		# bestModel<-optim(breaks, lkPartitionModel, dframe=dframe, intervals=intervals)
		# currentAIC<-bestModel$objective
	# }
}

# makeFooteogram<-function(thisMat, pca, intervals, scaler) {
	# famList<-read.csv("~/Dropbox/code/R/common_dat/taxonomy.csv")
	# bigList<-data.frame(taxon=as.character(famList$taxon[famList$taxon%in%rownames(pca$x)]), family=as.character(famList$occurrences.family_name[famList$taxon%in%rownames(pca$x)]), stringsAsFactors=FALSE)
	# bigList[bigList[,2] =="",2]<-as.character(famList$occurrences.parent_name[famList$taxon%in%rownames(pca$x)])[bigList[,2] ==""]
	# famLegend<-data.frame(family=sort(unique(bigList[,2])), symbol=array(data=c(15:17, 21:22,24), dim=length(sort(unique(bigList[,2])))), color=rainbow(length(sort(unique(bigList[,2])))), stringsAsFactors=FALSE)
	# bigList<-merge(bigList, famLegend)
	# for (int in seq_len(nrow(intervals))) {
		# setwd("~/Desktop/amandaFigs/")
		# png(filename=paste("pca_", regmatches(rownames(intervals)[int], regexpr(" ", rownames(intervals)[int]), invert=TRUE)[[1]], ".png", sep=""), width=scaler*(max(pca$x[,2])-min(pca$x[,2])), height=scaler*(max(pca$x[,3])-min(pca$x[,3])))
		# intSp<-which(!is.na(thisMat$FO) & !is.na(thisMat$LO) & thisMat$FO>intervals$ageTop[int] & thisMat$LO<intervals$ageBase[int] & !is.na(thisMat$bodyMass))
		# plot(pca$x[rownames(pca$x)%in%rownames(thisMat)[intSp],2], pca$x[rownames(pca$x)%in%rownames(thisMat)[intSp],3], xlim=c(min(pca$x[,2]),1.5*max(pca$x[,2])), ylim=c(min(pca$x[,3]),max(pca$x[,3])), xlab=paste("PC2 (",round(100*(pca$sdev[2]^2/sum(pca$sdev^2)), digits=1),"%)", sep=""), ylab=paste("PC3 (",round(100*(pca$sdev[3]^2/sum(pca$sdev^2)), digits=1),"%)", sep=""), type="n", main=rownames(intervals)[int])
		# # polygon(c(-10,10,10,-10), c(-10,-10,10,10), col="gray33")
		# lines(x=c(0,0), y=c(-100, 100), lty=3, col="gray50")
		# lines(x=c(-100,100), y=c(0, 0), lty=3, col="gray50")
		# # text(pca$x[rownames(pca$x)%in%rownames(thisMat)[intSp],2], pca$x[rownames(pca$x)%in%rownames(thisMat)[intSp],3], labels=rownames(rownames(pca$x)%in%rownames(thisMat)[intSp]), cex=0.5, col=famColors[match(bigList[match(rownames(pca$x), bigList[,1]),2], shortFam)])
		# points(pca$x[rownames(pca$x)%in%rownames(thisMat)[intSp],2], pca$x[rownames(pca$x)%in%rownames(thisMat)[intSp],3], cex=2.0, pch=bigList$symbol[match(rownames(pca$x)[rownames(pca$x)%in%rownames(thisMat)[intSp]], bigList$taxon)], col=bigList$color[match(rownames(pca$x)[rownames(pca$x)%in%rownames(thisMat)[intSp]], bigList$taxon)])
		# legend("bottomright", legend=famLegend$family, pch=famLegend$symbol, col=famLegend$color, box.col="gray50", bg="white", cex=1)
		# dev.off()
	# }
# }

# # lkOnePartition<-function(subVec, breaks) {
	# thisHist<-hist(subVec, breaks=breaks, plot=FALSE)
	# -sum(log((thisHist$counts/sum(thisHist$counts))^thisHist$counts))
# }
# lkp<-cmpfun(lkOnePartition)

# lkPartitionModel<-function(partitions, thisCol, dframe, intervals, breaks) {
	# if (any(diff(partitions)<=0.5) | any(partitions>max(intervals))) return (9999)
	# # lk<-vector(mode="numeric", length=length(partitions))
	# partitions[length(partitions)+1]<-max(intervals)
	# lk<- 0
	# top<-min(intervals)
	# for (i in seq_along(partitions)) {
		# # thisSample<-dframe[dframe$FO>top & dframe$LO<partitions[i], thisCol]
		# thisSample<-dframe[dframe$FO>top & dframe$FO<partitions[i], thisCol]
		# if (length(thisSample)<2) return(9999)
		# # lk[i]<- lkp(thisSample, breaks)
		# lk<- lk+lkp(thisSample, breaks)
		# top<-partitions[i]
	# }
	# # sum(lk,na.rm=TRUE)
	# lk
# }
# lkpm<-cmpfun(lkPartitionModel)

# lkPartitionModel_normal<-function(partitions, thisCol, dframe, intervals) {
	# if (any(diff(partitions)<=0.5) | any(partitions>max(intervals))) return (9999)
	# # lk<-vector(mode="numeric", length=length(partitions))
	# partitions[length(partitions)+1]<-max(intervals)
	# lk<- 0
	# top<-min(intervals)
	# for (i in seq_along(partitions)) {
		# # thisSample<-dframe[dframe$FO>top & dframe$LO<partitions[i], thisCol]
		# thisSample<-dframe[dframe$FO>top & dframe$FO<partitions[i], thisCol]
		# if (length(thisSample)<2) return(9999)
		# # lk[i] <- (-fitdistr(thisSample, "normal")$loglik)
		# # lk<- lk-fitdistr(thisSample, "normal")$loglik
		# lk<- lk-sum(sapply(thisSample, pnorm, mean=mean(thisSample), sd=sd(thisSample), log.p=TRUE))
		# top<-partitions[i]
	# }                                                                                                                                                                               
	# # sum(lk,na.rm=TRUE)
	# lk
# }
# lkpm_n<-cmpfun(lkPartitionModel_normal)

# # likelihoodOneSubset <- function(thisClass, thisSet, intervals, thisMat, traitCol, thisBreaks, thisBox) {
	# if (sum(thisSet==thisClass) > 1) { sum(apply(thisBox[which(thisSet==thisClass),], 1, FUN=dmultinom, prob=hist(thisMat[thisMat$FO > min(intervals[which(thisSet==thisClass),]) & thisMat$LO < max(intervals[which(thisSet==thisClass),]), traitCol], breaks=thisBreaks, plot=FALSE)$counts,log=TRUE))
	# } else dmultinom(thisBox[thisSet==thisClass,], prob=hist(thisMat[thisMat$FO > min(intervals[which(thisSet==thisClass),]) & thisMat$LO < max(intervals[which(thisSet==thisClass),]), traitCol], breaks=thisBreaks, plot=FALSE)$counts,log=TRUE)
# }

# likelihoodOneSetIntervals <- function(thisSet, intervals, thisMat, traitCol, thisBox) {
	# sum(sapply(seq_len(max(thisSet)), likelihoodOneSubset, thisSet, intervals, thisMat, traitCol, thisBreaks, thisBox))
# }

# getLikelihoodsForOneNumberOfRates <- function(nRates, intervals, thisMat, traitCol, thisBox) {
	# simplify2array(mclapply(modelList, likelihoodOneSetIntervals, intervals, thisMat, traitCol, thisBox, mc.cores=(detectCores()-2)))
# }

# lapply(seq_len(nrow(intervals)), getLikelihoodsForOneNumberOfRates, intervals, thisMat, traitCol, thisBox)
# lnl <- lapply(seq_len(5), getLikelihoodsForOneNumberOfRates, intervals, thisMat, traitCol, thisBox)

# getLikelihoodOneSetOfIntervals <- function(thisIntSet, thisCol, breaks) {
	# # if (length(thisIntSet) > 1) { thisRange <- range(unlist(intList[seq(from=thisIntSet[1], to=thisIntSet[2])])) 
	# # } else thisRange <- intList[[thisIntSet]] 
	# # thisFrame <- dframe[is.finite(dframe[,thisCol]) & is.finite(dframe$FO) & dframe$FO > thisRange[1] & dframe$FO <= thisRange[2],] #& is.finite(dframe$LO) 
	# sum(rowSums(thisCounts) * log(rowSums(thisCounts)/sum(thisCounts)))
	# # if (nrow(thisFrame) > 1) {
		# # thisHist <- hist(thisFrame[,thisCol], breaks=breaks, plot=FALSE)
		# # thisCounts <- thisHist$counts[thisHist$counts > 0]
		# # # dmultinom(x=thisCounts, prob=thisCounts, log=TRUE)
		# # sum(log(thisCounts/sum(thisCounts)) * thisCounts, na.rm=TRUE) # + lfactorial(sum(thisCounts)) - sum(lfactorial(thisCounts))
	# # } else NA
# }

getLikelihoodOneSetOfIntervals <- function(this.count.set) {
	# if (!is.null(dim(this.count.set))) { 
		# rs <- rowSums(this.count.set)
		# sum(rs[rs > 0] * log(rs[rs > 0]/sum(this.count.set)))
	# } else sum(this.count.set[this.count.set > 0] * log(this.count.set[this.count.set > 0]/sum(this.count.set)))
	if (!is.null(dim(this.count.set))) { 
		this.count.set <- this.count.set[,colSums(this.count.set)>0]		# remove any columns with all zeros
		if (is.null(dim(this.count.set))) return(0)
		if (ncol(this.count.set)==0) return(0)						# if all columnns removed, then zplit
		this.prob <- rowSums(this.count.set)/sum(this.count.set)
		sum(apply(this.count.set, 2, dmultinom, prob=this.prob, log=TRUE))
	} else {
		if (all(this.count.set==0)) { return(0) 
		} else dmultinom(x=this.count.set, prob=this.count.set/sum(this.count.set), log=TRUE)
	}
	# for (i in seq_len(ncol(this.count.set))) dmultinom(this.count.set[,1], prob=this.prob, log=TRUE)
}

cmp_getLikelihoodOneSetOfIntervals <- cmpfun(getLikelihoodOneSetOfIntervals)

getLikelihoodOneSetIntBreaks_core <- function(this.count.set, intBreaks) {
	lnL <- vector(mode="numeric", length=length(intBreaks)-1)
	for (i in seq_len(length(intBreaks) - 2)) {
		if (intBreaks[i] != intBreaks[i+1]) { lnL[i] <- cmp_getLikelihoodOneSetOfIntervals(this.count.set = this.count.set[, seq(from=intBreaks[i], to=(intBreaks[i + 1] + 1) ) ])
		} else lnL[i] <- cmp_getLikelihoodOneSetOfIntervals(this.count.set = this.count.set[, intBreaks[i]])
	}
	lnL[i+1] <- cmp_getLikelihoodOneSetOfIntervals(this.count.set = this.count.set[, seq(from=intBreaks[i+1], to=(intBreaks[i + 2]))])
	sum(lnL)
}

# getLikelihoodOneSetIntBreaks_heuristic <- function(newBreak, oldIntBreaks=NULL, thisCounts) {
	# intBreaks <- sort(c(1, newBreak, oldIntBreaks, dim(thisCounts[[1]])[2]), decreasing=TRUE)
	# sum(sapply(thisCounts, FUN=getLikelihoodOneSetIntBreaks_core, intBreaks=intBreaks))
	# # for (this.diet in seq_along(thisCounts)) getLikelihoodOneSetIntBreaks_core(this.count.set=thisCounts[[this.diet]], intBreaks)
# }

getLikelihoodOneSetIntBreaks <- function(intBreaks=NULL, thisCounts) {
	intBreaks <- sort(c(1, intBreaks, dim(thisCounts[[1]])[2]), decreasing=TRUE)
	sum(sapply(thisCounts, FUN=getLikelihoodOneSetIntBreaks_core, intBreaks=intBreaks))
}

getOldIntBreaksStarterList_recursive <- function(master.break.list, index, oldIntBreaks, extra.intvs) {
	this.break.list <- oldIntBreaks[index] + (-extra.intvs:extra.intvs)
	master.break.list <- unlist(lapply(master.break.list, function(x) lapply(this.break.list, function(y, x) sort(c(y, x), decreasing=TRUE), x=x)), recursive=FALSE)
	if (index<length(oldIntBreaks)) master.break.list <- getOldIntBreaksStarterList_recursive(master.break.list, index=index+1, oldIntBreaks, extra.intvs= extra.intvs) else return(master.break.list)
}

getHeuristicBreaklist <- function(n.intv, oldIntBreaks, extra.intvs=1) {
	master.break.list <- as.list(oldIntBreaks[1] + (-extra.intvs:extra.intvs))
	master.break.list <- unique(getOldIntBreaksStarterList_recursive(master.break.list=master.break.list, index=2, oldIntBreaks, extra.intvs=extra.intvs))
	master.break.list <- unlist(lapply(X=master.break.list, FUN=function(x) lapply(X=seq_len(n.intv - 1)[-x], FUN=function(y, x) sort(c(y, x), decreasing=TRUE), x=x)), recursive=FALSE)
	master.break.list <- master.break.list[sapply(master.break.list, function(x) all(x>0 & x<n.intv))]	### get rid of any elements with breaks beyond the intervals
	unique(master.break.list)
}

doHandleyTest <- function(thisCounts, n, use.LRT=FALSE, sig=0.01, do.heuristic=TRUE, extra.intvs=0, do.parallel=FALSE, this.cores=NULL) {
	if (do.parallel) {
		require(parallel)
		if (is.null(this.cores)) this.cores <- detectCores() - 2
	}

	if (!is.list(thisCounts)) thisCounts <- list(thisCounts)

	n.hist.classes <- dim(thisCounts[[1]])[1]
	n.intv <- dim(thisCounts[[1]])[2]

	optList=list()
	optList[[1]] <- list(nrates=1, optBreaks=NULL, optlnL= sum(sapply(thisCounts, cmp_getLikelihoodOneSetOfIntervals)))	
	k <- n.hist.classes * 1 # k is the number of size classes times the number of categories(=splits +1)
	optList[[1]]$AIC <- (-2 * optList[[1]]$optlnL) + (2 * k)
	optList[[1]]$AICc <- optList[[1]]$AIC + (2 * k * (k + 1) / (n - k - 1))
	optList[[1]]$k <- k
	optList[[1]]$n <- n

	flag=FALSE
	nrates <- 2
	oldIntBreaks <- NULL
	
	while(!flag & nrates <= (n.intv - 1)) {
	cat("Beginning Handley analysis for ", nrates, " distributions...\r")
		if (do.heuristic) {
			# if (nrates > 2) breakList <- seq_len(n.intv - 1)[-oldIntBreaks] else breakList <- seq_len(n.intv - 1) # breaks are at the base of intervals
			if (nrates > 2) breakList <- getHeuristicBreaklist(n.intv, oldIntBreaks, extra.intvs=extra.intvs) else breakList <- seq_len(n.intv - 1) # breaks are at the base of intervals
			# breakList <- sample(breakList, size=length(breakList))	#randomly reorders breakList to avoid solidifying the "early" breaks
		} else breakList <- combn(x=nrow(intervals), m=nrates-1, simplify=FALSE)

		# if (do.parallel) { lnL <- simplify2array(mclapply(X=breakList, FUN=getLikelihoodOneSetIntBreaks_heuristic, oldIntBreaks=oldIntBreaks, thisCounts=thisCounts, mc.cores=detectCores()-2))
		# } else lnL <- sapply(X=breakList, FUN=getLikelihoodOneSetIntBreaks_heuristic, oldIntBreaks=oldIntBreaks, thisCounts=thisCounts)
		if (do.parallel) { lnL <- simplify2array(mclapply(X=breakList, FUN=getLikelihoodOneSetIntBreaks, thisCounts=thisCounts, mc.cores=this.cores))
		} else lnL <- sapply(X=breakList, FUN=getLikelihoodOneSetIntBreaks, thisCounts=thisCounts)
		# for (this.break in seq_along(breakList)) print(getLikelihoodOneSetIntBreaks_heuristic(newBreak=breakList[this.break], thisCounts=thisCounts))
		# for (this.break in seq_along(breakList)) print(getLikelihoodOneSetIntBreaks(intBreaks=breakList[[this.break]], thisCounts=thisCounts))

		k <- n.hist.classes * nrates # k is the number of size classes x the number of "regimes"
		# if (do.heuristic) { this.optBreaks <- sort(unlist(c(oldIntBreaks, breakList[which(lnL == max(lnL, na.rm=TRUE))])), decreasing=TRUE)
		# } else this.optBreaks <- sort(breakList[[which.max(lnL)]], decreasing=TRUE)
		optList[[nrates]] <- list(nrates=nrates, 
								  optBreaks=sort(breakList[[which.max(lnL)]], decreasing=TRUE),
								  optlnL=max(lnL, na.rm=TRUE), 
								  AIC=(-2 * max(lnL, na.rm=TRUE)) + (2 * k), 
								  AICc=(-2 * max(lnL, na.rm=TRUE)) + (2 * k)  + (((2 * k) * (k + 1)) / (n - k - 1)),
								  k=k,
								  n=n)

		if (use.LRT) { if (pchisq(q=(2 * (optList[[nrates]]$optlnL - optList[[nrates-1]]$optlnL)), df=n.hist.classes, lower.tail=FALSE) >= sig) flag=TRUE
		} else if (optList[[nrates]]$AICc > optList[[nrates - 1]]$AICc) flag=TRUE
					
		oldIntBreaks <- optList[[nrates]]$optBreaks
		nrates <- nrates + 1
	}
	optList
}	

getTaxaInOneLoc <- function(thisLoc) {
	sort(occs$taxon[occs$collection_no %in% locs$collection_no[locs$name==thisLoc]])
}

getTaxaInJanisLocs <- function() {
	locs <- read.csv("~/Dropbox/code/R/amandaTeeth/dat/janis_et_al_2000_Localities.csv")
	locNames <- unique(locs$name)
	locList <- lapply(locNames, getTaxaInOneLoc)
	names (locList) <- locNames
	locList
}	

testNewTaxaOneRep <- function (intSp, thisVec, test=c("ks", "mw")) {
	test <- match.arg(test)
	thisTest <- list()
	# thisTest[[length(intSp)]] <- list(prevDist=NA, newDist=NA, ks=NA)
	for (i in seq(from=length(intSp) - 1, to=1)) {
		prevDist = thisVec[intSp[[i]][intSp[[i]] %in% intSp[[i+1]]]]
		newDist = thisVec[intSp[[i]][!intSp[[i]] %in% intSp[[i+1]]]]
		if (length(prevDist) > 0 & length(newDist) > 0) { 
			if(test=="ks") { thisTest[[i]] <- list( prevDist = prevDist, newDist = newDist, ks = ks.test(prevDist, newDist))
			} else thisTest[[i]] <- list( prevDist = prevDist, newDist = newDist, ks = wilcox.test(prevDist, newDist))
		} else thisTest[[i]] <- list( prevDist = prevDist, newDist = newDist, ks = NA)
	}
	sigVec <- c(sapply(thisTest, function(x) if (!is.na(x$ks)) x$ks$p.value else NA), NA)
	# cbind(intervals, sigVec)
	# plot(rowMeans(intervals), sigVec, xlim=c(65, 34), type="o", lty=3)
	# text(rowMeans(intervals), sigVec, labels=rowMeans(intervals), pos=4, cex=0.3)
	# abline(h=0.05, lty=2, lwd=0.5)
	sigVec
}

getHypsodontyFromSample <- function(thisRep) {
	hyps <- read.csv("~/Dropbox/code/common_dat/hypsodonty.csv", header=TRUE)
	hyps$hypsodonty[hyps$hypsodonty=="unknown"] <- NA
	hyps$hypsodonty <- ordered(hyps$hypsodonty, levels=c("brachydont", "mesodont", "hypsodont"))
	thisHyps <- sapply(thisRep, FUN=function(x) hyps$hypsodonty[hyps$taxon %in% x])
	sapply(thisHyps, table)
}

getOneIntervalHistogram <- function(this.intv.sample, breaks) {
	hist(bm.mat$bodyMass[this.intv.sample], breaks=breaks, plot=FALSE)$counts
}

getOneRepCounts <- function(this.rep.sample, breaks) {
	sapply(this.rep.sample, getOneIntervalHistogram, breaks=breaks, simplify = "array")
}

getCountCube <- function(rep.intv.sp, breaks) {
	sapply(rep.intv.sp, getOneRepCounts, breaks, simplify="array")
}

getOneCategoryCountsFromIntvSample <- function(this.intv.sample, this.category=NULL, split.by=NULL) {
	this.intv.sample[split.by[this.intv.sample]==this.category]
}

getOneCategoryCountsFromRep <- function(this.rep.sample, this.category=NULL, split.by=NULL) {
	lapply(this.rep.sample, FUN=getOneCategoryCountsFromIntvSample, this.category=this.category, split.by=split.by)
}

getAllRepCountsForOneCategory <- function(this.category, split.by, repIntSp) {
	lapply(repIntSp, FUN=getOneCategoryCountsFromRep, this.category=this.category, split.by=split.by)
}

getAllRepCountsSplit <- function(repIntSp, split.by) {
	rep.intv.sp.split <- lapply(unique(split.by), FUN=getAllRepCountsForOneCategory, split.by=split.by, repIntSp=repIntSp)
	names(rep.intv.sp.split) <- unique(split.by)
	rep.intv.sp.split
}

extractOneRepCountSplitCube <- function(this.rep, this.cube) {
	lapply(this.cube, function(x, this.rep) x[,,this.rep], this.rep=this.rep)
}


#############################################################################################################
#Get frequency of combinations for breaks ecology analysis
BreakComboFreq <- function(optList = NULL) {
	optFrame <- matrix(nrow = 1000, ncol = 15)
	for(ii in seq(1, length(optList),1))
	{
		for(mm in seq(1,length(optList[[ii]][[length(optList[[ii]])]]$optBreaks),1))
		{
			optFrame[[ii,mm]] <- optList[[ii]][[length(optList[[ii]])]]$optBreaks[[mm]]
		}
	}
	return(count(optFrame))
}

#regimeHist is meant to parse through the repIntSp list object to print the distribution of traits
#within the regimes demarcated by breaks or designated by optList (if breaks are not provided).
#netFreq and regime Freq sets whether to use proportions or actually counts for the net gain/loss and regime 
#distribution histograms, respectively.
regimeHist <- function(repIntSp = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq = FALSE, netPlotType = "absolute", plot.together = FALSE) {
	#function to generate histograms of species within each regime
	#get list of intervals that comprise a regime
repIntSp_regimes <- repIntSp
	if (is.null(breaks)) 
		breaks <- optList_bm_allReps[[1]][[10]]$optBreaks

	regimeBM.list <- list()
	regimeSp.List <- list()

	for (ii in seq(1, length(repIntSp_regimes), 1)) {
		names(repIntSp_regimes[[ii]]) <- str_remove(names(repIntSp_regimes[[ii]]), " Ma")
	}

	regimeSp <- unique(unlist(repIntSp_regimes[[1]][which(as.double(names(repIntSp_regimes[[ii]])) > breaks[1])]))
	regimeSp <- regimeSp[regimeSp %in% thisMat$species]
	regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]

	regimeSp.List[[1]] <- regimeSp
	regimeBM.list[[1]] <- regimeBM

	rm(regimeSp, regimeBM)

	for (mm in seq(2, length(breaks), 1)) {
		#get regimes for remaining sections
		
		#seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
		#max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
regimeSp <- unique(unlist(repIntSp_regimes[[1]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm - 1] & as.double(names(repIntSp_regimes[[ii]])) > 
			breaks[mm])]))
		regimeSp <- regimeSp[regimeSp %in% thisMat$species]
		regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]

		regimeSp.List[[mm]] <- regimeSp
		regimeBM.list[[mm]] <- regimeBM

		rm(regimeSp, regimeBM)
	}

	regimeSp <- unique(unlist(repIntSp_regimes[[1]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[length(breaks)])]))
	regimeSp <- regimeSp[regimeSp %in% thisMat$species]
	regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]

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
	bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3, 4)
	breakCol <- rainbow(length(bmBreaks))

	hist.list <- list()
	quartz()
	par(mfrow = c(2, length(regimeBM.list)/2 + 1))
	for (jj in seq(1, length(regimeBM.list), 1)) {
		hist.list[[jj]] <- hist(regimeBM.list[[jj]], main = names(regimeBM.list[[jj]]), breaks = bmBreaks, las = 1, col = breakCol, freq = regimeFreq, 
			xlab = "log Body Mass (kg)", ylim = c(0, 1))
	}

	regimeNetChange <- hist.list
	regimeNetChange[[length(regimeNetChange)]] <- NULL
	names(regimeNetChange) <- breaks

	for (tt in seq(2, length(hist.list), 1)) {
		regimeNetChange[[tt - 1]]$counts <- hist.list[[tt]]$counts - hist.list[[tt - 1]]$counts
		tt <- tt + 1
	}

	if (netPlotType == "pos/neg") {
		#net change between regimes
		quartz()
		par(mfrow = c(2, length(hist.list)/2))
		for (hh in seq(1, length(regimeNetChange), 1)) {
			plot(regimeNetChange[[hh]], main = names(regimeNetChange[hh]), freq = netFreq, las = 1, col = breakCol, xlab = "log Body Mass (kg)") #####for some reason is identical to regular hist
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

	####Get a hsitogram for each entry of repIntSp to avoid missing unique species in each regime
	####and then take median across all of the histograms.	

	return()
}

#regimeHist_countbox is meant to parse through the taxcube or countcube objects (matrices) to print the distribution 
#of traits within the regimes demarcated by breaks or designated by optList (if breaks are not provided).
#netFreq and regime Freq sets whether to use proportions or actually counts for the net gain/loss and regime 
#distribution histograms, respectively.
regimeHist_countBox <- function(countBox = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
																netPlotType = "absolute", plot.together = FALSE, grayscale = TRUE, plot.axes = FALSE) {
	if(class(countBox) == "matrix"){
		if(is.null(breaks)) breaks <- optList_bm_allReps[[1]][[length(optList_bm_allReps)]]$optBreaks
		
		regimeBM.list <- list()
		regimeSp.List <- list()
		
		colnames(countBox) <- str_remove(colnames(countBox), " Ma")
		
		regimeBM <- countBox[,which(as.double(colnames(countBox)) > breaks[1])]
		regimeBM.list[[1]] <- apply(regimeBM, c(1), sum)
		rm(regimeBM)
		
		for(mm in seq(2,length(breaks),1))
		{
			#get regimes for remaining sections
			
			#seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
			#max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
			regimeBM <- countBox[,which(as.double(colnames(countBox)) < breaks[mm-1] 
																	& as.double(colnames(countBox)) > breaks[mm])]
			regimeBM.list[[mm]] <- apply(regimeBM, c(1), sum)
			
			rm(regimeBM)
		}  
		
		regimeBM <- countBox[,which(as.double(colnames(countBox)) < breaks[length(breaks)])]
		regimeBM.list[[length(breaks)+1]] <- apply(regimeBM, c(1),sum)
		
		rm(regimeBM)
		
		regimeNames <- vector()
		regimeNames[1] <- paste(">",max(breaks,sep=""))
		
		for(ii in seq(2,length(regimeBM.list)-1,1))
		{
			regimeNames[ii] <- paste(paste(breaks[ii-1]," to ",sep=""),breaks[ii],sep="")
		}
		regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
		names(regimeBM.list) <- regimeNames 
		bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
		
		btwnBreaks <- vector()
		for(ii in seq(1, length(bmBreaks)-1,1)) btwnBreaks[ii] <- (bmBreaks[ii]+bmBreaks[ii+1])/2 
		
		if(grayscale == TRUE) breakCol <- gray.colors(length(bmBreaks), start=0, end=1)
		else breakCol <- rainbow(length(bmBreaks))
		
		#make vector that conbtains the btwBreak values at qunatity listed by countBox
		for(hh in seq(1, length(regimeBM.list),1))
		{
			
			regimeBM.list[[hh]] <- c(rep(btwnBreaks[1],regimeBM.list[[hh]][1]),rep(btwnBreaks[2],regimeBM.list[[hh]][2]),
															 rep(btwnBreaks[3],regimeBM.list[[hh]][3]),rep(btwnBreaks[4],regimeBM.list[[hh]][4]),
															 rep(btwnBreaks[5],regimeBM.list[[hh]][5]),rep(btwnBreaks[6],regimeBM.list[[hh]][6]))
		}
		
		hist.list <- list()
		for(jj in seq(1, length(regimeBM.list),1))
		{
			
			hist.list[[jj]] <- hist(regimeBM.list[[jj]], breaks = bmBreaks, plot = FALSE)
			hist.list[[jj]]$density <- hist.list[[jj]]$counts/sum(hist.list[[jj]]$counts)
		}
		
		net.plot <- hist.list
		net.plot[[length(regimeNetChange)]] <- NULL
		names(net.plot) <- breaks
		
		for(tt in seq(2,length(hist.list),1))
		{
			net.plot[[tt-1]]$counts <- hist.list[[tt]]$counts - hist.list[[tt-1]]$counts
			#tt <- tt+1
		}
		
		max(sapply(net.plot, function(x) max(x$counts)))
		min(sapply(net.plot, function(x) min(x$counts)))
		#get limits for plot dimensions and scale
		hist.plot_ylimMax <- max(sapply(hist.list, function(x) max(x$density)))
		net.plot_ylimMax <- max(sapply(net.plot, function(x) max(x$counts)))
		net.plot_ylimMin <- min(sapply(net.plot, function(x) min(x$counts)))
		
		if(plot.axes == TRUE) marginPlot <- c(2,3,0,0)
		else marginPlot <- c(0,0.5,0,0)
		
		quartz(width = 11, height = 4)
		par(mfrow = c(1, length(hist.list)),mar=marginPlot)
		axesCheck <- FALSE
		for(ii in seq(1, length(hist.list),1))
		{
			if(plot.axes == TRUE)
			{
				if(ii == 1) axesCheck <- TRUE 
				else axesCheck <- FALSE
			}
			plot(hist.list[[ii]], col = breakCol, axes = axesCheck, ylim = c(0,1), main = NULL, freq = FALSE)
		}
		quartz(width = 11, height = 4)
		par(mfrow = c(1, length(net.plot)+1),mar=marginPlot) #have +1 to plot lenght so plots are same size as those for hist.plot
		for(ii in seq(1, length(net.plot),1))
		{
			if(plot.axes == TRUE)
			{
				if(ii == 1) axesCheck <- TRUE 
				else axesCheck <- FALSE
			}
			plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = axesCheck, ylim = c(net.plot_ylimMin,net.plot_ylimMax), 
					 main = NULL)
		}
	}
	
	return()
}
#############################################################################################################################
regimeHist_HistMedian<- function(repIntSp = NULL, breaks = NULL, optList, thisMat, netFreq = TRUE, regimeFreq=FALSE,
																 netPlotType = "absolute", plot.together = FALSE) {
	repIntSp_regimes <- repIntSp
	if(is.null(breaks)) print("Error: Input a vector of break dates")
	
	regimeBM.list <- list()
	regimeSp.List <- list()
	all.hist <- list()
	
	for(ii in seq(1, length(repIntSp_regimes),1))
	{
		names(repIntSp_regimes[[ii]]) <- str_remove(names(repIntSp_regimes[[ii]]), " Ma")
	}
	
	for (ii in seq(1, length(repIntSp_regimes),1)) {
		regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) > breaks[1])]))
		regimeSp <- regimeSp[regimeSp %in% thisMat$species]
		regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
		
		regimeSp.List[[1]] <- regimeSp
		regimeBM.list[[1]] <- regimeBM
		
		rm(regimeSp, regimeBM)
		
		for(mm in seq(2,length(breaks),1))
		{
			#get regimes for remaining sections
			
			#seq(min(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),
			#max(which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])),1))
			regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[mm-1] 
																														 & as.double(names(repIntSp_regimes[[ii]])) > breaks[mm])]))
			regimeSp <- regimeSp[regimeSp %in% thisMat$species]
			regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
			
			regimeSp.List[[mm]] <- regimeSp
			regimeBM.list[[mm]] <- regimeBM
			
			rm(regimeSp, regimeBM)
		}
		
		regimeSp <- unique(unlist(repIntSp_regimes[[ii]][which(as.double(names(repIntSp_regimes[[ii]])) < breaks[length(breaks)])]))
		regimeSp <- regimeSp[regimeSp %in% thisMat$species]
		regimeBM <- thisMat[thisMat$species %in% regimeSp, "bodyMass"]
		
		regimeSp.List[[length(breaks)+1]] <- regimeSp
		regimeBM.list[[length(breaks)+1]] <- regimeBM
		
		rm(regimeSp, regimeBM)
		
		regimeNames <- vector()
		regimeNames[1] <- paste(">",max(breaks,sep=""))
		for(nn in seq(2,length(regimeBM.list)-1,1))
		{
			regimeNames[nn] <- paste(paste(breaks[nn-1]," to ",sep=""),breaks[nn],sep="")
		}
		regimeNames[length(regimeBM.list)] <- paste(breaks[length(regimeBM.list)-1],">",sep="")
		names(regimeBM.list) <- regimeNames 
		bmBreaks <- c(0, 0.69897, 1.39794, 2.176091, 2.69897, 3.0, 4.0)
		breakCol <- rainbow(length(bmBreaks))
		
		hist.list <- list()
		#quartz()
		#par(mfrow=c(2,length(regimeBM.list)/2+1))
		for(jj in seq(1, length(regimeBM.list),1))
		{
			hist.list[[jj]] <- hist(regimeBM.list[[jj]], main=names(regimeBM.list[[jj]]), breaks = bmBreaks, 
															las=1, col = breakCol, freq = regimeFreq, xlab="log Body Mass (kg)", plot = FALSE)
		}
		all.hist[[ii]] <- hist.list 
	}
	
	#get median across all hist counts and test to get median density
	median.density <- list()
	median.counts <- list()
	
	for(rr in seq(1, length(all.hist[[1]]),1)) #to go through each regime
	{
		regimeDensityMat <- matrix(nrow = length(all.hist),ncol= length(bmBreaks)-1)
		regimeCountMat <- matrix(nrow = length(all.hist),ncol= length(bmBreaks)-1)
		
		for(gg in seq(1, length(all.hist),1))
		{
			regimeDensityMat[gg,] <- all.hist[[gg]][[rr]]$density 
			regimeCountMat[gg,] <- all.hist[[gg]][[rr]]$counts
		}
		
		median.density[[rr]] <- apply(regimeDensityMat,2,median)
		median.counts[[rr]] <- apply(regimeCountMat,2,median)
	}
	hist.plot <- hist.list
	
	#overwite with median density and counts
	quartz(width = 11, height = 4)
	par(mfrow = c(1, length(hist.plot)),mar=c(0,0.5,0,0))
	for(pp in seq(1, length(hist.plot),1))
	{
		hist.plot[[pp]]$density <- median.density[[pp]]
		hist.plot[[pp]]$counts <- median.counts[[pp]]
		#plot(hist.plot[[pp]], col = breakCol, axes = FALSE)
	}
	
	#get net change between regimes
	net.plot <- hist.plot; net.plot[[length(hist.plot)]] <- NULL
	for(ii in seq(1, length(net.plot),1))
	{
		net.plot[[ii]]$counts <- hist.plot[[ii+1]]$counts-hist.plot[[ii]]$counts
		#plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = FALSE)
	}
	max(sapply(net.plot, function(x) max(x$counts)))
	min(sapply(net.plot, function(x) min(x$counts)))
	#get limits for plot dimensions and scale
	hist.plot_ylimMax <- max(sapply(hist.plot, function(x) max(x$density)))
	net.plot_ylimMax <- max(sapply(net.plot, function(x) max(x$counts)))
	net.plot_ylimMin <- min(sapply(net.plot, function(x) min(x$counts)))
	
	#plot both the hist and net changes
	quartz(width = 11, height = 4)
	par(mfrow = c(1, length(hist.plot)),mar=c(0,0.5,0,0))
	for(ii in seq(1, length(hist.plot),1))
	{
		plot(hist.plot[[ii]], col = breakCol, axes = TRUE, ylim = c(0,hist.plot_ylimMax), main = NULL)
	}
	quartz(width = 11, height = 4)
	par(mfrow = c(1, length(net.plot)),mar=c(0,0.5,0,0))
	for(ii in seq(1, length(net.plot),1))
	{
		plot(net.plot[[ii]], col=breakCol, freq=TRUE, axes = TRUE, ylim = c(net.plot_ylimMin,net.plot_ylimMax), 
				 main = NULL)
	}
	return() 
}

################################################################################################################
###Determine if range-through is working properly and then fill in the gaps if it isnt 
rangeThrough_alt <- function(IntMat) {
	#IntMat <- repIntCheck
	#IntMat.rangeThrough <- IntMat
	for(xx in seq_len(nrow(IntMat)))
	{
		int.index <- which(IntMat[xx,] == 1)
		int.new <- min(int.index):max(int.index)
		for (yy in int.new)
		{
			IntMat[xx,yy] <- 1
		}
	}
	
	repIntRangeThrough <- list()
	
	for(gg in seq_len(ncol(IntMat)))
	{
		specVec <- vector()
		for(hh in seq_len(nrow(IntMat)))
		{
			if(IntMat[hh,gg] == 1) specVec <- append(specVec, rownames(IntMat)[hh])
		}
		repIntRangeThrough[[gg]] <- specVec
	}
	names(repIntRangeThrough) <- colnames(IntMat)
	return(repIntRangeThrough)
}

checkRangeThrough <- function(repIntTest, ints) {
	#function to find if rangethrough is working
	int.label<- vector()
	species <- list()
	for(xx in seq(1, length(ints)-1,1)) 
	{
		species[[xx]] <- repIntTest[[xx]] #species[[xx-(min(ints)-1)]] 
		int.label[xx] <- (ints[xx]+ints[xx+1])/2
	}
	uni.sp <- unique(unlist(species))
	
	checkMat <- matrix(nrow = length(uni.sp), ncol = length(ints)-1)
	rownames(checkMat) <- uni.sp
	colnames(checkMat) <- int.label
	
	# find way to parse through each species and mark if they are in a given interval
	for(xx in seq(1, nrow(checkMat),1))
	{
		for(yy in seq(1, ncol(checkMat),1))
		{
			#for(zz in repIntTest[[yy+49]])
			##{
			test.true <- which(rownames(checkMat)[xx] == repIntTest[[yy]]) #repIntTest[[yy+(min(ints)-1)]])
			if(length(test.true) == 0)
			{
				checkMat[xx,yy] <- 0
			}
			if(length(test.true) >= 1) checkMat[xx,yy] <- 1
			#}
		}
	}
	return(checkMat)
}
################################################################################################################

