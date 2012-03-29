PeakML.Isotope.plotSamples <- function(isotopeChroms, metName, metFormula, metMass, stdRT, sampleType, sampleGroups, plotOrder, useArea, followCarbon){

	processRatioMtx <- function(ratioMtx){
		plotMtx <- matrix(nrow=nrow(ratioMtx), ncol=ncol(ratioMtx))
		rownames(plotMtx) <- row.names(ratioMtx)
		for (r in 1:nrow(ratioMtx)){
			sumr <- sum(ratioMtx[r,])
			for (c in 1:ncol(ratioMtx)){
				if (!is.na(ratioMtx[r,c])) {
					plotMtx[r,c] <- ratioMtx[r,c]/sumr
				} else {
					plotMtx[r,c] <- 0
				}
			}
		}
		plotMtx
	}
	
	mzList <- isotopeChroms[[1]]
	intList <- isotopeChroms[[2]]
	rtList <- isotopeChroms[[3]]
	
	numPeakGroups <- length(mzList)
	numSampleGroups <- length(mzList[[1]])
	numCarbons <- length(mzList[[1]][[1]])
	#numReplicates <- length(mzList[[1]][[1]][[1]])
	
	metName <- unlist(strsplit(as.character(metName), ", "))
	
	fillLabels <- c("UL" , paste(as.character("+"),c(1:numCarbons),sep=""))# create a list of isotop names
	fillColor <- c(1:numCarbons)
	if (numCarbons > 7) fillColor <- c("black", rainbow(numCarbons, alpha=0.5))
	
	for (peakGroup in 1:numPeakGroups){
		trendList <- PeakML.Isotope.getTrendList (intList, sampleGroups, useArea) [[peakGroup]]
		ratioMtx <- processRatioMtx (t(PeakML.Isotope.getRatioMtxList (intList, sampleGroups, useArea, metName[1]) [[peakGroup]]))
		
		par (mar=c(0,0,0,0))
		plot (c(1:10), c(1:10), xlab="", ylab="", pch="", axes=F) 
		
		if(length(metName)==1){
			text(1.5, 6, metName[1], cex = 2, pos=4)
			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, "Ion:", sampleType, sep = "  "), pos=4, cex = 1.5)
#			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, sep = "  "), pos=4, cex = 1.5)
		} else {
			text(1.5, 7, metName[1], cex = 2, pos=4)
			text(1.5, 5, paste(metName[2:length(metName)], collapse=","), cex = .5, pos=4)
			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, "Ion:", sampleType, sep = "  "), pos=4, cex = 1)
#			text(1.5, 4, paste("Formula:", metFormula, "Mass:", round(metMass,3), "Std.RT:", stdRT, sep = "  "), pos=4, cex = 1)
		}
		legend(1.5, 3,fill=fillColor, fillLabels[1:numCarbons], bty="n", horiz=TRUE)
		
		plot (1, 1, xlab="", ylab="", pch="", axes=F)
		text(1, 1, paste("G", as.character(peakGroup), sep=""), cex = 3)
	
		for (item in 1:length(plotOrder)){
			if(plotOrder[item] == "LEGEND"){

				par(mar=c(0.5,0.5,0.5,0.5))
				plot (1, 1, xlab="", ylab="", pch="", axes=F)
				legend("topleft",fill=fillColor, fillLabels[1:numCarbons], bty="n")

			} else if (plotOrder[item] == "TREND"){
			
				errBar <- function(x, y, upper, lower=upper, length=0.1,...){
					if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
					stop("vectors must be same length")
					arrows(x, y+upper, x, y-lower, angle=90, code=3, length=length, ...)
				}
				
				trendMtx <- PeakML.Isotope.getTrendMtx(trendList, sampleGroups)
				trendMtx.sum <- apply(trendMtx, 2, sum)
				trendMtx.sd <- apply(trendMtx, 2, sd)
				ylimit <- sum(max(trendMtx.sum), max(trendMtx.sd))
				# calculate std err
				#errMtx <- replace(trendMtx, trendMtx==0, NA)
				#stdErr <- apply(errMtx, 2, PeakML.Methods.getStdErr)

				par(mar=c(4,2.5,0.5,3))
				if (useArea == FALSE){
					ylabel <- paste(metName[1], "(mean peak height)", sep= " ")
				} else {
					ylabel <- paste(metName[1], "(mean peak area)", sep= " ")
				}

				barx <- barplot(trendMtx, beside=FALSE, col=fillColor, ylab=ylabel,ylim=c(0,ylimit), border=NA, axisnames=FALSE)
				axis(1, las=2, at=c(1:length(sampleGroups)), lab=sampleGroups, lwd=0, cex.axis=.8)
				errBar(barx, trendMtx.sum, trendMtx.sd/numCarbons, lwd=.3)

			} else if (plotOrder[item] == "LABELLED"){
			
				if(is.na(followCarbon)) followCarbon <- 2

				if (followCarbon == 2){
					tMtx <-PeakML.Isotope.getFCMtxAbun(trendList, sampleGroups, followCarbon)
					fcMtx <- tMtx[[1]]
					colvector <- tMtx[[2]]
					#print (fcMtx)
					#print(colvector)
				} else {
					fcMtx <- PeakML.Isotope.getFCMtx(trendList, sampleGroups, followCarbon)
				}
				
				if (useArea == FALSE){
					ylabel <- paste(metName[1], "(mean peak height)", sep= " ")
				} else {
					ylabel <- paste(metName[1], "(mean peak area)", sep= " ")
				}

				ylimit <- max(apply(fcMtx,2,sum))
				#cat("------------------")
				#print(ylimit)
				if (!ylimit==0){
					par(mar=c(4,2.5,0.5,3))

					if (followCarbon==2){
						overlap <- "#993399" # purple
						measuredLarger <- "#FF9988" # red > expected
						expectedLarger <- "#8877FF" # blue
						
						
						COLORS <- c(overlap, measuredLarger, expectedLarger)
						colvector[colvector==1] <- COLORS[1]
						colvector[colvector==2] <- COLORS[2]
						colvector[colvector==3] <- COLORS[3]
						#print(colvector)
					} else {
						colvector <- c(followCarbon)
					}
					
					barplot(fcMtx, beside=FALSE, col=colvector, axisnames=FALSE, ylab=ylabel, ylim=c(0,ylimit), border=NA)
					axis(1, las=3, at=c(1:length(sampleGroups)), lab=sampleGroups, lwd=0, cex.axis=.8)
					
					if (followCarbon==2){
						legend("topright",fill=COLORS, c("Overlap", "<Expected", ">Expected"), bty="n")
					}
				} else {
					par(mar=c(0.5,0.5,0.5,0.5))
					plot (1, 1, xlab="", ylab="", pch="", axes=F)
				}

			} else if (plotOrder[item] == "RATIO") {

				par(mar=c(4,4,3,2))
				barplot(t(ratioMtx), beside=FALSE, col=fillColor, ylab="% area under peak", ylim = c(0,1), border=NA, axisnames=FALSE)
				axis(1, las=3, at=c(1:length(row.names(ratioMtx))), lab=row.names(ratioMtx), lwd=0, cex.axis=.5)
				
			} else if (plotOrder[item] == "EMPTY"){

				# This is a an empty plot for flexibility e.g. define empty if a plot is not needed
				plot (1, 1, xlab="", ylab="", pch="", axes=F)
			} else{
				if (! length(unlist(rtList[[peakGroup]][[item]])) == 0){
					maxRT <- max(unlist(rtList[[peakGroup]][[item]]),na.rm=TRUE)
					minRT <- min(unlist(rtList[[peakGroup]][[item]]),na.rm=TRUE)
					maxIN <- max(unlist(intList[[peakGroup]][[item]]),na.rm=TRUE)
					minIN <- 0
				}else{
					maxRT <- 0
					minRT <- 0
					maxIN <- 0
					minIN <- 0
				}

				par(mar=c(4,4,0.5,0.5), mgp=c(1.5,0.5,0))
				plot (1, 1, pch="", xlab="", ylab="", xlim=c(minRT,maxRT), ylim=c(minIN,maxIN), cex.axis=0.75)
				for (isotop in 1:numCarbons){
					numReplicates <- length(rtList[[peakGroup]][[item]][[isotop]])
					for (rep in 1:numReplicates){
						points(rtList[[peakGroup]][[item]][[isotop]][[rep]], intList[[peakGroup]][[item]][[isotop]][[rep]], type="l",col=fillColor[[isotop]])
					}
				}
				legend("topright", plotOrder[item], bty="n", cex=.5)
			}
		}
	}
}
