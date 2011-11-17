PeakML.Isotope.getIsotopes <- function(peakDataMtx, mzXMLSrc, sampleNames, numCarbons, metMass, ppm, massCorrection,  stdRT=NA, stdRTWindow=NULL, fillGaps="ALLPEAKS"){
	#
	# Returns a list that has all isotopes of the given mass in the form isotope[[peak_group]][[isotope]][[sample]] <- peak id from peak data matrix
	# PRE: 
	#	peakDataMtx <- the original peakDataMtx
	# 	metData <- metName, metFormula, numCarbons, metMass, ppm, massWindow, stdRT, stdRTWindow, sampleType
	# POST: 
	#	returns a list all isotops like for each gid all isotops

	massWindow <- PeakML.Methods.getPPMWindow(metMass, ppm)

	stdRTWin <- NULL
	massFilterHits <- c()
	finalList<-NULL

	# get the mass filter hits from the peakML data within the range of the mass specified
	if (is.null(stdRTWindow)){
		massFilterHits <- which(peakDataMtx[,1]>=massWindow[[1]] & peakDataMtx[,1]<=massWindow[[2]])
	} else{
		if (!is.na(stdRT)){
			stdRTWindow <- PeakML.Methods.getRTWindow(stdRT, stdRTWindow)
			massFilterHits <- which(peakDataMtx[,1]>=massWindow[[1]] & peakDataMtx[,1]<=massWindow[[2]] & peakDataMtx[,4]>=stdRTWindow[[1]] & peakDataMtx[,4]<=stdRTWindow[[2]])
		} else {
			cat ("\tYou have specified a standard RT window, however no standard RT was provided. Defaulting to no standard RT filter mode.\n")
			stdRTWindow <- NULL
			massFilterHits <- which(peakDataMtx[,1]>=massWindow[[1]] & peakDataMtx[,1]<=massWindow[[2]])
		}
	}
	if (is.null(stdRTWindow) & length(massFilterHits==0)){
		stdRTWindow <- NULL
		massFilterHits <- which(peakDataMtx[,1]>=massWindow[[1]] & peakDataMtx[,1]<=massWindow[[2]])
	}
	
	pdMtxMassFiltered <- rbind(peakDataMtx[massFilterHits,])
	pdMtxMassFiltered <- cbind(pdMtxMassFiltered, massFilterHits) # append the original massFilterHits to keep track of the id
	colnames(pdMtxMassFiltered)[12] <- "" # remove the colname from the matrix to keep it all similar
	
	uniqueGroups <- unique(pdMtxMassFiltered[,10]) 		# to arrage the final list in the form of groups -> samples -> peakdata

	if (is.null(stdRTWindow)){
		if (length(uniqueGroups)==0){
			return (NULL)
		} else {
			finalList <- vector("list",length(uniqueGroups))	# create a list that has the length of the uniqueGroups	
		}
	} else {
		if (length(uniqueGroups)==0){
			finalList <- vector("list",1)
		} else {
			finalList <- vector("list",length(uniqueGroups))	# create a list that has the length of the uniqueGroups	
		}
	}

	carbonsList <- vector("list", numCarbons)		# numCarbons+1 as we need to store the unlabelled as well
	for (a in 1:length(finalList)){
		for (b in 1:length(carbonsList)){
			carbonsList[[b]] <- vector("list", length(unique(peakDataMtx[,9])))   # For each group allocate a list whose size is of uni. samples
		}
		finalList[[a]] <- carbonsList
	}

	for (gid in 1:length(finalList)){
		cat(".")
		
		possibleSamples <- unique(peakDataMtx[,9])
		
		pdMtxMassFilteredGroup <- NULL
		selGroup <- which(pdMtxMassFiltered[,10]==uniqueGroups[gid])	# select all hits in the same group
		pdMtxMassFilteredGroup <- rbind(pdMtxMassFilteredGroup, pdMtxMassFiltered[selGroup,]) # This matrix now contains the unlabelled peakset

		if (is.null(stdRTWindow)){
			if (length(pdMtxMassFilteredGroup)!=0){
				massAve <- mean(pdMtxMassFilteredGroup[,1]) # Uses the mean of average mass
				rtWindow <- c(min(pdMtxMassFilteredGroup[,5]), max(pdMtxMassFilteredGroup[,6]))			# Uses the min and max RT
			}
		} else {
			if (length(pdMtxMassFilteredGroup)!=0){
				massAve <- mean(pdMtxMassFilteredGroup[,1]) # Uses the mean of average mass
				rtWindow <- c(min(pdMtxMassFilteredGroup[,5]), max(pdMtxMassFilteredGroup[,6]))			# Uses the min and max RT
				rtWindow <- stdRTWindow			# Uses the min and max RT
			} else {
				massAve <- metMass
				rtWindow <- stdRTWindow
			}
		}

		# STORING THE UNLABELLED
		unlabledList <- vector("list", length(unique(peakDataMtx[,9])))
		if (length(pdMtxMassFilteredGroup)!=0){
			for (rown in 1:nrow(pdMtxMassFilteredGroup)){
				if (fillGaps=="ALLPEAKS") {
					sampleName <- sampleNames[pdMtxMassFilteredGroup[rown,9]]
					aveMass <- pdMtxMassFilteredGroup[rown,1]
					mWindow <- PeakML.Methods.getPPMWindow(aveMass, ppm)
					gapFillHits <- PeakML.Methods.getRawSignals(mzXMLSrc, sampleName, rtWindow, mWindow, massCorrection)
					if (max(gapFillHits[2,] != -1)){
						unlabledList[[pdMtxMassFilteredGroup[rown,9]]] <- list("gapfilled", gapFillHits)
					}
				} else {
					unlabledList[[pdMtxMassFilteredGroup[rown,9]]] <- pdMtxMassFilteredGroup[rown,12]  
				}
			}
		}
		finalList[[gid]][[1]] <- unlabledList # stores peakid of unlabelled, in each selected group sample
	
		# STORING THE LABELLED
		for (carbon in 1:(length(carbonsList)-1)){
			massWindow <- PeakML.Isotope.getMassWindow(massAve, carbon, ppm)		# Window of isotope masses

			filterHits <- c()
			if (!is.null(stdRTWindow)){
				filterHits <- PeakML.Isotope.filterPeaks(peakDataMtx, stdRTWindow, massWindow)
				if (length(filterHits)!=0){
					rtWindow <- stdRTWindow # for gap filler to use base rt window if stdRT window could not find any hits
				}
			}
			if (length(filterHits)==0){
				filterHits <- PeakML.Isotope.filterPeaks(peakDataMtx, rtWindow, massWindow)
			}
			
			if (fillGaps=="NONE"){
				if (length(filterHits) != 0) {
					for (fh in 1:length(filterHits)){
						finalList[[gid]][[carbon+1]][[peakDataMtx[filterHits[[fh]],9]]] <- filterHits[[fh]] # peakData[[9]] is the measurement id
					}
				}
			}else if(fillGaps=="MISSINGPEAKS"){
				if (length(filterHits != 0)){
					storedSamples <- c()
					for (fh in 1:length(filterHits)) {
						storedSamples <- c(storedSamples, peakDataMtx[filterHits[[fh]],9]) # samples for which peaks are exisiting
						finalList[[gid]][[carbon+1]][[peakDataMtx[filterHits[[fh]],9]]] <- filterHits[[fh]] # storing index of peaks that are existing
					}
					emptySamples <- possibleSamples[-c(storedSamples)]
				} else {
					emptySamples <- possibleSamples
				}

				for (es in emptySamples){
					sampleName <- sampleNames[es]
					gapFillHits <- PeakML.Methods.getRawSignals(mzXMLSrc, sampleName, rtWindow, massWindow, massCorrection)
					if (max(gapFillHits[2,] != -1)){
						finalList[[gid]][[carbon+1]][[es]] <- list("gapfilled", gapFillHits)
					}
				}
			} else if(fillGaps=='ALLPEAKS'){
				for (es in possibleSamples){
					sampleName <- sampleNames[es]
					gapFillHits <- PeakML.Methods.getRawSignals(mzXMLSrc, sampleName, rtWindow, massWindow, massCorrection)
					if (max(gapFillHits[2,] != -1)){
						finalList[[gid]][[carbon+1]][[es]] <- list("gapfilled", gapFillHits)
					}
				}
			}
		}
		
		if (fillGaps=="ALLPEAKS" & !is.null(stdRTWindow)){
			if (is.infinite(max(unlist(unlabledList)))){
				massWindow <- PeakML.Methods.getPPMWindow(massAve, ppm)
				for (es in possibleSamples){
					sampleName <- sampleNames[es]
					gapFillHits <- PeakML.Methods.getRawSignals(mzXMLSrc, sampleName, rtWindow, massWindow, massCorrection)
					if (max(gapFillHits[2,] != -1)){
						gapFillHits[2,] <- PeakML.Methods.baseCorrection (gapFillHits[2,])
						unlabledList[[es]] <- list("gapfilled", gapFillHits)
					}
				}
				finalList[[gid]][[1]] <- unlabledList
			}
		}
	}
	finalList
}
