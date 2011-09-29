PeakML.Read <- function(filename, ionisation = "detect", Rawpath=NULL){
	# Return the peak and chrom data from the file specified
	# PRE: peakML filename
	# POST:
	#	list(peakDataMtx, chromDataList)
	
	
	# Create a new java javaProject and allocate memory to it and load the PeakML file into it.
	sampleLookUP <- function (x)
	{
		phenoData <- sampleClasses[sampleGroups[peakDataMtx[,9]==x][1]]
		phenoData
	}
	
	getRTScans <- function(project, sampleNames){
		corRT <- vector ("list",length(sampleNames))
		rawRT <- vector ("list",length(sampleNames))
		for (measID in 1:length(sampleNames)){
			corRT[[measID]] <- .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes",as.integer(measID-1))
			temp <- rep(NA, length(corRT[[measID]]))
			for (scan in 1:length(temp)){
				temp[scan] <- as.numeric(.jcall(project, returnSig="S", method="getScanAnnotation", as.integer(measID-1), as.integer(scan-1), as.character("RT_raw")))
			}
			rawRT[[measID]] <- temp
		}
		list(corRT, rawRT)
	}

	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
	st <- system.time(.jcall(project, returnSig="V", method="load",filename))
	cat ("Done in:",st[3],"s \n")
	
	massCorrection <- PeakML.Methods.getMassCorrection(project, ionisation)
	samPath <- PeakML.Methods.getRawDataPaths(project, Rawpath)
	sampleNames <- samPath[[1]]
	rawDataFullPaths <- samPath[[2]]
	peakDataMtx <- PeakML.Methods.getPeakData(project, massCorrection)
	chromDataList <- PeakML.Methods.getChromData(project, peakDataMtx, massCorrection)
	rtScanList <- getRTScans(project, sampleNames)
	sampleGroups <- peakDataMtx[,11]
	sampleClasses <- .jcall(project,returnSig="[S", method="getSetNames")
	phenoData <- unlist(lapply(1:length(sampleNames),sampleLookUP))

	rv = list()
	rv$peakDataMtx <- peakDataMtx
	rv$chromDataList <- chromDataList
	rv$sampleClasses <- sampleClasses
	rv$sampleNames <- sampleNames
	rv$massCorrection <- massCorrection
	rv$sampleGroups <- sampleGroups
	rv$phenoData <- phenoData
	rv$correctedRTList <- rtScanList[[1]]
	rv$rawRTList <- rtScanList[[2]]
	rv$fileName <- filename
	rv$rawDataFullPaths <- rawDataFullPaths
	rv
}
