PeakML.Write <- function(peakDataMtx, chromDataList, sampleNames, rawdataFullPaths, phenoData, corRT, rawRT, outFileName, ionisation="neutral",GroupAnnotations=NULL){
	# PRE:
	#	peakDataMtx <- PeakML.Methods.getPeakData$peakDataMtx
	#	chromDataList <- PeaKML.Methods.getChromData or PeakML.Methods.getPeakData$chromDataList
	#	sampleName <- sample names; PeakML.Methods.getPeakData$sampleNames
	#	rawDataFullPaths <- full paths to the raw data; PeakML.Methods.getPeakData$rawDataFullPaths
	#	phenoData <- the pheno data PeakML.Methods.getPeakData$phenoData
	#	corRT <- the corrected retention time PeakML.Methods.getPeakData$correctedRTList
	#	rawRT <- the raw retention time; PeakML.Methods.getPeakData$rawRTList
	#	outFileName <- the output file name 
	#	ionisation <- the ionisation, leave it as default
	# POST:
	#	write the peakml the the file specified
	
	project <- .jnew("peakml/util/rjava/Project", sampleNames, rawdataFullPaths, phenoData)

	for (measID in 1:length(sampleNames)){
		for (scanID in 1:length(rawRT[[measID]])){
			.jcall(project, returnSig="V", method="addScanInfo", as.integer(measID-1), as.numeric(corRT[[measID]][scanID]), as.character(ionisation))
			.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(measID-1), as.integer(scanID-1), as.character("RT_raw"), as.character(rawRT[[measID]][scanID]))
		}
	}

	for (peakID in 1:length(chromDataList)){
		chrom <- chromDataList[[peakID]]
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(peakDataMtx[peakID,9]-1), as.integer(chrom[4,]), as.numeric(chrom[3,]),as.numeric(chrom[1,]), as.numeric(chrom[2,]), as.character(ionisation))
	}
	
	setIndexes <- vector("list",length(unique(peakDataMtx[,10])))
	for (sid in 1:length(setIndexes)){
		setIndexes[[sid]] <- which(peakDataMtx[,10]==sid)
	}
	
	for (sid in 1:length(setIndexes)){
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(setIndexes[[sid]]-1))
	}

	if (!is.null(GroupAnnotations))
	{
		PeakML.Methods.writeGroupAnnotations (project, GroupAnnotations)
	}
	.jcall(project, returnSig="V", method="write", outFileName)
}
