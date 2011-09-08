PeakML.Methods.getChromData <- function(javaProject, peakDataMtx, massCorrection)
{	
	# Generates a chromatorgram data list from peakDataMtx
	# PRE:	javaProject: the javaProject objects
	# 	peakDataMtx: the peak data matrix
	# POST:	Chromatogram data list
	
	## Create an empty list of length nrow(peakDataMtx)

	chromDataList <- vector("list",nrow(peakDataMtx))
	
	for (chrnum in 1:nrow(peakDataMtx))
	{
		retentiontimes <- .jcall(javaProject, returnSig="[D", method="getRetentionTimes", as.integer(chrnum-1))
		intensities <- .jcall(javaProject, returnSig="[D", method="getIntensities", as.integer(chrnum-1))
		masses <- .jcall(javaProject, returnSig="[D", method="getMasses", as.integer(chrnum-1))+(massCorrection)  ## We have to correct masses for ionisation mode
		scanids <- .jcall(javaProject, returnSig="[I", method="getScanIDs", as.integer(chrnum-1))
		chromDataList[[chrnum]] <- rbind(masses,intensities,retentiontimes,scanids)
	}

	chromDataList
}
