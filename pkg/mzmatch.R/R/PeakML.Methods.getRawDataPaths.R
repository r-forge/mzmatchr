PeakML.Methods.getRawDataPaths <- function(javaProject, Rawpath = NULL){
	# PRE: 
	#	javaProject, path of raw files if different from the ones referened in the peakml files
	# POST:
	#	Return the samples names and the full path to where the .mzXML files are located
	sampleNames <- .jcall(javaProject, returnSig="[S", method="getMeasurementNames")
	
	if (!is.null(Rawpath)){
		dirContent <- dir(Rawpath, recursive = TRUE, full.names=TRUE)
		fileid <- rep(NA,length(sampleNames))
		for (filenum in 1:length(sampleNames))
		{
			fileid[filenum] <- grep (sampleNames[filenum],dirContent)
		}
	  	if(length(which(is.na(fileid)))==0){
			rawDataPaths <- dirContent[fileid]
			cat(paste("Raw data file located at: ", rawDataPaths, "\n", sep=""))
		} else{
			cat(paste("Raw data file: ",sampleNames[which(is.na(fileid))],"cant't be found in folder: ", Rawpath, "\n", sep=""))
			rawDataPaths <- NULL
		}
	} else{
		rawDataPaths <- .jcall(javaProject, returnSig="[S", method="getFileNames")
		mzFiles <- which(file.exists(rawDataPaths)==TRUE)
		if (length(mzFiles) == length(rawDataPaths)){
			cat(paste("Raw data file located at: ", rawDataPaths, "\n", sep=""))
		} else {
			cat(paste("Raw data file can't be read from location': ", rawDataPaths, "\n", sep=""))
			rawDataPaths <- NULL
		}
	}
	list(sampleNames,rawDataPaths)
}
