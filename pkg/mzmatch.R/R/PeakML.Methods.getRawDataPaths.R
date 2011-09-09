PeakML.Methods.getRawDataPaths <- function(javaProject, rawPath = NULL){
	# PRE: 
	#	javaProject, path of raw files if different from the ones referened in the peakml files
	# POST:
	#	Return the samples names and the full path to where the .mzXML files are located
	sampleNames <- .jcall(javaProject, returnSig="[S", method="getMeasurementNames")
	
	if (!is.null(rawPath)){
		dirContent <- dir(rawPath, recursive = TRUE)
		dirContent <- sub(".mzXML", "", dirContent)
		mzFiles <- which(dirContent%in%sampleNames)
	  	if(length(dirContent[dirContent%in%sampleNames])==length(sampleNames)){
			rawDataPaths <- dir(rawPath, full.names=TRUE, recursive=TRUE)[mzFiles]
			cat(paste("Raw data file located at: ", rawDataPaths, "\n", sep=""))
		} else{
			cat(paste("Raw data file: ",sampleNames,"cant't be found in folder: ", rawPath, "\n", sep=""))
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
