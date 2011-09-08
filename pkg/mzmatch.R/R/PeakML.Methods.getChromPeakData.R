PeakML.Methods.getChromPeakData <- function(filename, ionisation = "detect"){
	# Return the peak and chrom data from the file specified
	# PRE: peakML filename
	# POST:
	#	list(peakDataMtx, chromDataList)
	
	
	# Create a new java javaProject and allocate memory to it and load the PeakML file into it.
	javaProject <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory. Please wait. \n")
	st <- system.time(.jcall(javaProject, returnSig="V", method="load",filename))
	cat ("Loaded peakML file into memory in:",st[3],"s \n")
	
	massCorrection <- PeakML.Methods.getMassCorrection(javaProject, ionisation)
	peakDataMtx <- PeakML.Methods.getPeakData(javaProject, massCorrection)
	chromDataList <- PeakML.Methods.getChromData(javaProject, peakDataMtx, massCorrection)
	
	rv <- vector("list",4)
	rv[[1]] <- peakDataMtx
	rv[[2]] <- chromDataList
	rv[[3]] <- .jcall(javaProject,returnSig="[S", method="getSetNames") # sample classes 
	rv[[4]] <- .jcall(javaProject,returnSig="[S", method="getMeasurementNames") # measurement names
	rv[[5]] <- massCorrection
	
	rv
}
