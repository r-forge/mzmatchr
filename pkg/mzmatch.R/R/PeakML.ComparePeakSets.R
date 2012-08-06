PeakML.ComparePeakSets <- function(standard_filename, filename, stdionisation="detect",ionisation="detect", Rawpath=NULL, outputfile,ppm=5,rtwin=20)
{
	st <- system.time (stdPeakMLdata <- PeakML.Read (standard_filename,stdionisation,Rawpath=NULL))
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))

	rawdatafullpaths <- PeakMLdata$rawDataFullPaths
	
	if (is.null(rawdatafullpaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	stdPeakMLDataTable <- PeakML.Methods.getCompleteTable (stdPeakMLdata,sumintensity=FALSE)
	stdpeakMasses <- apply(stdPeakMLDataTable[[2]],2,median,na.rm=TRUE)-stdPeakMLdata$massCorrection[[1]]
	stdpeakRTs <- apply(stdPeakMLDataTable[[3]],2,median,na.rm=TRUE)

	PeakMLDataTable <- PeakML.Methods.getCompleteTable (PeakMLdata,sumintensity=FALSE)
	peakMasses <- apply(PeakMLDataTable[[2]],2,median,na.rm=TRUE)-PeakMLdata$massCorrection[[1]]
	peakRTs <- apply(PeakMLDataTable[[3]],2,median,na.rm=TRUE)

	selectedSets <- NULL
	stdMatchedSets <- NULL
	for (setnum in 1:length(stdpeakMasses))
	{
		setMass <- unlist(PeakML.Methods.getPPMWindow(stdpeakMasses[setnum],ppm))
		setRT <- c(stdpeakRTs[setnum]-rtwin,stdpeakRTs[setnum]+rtwin)

		hits <- which(peakMasses>=setMass[1] & peakMasses<=setMass[2] & peakRTs>=setRT[1] & peakRTs<=setRT[2])

		if (length(hits)!=0)
		{
			selectedSets <- append(selectedSets,hits)
			stdMatchedSets <- append(stdMatchedSets,setnum)
		}
	}
	selectedSets <- unique(selectedSets)
	stdMatchedSets <- unique(stdMatchedSets)

	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=selectedSets)

	## Write no matching
	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("Not_matched_",outputfile,sep=""), sets=c(1:length(peakRTs))[-c(selectedSets)])

	# Not macthed from std file
	PeakML.Methods.extractPeakGroups (PeakMLData=stdPeakMLdata, outputfile=paste("STD_not_matched_",standard_filename,sep=""), sets=c(1:length(stdpeakRTs))[-c(stdMatchedSets)])
}
