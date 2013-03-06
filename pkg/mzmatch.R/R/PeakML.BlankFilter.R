PeakML.BlankFilter <- function (filename, ionisation="detect", Rawpath=NULL,outputfile,BlankSample=NULL)
{
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata)
	Masses <- apply(PeakTable[[2]],2,median,na.rm=TRUE)
	RTs <- apply(PeakTable[[3]],2,median,na.rm=TRUE)

	## Check which samples are blanks
	blanksamples <- which (PeakMLdata$phenoData==BlankSample)
	if (length(blanksamples)<1)
	{
		stop ("No samples matching a given names for \"blanks\" are found. Please check if function parameter `BlankSample` is set correctly.")
	} 
	Intensities.blank <- Intensities <- PeakTable[[1]][blanksamples,]
	Intensities.blank[is.na(Intensities.blank)] <- 0
	Intensities.blank <- apply(Intensities.blank,2,max)
	Intensities.rest <- apply(PeakTable[[1]][-c(blanksamples),],2,max,na.rm=TRUE)
	to.remove <- which(Intensities.blank > Intensities.rest)

	## Write no matching
	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("blank_removed",outputfile,sep=""), sets=c(1:length(RTs))[to.remove])

	# Not macthed from std file
	PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=c(1:length(RTs))[-c(to.remove)])
}
