PeakML.BlankFilter <- function (filename, ionisation="detect", Rawpath=NULL,outputfile,BlankSample=NULL)
{
	st <- system.time (PeakMLdata <- PeakML.Read (filename,ionisation,Rawpath))
	PeakTable <- PeakML.Methods.getCompleteTable (PeakMLdata)
	numberOfpeakSets <- ncol(PeakTable[[1]])

	## Check which samples are blanks
	blanksamples <- which (PeakMLdata$phenoData==BlankSample)
	if (length(blanksamples)<1)
	{
		stop ("No samples matching a given names for \"blanks\" are found. Please check if function parameter `BlankSample` is set correctly.")
	} 
	Intensities.blank <- PeakTable[[1]][blanksamples,]
	Intensities.blank[is.na(Intensities.blank)] <- 0
	Intensities.blank <- rbind(Intensities.blank,NULL)
	if (nrow(Intensities.blank)>1)
	{
		Intensities.blank <- apply(Intensities.blank,2,max)
	}
	Intensities.rest <- PeakTable[[1]][-c(blanksamples),]
	Intensities.rest[is.na(Intensities.rest)] <- 0
	Intensities.rest <- apply(Intensities.rest,2,max,na.rm=TRUE)
	to.remove <- which(Intensities.blank >= Intensities.rest)

	## Write no matching
	if (length(to.remove)>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=paste("blank_removed_",outputfile,sep=""), sets=c(1:numberOfpeakSets)[to.remove])
	}
	# Not macthed from std file
	if (length(to.remove)>0)
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=c(1:numberOfpeakSets)[-c(to.remove)])
	} else
	{
		PeakML.Methods.extractPeakGroups (PeakMLData=PeakMLdata, outputfile=outputfile, sets=c(1:numberOfpeakSets))
	}
}
