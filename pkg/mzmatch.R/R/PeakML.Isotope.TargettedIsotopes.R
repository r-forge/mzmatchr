PeakML.Isotope.TargettedIsotopes <- function(baseDir, molFormulaFile, outFileName, mzXMLSrc=NULL, outDirectory = "targettedIsotops",  peakMLFile="final_combined_related_identified.peakml", sampleGroups = NULL, layoutMtx = NULL, ppm = 3, trendPlots = NULL, fillGaps = "ALLPEAKS", useArea = FALSE, stdRTWindow = NULL, baseCorrection=TRUE){
	# PRE:
	#	peakMLFiles: the complete peakml dataset
	#	molFormulaFile: file containing the list of molecules whoes isotops has to be found out
	#	outDirectory: is the outDirectory where the output has to be saved
	#	mzXMLSrc: is the source of the mzXML file is not in the current working outDirectory
	#	include_ionisation: not necessary in this cases
	#	ionisation: set this if include_ionisation=TRUE
	#	loadSavedData: load from the saved peakml file
	# 	sampleType: the sample type eg. NEG, POS etc
	# POST:
	#	vector containing the list of isotops
	
	## Reads the peakml file & prepare the parameters to scan for isotops
	## --------------------------------------------------------------------
	cat("Indentifying isotopes in sample\n")
	#setwd (paste(baseDir, sampleType, sep="/"))
	setwd (baseDir)
	
	if (is.null(mzXMLSrc)){
	#if (!is.null(mzXMLSrc)){
	#	mzXMLSrc <- paste(mzXMLSrc, sampleType, sep="/")
	#} else {
		stop ("Please provide the location of the raw data (mzXML) files ")
	}
	
	if (file.exists("cpData.Rdata") == TRUE){
		load("cpData.Rdata")
	} else{
		chromPeakData <- PeakML.Read(peakMLFile, ionisation = "neutral", mzXMLSrc)
		save("chromPeakData", file="cpData.Rdata")
	}
	
	peakDataMtx <- chromPeakData$peakDataMtx
	chromDataList <- chromPeakData$chromDataList
	sampleClasses <- chromPeakData$sampleClasses
	sampleNames <- chromPeakData$sampleNames
	massCorrection <- PeakML.Methods.getMassCorrection(filename=peakMLFile)
	phenoData <- PeakML.Methods.getPhenoData(sampleClasses, sampleNames, peakDataMtx)


	if (is.null(sampleGroups)) sampleGroups <- unique(phenoData)		# To enable the user to change the order of the samples

	if (is.null(trendPlots)) trendPlots <- c("RATIO","TREND", "LABELLED", "TOTRATIO")
	
	if (length(sampleGroups)>22) {
		if (is.null(layoutMtx)) stop("You have more than 22 samples to plot. Please specify an appropriate layout matrix.\n") else plotOrder <- c(sampleGroups, trendPlots)
	} else {
		numSG <- length(sampleGroups)
		if (numSG <= 7){
			if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,2, 3,4,5,6,7,8,9, 10,11,11,12,12,13,13),3,7, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 7-numSG), trendPlots)
		} else if (numSG > 7 & numSG <=14){
			if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,2, 3,4,5,6,7,8,9, 10,11,12,13,14,15,16, 17,18,18,19,19,20,20),4,7, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 14-numSG), trendPlots)
		} else if (numSG > 14 & numSG <=18){
			if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,1,2,2, 3,4,5,6,7,8,9,10,11, 12,13,14,15,16,17,18,19,20, 21,21, 22,22,22, 23,23, 24,24),4,9, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 18-numSG), trendPlots)
		} else if (numSG > 18 & numSG <=22){
			 if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,1,1,1,1,1,2,2, 3,4,5,6,7,8,9,10,11,12,13, 14,15,16,17,18,19,20,21,22,23,24, 25,25, 26,26,26, 27,27,27, 28,28,28),4,11, byrow=TRUE)
			plotOrder <- c(sampleGroups, rep("EMPTY", 22-numSG), trendPlots)
		}

		#if (is.null(layoutMtx)) layoutMtx <- matrix(c(1,1,1,1,2, 3,4,5,6,7, 8,9,10,11,12, 13 ,14,14, 15,15),4,5, byrow=TRUE)
		#if (length(sampleGroups)<10) plotOrder <- c(sampleGroups, rep("EMPTY", 10-length(sampleGroups)), trendPlots)
		#if (length(sampleGroups)==10) plotOrder <- c(sampleGroups, trendPlots)
	}
	
	PeakML.Isotope.processTargettedIsotopes(molFormulaFile, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, 
	fillGaps, massCorrection, useArea, baseCorrection)

}

