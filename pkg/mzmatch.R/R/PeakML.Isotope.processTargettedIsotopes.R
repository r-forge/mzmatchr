PeakML.Isotope.processTargettedIsotopes <- function (molFormulaFile, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleType,
	sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, 
	fillGaps, massCorrection, useArea){

	readTargetsFromFile<- function(inputFile){
		# PRE: 
		#	inputFile: a tab separated csv file that conforms to RCreateXMLDB format e.i. "id", "name", "formula" as column headings
		# POST:
		#	Contents of the input file as a dataFrame that has masses added in column "mass"
	
		# Load the java project where the java class is located with dummy parameters

		molFrame <- read.csv(inputFile, sep="\t") # read the file as a data frame
		molMasses <- NULL
		for (imol in 1:length(molFrame$formula)){
			molMasses <- c (molMasses, PeakML.Methods.formula2mass(as.character(molFrame$formula)[imol]))
		}
		molFrame$mass <- molMasses
		molFrame
	}

	dir.create (outDirectory, showWarnings = FALSE)
	# To generate the tab delimited file
	csvFile <- paste(outDirectory, "/", outFileName, ".csv", sep ="")
	cat("\n", file=csvFile)
	# To generate the pdf plots
	pdfFile <- paste(outDirectory, "/", outFileName, ".pdf", sep="")
	pdf (file=pdfFile, page="a4", height=10, width=7)
	# Create the layout for the pdf
	layout(layoutMtx, heights=c(0.4, rep(1, nrow(layoutMtx)-1)),TRUE)
	# Reading the list of targets in the mol formula file
	molFrame <- readTargetsFromFile(molFormulaFile) # reading the molformula file
	# To save the abundance matrix if needed for later processing .
	molAbunList <- vector("list", nrow(molFrame))

	for (i in 1:nrow(molFrame)){
		metName <- as.character(molFrame$name[i])
		metFormula <- as.character(molFrame$formula[i])
		numCarbons <- PeakML.Methods.getCarbon(metFormula)
		metMass <- as.numeric(molFrame$mass[i])
#		massWindow <- PeakML.Methods.getPPMWindow(metMass, ppm)
		stdRT <- as.numeric(molFrame$rt[i]) * 60
		followCarbon <- as.numeric(molFrame$follow[i])+1

		cat(metName, ":\n")
		if (numCarbons==0){
			cat("\tThere are no carbons in ", metName, ", hence skipping. \n")
			next()
		}
		numCarbons <- numCarbons + 1 					# This is to account for the basal peaks as well.
#		metData <- list(metName, metFormula, numCarbons, metMass, ppm, massWindow, stdRT, stdRTWindow, sampleType)

		cat ("\tIdentifying isotopes: ")
		# get the UID of isotops
#		isotopeList <- PeakML.Isotope.getIsotopes (peakDataMtx, metData, sampleNames, mzXMLSrc, fillGaps, massCorrection)
		isotopeList <- PeakML.Isotope.getIsotopes (peakDataMtx, mzXMLSrc, sampleNames, numCarbons, metMass, ppm, massCorrection, stdRT, stdRTWindow, fillGaps)
		
		if (!is.null(unlist(isotopeList))){
			cat ("\n\tGenerating the plots. \n")
			isotopeChroms <- PeakML.Isotope.getChromData (isotopeList, chromDataList, phenoData, sampleGroups)
#			plotSamples(isotopeChroms, plotOrder, sampleGroups, metData, useArea, followCarbon)
			PeakML.Isotope.plotSamples(isotopeChroms, metName, metFormula, metMass, stdRT, sampleType, sampleGroups, plotOrder, useArea, followCarbon)
			
			abunMtxList <- PeakML.Isotope.getAbunMtxList(isotopeChroms, sampleGroups, useArea)
			molAbunList[[metName]] <- abunMtxList
			cat("Metabolite: ", toupper(metName), "\t Formula: ", metFormula, "\tMass: ", metMass, "\n", file=csvFile, append=TRUE)
			cat("-----------------------------------------------------------------------------\n", file=csvFile, append=TRUE)
			for (pkgrp in 1:length(abunMtxList)){
				cat(paste(metName,"_", metFormula, "_G", pkgrp, "\t", paste(sampleNames, collapse="\t") , "\n"), file=csvFile, append=TRUE)
				write.table(abunMtxList[[pkgrp]] , sep="\t", na= " ", file=csvFile, quote=FALSE, col.names=FALSE, append=TRUE) 
				cat("\n", file=csvFile, append=TRUE)
			}
			cat("\n")

		} else{
			cat("\tNo peaks found with mass:", metMass ," and its isotopes\n")
		}
	}
	save("molAbunList", file="abunList.Rdata")
	dev.off()
}
