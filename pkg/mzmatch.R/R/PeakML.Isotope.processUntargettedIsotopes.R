#PeakML.Isotope.UntargettedIsotopes(baseDir="/home/cbaunni/Project/data_analysis/hct1162/trunk/NEG/", outFileName="utgt_test", mzXMLSrc="/home/cbaunni/Project/data_analysis/hct1162/mzxml/", fillGaps="NONE", useArea=TRUE)

PeakML.Isotope.processUntargettedIsotopes <- function (peakMLFile, databases, outDirectory, outFileName, layoutMtx, ppm, stdRTWindow, sampleType,
	sampleNames, peakDataMtx, chromDataList, phenoData, sampleGroups, plotOrder, mzXMLSrc, 
	fillGaps, massCorrection, useArea){

	dir.create (outDirectory, showWarnings = FALSE)
	# To generate the tab delimited file
	csvFile <- paste(outDirectory, "/", outFileName, ".csv", sep ="")
	cat("\n", file=csvFile)
	# To generate the pdf plots
	pdfFile <- paste(outDirectory, "/", outFileName, ".pdf", sep="")
	pdf (file=pdfFile, paper="a4", height=10, width=7)
	# Create the layout for the pdf
	layout(layoutMtx, heights=c(0.4, rep(1, nrow(layoutMtx)-1)),TRUE)
	# Reading the list of targets in the mol formula file
	molFrame <- PeakML.Isotope.DB2Text(peakMLFile, databases)
	# To save the abundance matrix if needed for later processing .
	molAbunList <- vector("list", nrow(molFrame))

	for (i in 1:nrow(molFrame)){
		metName <- as.character(molFrame$name[i])
		metFormula <- strsplit(as.character(molFrame$formula[i]), ",")[[1]][1]#as.character(molFrame$formula[i]) **********GLITCH**********
		numCarbons <- PeakML.Methods.getCarbon(metFormula)
		metMass <- as.numeric(molFrame$mass[i])
#		massWindow <- PeakML.Methods.getPPMWindow(metMass, ppm)
		stdRT <- as.numeric(molFrame$rt[i])# * 60
		
		if (is.null(molFrame$follow[i])){
			followCarbon <-  numCarbons + 1
		}else{
			followCarbon <- as.numeric(molFrame$follow[i])+1
		}
		
		
		cat(metName, ":\n")
		cat(metFormula, numCarbons, "\n")

		if (numCarbons==0){
			cat("\tThere are no carbons in ", metName, ", hence skipping. \n")
			next()
		}
		cat("\n")
		
		numCarbons <- numCarbons + 1 					# This is to account for the basal peaks as well.

		cat ("\tIdentifying isotopes: ")
		# get the UID of isotops
		isotopeList <- PeakML.Isotope.getIsotopes (peakDataMtx, mzXMLSrc, sampleNames, numCarbons, metMass, ppm, massCorrection, stdRT, stdRTWindow, fillGaps)
		
		
		# This is where to include the isotop level filter.
		isotopeList <- PeakML.Isotope.filteroutUnlabelled(isotopeList, numCarbons, fillGaps, sampleNames, stringency=30)
		
		
		if (!is.null(unlist(isotopeList))){
			cat ("\n\tGenerating the plots. \n")
			isotopeChroms <- PeakML.Isotope.getChromData (isotopeList, chromDataList, phenoData, sampleGroups)
			PeakML.Isotope.plotSamples(isotopeChroms, metName, metFormula, metMass, stdRT, sampleType, sampleGroups, plotOrder, useArea, followCarbon)
			
			ratioMtxList <- PeakML.Isotope.getRatioMtxList(isotopeChroms[[2]], sampleGroups, useArea, metName)
			
			molAbunList[[metName]] <- ratioMtxList
			
			cat("Metabolite: ", toupper(metName), "\t Formula: ", metFormula, "\tMass: ", metMass, "\n", file=csvFile, append=TRUE)
			cat("-----------------------------------------------------------------------------\n", file=csvFile, append=TRUE)
			for (pkgrp in 1:length(ratioMtxList)){ 
				cat("Group: ", pkgrp, "\n", file=csvFile, append = TRUE)
				sNames <- paste(sampleNames, collapse="\t")
				cat(paste(metName, metFormula, "\t", sNames , "\n"), file=csvFile, append=TRUE)
				write.table(ratioMtxList[[pkgrp]] , sep="\t", na= " ", file=csvFile, quote=FALSE, col.names=FALSE, append=TRUE) 
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
