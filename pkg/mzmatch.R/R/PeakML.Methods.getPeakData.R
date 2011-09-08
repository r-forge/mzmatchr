PeakML.Methods.getPeakData <- function(javaProject, massCorrection){
	# Returns peak data from mass chromatograms generated from PeakML files
	# PRE:  javaProject: The java project
	#	ionisation: "detect" will detect the ionisation mode from PeakML files
	# POST: peakDataMtx: A matrix containing matrix of 11 columns, column names : "AVGMZ","MINMZ","MAXMZ","AVGRT (this is recalculated at maximum 
	#	intensity of peak, warning -this approach differs from that one which is used in XCMS)", "MINRT",  "MAXRT", "SUMINTENSITY",
	#	"MAXINTENSITY", "measurement id", "group id", "peakcount" ,
	#	"chromData: A list of the form "[id] masses,intensities,retentiontimes,scanids"
	
	# Extract mass chromatograms
	cat ("Extracting mass chromatograms from peakML data file.\n")
	st <- system.time(masschromatograms <- .jcall(javaProject,returnSig="[[D", method="getMassChromatograms"))
	cat ("Mass chromatograms extracted in:",st[3],"s \n")

	# This table reads generic masschrograms information from PeakML file, which should be present in all PeakML files
	reconPeakTable <- function (chrom)
	{
		peakMLData <- .jevalArray(masschromatograms[[chrom]])
		peakDataMtx <- matrix(ncol=11,nrow=1)
		# AVGMZ
		peakDataMtx[1] <-  peakMLData[8]+(massCorrection)
		# MINMZ
		peakDataMtx[2] <-  peakMLData[6]+(massCorrection)
		# MAXMZ
		peakDataMtx[3] <-  peakMLData[7]+(massCorrection)
		## Calculate RT at the maximum intensity of the peak to avoid shifted RT's for peaks with long tails
		retentiontimes <- .jcall(javaProject, returnSig="[D", method="getRetentionTimes", as.integer(chrom-1))
		intensities <- .jcall(javaProject, returnSig="[D", method="getIntensities", as.integer(chrom-1))
		peakDataMtx[4] <- retentiontimes[which(intensities==max(intensities))[1]]
		# MINRT
		peakDataMtx[5] <-  peakMLData[1]
		# MAXRT
		peakDataMtx[6] <-  peakMLData[2]
		# SUMINTENSITY
		peakDataMtx[7] <-  peakMLData[10]
		# MAXINTENSITY
		peakDataMtx[8] <-  peakMLData[9]
		# measurement id, +1 because java index starts at 0
		peakDataMtx[9] <- peakMLData[11]+1
		# group indexes, in next loop converted to list as required by xcms structure
		peakDataMtx[10] <- peakMLData[13]+1
		## For this group table we need to count how many peaks from each sample class are grouped together.
		## So this vector will be used to genarate proper groups table structure
		peakDataMtx[11] <- peakMLData[12]+1
		peakDataMtx
	}

	cat ("Extracting peak data from PeakMl file,\n")
	st <- system.time(peakDataMtx <- do.call (rbind,lapply(1:length(masschromatograms),reconPeakTable)))
	cat ("Peak data created in ",st[3],"s \n")

	peakDataMtx
}

