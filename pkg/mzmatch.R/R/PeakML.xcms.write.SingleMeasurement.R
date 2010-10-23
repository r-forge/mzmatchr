PeakML.xcms.write.SingleMeasurement <- function(xset, filename,ionisation="detect",addscans=5,ppm=5,writeRejected=FALSE)
{
	#### Java calls	
	## public int getNrPeaksets() - get the number of PeakSets (groups) in memory 
	## public int getNrMassChromatograms() - get the number of masstraces in memory
	####

	## Global header annotations
	## public void addHeaderAnnotation(String label, String value)	
	## public String getHeaderAnnotation(String label)
	
	## Create scan info, for each ot the measurements
	## public void addScanInfo(int measurementid, double retentiontime, String polarity)
	
	## Add extra info to scans, for example RAW RTs
	## public void addScanAnnotation(int measurementid, int scanid, String label, String value)
	## public String getScanAnnotation(int measurementid, int scanid, String label)

	## Groups or sets		
	## public String getGroupAnnotation(int groupid, String name)
	## public void addGroupAnnotation(int groupid, String label, String value)

	## Masschromatogramms level
	## public String getAnnotation(int index, String name)
	## public void addAnnotation(int index, String label, String value)

	### ionisation: detect - default, from mzXML, negative, positive, unknown

	# create a new Project

	project <- .jnew("peakml/util/rjava/ProjectSingleMeasurement", as.character(rownames(xset@phenoData)), as.character(xset@filepaths))
	## Setting up counter for mass chromatograms with length 0	
	rejected <- NULL
	accepted <- NULL
	
	# load the raw data
	rawdata <- xcmsRaw(xset@filepaths[1])

	# correction for ionisation mode		
	if (ionisation=="detect") 
	{
		ionisation <- as.character(rawdata@polarity[1])
	}

	## Detect which peak pickging algorithm were used to generate XCMS object, for matchedFilter peaks table contain columns "intf", for centWave columns name ir "intb".
	intcolname <- colnames(xset@peaks)[8]
	if (intcolname=="intf")
	{
		peakpicker="matchedFilter"
	}

	if (intcolname=="intb")
	{
		peakpicker="centWave"
	}
	
	if (intcolname!="intb" & intcolname!="intf")
	{
		peakpicker="unknown"	
	}

	# Insert peakpicker method name ir peakML file header.
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peakproc"),as.character(peakpicker))

	# Insert column names of XCMS set peaktable in header
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peaktablecolnames"),as.character(paste(colnames(xset@peaks),collapse=";")))
	

	# Insert column names of XCMS set groups table in header
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("grouptablecolnames"),as.character(paste(colnames(xset@groups),collapse=";")))

	# Inserting Scan numbers, RT corrected and RT raw for every sample
	for (scannum in 1:length(xset@rt$corrected[[1]]))
	{
		.jcall(project, returnSig="V", method="addScanInfo", as.integer(0),xset@rt$corrected[[1]][scannum],as.character(ionisation))
	}
	
	## Appending RT raw values to PeakMl header
	for (scannum in 1:length(xset@rt$corrected[[1]]))
	{
		#public void addScanAnnotation(int measurementid, int scanid, String label, String value)
		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(0),as.integer(scannum-1),as.character("RT_raw"),as.character(xset@rt$raw[[1]][scannum]))
	}		
	
	# bollocks - now we need to load all the files again and retrieve the mass chromatogram data ... idiots; who thinks of these things
	cat("retrieving raw data\n")
	
	#print (measurementid)		
	cat("-", xset@filepaths[1], "\n")
	
	## Make a list of masschromatograms for every peak, so we can filter out peaks with crappy peakshape.
	

	chromatograms <- vector("list",nrow(xset@peaks))
	for (peakid in 1:nrow(xset@peaks))
	{
		#cat (peakid," ")
		##cat (nrpeaks[peakid]," ")		
		# retrieve the current peak and check whether it is part of the current measurement 
		# this assumes that all the peaks are sorted on xset@peaks["sample"]
		peak <- xset@peaks[peakid,]
		mz <- peak["mz"]
		mzdiff <- mz*ppm*10^-6

		# first locate the index of the rtmin/rtmax in the corrected retention time table
		rettime <- as.numeric(xset@rt$corrected[1][[1]])
		indx_start <- which(rettime == peak["rtmin"])[1]-addscans
		indx_finis <- which(rettime == peak["rtmax"])[1]+addscans
		
		if (indx_start < 1)
			indx_start <- 1
		if (indx_finis > length(rawdata@scantime))
			indx_finis <- length(rawdata@scantime)
		length <- indx_finis - indx_start

		# now we can locate the retention time as they will be in the raw data
		rt_start <- xset@rt$raw[1][[1]][indx_start]
		rt_finis <- xset@rt$raw[1][[1]][indx_finis]

		

		######
		#  Extract data form raw data files
		######
				 		
		## Extract extraplated peaks from raw data file
		C <- rawEIC (rawdata,mzrange=cbind(mz-mzdiff, mz+mzdiff),scanrange=cbind(indx_start,indx_finis))

		## To get the accurate mass measured at each scan, I use this extra call
		Ca <- rawMat (rawdata,mzrange=cbind(mz-mzdiff, mz+mzdiff),scanrange=cbind(indx_start,indx_finis))

		##plot (C[,3],type="l",axes=F)
		
		## At some cases, two or more identical RT's are extracted, getting rid of them			
		
		scans <- C$scan

		if (length(scans)<3 | all(C$intensity==0)) 
		{
			scans <- c(-1,-1,-1)
			retentiontimes <- c(-1,-1,-1)
			masses <- c(-1,-1,-1)
			intensities <- c(-1,-1,-1)
			rejected <- append(rejected,peakid)
		} else {
			retentiontimes <-  unlist(xset@rt$corrected[1])[scans]
			intensities <- C$intensity
			masses <- rep(peak["mz"],length(retentiontimes))
			# Fill in these correct masses, for which signal is detected, this is needed to plot complete peks
			detectedrts <- which(retentiontimes%in%Ca[,1])
			for (rtt in 1:length(detectedrts))
			{
				val <- retentiontimes[detectedrts[rtt]]
				val <- rbind(Ca[Ca[,1]==val,],NULL)
				if (nrow(val)>1)
				{
					val <- rbind(val[which(max(val[,3])==val[,3]),],NULL)
				}
				masses[detectedrts[rtt]] <- val[,2]
			}
			accepted <- append(accepted,peakid)
		}
		chromatograms[[peakid]] <- rbind(scans,retentiontimes,masses,intensities)
	}

	## Now we can filter on mass chromatogram as we want and how much we want. 
	for (chromnum in 1:length(accepted))
	{
		chrom <- chromatograms[[accepted[chromnum]]]
		if (chrom[4,1] >= max(chrom[4,]) | chrom[4,ncol(chrom)] >= max(chrom[4,]))
		{
			rejected <- append (rejected,accepted[[chromnum]])
		}
	}
	accepted <- accepted[-c(which(accepted%in%rejected))]
	
	# finally we can store the mass chromatogram in our project
		# subtract 1 from the measurementid to get in line with java
	for (chromnum in 1:length(accepted))
	{
		chrom <- chromatograms[[accepted[chromnum]]]
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(0), as.integer(chrom[1,]), chrom[2,],chrom[3,], chrom[4,], as.character(ionisation))
		# cat(peakid,length(scans),length(masses),"\n")
	}
	# now the mass chromatogram data has been collected the sets can be created - *sigh* memory consumption is such a bitch
	# -> this assumes that the sorting remains constant
	for (i in 1:length(accepted))
	{
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(i-1))
	}

	# and finally store the resulting data
	.jcall(project, returnSig="V", method="write", filename)
	
	## Write rejected chromatograms in saparate file
	if (writeRejected==TRUE)
	{
		project <- .jnew("peakml/util/rjava/ProjectSingleMeasurement", as.character(rownames(xset@phenoData)), as.character(xset@filepaths))
		for (scannum in 1:length(xset@rt$corrected[[1]]))
		{
			.jcall(project, returnSig="V", method="addScanInfo", as.integer(0),xset@rt$corrected[[1]][scannum],as.character(ionisation))
		}
		for (chromnum in 1:length(rejected))
		{
			chrom <- chromatograms[[rejected[chromnum]]]
			.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(0), as.integer(chrom[1,]), chrom[2,],chrom[3,], chrom[4,], as.character(ionisation))
		}
		for (i in 1:length(rejected))
		{
			.jcall(project, returnSig="V", method="addPeakSet", as.integer(i-1))
		}
		filename2 <- sub (".peakml","",filename)
		filename2 <- paste(filename2,"_rejected.peakml",sep="")
		.jcall(project, returnSig="V", method="write", filename2)
	}

}
