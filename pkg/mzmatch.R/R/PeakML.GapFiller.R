PeakML.GapFiller <- function(filename,ionisation="detect",Rawpath=NULL,outputfile,ppm=0,rtwin=0)
{
	#####
	# Methods
	#####

	# .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(0))
	# .jcall(project, returnSig="[D", method="getMasses", as.integer(0))
	# .jcall(project, returnSig="[D", method="getIntensities", as.integer(0))
	# .jcall(project, returnSig="[D", method="getScanIDs", as.integer(0))
	# .jcall(project, returnSig="I", method="getMeasurementID", 0) - for masschromatogram, by giving index of peak.
	# .jcall(project,returnSig="[[D", method="getMassChromatograms")	
	# .jcall(project,returnSig="[S", method="getSetNames")
	# .jcall(project,returnSig="[S", method="getMeasurementNames")
	# .jcall(project,returnSig="[S", method="getFileNames")	
	# .jcall(project,returnSig="S", method="getAnnotation",as.integer(0),"string")
	# .jcall(project,returnSig="S", method="getGroupAnnotation",as.integer(0),"string")
	# 				string = relation.id or identification
	# .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes","MeasurementName")
	# .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes",as.integer(MeasurementID))
	
	
	#####
	# values for GetMassChromatograms
	#####
	
	#final int MINRT   =  0;		1
  	#final int MAXRT   =  1;		2
  	#final int AVGRT   =  2;		3
  	#final int MINSCAN  =  3;		4
  	#final int MAXSCAN  =  4;		5
  	#final int MINMZ   =  5;		6
  	#final int MAXMZ   =  6;		7
  	#final int AVGMZ   =  7;		8
  	#final int INTENSITY  =  8;		9
  	#final int SUMINTENSITY =  9; 	10
  	#final int MEASUREMENTID = 10; 	11
  	#final int SETID   = 11; 		12
  	#final int GROUPID  = 12; 		13


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
	
	## public int getNrScans(int measurementid)
	## public String getIonisation(int measurementid, int scanid)
	
	##.jcall (project,return="D",method="formulaToMass",as.character("[M1];[C4H4]n[+]"))

	#Correction factor by ionisation mode
	hydrogen <- 1.00782503214
	electron <- 0.00054857990924
	protonmass <- hydrogen-electron
		
	
	# create a new Project
	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
	st <- system.time(.jcall(project, returnSig="V", method="load",filename))
	cat ("Done in:",st[3],"s \n")
	
	#Detect ionisation mode, only first scan from first sample is used, no ionisation switching is supported yet
	
	
	if (ionisation=="detect")
	{
		ionisation<- .jcall (project,return="S",method="getIonisation",as.integer(0),as.integer(0))
	}
	
	if (ionisation!="negative" & ionisation!="positive")
	{
		protonCoef <- 0
	} 
	if (ionisation=="positive") 
	{
		protonCoef <- 1
	}
	if (ionisation=="negative") 
	{
		protonCoef <- -1
	} 
	
	## Read samplenames from peakml file	
	samplenames <- .jcall(project,returnSig="[S", method="getMeasurementNames")
	
	# Check if data directory for raw files ir readable, and if files are present there, 
	# if RaeDataPath is not set, check for raw data in default folder which is defined in PeakML file	
	if (!is.null(Rawpath))
	{
		directorycontent <- dir(Rawpath,recursive=TRUE)
		directorycontent <- sub(".mzXML","",directorycontent)
		matchedfiles <- which(directorycontent%in%samplenames)
	  	if(length(directorycontent[directorycontent%in%samplenames])==length(samplenames))
		{
			fileslocated="Y"	
			rawdatafullpaths <- dir(Rawpath,full.names=TRUE,recursive=TRUE)[matchedfiles]	
			cat(paste ("Raw data file located at: ",rawdatafullpaths,"\n",sep=""),"\n")
		} else
		{
			fileslocated="N"
			errormessage = paste("Raw data file: ",samplenames,"cant't be found in folder: ",Rawpath,"\n",sep="")
		}
	} else
	{
		rawdatafullpaths <- .jcall(project,returnSig="[S", method="getFileNames")
		## Check if files are readable at location defined in peakml file		
		locatedfiles <- which(file.exists(rawdatafullpaths)==TRUE)	
		if (length(locatedfiles)!=length(rawdatafullpaths))	
		{
		fileslocated="N"
		if (length(locatedfiles)!=0)
		{
			rawdatafullpaths <- rawdatafullpaths[-c(locatedfiles)]			
		}	
		errormessage = paste("Raw data file can't be read from location': ",rawdatafullpaths,"\n",sep="")
		} else
		{
			fileslocated="Y"
			cat(paste ("Raw data file located at: ",rawdatafullpaths,"\n",sep=""),"\n")
		}
	}

	if (fileslocated=="N")
	{
		cat (errormessage)
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	# Extract mass chromatograms
	cat ("Extracting mass chromatograms from peakml data file")	
	st <- system.time(masschromatograms <- .jcall(project,returnSig="[[D", method="getMassChromatograms"))
	cat (", done in:",st[3],"s \n")

	## Sort masschromatograms, check out which peak for which sample are missing?

	## This table reads generic masschrograms information from PeakML file, which should be present in all PeakML files
	## Output: matrix of 11 columns, column names : "AVGMZ","MINMZ","MAXMZ","AVGRT (this is recalculated at maximum intensity of peak, warning - this approach differs from that one which is used in XCMS)", "MINRT","MAXRT","SUMINTENSITY","MAXINTENSITY","measurement id","group id","peakcount"
	reconstructpeaktable <- function (i)
	{
		peakmldata <- .jevalArray(masschromatograms[[i]])
		peakdata <- matrix(ncol=11,nrow=1)		
		# AVGMZ		
		peakdata[1] <-  peakmldata[8]+(protonmass*protonCoef)
		# MINMZ		
		peakdata[2] <-  peakmldata[6]+(protonmass*protonCoef)
		# MAXMZ
		peakdata[3] <-  peakmldata[7]+(protonmass*protonCoef)
		# AVGRT				
		#peakdata[4] <-  peakmldata[3]
		## Calculate RT at the maximum intensity of the peak to avoid shifted RT's for peaks with long tails
		retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(i-1))
		intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(i-1))	
		peakdata[4] <- retentiontimes[which(intensities==max(intensities))[1]]
		# MINRT		
		peakdata[5] <-  peakmldata[1]
		# MAXRT		
		peakdata[6] <-  peakmldata[2]
		# SUMINTENSITY		
		peakdata[7] <-  peakmldata[10]
		# MAXINTENSITY
		peakdata[8] <-  peakmldata[9]
		# measurement id, +1 because java index starts at 0		
		peakdata[9] <- peakmldata[11]+1
		# group indexes, in next loop converted to list as required by xcms structure		
		peakdata[10] <- peakmldata[13]+1
		## Fort this stupid group table we need to count how many peaks from each sample class are grouped together.
		## So this vector will be used to genarate proper groups table structure
		peakdata[11] <- peakmldata[12]+1	
		peakdata	
	}	
		
	cat ("Extracting peak data from PeakMl file,")
	ST <- system.time(peakdata <- do.call (rbind,lapply(1:length(masschromatograms),reconstructpeaktable)))
	cat ("done in",ST[3],"s \n")
	
	## Generate a matrix of chromatogram numbers and peaksets

	## Vector of all possible peaks for every sample in peakset
	numchromsexpected <- unlist(lapply(1:max(peakdata[,10]),function (x) rep(x,length(samplenames))))

	## This table is used to generate information for chromatograms extraction/filling. 1st column - peakset number. 2nd-column presence/absence of the chromatogram in the peakml input file. 3rd column - data file number

	numchromsexpected <- cbind(numchromsexpected,NA,c(1:length(samplenames)))

	for (setnum in 1:max(peakdata[,10]))
	{
		inset <- c(1:length(samplenames))
		hit <- peakdata[peakdata[,10]==setnum,9]
		numchromsexpected[which(numchromsexpected[,1]==setnum),2] <- as.numeric(inset%in%hit)
	}

	## List object with mass chromatograms which are already in file
	chromslist <- vector("list",nrow(numchromsexpected))
	
	## Insert chromatograms from peakml file in the list
	detected <- which(numchromsexpected[,2]==1)
	
	system.time(
	for (chrnum in 1:nrow(peakdata))
	{
		retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(chrnum-1))
		intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(chrnum-1))
		## We have to correct masses for ionisation mode, that why there is a protonCoef defined
		masses <- .jcall(project, returnSig="[D", method="getMasses", as.integer(chrnum-1))+(protonmass*protonCoef)
		scanids <- .jcall(project, returnSig="[I", method="getScanIDs", as.integer(chrnum-1))
		chromslist[[detected[chrnum]]] <- rbind(masses,intensities,retentiontimes,scanids)
	})

	## Now fill in missing peaks
	notdetected <- which(numchromsexpected[,2]==0)
	whichfiles <- numchromsexpected[notdetected,3]
	samplenums <- unique(whichfiles)
	zerocount <- 0
	nonzerocount <- 0
	fillednums <- NULL

	system.time(
	for (filenum in 1:length(samplenums))
	{
		file <- samplenums[filenum]
		rawfile <- xcmsRaw(rawdatafullpaths[file])
		## detect which peaks to fill in for current data file
		fillinnums <- notdetected[whichfiles==file]
		
		FillinPeaks <- function(peaknum)
		{
			whichpeakset <- numchromsexpected[fillinnums[peaknum],1]
			## get a peaset RT and mass window from peakdata for all detected samples
			subtable <- peakdata[peakdata[,10]==whichpeakset,]
			subtable <- rbind(subtable,NULL)
			subtable <- apply(subtable,2,mean,na.rm=TRUE)

			## Eextract chromatogram
			# now we can locate the retention time as they will be in the raw data
			rt_start <- subtable[5]-rtwin
			rt_finis <- subtable[6]+rtwin
			mz_start <- subtable[2]
			mz_finis <- subtable[3]
			mz_start<- mz_start-(mz_start*ppm/10^6)
			mz_finis<- mz_finis+(mz_finis*ppm/10^6)
			######
			#  Extract data form raw data files
			######
				 
			C <- rawMat (rawfile,mzrange=cbind(mz_start, mz_finis),rtrange=c(rt_start,rt_finis))
			C <- rbind (C,NULL)

			## At some cases, two or more identical RT's are extracted, getting rid of them			
			repeatingRTS <- as.numeric(names(which(table(C[,1])>=2)))
			
			## Removing values with repeating RT's
			## Row with largest intensity are selected				
			if (length(repeatingRTS!=0))
			{			
				for (z in 1:length(repeatingRTS))
				{	
					Csub <- which(round(C[,1],5)==round(repeatingRTS[z],5))
					Csub <- Csub[-c(which(C[Csub,3]==max(C[Csub,3]))[1])]
					C <- C[-c(Csub),]
					C <- rbind(C,NULL)
				}
			}
			
			scanids <- which(rawfile@scantime%in%C[,1])
			## if RT correction was applied, scans should be extracted from raw RT's
			
			if (length(scanids)<3) 
			{
				scanids <- c(-1,-1,-1)
				retentiontimes <- c(-1,-1,-1)
				masses <- c(-1,-1,-1)
				intensities <- c(-1,-1,-1)
				zerocount <- zerocount+1
			} else {
				retentiontimes <- C[,1]
				#masses <- C[,2]+(protonmass*protonCoef)
				masses <- C[,2]
				intensities <- C[,3]
				nonzerocount <- nonzerocount+1
				fillednums <- append(fillednums,fillinnums[peaknum])
			}
			assign ("zerocount",zerocount,envir=.GlobalEnv)
			assign ("nonzerocount",nonzerocount,envir=.GlobalEnv)
			assign ("fillednums",fillednums,envir=.GlobalEnv)
			OUT <- rbind(masses,intensities,retentiontimes,scanids)
		}
		filledlist <- lapply(1:length(fillinnums),FillinPeaks)
		for (i in 1:length(filledlist))
		{
			chromslist[[fillinnums[i]]] <- filledlist[[i]]
		}
		rm (rawfile,filledlist)
	}
	)

	##SetNames, if peakml file has a several peaksets, they will be restored.
	samplegroups=peakdata[,11]
	sampleclasses <- .jcall(project,returnSig="[S", method="getSetNames")

	samplelookfunction <- function (x)	
	{
		phenoData <- sampleclasses[samplegroups[peakdata[,9]==x][1]]
		phenoData
	}
	phenoData <- unlist(lapply(1:length(samplenames),samplelookfunction)) 

	## Raw and corrected RT times, do we need these numbers?
	RAWrt <- vector ("list",length(samplenames))
	Corrt <- vector ("list",length(samplenames))

	for (measurementid in 1:length(samplenames))
	{
		Corrt[[measurementid]] <- .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes",as.integer(measurementid-1))
		raw_rt <- rep(NA,length(Corrt[[measurementid]]))
		for (scannum in 1:length(raw_rt))
		{
			raw_rt[scannum] <- as.numeric(.jcall(project,returnSig="S",method="getScanAnnotation",as.integer(measurementid-1), as.integer(scannum-1), as.character("RT_raw")))
		}
		RAWrt[[measurementid]] <- raw_rt
	}

	## Wipe out some unused large objects before writing
	rm (masschromatograms,project)

	## Write all of this out
	
	project <- .jnew("peakml/util/rjava/Project", samplenames, rawdatafullpaths, as.character(phenoData))
	# Insert peakpicker method name ir peakML file header.
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peakproc"),as.character("unknown"))
	
	# Inserting Scan numbers, RT corrected and RT raw for every sample
	for (measurementid in 1:length(samplenames))
	{
		for (scannum in 1:length(RAWrt[[measurementid]]))
		{
		.jcall(project, returnSig="V", method="addScanInfo", as.integer(measurementid-1),Corrt[[measurementid]][scannum],as.character(ionisation))
		}
	}
	
	## Appending RT raw values to PeakMl header
	for (measurementid in 1:length(samplenames))
	{
		for (scannum in 1:length(RAWrt[[measurementid]]))
		{
		#public void addScanAnnotation(int measurementid, int scanid, String label, String value)
		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(measurementid-1),as.integer(scannum-1),as.character("RT_raw"),as.character(RAWrt[[measurementid]][scannum]))
		}		
	}

	# finally we can store the mass chromatogram in our project
	# subtract 1 from the measurementid to get in line with java
	# IONISATION is set to neutral, as the data sets which are in peakml file are already recalculated.
	for (i in 1:length(chromslist))
	{
		chrom <- chromslist[[i]]
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(numchromsexpected[i,3]-1), as.integer(chrom[4,]), as.numeric(chrom[3,]),as.numeric(chrom[1,]), as.numeric(chrom[2,]), as.character(ionisation))
	}

	# now the mass chromatogram data has been collected the sets can be created - *sigh* memory consumption is such a bitch
	# -> this assumes that the sorting remains constant
	## Make a list with peak numbers which should be grouped together in sets
	setindexes <- vector("list",length(unique(numchromsexpected[,1])))
	
	for (indexnumber in 1:length(setindexes))
	{	
		setindexes[[indexnumber]] <- which(numchromsexpected[,1]==indexnumber)
	}

	for (ind in 1:length(setindexes))
	{
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(setindexes[[ind]]-1))
	}

	## Add Annotations to the peaks which are filled in, and which not.
	## Now we can add extra atributes to masschromatogram sets (groups in XCMS)
	## public void addGroupAnnotation(int groupid, String label, String value)
	fillindex <- rep (0,length(setindexes))
	## Peak sets which were appended
	fillindex[unique(numchromsexpected[fillednums,1])] <- 1

	for (groupnumber in 1:length(setindexes))
	{	
		.jcall(project, returnSig="V",method="addGroupAnnotation",as.integer(groupnumber-1),
		as.character("gap.filled"),as.character(fillindex[groupnumber]))
	}

	## Finally write file
	#cat(zerocount," of fillled in mass chromatograms are removed (length <3 scans)\n")	
	#cat (nonzerocount,"mass chromatogramms filled in.\n")
	.jcall(project, returnSig="V", method="write", outputfile)
}

