PeakML.Normalization <- function(filename,ionisation="detect",Rawpath=NULL,outputfile,values=NULL)
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
	
	if (length(samplenames)!=length(values) & !is.null(values))
	{
		cat("Number of the samples in the peakml file does not match a number of normalization values given in function argument \"values\". Length of the normalization values must be the same as the sample number in peakml file.\n")
		stop ()
	}
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
		cat ("Some of the raw data files are not accessible. Please set \"Rawpath\" argument with location where files are located.\n")
		stop ()
	}

	# Extract mass chromatograms
	cat ("Extracting mass chromatograms from peakml data file")	
	st <- system.time(masschromatograms <- .jcall(project,returnSig="[[D", method="getMassChromatograms"))
	cat (", done in:",st[3],"s \n")

	## Sort masschromatograms, check out which peak for which sample are missing?

	## This table reads generic masschrograms information from PeakML file, which should be present in all PeakML files
	## Output: matrix of 2 columns, column names : "measurement id","group id".
	reconstructpeaktable <- function (i)
	{
		peakmldata <- .jevalArray(masschromatograms[[i]])
		peakdata <- matrix(ncol=3,nrow=1)		
		
		# measurement id, +1 because java index starts at 0		
		peakdata[1] <- peakmldata[11]+1
		# group indexes, in next loop converted to list as required by xcms structure		
		peakdata[2] <- peakmldata[13]+1
		## So this vector will be used to genarate proper groups table structure
		peakdata[3] <- peakmldata[12]+1
		peakdata
	}
		
	cat ("Extracting peak data from PeakMl file,")
	ST <- system.time(peakdata <- do.call (rbind,lapply(1:length(masschromatograms),reconstructpeaktable)))
	cat ("done in",ST[3],"s \n")

	## Append peakdata with a vector of normalisation values which should be used for every vector of intensities

	## if values=NULL use TIC's of the raw data files for normalisation
	if (is.null(values))
	{
		values <- rep(NA,length(samplenames))
		for (file in 1:length(samplenames))
		{
			rawfile <- xcmsRaw(rawdatafullpaths[file])
			values[file] <- sum(rawfile@tic)
		}
	}
	
	normvals <- rep(NA,nrow(peakdata))
	for (nval in 1:length(samplenames))
	{
		normvals[peakdata[,1]==nval] <- values[nval]
	}

	## Object: List, with mass chromatograms 
	chromslist <- vector("list",nrow(peakdata))
	
	system.time(
	for (chrnum in 1:nrow(peakdata))
	{
		retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(chrnum-1))
		intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(chrnum-1))/normvals[chrnum]
		## We have to correct masses for ionisation mode, that why there is a protonCoef defined
		masses <- .jcall(project, returnSig="[D", method="getMasses", as.integer(chrnum-1))+(protonmass*protonCoef)
		scanids <- .jcall(project, returnSig="[I", method="getScanIDs", as.integer(chrnum-1))
		chromslist[[chrnum]] <- rbind(masses,intensities,retentiontimes,scanids)
	})

	##SetNames, if peakml file has a several peaksets, they will be restored.
	samplegroups=peakdata[,3]
	sampleclasses <- .jcall(project,returnSig="[S", method="getSetNames")

	samplelookfunction <- function (x)	
	{
		phenoData <- sampleclasses[samplegroups[peakdata[,1]==x][1]]
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
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(peakdata[i,1]-1), as.integer(chrom[4,]), as.numeric(chrom[3,]),as.numeric(chrom[1,]), as.numeric(chrom[2,]), as.character(ionisation))
	}

	# now the mass chromatogram data has been collected the sets can be created - *sigh* memory consumption is such a bitch
	# -> this assumes that the sorting remains constant
	## Make a list with peak numbers which should be grouped together in sets
	setindexes <- vector("list",length(unique(peakdata[,2])))
	
	for (indexnumber in 1:length(setindexes))
	{	
		setindexes[[indexnumber]] <- which(peakdata[,2]==indexnumber)
	}

	for (ind in 1:length(setindexes))
	{
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(setindexes[[ind]]-1))
	}

	## Finally write file
	.jcall(project, returnSig="V", method="write", outputfile)
}
