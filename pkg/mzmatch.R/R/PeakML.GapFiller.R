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

	# create a new Project
	project <- .jnew("peakml/util/rjava/Project", rep("A",3), rep("A",3), rep("A",3))
	cat ("Loading peakML file in memory (it can take some time, sorry) \n")	
	st <- system.time(.jcall(project, returnSig="V", method="load",filename))
	cat ("Done in:",st[3],"s \n")
	
	protonMass <- PeakML.Methods.getProtonMass()
	protonCoef <- PeakML.Methods.getProtonCoef(project, ionisation)
	massCorrection <- PeakML.Methods.getMassCorrection(project, ionisation)
	samPath <- PeakML.Methods.getRawDataPaths(project, Rawpath)
	samplenames <- samPath[[1]]
	rawdatafullpaths <- samPath[[2]]
	
	
	if (is.null(rawdatafullpaths)){
		cat ("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
		stop ()
	}

	## Generate a matrix of chromatogram numbers and peaksets
	peakdata <- PeakML.Methods.getPeakData(project, massCorrection)

	## Vector of all possible peaks for every sample in peakset
	numchromsexpected <- unlist(lapply(1:max(peakdata[,10]),function (x) rep(x,length(samplenames))))

	## This table is used to generate information for chromatograms extraction/filling. 1st column - peakset number. 2nd-column presence/absence of the chromatogram in the peakml input file. 3rd column - data file number
	numchromsexpected <- cbind(numchromsexpected,NA,c(1:length(samplenames)))

	for (setnum in 1:max(peakdata[,10])){
		inset <- c(1:length(samplenames))
		hit <- peakdata[peakdata[,10]==setnum,9]
		numchromsexpected[which(numchromsexpected[,1]==setnum),2] <- as.numeric(inset%in%hit)
	}

	## List object with mass chromatograms which are already in file
	chromslist <- vector("list",nrow(numchromsexpected))
	
	## Insert chromatograms from peakml file in the list
	detected <- which(numchromsexpected[,2]==1)
	
	system.time(
	for (chrnum in 1:nrow(peakdata)){
		retentiontimes <- .jcall(project, returnSig="[D", method="getRetentionTimes", as.integer(chrnum-1))
		intensities <- .jcall(project, returnSig="[D", method="getIntensities", as.integer(chrnum-1))
		## We have to correct masses for ionisation mode, that why there is a protonCoef defined
		masses <- .jcall(project, returnSig="[D", method="getMasses", as.integer(chrnum-1))+(protonMass*protonCoef)
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

	if (length(samplenums!=0)){
		system.time(
		for (filenum in 1:length(samplenums)){
			file <- samplenums[filenum]
			rawfile <- xcmsRaw(rawdatafullpaths[file])
			## detect which peaks to fill in for current data file
			fillinnums <- notdetected[whichfiles==file]
			
			FillinPeaks <- function(peaknum){
				whichpeakset <- numchromsexpected[fillinnums[peaknum],1]
				## get a peakset RT and mass window from peakdata for all detected samples
				subtable <- peakdata[peakdata[,10]==whichpeakset,]
				subtable <- rbind(subtable,NULL)
				#subtable <- apply(subtable,2,mean,na.rm=TRUE)

				## Extract chromatogram
				# now we can locate the retention time as they will be in the raw data
				rt_start <- min(subtable[,5])-rtwin
				rt_finis <- max(subtable[,6])+rtwin
				if (rt_finis > max(rawfile@scantime)){
					rt_finis <- max(rawfile@scantime)
				}
				if (rt_start > max(rawfile@scantime)){
					rt_start <- max(rawfile@scantime)
				}
				if (rt_finis < min(rawfile@scantime)){
					rt_finis <- min(rawfile@scantime)
				}
				if (rt_start < min(rawfile@scantime)){
					rt_start <- min(rawfile@scantime)
				}
				
				mz_start <- min(subtable[,2])
				mz_finis <- max(subtable[,3])
				mz_start<- mz_start-(mz_start*ppm/10^6)
				mz_finis<- mz_finis+(mz_finis*ppm/10^6)
				######
				#  Extract data form raw data files
				######
				if ((rt_finis-rt_start)<5){
					C <- c(1,1,1)
				} else{
					C <- rawMat (rawfile,mzrange=cbind(mz_start, mz_finis),rtrange=c(rt_start,rt_finis))
				}
				C <- rbind (C,NULL)
	
				## At some cases, two or more identical RT's are extracted, getting rid of them
				repeatingRTS <- as.numeric(names(which(table(C[,1])>=2)))
				
				## Removing values with repeating RT's
				## Row with largest intensity are selected
				if (length(repeatingRTS!=0)){
					for (z in 1:length(repeatingRTS)){
						Csub <- which(round(C[,1],5)==round(repeatingRTS[z],5))
						Csub <- Csub[-c(which(C[Csub,3]==max(C[Csub,3]))[1])]
						C <- C[-c(Csub),]
						C <- rbind(C,NULL)
					}
				}
			
				scanids <- which(rawfile@scantime%in%C[,1])
				## if RT correction was applied, scans should be extracted from raw RT's
			
				## Remove artefacts
				if (length(scanids)<=3 | length(unique(C[,3]))<=3){
					scanids <- c(-1,-1,-1)
					retentiontimes <- c(-1,-1,-1)
					masses <- c(-1,-1,-1)
					intensities <- c(-1,-1,-1)
					zerocount <- zerocount+1
				} else {
					retentiontimes <- C[,1]
					#masses <- C[,2]+(protonMass*protonCoef)
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
			
			## For debug purpose only
			#for (bd in 1: length(fillinnums))
			#{
			#	cat (bd,"\n")
			#	FillinPeaks(bd)
			#}

			assign ("zerocount",zerocount,envir=.GlobalEnv)
			assign ("nonzerocount",nonzerocount,envir=.GlobalEnv)
			assign ("fillednums",fillednums,envir=.GlobalEnv)
			
			for (i in 1:length(filledlist)){
				chromslist[[fillinnums[i]]] <- filledlist[[i]]
			}
			rm (rawfile,filledlist)
		}
		)
	}
	
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

	for (measurementid in 1:length(samplenames)){
		Corrt[[measurementid]] <- .jcall(project,returnSig="[D", method="getMeasurementRetentionTimes",as.integer(measurementid-1))
		raw_rt <- rep(NA,length(Corrt[[measurementid]]))
		for (scannum in 1:length(raw_rt)){
			raw_rt[scannum] <- as.numeric(.jcall(project,returnSig="S",method="getScanAnnotation",as.integer(measurementid-1), as.integer(scannum-1), as.character("RT_raw")))
		}
		RAWrt[[measurementid]] <- raw_rt
	}

	## Wipe out some unused large objects before writing
#	rm (masschromatograms,project)
	rm (project) # as masschromatograms is inside a function now
	## Write all of this out
	
	project <- .jnew("peakml/util/rjava/Project", samplenames, rawdatafullpaths, as.character(phenoData))
	# Insert peakpicker method name ir peakML file header.
	## public void addHeaderAnnotation(String label, String value)	
	.jcall(project, returnSig="V", method="addHeaderAnnotation",as.character("peakproc"),as.character("unknown"))
	
	# Inserting Scan numbers, RT corrected and RT raw for every sample
	for (measurementid in 1:length(samplenames)){
		for (scannum in 1:length(RAWrt[[measurementid]])){
		.jcall(project, returnSig="V", method="addScanInfo", as.integer(measurementid-1),Corrt[[measurementid]][scannum],as.character(ionisation))
		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(measurementid-1),as.integer(scannum-1),as.character("RT_raw"),as.character(RAWrt[[measurementid]][scannum]))
		}
	}
	
	## Appending RT raw values to PeakMl header
#	for (measurementid in 1:length(samplenames)){
#		for (scannum in 1:length(RAWrt[[measurementid]])){
#		#public void addScanAnnotation(int measurementid, int scanid, String label, String value)
#		.jcall(project, returnSig="V", method="addScanAnnotation", as.integer(measurementid-1),as.integer(scannum-1),as.character("RT_raw"),as.character(RAWrt[[measurementid]][scannum]))
#		}
#	}

	# finally we can store the mass chromatogram in our project
	# subtract 1 from the measurementid to get in line with java
	# IONISATION is set to neutral, as the data sets which are in peakml file are already recalculated.
	for (i in 1:length(chromslist)){
		chrom <- chromslist[[i]]
		.jcall(project, returnSig="V", method="addMassChromatogram", as.integer(numchromsexpected[i,3]-1), as.integer(chrom[4,]), as.numeric(chrom[3,]),as.numeric(chrom[1,]), as.numeric(chrom[2,]), as.character(ionisation))
	}

	# now the mass chromatogram data has been collected the sets can be created - *sigh* memory consumption is such a bitch
	# -> this assumes that the sorting remains constant
	## Make a list with peak numbers which should be grouped together in sets
	setindexes <- vector("list",length(unique(numchromsexpected[,1])))
	
	for (indexnumber in 1:length(setindexes)){
		setindexes[[indexnumber]] <- which(numchromsexpected[,1]==indexnumber)
	}
	for (ind in 1:length(setindexes)){
		.jcall(project, returnSig="V", method="addPeakSet", as.integer(setindexes[[ind]]-1))
	}

	## Add Annotations to the peaks which are filled in, and which not.
	## Now we can add extra atributes to masschromatogram sets (groups in XCMS)
	## public void addGroupAnnotation(int groupid, String label, String value)
	fillindex <- rep (0,length(setindexes))
#	## Peak sets which were appended
	fillindex[unique(numchromsexpected[fillednums,1])] <- 1
	for (groupnumber in 1:length(setindexes)){
		.jcall(project, returnSig="V",method="addGroupAnnotation",as.integer(groupnumber-1), as.character("gap.filled"),as.character(fillindex[groupnumber]))
	}

	## Finally write file
	#cat(zerocount," of fillled in mass chromatograms are removed (length <3 scans)\n")
	#cat (nonzerocount,"mass chromatogramms filled in.\n")
	.jcall(project, returnSig="V", method="write", outputfile)
}
