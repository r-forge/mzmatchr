### This needs to be fixed again for the data with MS2


"PeakML.GapFiller.fixed.UNSTABLE" <- function (filename, ionisation = "detect", Rawpath = NULL, outputfile, 
          ppm = 0, rtwin = 0, nSlaves = 1, fillAll = FALSE) 
{
	require (rJava)
	version.1 <- get("version.1", envir = .GlobalEnv)
    ## this variable is not used anywhere else - doesn't make sense!
    
    ### this is an internal function
    ### I will consider it when is used in the main function
    FillinPeaks <- function(peaknum) {
        whichpeakset <- numchromsexpected[fillinnums[peaknum], 1]
        subtable <- PeakMLdata$peakDataMtx[PeakMLdata$peakDataMtx[, 
                                                                  10] == whichpeakset, ]
        subtable <- rbind(subtable, NULL)
        rt_start <- min(subtable[, 5]) - rtwin
        rt_finis <- max(subtable[, 6]) + rtwin
        if (rt_finis > max(correctedRT)) {
            rt_finis <- max(correctedRT)
        }
        if (rt_start > max(correctedRT)) {
            rt_start <- max(correctedRT)
        }
        if (rt_finis < min(correctedRT)) {
            rt_finis <- min(correctedRT)
        }
        if (rt_start < min(correctedRT)) {
            rt_start <- min(correctedRT)
        }
        mz_start <- min(subtable[, 2])
        mz_finis <- max(subtable[, 3])
        mz_start <- mz_start - (mz_start * ppm/10^6)
        mz_finis <- mz_finis + (mz_finis * ppm/10^6)
        scan_start <- which(correctedRT >= rt_start)[1] - 1
        if (scan_start == 0) 
            scan_start = 1
        scan_finis <- which(correctedRT >= rt_finis)[1]
        C <- try(PeakML.Methods.getRawMat_2(allRawPeaks, scan_start, 
                                          scan_finis, mz_start, mz_finis, correctedRT, uncorrectedRT))
        if (class(C) == "try-error") {
            C <- c(1, 1, 1, 1, 1)
        }
        C <- rbind(C, NULL)
        if (nrow(C) <= 3 | length(unique(C[, 5])) <= 3) {
            scanids <- c(-1, -1, -1)
            retentiontimes <- c(-1, -1, -1)
            masses <- c(-1, -1, -1)
            intensities <- c(-1, -1, -1)
        }
        else {
            scanids <- C[, 3]
            retentiontimes <- C[, 2]
            masses <- C[, 4]
            intensities <- C[, 5]
        }
        OUT <- rbind(masses, intensities, retentiontimes, scanids - 
                         1)
    }
    
    st <- system.time(PeakMLdata <- PeakML.Read(filename, ionisation, Rawpath))
    # here we are just reading the PeakML input file while recording the time needed
    
    ionisation <- PeakMLdata$massCorrection[[2]]
    ### here I am reading the ionisation from the file... making the input parameter completely useless
    
    massCorrection <- PeakMLdata$massCorrection[[1]]
    ### getting the flag about the mass correction... is not used after this...non sense
    
    samplenames <- PeakMLdata$sampleNames
    ## getting the sample names
    
    rawdatafullpaths <- PeakMLdata$rawDataFullPaths
    #### getting the patha of the full data
    
    if (is.null(rawdatafullpaths)) {
        cat("Some of the raw data files are not accessible, we will not be able to fill in missing peaks. Please set \"Rawpath\" argument with location where files can be located\n")
        stop()
    }
    # this is self explanatory
    
    numchromsexpected <- unlist(lapply(1:max(PeakMLdata$peakDataMtx[,10]), 
                                       function(x) rep(x, length(samplenames))))
    # PeakMLdata$peakDataMtx[,10] contains the "GROUPIDs"....
    # we are basically saying that we are expecting to detect the same groupid
    # in each sample. Instead this is not the case (dim(PeakMLdata$peakDataMtx)[1] <= length(numchromsexpected))
    
    numchromsexpected <- cbind(numchromsexpected, NA, NA, NA)
    ### transforming numchromsexpected into a matrix where the first column is
    ### the old numchromsexpected and the next three are filled with NAs
    
    
    ######
    ### This for loops on the groupids
    inset <- c(1:length(samplenames)) ## this was inside the for, for no reason
    
    for (setnum in 1:max(PeakMLdata$peakDataMtx[, 10])) {
        rownums <- which(PeakMLdata$peakDataMtx[, 10] == setnum)
        ## identifing the rows in the peakDataMtx that share the same
        ## set num groupid 
        
        hit <- PeakMLdata$peakDataMtx[rownums, 9]
        ### in this vector we are extracting the measurement ids of the above
        ### mentioned rows
        
        oneEqualsSetnum <- which(numchromsexpected[, 1] == setnum)
        ### this vector is identifing which numchromsexpected has same current
        ### groupid (it will be a vector of length N samples)
        
        numchromsexpected[oneEqualsSetnum, 2] <- as.numeric(inset %in% hit)
        ### In the second column of the numchromsexpected where are putting
        ### 1 if there is a measure there, zero otherwise.
        
        missed <- which(inset %in% hit == FALSE)
        ### in this vector we are saving the positions of which of the nsamples
        ### expected measurements are not there
        
        
        if (length(missed) > 0) {
            detectedpeaks <- c(rep(1, length(hit)), rep(0, length(missed)))
            hit <- append(hit, missed)
            rownums <- append(rownums, rep(0, length(missed)))
        }else {
            detectedpeaks <- rep(1, length(hit))
        }
        ## if there are no missing peaks detectedpeaks is a vector of ones of the
        ## same length of hit.
        ## otherwise, detectedpeaks is c(rep(1, length(hit)), rep(0, length(missed)))
        ## hit <- append(hit, missed)
        ## rownums <- append(rownums, rep(0, length(missed)))
        
        
        ## this last 3 vectors are stored in the numchromsexpected matrix
        numchromsexpected[oneEqualsSetnum, 3] <- hit
        numchromsexpected[oneEqualsSetnum, 2] <- detectedpeaks
        numchromsexpected[oneEqualsSetnum, 4] <- rownums
    }
    #####
    ####  not sure yet why numchromsexpected is orgnized this way, but it looks
    ### like we are mapping which measuraments we have to retrive from the raw data
    
    
    colnames(numchromsexpected) <- NULL
    ### removing column names of this matrix...no reason for this
    
    chromslist <- vector("list", nrow(numchromsexpected))
    ### creating a list with as many elements as numchromsexpected (Nsamples*NgroupIDs)
    
    if (fillAll == TRUE) {
        detectedchromatograms <- numchromsexpected[, 2]
        numchromsexpected[, 2] <- 0
    }
    ### if the fillAll flag is set to true, we store the second column of
    ### numchromsexpected in detectedchromatograms befor filling it with zeros 
    
    notdetected <- which(numchromsexpected[, 2] == 0)
    ### finding which chroms are not detected, or better the ones that I have to retrive from the raw data
    
    whichfiles <- numchromsexpected[notdetected, 3]
    ### for each missing chroms I am storing in this vector the which samples are
    ### missing that group id...
    
    samplenums <- unique(whichfiles)
    #### how many samples have at least one thing missing?
    
    #if (length(samplenums != 0)) { ## I don't think this made sense at all.
    if (length(samplenums) != 0) {
        ##########################
        ######### In this for is going over the files one by one and working on filling it
        #########################
        for (filenum in 1:length(samplenums)) {
            samplefile <- samplenums[filenum]
            ### selecting the file
            cat("Working on file: ", rawdatafullpaths[samplefile], 
                "\n")
            
            rawfile <- openMSfile(rawdatafullpaths[samplefile], 
                                  verbose = FALSE)
            ## opening the raw data as a mzR object
            
            allRawPeaks <- peaks(rawfile)
            #accessing the MS raw data ....mzR function
            ### it is a list... each element is a scanid and contains a matrix
            ### with masses (column1) and intensities (column2)
            
            correctedRT <- as.numeric(PeakMLdata$correctedRTList[[samplefile]])
            ### extracting corrected RTs from peakml file
            uncorrectedRT <- header(rawfile)$retentionTime
            ### extracting RTs from raw file
            
            if (all(correctedRT == uncorrectedRT)) {
                rtCorrection <- FALSE
            }else {
                rtCorrection <- TRUE
            }
            ### checking if rt correction were applied or not
            
            fillinnums <- notdetected[whichfiles == samplefile]
            #### here I am checking which notdectected chroms are related to the
            #### current sample
            
            filledlist <- lapply(1:length(fillinnums), FillinPeaks)
            ### here I am applying to each of this not detected chroms the
            ### FillinPeaks function defined above (which is what is not working properly)
            
            
            for (i in 1:length(filledlist)) {
                chromslist[[fillinnums[i]]] <- filledlist[[i]]
            }
            
            mzR::close(rawfile)
            rm(filledlist)
        }
    ##########################
    ######### In this for is going over the files one by one and working on filling it
    #########################
    }
    project <- .jnew("peakml/util/rjava/Project", samplenames, 
                     rawdatafullpaths, as.character(PeakMLdata$phenoData))
    .jcall(project, returnSig = "V", method = "addHeaderAnnotation", 
           as.character("peakproc"), as.character("XCMS_Gapfilled"))
    for (measurementid in 1:length(samplenames)) {
        for (scannum in 1:length(PeakMLdata$correctedRTList[[measurementid]])) {
            .jcall(project, returnSig = "V", method = "addScanInfo", 
                   as.integer(measurementid - 1), as.numeric(PeakMLdata$correctedRTList[[measurementid]][scannum]), 
                   as.character(ionisation))
            .jcall(project, returnSig = "V", method = "addScanAnnotation", 
                   as.integer(measurementid - 1), as.integer(scannum - 
                                                                 1), as.character("RT_raw"), as.character(PeakMLdata$rawRTList[[measurementid]][scannum]))
        }
    }
    for (i in 1:length(chromslist)) {
        if (numchromsexpected[i, 2] == 0) {
            chrom <- chromslist[[i]]
        }
        else {
            ind <- numchromsexpected[i, 4]
            chrom <- PeakMLdata$chromDataList[[ind]]
        }
        .jcall(project, returnSig = "V", method = "addMassChromatogram", 
               as.integer(numchromsexpected[i, 3] - 1), as.integer(chrom[4, 
                                                                         ]), as.numeric(chrom[3, ]), as.numeric(chrom[1, 
                                                                                                                      ]), as.numeric(chrom[2, ]), as.character(ionisation))
    }
    setindexes <- vector("list", length(unique(numchromsexpected[, 
                                                                 1])))
    for (indexnumber in 1:length(setindexes)) {
        setindexes[[indexnumber]] <- which(numchromsexpected[, 
                                                             1] == indexnumber)
    }
    for (ind in 1:length(setindexes)) {
        .jcall(project, returnSig = "V", method = "addPeakSet", 
               as.integer(setindexes[[ind]] - 1))
    }
    if (!is.null(PeakMLdata$GroupAnnotations)) {
        PeakML.Methods.writeGroupAnnotations(project, PeakMLdata$GroupAnnotations)
    }
    .jcall(project, returnSig = "V", method = "write", outputfile)
}