mzmatch.R.Setup <- function (projectFolder=NULL, samplelist=NULL, outputfolder="peakml")
{
	if (is.null(projectFolder))
	{
		projectFolder <- tclvalue (tkchooseDirectory())
	}
	setwd (projectFolder)

	if (is.null(samplelist))
	{
		mzXMLfiles.fullnames <- dir(full.names=TRUE,pattern="\\.mzXML$",recursive=TRUE)
		mzXMLfiles.fullnames <- append(mzXMLfiles.fullnames,dir(full.names=TRUE,pattern="\\.mzData$",recursive=TRUE))
		mzXMLfiles.fullnames <- append(mzXMLfiles.fullnames,dir(full.names=TRUE,pattern="\\.mzML$",recursive=TRUE))

		mzXMLfiles.shortnames <- dir(full.names=FALSE,pattern="\\.mzXML$",recursive=TRUE)
		mzXMLfiles.shortnames <- append(mzXMLfiles.shortnames,dir(full.names=FALSE, pattern="\\.mzData$",recursive=TRUE))
		mzXMLfiles.shortnames <- append(mzXMLfiles.shortnames,dir(full.names=FALSE, pattern="\\.mzML$",recursive=TRUE))
		
		outputfilenames <- paste("peakml/",sub(".mzXML", "", mzXMLfiles.fullnames), ".peakml", sep="")
		outputfilenames <- sub(".mzData", ".peakml", mzXMLfiles.fullnames)
		outputfilenames <- sub(".mzML", ".peakml", mzXMLfiles.fullnames)

		sampleList <- data.frame (filenames=mzXMLfiles.fullnames, sampleClass=rep("",length(mzXMLfiles.fullnames)), globalClass=rep("",length(mzXMLfiles.fullnames)))
		fix (sampleList)
		## fix writed edited object to the .GlobalEnve, we have to read it back into fucntion environmnet.
		sampleList <- get ("sampleList", envir=.GlobalEnv)
		sampleList[,1] <- as.character(sampleList[,1])
		write.table (sampleList,file="sample_setup.tsv",sep="\t",row.names=FALSE)
	}
	else
	{
		sampleList <- read.table (samplelist, sep="\t", header=TRUE)
		mzXMLfiles.shortnames <- rep(NA,nrow(sampleList))
		sampleList[,1] <- as.character(sampleList[,1])
		for (i in 1:nrow(sampleList))
		{
			mzXMLfiles.shortnames[i] <- PeakML.Methods.extractFileName(sampleList$filenames[i])
		}
		outputfilenames <- paste(outputfolder,"/",sub(".mzXML", "", mzXMLfiles.shortnames), ".peakml", sep="")
	}
	if (!file.exists(outputfolder))
	{
		dir.create (outputfolder)
	}
	sampleList$outputfilenames <- outputfilenames
	assign("sampleList", sampleList, envir=.GlobalEnv)
}
