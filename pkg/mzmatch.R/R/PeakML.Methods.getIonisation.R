PeakML.Methods.getIonisation <- function(filename){
	
	# Binary file like cdf format, can't be parsed by XML.
	tryCatch(xmlfile <- xmlParse(filename), 
		error= function(e) {cat ("Measurement date and time information can be extracted only from mzML data files. \n")
		mzMLtree <- NULL}
	)

	if (!is.null(mzMLtree))
	{
		
		#x reffers there to the default namespace of the document				
		#ionisation <- unlist(getNodeSet(xmlRoot(mzMLtree), "//x:run/@startTimeStamp", "x"))
		getNodeSet(xmlRoot(mzMLtree), "//x:cvParam[@name='negative scan']", "x")

		startTimeStamp <- sub ("T"," ",startTimeStamp)
		rm (mzMLtree)
		gc ()
	} else
	{
		startTimeStamp <- NA
	}
	if (length(startTimeStamp)==0)
	{
		# mzXML is XML file (so will be parsed by xmlParse), but does not contain time stamp.
		cat ("Measurement date and time information can be extracted only from mzML data files. \n")
		startTimeStamp <- NA
	}
	startTimeStamp
}
