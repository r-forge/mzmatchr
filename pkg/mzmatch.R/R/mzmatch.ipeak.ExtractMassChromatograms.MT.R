mzmatch.ipeak.ExtractMassChromatograms <- function(i=NULL, o=NULL, masstraces=NULL, label=NULL, threshold=NULL, ppm=NULL, h=NULL, v=NULL, nproc=3)
{
	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xms1425m -Xmx1425m -Xss160k -XX:+UseParallelGC -XX:ParallelGCThreads=10 -cp"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	mzmatch <- paste(java, " ", .find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
	if (!file.exists(o))
		dir.create (o)
	## setup the tool
	tool <- paste(mzmatch, "mzmatch.ipeak.ExtractMassChromatograms")
	if (!is.null(masstraces))
		tool <- paste(tool, "-masstraces", masstraces)
	if (!is.null(label))
		tool <- paste(tool, "-label", label)
	if (!is.null(threshold))
		tool <- paste(tool, "-threshold", threshold)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")
	while (length(i)!=0)
	
	{
		if (!is.null(i))
			tool <- paste(tool, "-i", i)
		if (!is.null(o))
			tool <- paste(tool, "-o", o)
		
	}	


	system(tool)
}
