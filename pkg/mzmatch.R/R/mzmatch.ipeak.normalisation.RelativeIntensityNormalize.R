mzmatch.ipeak.normalisation.RelativeIntensityNormalize <- function(JHeapSize=1425,i=NULL, o=NULL, sets=NULL, v=NULL, h=NULL)
{
	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xss128k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	## JHeapSize - define amount of RAM availiable for java VM
	JHeapSize <- paste(JHeapSize,"m",sep="")
	java <- paste(java," -Xms",JHeapSize," -Xmx",JHeapSize," -cp",sep="")
	mzmatch <- paste(java, " ", .find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
	
	## setup the tool
	tool <- paste(mzmatch, "mzmatch.ipeak.normalisation.RelativeIntensityNormalize")
	if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(sets))
		tool <- paste(tool, "-sets", sets)
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")

	system(tool)
}
