mzmatch.ipeak.util.Identify <- function(JHeapSize=1425,i=NULL, o=NULL, ppm=NULL, databases=NULL, minrt=NULL, maxrt=NULL, rtwindow=NULL, h=NULL, v=NULL)
{
	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xss160k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	## JHeapSize - define amount of RAM availiable for java VM
	JHeapSize <- paste(JHeapSize,"m",sep="")
	java <- paste(java," -Xms",JHeapSize," -Xmx",JHeapSize," -cp",sep="")
	mzmatch <- paste(java, " ", .find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
	
	## setup the tool
	tool <- paste(mzmatch, "mzmatch.ipeak.util.Identify")
		if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(databases))
		tool <- paste(tool, "-databases", databases)
	if (!is.null(minrt))
		tool <- paste(tool, "-minrt", minrt)
	if (!is.null(maxrt))
		tool <- paste(tool, "-maxrt", maxrt)
	if (!is.null(rtwindow))
		tool <- paste(tool, "-rtwindow", rtwindow)
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

	system(tool)
}
