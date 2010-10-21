mzmatch.ipeak.Combine <- function(JHeapSize=1425,i=NULL, o=NULL, label=NULL, labels=NULL, ppm=NULL, rtwindow=NULL, combination=NULL, h=NULL, v=NULL)
{
	## define the java runtime parameters
	java <- "java -da -dsa -Xmn1g -Xss128k -XX:+UseParallelGC -XX:ParallelGCThreads=10"
	## locate the mzmatch.jar file (it's included in the peakmlR package)
	## JHeapSize - define amount of RAM availiable for java VM
	JHeapSize <- paste(JHeapSize,"m",sep="")
	java <- paste(java," -Xms",JHeapSize," -Xmx",JHeapSize," -cp",sep="")
	mzmatch <- paste(java, " ", .find.package("mzmatch.R"), "/java/mzmatch.jar", sep="")
	
	## setup the tool
	tool <- paste(mzmatch, "mzmatch.ipeak.Combine")
	if (!is.null(i))
		tool <- paste(tool, "-i", i)
	if (!is.null(o))
		tool <- paste(tool, "-o", o)
	if (!is.null(label))
		tool <- paste(tool, "-label", label)
	if (!is.null(labels))
		tool <- paste(tool, "-labels", labels)
	if (!is.null(ppm))
		tool <- paste(tool, "-ppm", ppm)
	if (!is.null(rtwindow))
		tool <- paste(tool, "-rtwindow", rtwindow)
	if (!is.null(combination))
		tool <- paste(tool, "-combination", combination)
	if (!is.null(h) && h==T)
		tool <- paste(tool, "-h")
	if (!is.null(v) && v==T)
		tool <- paste(tool, "-v")

	system(tool)
}
