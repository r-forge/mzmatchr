mzmatch.init <- function (memorysize=1024)
{	
	lib <- paste(.find.package("mzmatch.R"),"/java/mzmatch.jar",sep="")
	params <- paste("-Xmx",memorysize,"m",sep="")
	.jinit(classpath=lib, parameters=params, force.init=TRUE)
}
