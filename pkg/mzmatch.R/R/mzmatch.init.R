mzmatch.init <- function (memorysize=1024)
{	
	DBS.path <- paste(.find.package("mzmatch.R"),"/dbs",sep="")
	mzmatch.path <- paste(.find.package("mzmatch.R"),"/java",sep="")

	# Install required JAR files and XML databeses on the first run
	if (!file.exists(DBS.path))
	{
		current.dir <- getwd ()
		setwd (.find.package("mzmatch.R"))
		download.file("http://sourceforge.net/projects/mzmatch/files/mzmatch.R/dbs.zip/download","dbs.zip",cacheOK=FALSE,mode="wb")
		unzip ("dbs.zip")
		file.remove ("dbs.zip")
		setwd (current.dir)
	}
	
	if (!file.exists(mzmatch.path))
	{
		current.dir <- getwd ()
		dir.create (mzmatch.path)
		setwd (mzmatch.path)
		download.file("http://sourceforge.net/projects/mzmatch/files/mzmatch.R/mzmatch.jar/download","mzmatch.jar",cacheOK=FALSE,mode="wb")
		setwd (current.dir)
	}

	lib <- paste(.find.package("mzmatch.R"),"/java/mzmatch.jar",sep="")
	params <- paste("-Xmx",memorysize,"m",sep="")
	.jinit(classpath=lib, parameters=params, force.init=TRUE)
}
