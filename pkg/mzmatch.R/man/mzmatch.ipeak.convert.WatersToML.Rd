\name{mzmatch.ipeak.convert.WatersToML}
\alias{mzmatch.ipeak.convert.WatersToML}
\title{mzmatch.ipeak.convert.WatersToML}

\description{
Convert Waters raw data files to mzXML/mzDATA files.
}

\usage{
mzmatch.ipeak.convert.WatersToML(JHeapSize=1425, i, o, format, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
A list of valid Waters RAW files to be converted to the chosen format.
}
\item{o}{
The directory where to output the files. For each filterline encountered in the RAW file, a new sub-directory will be created. In the case of multiple controllers, the controller-number will also be inserted in the path.
}
\item{format}{
The output format to use. The current choices are limited to mzXML and mzML as these are the most used.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Convert Waters raw data files to mzXML/mzDATA files.
}

\value{
This function returns no value.
}

\references{
0. PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.


}

\author{
RA Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
\code{\link{PeakML.xcms.read}}
}

\keyword{mzMatch}
\keyword{peakML}
\keyword{java tools}
