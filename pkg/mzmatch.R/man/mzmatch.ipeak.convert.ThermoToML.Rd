\name{mzmatch.ipeak.convert.ThermoToML}
\alias{mzmatch.ipeak.convert.ThermoToML}
\title{mzmatch.ipeak.convert.ThermoToML}

\description{
Converts the contents of a Thermo RAW file to one or more open standard files (currently only mzXML and mzML are supported).
}

\usage{
mzmatch.ipeak.convert.ThermoToML(JHeapSize=1425, i, o, format, centroid, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
A list of valid Thermo RAW files to be converted to the chosen format.
}
\item{o}{
The directory where to output the files. For each filterline encountered in the RAW file, a new sub-directory will be created. In the case of multiple controllers, the controller-number will also be inserted in the path.
}
\item{format}{
The output format to use. The current choices are limited to mzXML and mzML as these are the most used.
}
\item{centroid}{
When this is set, the data will be centroided. When the data is already in centroid mode this will have no effect.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Converts the contents of a Thermo RAW file to one or more open standard files (currently only mzXML and mzML are supported). Thermo RAW files can contain multiple controllers and multiple methods (separated by filterline). These are all split into separate files, distinguished by directories with the controller-number and the filterline (i.e. the original filename is maintained, so tools working with a filename database will continue to work correctly). When only one controller is present in the file, the controller-number is omited in the directory naming.

The data can be centroided by this application, for which the functionality in the Thermo RAWfile DLL is used. To do so, specify the commandline parameter centroid.

WARNINGS:
* This tool makes direct use of the Thermo Rawfile DLL and can only be used on Microsoft Windows platforms.
* Due to a bug in the Thermo DLL no more than 8 files can be processed in a single run
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
