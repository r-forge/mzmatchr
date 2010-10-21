\name{mzmatch.ipeak.util.ThermoStatusLog}
\alias{mzmatch.ipeak.util.ThermoStatusLog}
\title{mzmatch.ipeak.util.ThermoStatusLog}

\description{
Extracts the status-logs from the given Thermo RAW files (for convenience it also includes the trailer-extra logs) and outputs the logs of each controller in a separate file stored in the RAW file into a separate tab-delimited text-file.
}

\usage{
mzmatch.ipeak.util.ThermoStatusLog(JHeapSize=1425, i, o, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file(s). These should be in the Thermo Raw file format produced by Xcalibur.
}
\item{o}{
Option for the output destination. This is required to be a directory where all of the generated  files are written. When the directory (or path leading to the directory does not exist, it is created). For the filename convention used for the ouput files, see the general description.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Extracts the status-logs from the given Thermo RAW files (for convenience it also includes the trailer-extra logs) and outputs the logs of each controller in a separate file stored in the RAW file into a separate tab-delimited text-file. The columns contain the different types of logs stored and the rows the values for each scan.  The fact that a single RAW file can contain multiple controllers, necessitates that the output is a directory where multiple files can be written.

No warrenties are given about the contents, as this seems to be largely dependent on which machine was used.
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
