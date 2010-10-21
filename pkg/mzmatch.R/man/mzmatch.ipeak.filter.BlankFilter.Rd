\name{mzmatch.ipeak.filter.BlankFilter}
\alias{mzmatch.ipeak.filter.BlankFilter}
\title{mzmatch.ipeak.filter.BlankFilter}

\description{
With this tool a blank measurement can be subtracted from a normal measurement.
}

\usage{
mzmatch.ipeak.filter.BlankFilter(JHeapSize=1425, i, blank, o, rejected, ppm, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
The input file (in PeakML format) containing the peaks from the normal measurement. When this option is not set the input is read from the standard input.
}
\item{blank}{
The blank file (in PeakML format) containing the peaks from the blank measurement. This option is required.
}
\item{o}{
The output file (in PeakML format) containing all the peaks from the normal measurement, which did not have a match in the blank measurement. When this option is not set the output is written to the standard output.
}
\item{rejected}{
Optional file (in PeakML format) containing all the rejected peaks from the normal measurement plus their match from the blank measurement. This file can be used for inspection and scavenging purposes.
}
\item{ppm}{
The precision of the data in parts-per-million. This accuracy value is used for matching the peaks from the blank to the peaks in the measurement.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
With this tool a blank measurement can be subtracted from a normal measurement. Both measurements should be from the same batch and preferably the last blank run before the normal measurement should be subtracted. Peaks from the measurements are automatically matched together with the ppm-setting. When the difference in intensity level is less than 3-fold two peaks are considered a match and the peak from the normal measurement is rejected.

With the rejected option all the rejected peaks can be stored in a separate file. Both the normal measurement peaks as the matched blank peaks are stored in a set in this file for inspection.
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
