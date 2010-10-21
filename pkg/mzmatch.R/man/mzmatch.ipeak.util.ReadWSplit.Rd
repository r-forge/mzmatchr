\name{mzmatch.ipeak.util.ReadWSplit}
\alias{mzmatch.ipeak.util.ReadWSplit}
\title{mzmatch.ipeak.util.ReadWSplit}

\description{
This tool handles the mess ReadW makes of the new RAW-file versions (as of Xcalibur 4.
}

\usage{
mzmatch.ipeak.util.ReadWSplit(JHeapSize=1425, i, o, fullms, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input files. These files should be in the mzXML file format produced by ReadW together with Xcalibur >=4.2.
}
\item{o}{
Option for the output destination. This is required to be a directory where all of the generated  files are written. When the directory (or path leading to the directory does not exist, it is created). For the filename convention used for the ouput files, see the general description.
}
\item{fullms}{
When this is set, the 'full lock ms' and 'full ms' scans are combined in a single file. This is specific for exactive platforms.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
This tool handles the mess ReadW makes of the new RAW-file versions (as of Xcalibur 4.2). With these new version Thermo is giving access to the both the lockmass corrected and the original data. With the advent of the exactive there can now possibly also be both negative as well as positive scans in a single file. In order to circumvent problems downstream of current pipeline software, each file is split into all of its different contents based on the FilterLine (used in Xcalibur to differentiate between the different data streams).

The input file is split up into all of the unique filter-lines (Xcalibur specific). In the case of the exactive there are 4 different lines: Full/neg, Full/pos, Lock/neg and Lock/pos, which represent the different datastreams (Full=collected data; Lock=lockmass corrected data). These are then written to separate files (Lock_neg_'input filename', etc). For analysis it is advised to only use the Lock-files, as these will provide the best mass accuracy.

Remark
This tool was specifically made for files produced with the Orbitrap Exactive, but may have application beyond this platform.
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
