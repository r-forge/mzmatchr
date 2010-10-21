\name{mzmatch.ipeak.db.CreateKEGGDB}
\alias{mzmatch.ipeak.db.CreateKEGGDB}
\title{mzmatch.ipeak.db.CreateKEGGDB}

\description{
General description of the function of the tool.
}

\usage{
mzmatch.ipeak.db.CreateKEGGDB(JHeapSize=1425, i, inchi, o, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
KEGG compounds data file in plain text format. Can be downloaded from: ftp://ftp.genome.jp/pub/kegg/ligand/compound/
}
\item{inchi}{
InCHI strings files. Can be downloaded from: ftp://ftp.genome.jp/pub/kegg/ligand/compound/
}
\item{o}{
InCHI strings files. Can be downloaded from: ftp://ftp.genome.jp/pub/kegg/ligand/compound/
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
General description of the function of the tool.
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
