\name{mzmatch.ipeak.util.Identify}
\alias{mzmatch.ipeak.util.Identify}
\title{mzmatch.ipeak.util.Identify}

\description{
Identifies the contents of the given PeakML file with the given databases. The databases are expected to have the format of the example below (a standard file format for this would be preferable) and contain all of the compounds to be tested for. Within the tool-chain, in the package 'mzmatch.ipeak.db', several tools are provided for converting downloadable files (usually on an ftp-server) from the major metabolite databases. Before identifying your files, please make sure you have the most recent version.
Matching is performed only on mass, which is taken from the topmost structure. In other words, if the PeakML file contains a list of mass chromatograms the mass of each individual mass chromatogram is matched to the database. However, if the PeakML file contains a list of matched mass chromatograms, the mean mass of the matched mass chromatograms is used for identification.
}

\usage{
mzmatch.ipeak.util.Identify(JHeapSize=1425, i, o, ppm, databases, minrt, maxrt, rtwindow, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input file.
}
\item{o}{
Option for the output file.
}
\item{rtwindow}{
Optional retention time window for finding matches from the databases. If this is set and the database contains entries with previously registered retention times for molecules. When this value is not set, the stored retention times are ignored.
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for matching the masses to those found in the supplied databases. This value is obligitory.
}
\item{databases}{
Option for the molecule databases to match the contents of the input file to. These files should adhere to the compound-xml format.
}
\item{minrt}{
Optional minimum retention time for excluding signals from the input-file. This is for example convenient for excluding lipids on an LC/MS setup with a HILIC column.
}
\item{maxrt}{
Optional maximum retention time for excluding signals from the input-file. This is for example convenient for including lipids on an LC/MS setup with a HILIC column.
}
}

\details{
When a peak is positively identified the annotation 'identification' is extended with the unique database ID corresponding to the match. This will keep the clutter in the PeakML file to a minimum and the information associated to the tag up-to-date when the database is updated with a new version. Additionally, the tag provides a convenient wayof removing false identifications from the PeakML file from a UI environment.
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

\keyword{mzmatch}
\keyword{peakml}
\keyword{java tools}
