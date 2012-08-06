\name{PeakML.ComparePeakSets}
\Rdversion{1.1}
\alias{PeakML.ComparePeakSets}
\title{PeakML.ComparePeakSets}
\description{
 This function will filter peak sets from one peakml file matching similar peak sets (based on RT and mass window) in the second file.
}
\usage{
PeakML.ComparePeakSets (standard_filename, filename, stdionisation="detect",ionisation="detect", Rawpath=NULL, outputfile,ppm=5,rtwin=20)
}
\arguments{
\item{standard_filename}{
	Name ot the standard input file in peakML format.
}
\item{filename}{
	Name ot the file which has to be processed in peakML format.
}
\item{stdionisation}{
	Defaulted to "detect". See documentation for PeakML.Methods.getProtonCoef.
}
\item{ionisation}{
	Defaulted to "detect". See documentation for PeakML.Methods.getProtonCoef.
}
\item{Rawpath}{
	Path to the raw data if different from the one stored in the PeakML file.
}
\item{outputfile}{
	Name of the resulting peakML file.
}
\item{ppm}{
	Mass window to much peak grops in ppm.
}
\item{rtwin}{
	Retention time window to match peak groups in seconds.
}
}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.Methods.getGroupAnnotations function allows to access attributes defined in PeakML file.
}
\value{
	List of 3 PeakML data files: matched peak sets, not matched from processed file (file prefix Not_matched) and not matched peak sets from standards file (file prefix STD_not_matched).
}
\references{
  PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.Read}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
