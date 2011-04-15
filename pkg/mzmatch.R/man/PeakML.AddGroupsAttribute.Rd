\name{PeakML.AddGroupsAttribute}
\Rdversion{1.1}
\alias{PeakML.AddGroupsAttribute}
\title{PeakML.AddGroupsAttribute}
\description{
  This function appends attributes to already existing peakml data file.
}
\usage{
  PeakML.AddGroupsAttribute (filename,outputfile,GroupAttributes=NULL,ionisation="detect",GrouptAttrName,GroupAttValues)
}
\arguments{
\item{filename}{
  Input file.
}
\item{outputfile}{
  The filename for the PeakML file. The path leading up to this file needs to exist, otherwise the function
  will fail.
}
\item{GroupAttributes}{
  Character vector indicating what existing attributes in peakml file must be stored. Input is a character vector containing attribute names. 
}
\item{ionisation}{
  Character value to be set on to "positive" if positive ionisation mode data are processed. And "negative" if negative ionisaton mode was used. Default value is seto to "detect", ionisation mode will be detected from the RAW data file.
}
\item{GrouptAttrName}{
  Names of the new group attributes. Character vector.
}
\item{GroupAttValues}{
  Matrix or vector of attribute values which have to be added to the peakml file. Number of columns must match the number of the peaksets in the file.
}


}
\details{
  The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry
  data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly,
  to exchange data between different data analysis software. This will give complete flexible control to the users of
  data analysis software to cherry pick those components best suited to their needs.
  
  The PeakML.xcms.write function extends the XCMS package with write functionality for exporting data to other applications.
  For example, the results of an XCMS step can then be loaded in the user interface environment PeakML Viewer and manipulated
  there.
}
\value{
  This function returns no value.
}
\references{
  PeakML/mzMatch - a file format, Java library, R library and tool-chain for mass spectrometry data analysis. In preparation.
}
\author{
  Andris Jankevics (a.jankevics@rug.nl)
  Richard Scheltema (r.a.scheltema@rug.nl)
}

\seealso{
  \code{\link{PeakML.xcms.read}}
}

\keyword{mzmatch.R}
\keyword{ XCMS }
