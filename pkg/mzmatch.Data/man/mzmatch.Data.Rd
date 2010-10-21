\name{mzmatch.Data}
\alias{mzmatch.Data}
\docType{data}
\title{Sample Orbitrap LC/MS data}
\description{
Sample Orbitrap LC/MS data files acquired in Positive ionisation mode on Zic HILIC column. Standart amino acids mixtures were analysed at three different concentrations. Data files are organised in three folders where each folder represents a sample concentration. Every folder contains three files which are three technical replicates of the single analytical sample.
}

\details{
  The corresponding raw mzXML files are located in the \code{orbitrap} subdirectory of this package.
}

\examples{
## The directory with the mzXML LC/MS files
mzdatapath <- file.path(.find.package("mzmatch.Data"), "orbitrap")
mzdatapath
files <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE)

if (require(xcms)) {

## Access raw data file
x <- xcmsRaw(files[1])
x

}
}

\keyword{datasets}
