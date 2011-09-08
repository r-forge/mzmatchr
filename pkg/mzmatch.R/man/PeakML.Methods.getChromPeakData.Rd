\name{PeakML.Methods.getChromPeakData}
\Rdversion{1.1}
\alias{PeakML.Methods.getChromPeakData}
\title{PeakML.Methods.getChromPeakData}
\description{
	This function returns a list containing peak data matrix, list of chromatograms for each peak in the peak data matrix, sample classes, measurement names and mass correction applied to the data set.
}
\usage{
	PeakML.Methods.getChromPeakData(filename, ionisation="detect")
}

\arguments{
	\item{filename}{
		The PeakML file containing peak sets.
	}
	\item{ionisation}{
		Defaulted to "detect". See documentation for PeakML.Methods.getProtonCoef.
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
	R list containing peak data matrix, list of chromatograms for each peak in the peak data matrix, sample classes, measurement names and mass correction applied to the data set.
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getMassCorrection}, \link{PeakML.Methods.getPeakData}, \link{PeakML.Methods.getChromData}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
