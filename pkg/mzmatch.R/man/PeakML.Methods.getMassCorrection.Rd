\name{PeakML.Methods.getMassCorrection}
\Rdversion{1.1}
\alias{PeakML.Methods.getMassCorrection}
\title{PeakML.Methods.getMassCorrection}
\description{

	Returns the mass correction based on the ionisation mode. If ionisation argument is set to "detect", then the ionisation mode is detected from the first scan in the first sample. No ionisation switching is supported in mzmatch.R yet.
}
\usage{
	PeakML.Methods.getMassCorrection(javaProject, ionisation="detect")
}

\arguments{
	\item{javaProject}{
		The mzmatch java project
	}
	\item{ionisation}{
		Defaults to "detect", to detect the ionisation from the scan data. Other accepted values are "positive" (proton coef = 1), "negative" (proton coef. = -1) and empty string (proton coef=0).
	}
}

\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}
\value{
  getMassCorrection returns mass correction based on the ionisation mode, a numeric value.
}
\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
  Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
  Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
  \code{\link{PeakML.Methods.getProtonMass}, \link{PeakML.Methods.getProtonCoef}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }