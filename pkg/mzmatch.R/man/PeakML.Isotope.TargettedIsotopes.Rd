\name{PeakML.Isotope.TargettedIsotopes}
\Rdversion{1.1}
\alias{PeakML.Isotope.TargettedIsotopes}
\title{PeakML.Isotope.TargettedIsotopes}
\description{
	Plot isotopes for all compounds listed in the molFormulaFile.
}
\usage{
	PeakML.Isotope.TargettedIsotopes(baseDir, sampleType = "NEG", molFormulaFile, outDirectory = "targettedIsotops", outFileName, mzXMLSrc=NULL, peakMLFile, sampleGroups = NULL, layoutMtx = NULL,massCorrection, ppm = 2, trendPlots = NULL, loadSavedData = TRUE, fillGaps, useArea=FALSE, stdRTWindow=NULL)
}

\arguments{
	\item{baseDir}{
		The parent directory where PeakML files are located.
	}
	\item{sampleType}{
		The sample type, eg., negative ionisation.
	}
	\item{molFormulaFile}{
		A tab separated file containing compound id, compound name, compound formula, followCarbon (the labelling of the carbon to be tracked) along the rows.
	}
	\item{outDirectory}{
		The output directory.
	}
	\item{outFileName}{
		The output filename
	}
	\item{mzXMLSrc}{
		The path to the raw data mzXML files.
	}
	\item{peakMLFile}{
		The name of the PeakML file which contains the combined peaks of all samples that are filtered, gap filled and identified. See the advanced metabolomics data processing tutorial at http://mzmatch.sourceforge.net/tutorial.mzmatch.r.advanced.html for more information.
	}
	\item{sampleGroups}{
		The sample groups.
	}
	\item{layoutMtx}{
		The layout matrix that specifies the layout of the plot to be generated. Leave this as NULL if the number of sample to be plotted is less than or equal to 10. See layout in R more details.
	}
	\item{massCorrection}{
		Requied if gap fill is TRUE. This is to correct the natural mass with the ionised mass. See PeakML.Methods.getMassCorrection(). see PeakML.Methods.getMassCorrection. 
	}
	\item{ppm}{
		The ppm to specify the mass window for searching the isotopes.
	}
	\item{trendPlots}{
		Leave it to NULL if the default trend plots have to be plots. Else specify the name of the custom trend plots as a vector. Corresponding plot function has to be incorporated into the PeakML.Isotope.plotSamples function.
	}
	\item{loadSavedData}{
		Set it TRUE (normally when the script is run for the first time) if the peak data matrix and chrom data list are to be directly loaded from the working directory.
	}
	\item{fillGaps}{
		"ALLPEAKS" if all the peaks are to gap filled, requied to fill missing signals in some part of the curve. "MISSINGPEAKS", if only the missing peaks have to be filled Or "NONE" if no peaks have to be gap filled. Defaults to "ALLPEAKS".
	}
	\item{useArea}{
		TRUE if sum of intensities has to be used instead of the maximum of intensities for peak intensity calculations.
	}
	\item{stdRTWindow}{
		The standard retention time  window. Leave it to NULL if no standard retention time is available for the compound.
	}
}


\details{
The PeakML file format defines an open and extensible file format for storing extracted hyphenated mass spectrometry data. As such, the format will enable users to store intermediate data in a structured manner, but more importantly, to exchange data between different data analysis software. This will give complete flexible control to the users of data analysis software to cherry pick those components best suited to their needs. The PeakML.Methods contains a set of R functions that will enable end users to rapidly write scripts to manipulate PeakML files.
}

\value{
	PeakML.Isotope.TargettedIsotopes generates a pdf document containing the plots of isotopes of each compound listed in the file (molFormulaFile) specifed.
}

\references{
	PeakML/mzMatch: A File Format, Java Library, R Library, and Tool-Chain for Mass Spectrometry Data Analysis, Richard A. Scheltema, Andris Jankevics, Ritsert C. Jansen, Morris A. Swertz, Rainer Breitling, Analytical Chemistry, 2011 83 (7), 2786-2793.

}
\author{
	Achuthanuni Chokkathukalam (unni.chokkathukalam@glasgow.ac.uk)
	Andris Jankevics (a.jankevics@rug.nl)
}

\seealso{
	\code{\link{PeakML.Isotope.getIsotopes, PeakML.Isotope.getChromData, PeakML.Isotope.plotSamples, PeakML.Isotope.processTargettedIsotopes}}
}

\keyword{ mzmatch.R }
\keyword{ peakML }
