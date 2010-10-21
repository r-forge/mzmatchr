\name{mzmatch.ipeak.ExtractMassChromatograms}
\alias{mzmatch.ipeak.ExtractMassChromatograms}
\title{mzmatch.ipeak.ExtractMassChromatograms}

\description{
Extracts mass chromatograms (x-axis: RT; y-axis: Intensity) from 2D mass spectrometry data (LC/MS or GC/MS).
}

\usage{
mzmatch.ipeak.ExtractMassChromatograms(JHeapSize=1425, i, o, masstraces, label, threshold, ppm, h, v)
}

\arguments{
\item{JHeapSize}{
	Amount of RAM memory in megabytes, which should be reserved by the Java virtual machine. The default value is 1425 megabytes.}
\item{i}{
Option for the input files, which should be in one of the open standard file formats (mzML, mzXML or mzData) and contain data from a 2D mass spectrometry setup (LC/MS or GC/MS).
When this option has not been set, the input is read from the stdin (allowing for pipeline building). When a single input file is defined, the output '-o' should contain the output filename. When multiple input files are defined, the output '-o' should define an output directory.
For now only centroid input data is supported.
}
\item{o}{
Option for the ouput file(s); refer to the option input '-i' for a description of behaviours with regards to multiple input files. The extracted mass chromatograms are written here in the PeakML format.
When this option has not been set the output is written to the standard output (works only when there is a single input file).Be sure to unset the verbose option when setting up a pipeline reading and writing from the standard in- and outputs.
}
\item{masstraces}{
Optional output file where the mass traces are written (only useful when the option 'threshold' has been defined), which can be used to debug the mass trace breakup approach.
}
\item{label}{
Optional label for the file, which will be stored in the header of the resulting file. The label is used for display purposes in UI environments.
}
\item{threshold}{
The percentage threshold value for breaking the mass traces up, as a percentage of the most intense portion of a mass trace. The threshold value is a percentage and required to be between 0 and 1.
}
\item{ppm}{
The accuracy of the measurement in parts-per-milion. This value is used for the collection of the data-points belonging to a mass trace and needs to be reasonable for the equipment used to make the measurement (the LTQ-Orbitrap manages approximatetly 3 ppm).
}
\item{h}{
When this is set, the help is shown.
}
\item{v}{
When this is set, the progress is shown on the standard output.
}

}

\details{
Extracts mass chromatograms (x-axis: RT; y-axis: Intensity) from 2D mass spectrometry data (LC/MS or GC/MS). The raw data is loaded from the open standard file formats (mzML, mzXML or mzData) and all of the individidual mass traces (M/Z +/- ppm over the whole scan range) are retrieved. When the option 'threshold' is defined, the individual mass traces are broken up into individual mass chromatograms (ie the isomers are separated). This is achieved by cutting peaks out of the mass trace where the threshold is reached (as a percentage of the most intense portion of the mass trace). This is an iterative process, where the sides are then analyzed in the same fashion.

The method employed here for retrieving mass chromatograms is greedy and extracts everything (although a modicum of noise reduction is applied to reduce the amount of fragments from broken up mass chromatograms). In order to reduce the resulting noise patterns, tools like 'mzmatch.filter.NoiseFilter' and 'mzmatch.filter.RSDFilter' can be employed.

The resulting output file is in PeakML-format, containing a list of all the extracted mass chromatograms. When the cutoff value has been selected, one can also specify a file (option 'masstraces') where the extracted mass chromatograms are overlaid on the mass traces they have been cut from.

Remarks
1. At this time only centroid data is supported.
2. NetCDF is not supported as it misses necessary meta-information
3. Direct injection data will not yield correct results

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
