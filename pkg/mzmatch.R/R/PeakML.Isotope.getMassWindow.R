PeakML.Isotope.getMassWindow <- function(mass, carbons, ppm){
	# PRE: 
	# 	mass : unlabelled mass of the compound
	#	carbons: number of labelled carbons
	#	ppm: the ppm window required.
	# POST: 
	# 	Mass window of isotop
	
	c13mass <- 1.0033 * carbons
	rv <- PeakML.Methods.getPPMWindow(mass + c13mass, ppm)
	rv
}

