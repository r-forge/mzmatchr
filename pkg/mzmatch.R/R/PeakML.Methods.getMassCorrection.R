PeakML.Methods.getMassCorrection <- function(javaProject, ionisation="detect"){
	# PRE:
	# 	the jave project, ionisation mode (see PeakML.Methods.getProtonCoef)
	# POST:
	#	returns mass correction as numeric
	
	protonMass <- PeakML.Methods.getProtonMass()
	protonCoef <- PeakML.Methods.getProtonCoef(javaProject, ionisation)
	rv <- protonMass * protonCoef[[1]]
	rv
}
