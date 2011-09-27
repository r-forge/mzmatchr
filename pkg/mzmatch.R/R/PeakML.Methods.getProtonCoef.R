PeakML.Methods.getProtonCoef <- function(javaProject, ionisation="detect"){
	# PRE:
	# 	Java project, inonisation (NULL, detect, negative, positive)
	# POST:
	#	protonCoef, mass correction factor
	
	#Detect ionisation mode, only first scan from first sample is used, no ionisation switching is supported yet
	if (ionisation=="detect") {
		ionisation<- .jcall (javaProject,return="S",method="getIonisation",as.integer(0),as.integer(0))
	}

	if (ionisation=="positive") {
		protonCoef <- 1
	} else if (ionisation=="negative") {
		protonCoef <- -1
	} else	protonCoef <- 0
	
	list(protonCoef,ionisation)
}
