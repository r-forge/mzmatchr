PeakML.Methods.getCarbon <- function(formula){
	# Returns the number of carbons in a compounds from the formula
	
	carbons <- 0
	formula <- as.character(formula)
	index <- grep("C", strsplit(formula,"")[[1]])
	if (!length(index)==0){
		subStr <- substr(formula, index, index + 2)
		carbons <- as.numeric(gsub("\\D", "", subStr))
	}
	carbons
}
