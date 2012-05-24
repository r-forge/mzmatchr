PeakML.Isotope.getFCMtxAbun <- function (trendList, sampleGroups, followCarbon){

	numCarbons <-  length(trendList[[1]])
	plotMtx <- matrix(nrow = 3, ncol = length(sampleGroups))
	dimnames(plotMtx) <- list(c("bottom", "top1", "top2"), sampleGroups)
	#colVector <- plotMtx
	
	N <- numCarbons
	M <- followCarbon - 1
	P <- 0.0112 # natural abundance of Carbon 1.12%
	
	NatAbuC <- choose(N,M)*(P)^M*((1-P)^(N-M))
	
	for (sam in 1:length(sampleGroups)){
		for (row in 1:length(rownames(plotMtx))){			# This is to plot the natural abundance

			x <- trendList[[sam]][[followCarbon]]
			y <- trendList[[sam]][[1]]
			
			if (is.null(x)){
				x <- 0
			}
			if (is.null(y)){
				y <- 0
			}
			
			z <- y * NatAbuC

			bottom <- min(x,z)
			top <- max(x,z) - min(x,z)
			
			VAL <- 0.000
			if (row == 1){
				VAL <- bottom
			} else if(row == 2){
				if ((x-z)<=0){
					VAL <- top
				}
			} else if(row ==3){
				if ((x-z)>=0){
					VAL <- top
				}
			}

			plotMtx[row, sam] <- VAL
#			colVector[row,sam] <- colr
		}
	}
#	colVector <- as.vector(colVector)
#	list(plotMtx, colVector)
	plotMtx
}






#PeakML.Isotope.getFCMtxAbun <- function (trendList, sampleGroups, followCarbon){

#	numCarbons <-  length(trendList[[1]])
#	plotMtx <- matrix(nrow = 2, ncol = length(sampleGroups))
#	dimnames(plotMtx) <- list(c("bottom", "top"), sampleGroups)
#	colVector <- plotMtx
#	
#	N <- numCarbons
#	M <- followCarbon - 1
#	P <- 0.0112 # natural abundance of Carbon 1.12%
#	
#	NatAbuC <- choose(N,M)*(P)^M*((1-P)^(N-M))
#	
#	for (sam in 1:length(sampleGroups)){
#		for (row in 1:length(rownames(plotMtx))){			# This is to plot the natural abundance

#			x <- trendList[[sam]][[followCarbon]]
#			y <- trendList[[sam]][[1]]
#			
#			if (is.null(x)){
#				x <- 0
#			}
#			if (is.null(y)){
#				y <- 0
#			}
#			
#			z <- y * NatAbuC

#			bottom <- min(x,z)
#			top <- max(x,z) - min(x,z)
#			
#			if (row == 1){
#				VAL <- bottom
#				colr <- 1
#			} else if(row == 2){
#				VAL <- top
#				if ((x-z)<=0){
#					colr <- 2
#				} else {
#					colr <- 3
#				}
#			}
#			plotMtx[row, sam] <- VAL
#			colVector[row,sam] <- colr
#		}
#	}
#	colVector <- as.vector(colVector)
#	list(plotMtx, colVector)
#}
