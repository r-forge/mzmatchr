PeakML.RankProducts <- function (sampleGroups, inputTable, groupA, groupB, numberOfPermutations=1000, calculateProduct=TRUE, MinNumOfValidPairs=NA, replaceNA=FALSE, RandomPairs=NA)
{
	## Class labels
	HITS <- c(which(sampleGroups==groupA), which(sampleGroups==groupB))


	# number of pairs equal to the number of the smaples in the smallest class
	npairs <- min(length(which(sampleGroups==groupA)),length(which(sampleGroups==groupB)))


	## Number of RandomPairs to generate, if set to NA, n*n is used
	if (is.na(RandomPairs))
	{
		RandomPairs <- npairs^2
	}

	if (is.na(MinNumOfValidPairs))
	{
		MinNumOfValidPairs <- npairs/2
	}

	sampleGroups <- sampleGroups[c(HITS)]
	IntensitiesData <- inputTable[c(HITS),]

	## Estimate detection limit value, to use as a replacement for NA's, if such option is set.
	detection.limit <- min(IntensitiesData[IntensitiesData!=0],na.rm=TRUE)

	## Prepare ratios for paired samples
	gr1 <- IntensitiesData[sampleGroups==groupA,]
	gr2 <- IntensitiesData[sampleGroups==groupB,]

	## If in one group there was not detection in all samples, but there are signals detected in group B, values are replaced with detection.limit value
	## Also check which metabolites doesn not match MinNumOfValidPairs criteria. 
	## Variables gr1NA and gr2NA will store information about which measurements don't contain values and those will be ignored during random pairing.

	badrows <- NULL
	gr1NA <- vector ("list", ncol(gr1))
	gr2NA <- vector ("list", ncol(gr2))

	for (coln in 1:ncol(gr1)) 
	{
		toreplacegr1 <- which(is.na(gr1[,coln]))
		toreplacegr2 <- which(is.na(gr2[,coln]))

		gr1NA[[coln]] <- toreplacegr1
		gr2NA[[coln]] <- toreplacegr2

		toreplacegr1 <- length(toreplacegr1)
		toreplacegr2 <- length(toreplacegr2)

		if (MinNumOfValidPairs > min((c(nrow(gr1)-toreplacegr1,nrow(gr2)-toreplacegr2))))
		{
			badrows <- append(badrows,coln)
			gr1NA[[coln]] <- NULL
			gr2NA[[coln]] <- NA
		}

		if (toreplacegr1==nrow(gr1)) toreplacegr1 <- 1 else toreplacegr1 <- 0
		if (toreplacegr2==nrow(gr2)) toreplacegr2 <- 1 else toreplacegr2 <- 0
		if (toreplacegr1==1 & toreplacegr2==0)
		{
			gr1[,coln] <- detection.limit
		}
		if (toreplacegr2==1 & toreplacegr1==0)
		{
			gr2[,coln] <- detection.limit
		}
	}

	# If otion is set all NA left are replaced with detection.limit value.
	if (replaceNA==TRUE)
	{
		hh <- which(is.na(gr1))
		if (length(hh)!=0)
		{
			gr1[hh] <- detection.limit
		}
		hh <- which(is.na(gr2))
		if (length(hh)!=0)
		{
			gr2[hh] <- detection.limit
		}
	}

	## Create list of random sample matching, and ratios. If RandomPairs is set to 1, original sample order is used. The first ratios are calculated using initial sample order. If 'RandomPairs' is set to 1, only initial sample order is used for pairing.

	ratios.list <- vector("list",RandomPairs)
	if (calculateProduct==TRUE)
	{
		ratios.list.Class2 <- vector("list",RandomPairs)
		for(i in 1:RandomPairs)
		{
			ratios.list.Class2[[i]] <- matrix(data = NA, nrow = npairs, ncol = ncol(gr1))
		}
	}

	for(i in 1:RandomPairs)
	{
		ratios.list[[i]] <- matrix(data = NA, nrow = npairs, ncol = ncol(gr1))
	}

	ratios.list[[1]] <- gr1[1:npairs,]/gr2[1:npairs,]

	if (calculateProduct==TRUE)
	{
		ratios.list.Class2[[1]] <- gr2[1:npairs,]/gr1[1:npairs,]
	}

	for (coln in 1:ncol(gr1))
	{
		val1 <- gr1[,coln]
		val2 <- gr2[,coln]
		## get rid of NA's
		val1 <- sort(val1)
		val2 <- sort(val2)

		# calculate ratios of ranodmly paired samples
		for (rp in 2:RandomPairs)
		{
			val1gr <- sample(val1)[1:npairs]
			val2gr <- sample(val2)[1:npairs]
			ratios.list[[rp]][,coln] <- val1gr/val2gr
			if (calculateProduct==TRUE)
			{
				ratios.list.Class2[[rp]][,coln] <- val2gr/val1gr
			}
		}
	}


	## Rank products calculation
	prod_out <- list ()
	if (calculateProduct==TRUE)
	{
		prod_out.Class2 <- list()
	}
	
	for (randpair in 1:RandomPairs)
	{
		inputmatrix <- log (t(ratios.list[[randpair]]))

		## replace values in badrows with NA, after "order" function rank to those rows will be set to NA, and ignored in sum or prod calculation
		if (length(badrows)!=0)
		{
			inputmatrix[badrows,] <- NA
		}

		output <- matrix(ncol=ncol(inputmatrix),nrow=nrow(inputmatrix))

		if (calculateProduct==TRUE)
		{
			inputmatrix.Class2 <- log (t(ratios.list.Class2[[randpair]]))
			if (length(badrows)!=0)
			{
				inputmatrix.Class2[badrows,] <- NA
			}
			output.Class2 <- matrix(ncol=ncol(inputmatrix.Class2),nrow=nrow(inputmatrix.Class2))
		}

		for (coln in 1:ncol(inputmatrix))
		{
			# After assigning rank with "order" function, ranks for NA's are set back to NA
			NAs <- which(is.na(inputmatrix[,coln]))
			if (length(NAs)>0)
			{
				output[order(inputmatrix[,coln],na.last=NA),coln] <- c(1:(nrow(inputmatrix)-length(NAs)))
			} else
			{
				output[order(inputmatrix[,coln]),coln] <- c(1:nrow(inputmatrix))
			}
		}

		if (calculateProduct==TRUE)
		{
			for (coln in 1:ncol(inputmatrix.Class2))
			{
				NAs <- which(is.na(inputmatrix.Class2[,coln]))
				if (length(NAs)>0)
				{
					output.Class2[order(inputmatrix.Class2[,coln],na.last=NA),coln] <- c(1:(nrow(inputmatrix.Class2)-length(NAs)))
				} else
				{
					output.Class2[order(inputmatrix.Class2[,coln]),coln] <- c(1:nrow(inputmatrix.Class2))
				}
			}
		}

	
		if (calculateProduct==TRUE)
		{
			ProdFunc <- function(x, DATA) exp(mean(log(DATA[x,]),na.rm=TRUE))
			
			products <- sapply (1:nrow(output),ProdFunc, DATA=output)
			products.Class2 <- sapply (1:nrow(output.Class2),ProdFunc, DATA=output.Class2)
			prod_out.Class2[[randpair]] <- products.Class2
		} else
		{
			products <- apply(output,1,mean,na.rm=TRUE)
		}
		prod_out[[randpair]] <- products
	}

	prod_out <- do.call (rbind,prod_out)
	ranks_out <- apply(prod_out,2,median,na.rm=TRUE)
	
	if (calculateProduct==TRUE)
	{
		prod_out.Class2 <- do.call (rbind,prod_out.Class2)
		ranks_out.Class2 <- apply(prod_out.Class2,2,median,na.rm=TRUE)
	}

	## pfp permutations
	permutationsList <- vector("list",numberOfPermutations)
	if (calculateProduct==TRUE)
	{
		permutationsList.Class2 <- vector("list",numberOfPermutations)
	}

	for (permN in 1:numberOfPermutations)
	{
		
		cat (permN,"\n")
		permuteddata <- do.call(cbind,lapply(1:ncol(output),function (x) sample(output[,x])))
		if (calculateProduct==TRUE)
		{
			permuteddata.Class2 <- do.call(cbind,lapply(1:ncol(output.Class2),function (x) sample(output.Class2[,x])))
			if (length(badrows)>0)
			{
				permutationsList[[permN]] <- sapply (1:nrow(permuteddata),ProdFunc,DATA=permuteddata)[-c(badrows)]
			} else
			{
				permutationsList[[permN]] <- sapply (1:nrow(permuteddata),ProdFunc,DATA=permuteddata)
			}
			if (length(badrows)>0)
			{
				permutationsList.Class2[[permN]] <- sapply (1:nrow(permuteddata.Class2),ProdFunc,DATA=permuteddata.Class2)[-c(badrows)]
			} else
			{
				permutationsList.Class2[[permN]] <- sapply (1:nrow(permuteddata.Class2),ProdFunc,DATA=permuteddata.Class2)
			}
		} else
		{
			if (length(badrows)>0)
			{
				permutationsList[[permN]] <- apply(permuteddata,1,mean,na.rm=TRUE)[-c(badrows)]
			} else
			{
				permutationsList[[permN]] <- apply(permuteddata,1,mean,na.rm=TRUE)
			}
		}
	}

	permutationsList <- unlist(permutationsList)

	# erp - average expected value for the rank product. c/p. c - count how many times the rank products of the genes in the permutations are smaller or equal to the observed rank product. p - number of permutations
	# pfp - percentage of false positives. erp(g)/rank(g) where rank(g) is the rank of gene g in a list of all n genes sorted by increasing RP
	erps <- matrix(nrow=length(ranks_out),ncol=2)
	pfps <- matrix(nrow=length(ranks_out),ncol=2)

	if (calculateProduct==TRUE)
	{
		permutationsList.Class2 <- unlist(permutationsList.Class2)
		for (rn in 1:length(ranks_out))
		{
			erps[rn,1] <- length(which(permutationsList <= ranks_out[rn]))/numberOfPermutations
			erps[rn,2] <- length(which(permutationsList.Class2 <= ranks_out.Class2[rn]))/numberOfPermutations
		}
	} else
	{

		for (rn in 1:length(ranks_out))
		{
			erps[rn,1] <- length(which(permutationsList <= ranks_out[rn]))/numberOfPermutations
			erps[rn,2] <- length(which(permutationsList >= ranks_out[rn]))/numberOfPermutations
		}
	}

	if (length(badrows)!=0)
	{
		erps[badrows,] <- NA
	}
	colnames (erps) <- c("Class1","Class2")


	rankorder <- rep(NA,nrow(erps))
	rankorder[order(erps[,1],na.last=NA)] <- 1:length(order(erps[,1],na.last=NA))
	pfps[,1] <- erps[,1]/rankorder


	rankorder <- rep(NA,nrow(erps))
	rankorder[order(erps[,2],na.last=NA)] <- 1:length(order(erps[,2],na.last=NA))
	pfps[,2] <- erps[,2]/rankorder

	colnames (pfps) <- c("Class1","Class2")

	cat ("done  ")
	## Output results
	out <- list()
	if (calculateProduct==TRUE)
	{
		out$ranks <- cbind(ranks_out,ranks_out.Class2)
		colnames (pfps) <- c("Class1","Class2")
		out$RandomPairs_ranks.Class2 <- prod_out.Class2
		out$permutationsList.Class2 <- permutationsList.Class2
	} else
	{
		out$ranks <- ranks_out
	}
	out$groups <- sampleGroups
	out$sampleindex <- HITS
	out$erp <- erps
	out$pfp <- pfps
	out$RandomPairs_ranks <- prod_out
	out$permutationsList <- permutationsList
	out
}
