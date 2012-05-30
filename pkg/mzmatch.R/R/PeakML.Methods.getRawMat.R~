PeakML.Methods.getRawMat <- function (rawfile,scan_start,scan_finis, mz_start, mz_finis,correctedRT,uncorrectedRT)
{
	Rawpeaks <- peaks(rawfile,scans=c(scan_start:scan_finis))
	rawRT <- uncorrectedRT[scan_start:scan_finis]
	cRT <- correctedRT[scan_start:scan_finis]
	scans <- scan_start:scan_finis

	peakExtract <- function (ind)
	{
		MZpeakset <- Rawpeaks[[ind]]
		hit <- which(MZpeakset[,1] >= mz_start & MZpeakset[,1] <= mz_finis)
		# if more than one peak match the criteria, mass for highest intense peak is used as output.
		if (length(hit)>0)
		{
			dout <- MZpeakset[hit,]
			dout <- rbind(dout,NULL)
			out <- c(rawRT[ind],cRT[ind],scans[ind],mean(dout[which(dout[,2]==max(dout[,2]))[1],1]),max(dout[,2]))
		} else
		{
			out <- NULL
		}
		out
	}

	RawMat <- do.call(rbind,lapply (1:length(Rawpeaks),peakExtract))
	if (is.null(RawMat))
	{
		RawMat <- c(1,1,1,1,1)
	}
	RawMat <- rbind(RawMat,NULL)
	if (!is.null(RawMat))
	{
		colnames (RawMat) <- c("rawRT","correctedRT","scan_id","mz","intensity")
	}
	RawMat
}
