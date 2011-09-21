RCreateXMLDB <- function (data,outputfile)
{
	require (XML)
	## These are default node names in database
	defnames <- c("id","name","formula","inchi","description","synonyms")
	## detected in DB
	childnames <- colnames (data)
	## Assign column numbers to correct nodes, unset nodes set to NULL
	id <- NULL
	name <- NULL
	defined <- which(defnames%in%childnames)
	if (length(defined!=0))
	{
		for (i in 1:length(defined))
		{
			assign(defnames[defined[i]],which(defnames[defined[i]]==childnames))
		}
	}

	undefined <- which(!defnames%in%childnames)
	if (length(undefined!=0))
	{
		for (i in 1:length(undefined))
		{
			assign(defnames[undefined[i]],NULL)
		}
	}

	if (is.null(id) | is.null(name) | is.null(formula))
	{
		cat ("Values for fields: id,name and formula must be set.\n")
	} else
	{
		## For every row in input data table, make child nodes with subnodes
		##XMLfile <- newXMLDoc ("treedoc")
		compounds <- xmlNode("compounds")
		nodeGen <- function(nodenum)
		{
			## These silly assignments are needed to avoid "NOTE" messages during R CMD check routine
			inchi <- inchi
			description <- description
			synonyms <- synonyms

			compound <- xmlNode ("compound",xmlNode ("id", data[nodenum,id]))
			compound[[2]] <- xmlNode ("name", data[nodenum,name])
			compound[[3]] <- xmlNode ("formula", data[nodenum,formula])
			if (is.null(inchi))
			{
				compound[[4]] <- xmlNode ("inchi")
			} else
			{
				compound[[4]] <- xmlNode ("formula", data[nodenum,inchi])
			}
			if (is.null(description))
			{
				compound[[5]] <- xmlNode ("description")
			} else
			{
				compound[[5]] <- xmlNode ("description", data[nodenum,description])
			}
			if (is.null(synonyms))
			{
				compound[[6]] <- xmlNode ("synonyms")
			} else
			{
				compound[[6]] <- xmlNode ("synonyms", data[nodenum,synonyms])
			}
			compound
		}

	system.time(compoundlist <- lapply(1:nrow(data),nodeGen))
	compounds <- addChildren(compounds,kids=compoundlist)
	saveXML(doc=compounds,file=outputfile,prefix="<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n \n")
	}
}

