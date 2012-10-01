###############################################################
###############################################################
###maxORmin
maxORminM <- function(matM, idCol, byCol, absolute=TRUE, decreasing=TRUE) {
	if (absolute) {
		matM <- matM[order(abs(matM[, byCol]), decreasing=decreasing),]
		matM <- matM[!duplicated(matM[,idCol]),]
	} else {
		matM <- matM[order(matM[, byCol], decreasing=decreasing),]
		matM <- matM[!duplicated(matM[,idCol]),]
	}
}



### ###############################################################
### ###############################################################
### ###max Variance
### maxVar <- function(matM, idCol, ...) {
### 	num <- sapply(matM,mode)%in%"numeric"
### 	data <- matM[, num, FALSE]
### 	dataVar <- apply(data, 1, var, ...)
### 	matM <- matM[order(dataVar, decreasing=TRUE),]
### 	matM <- matM[!duplicated(matM[, idCol]),]
### }


###############################################################
###############################################################
###geometric mean
geoMeanM <- function(x, ...) {
	if (any(x<0)) {
		stop("There are negative values, impossible to compute geometric Mean")
	} else {
		out <- exp(mean(log(x), ...))
	}
	return(as.numeric(out))
}


###############################################################
###############################################################
###random
randomM <- function(x, ...) {
	out <- sample(x, size=1, ...)
	return(as.numeric(out))
}


###############################################################
###############################################################
###mean
meanM <- function(x, ...) {
	out <- mean(x, ...)
	return(as.numeric(out))
}


###############################################################
###############################################################
###median
medianM <- function(x, ...) {
	out <- median(x, ...)
	return(as.numeric(out))
}


###############################################################
###############################################################
#####noRedundancy
filterRedundant <- function(object,
			    method=c("maxORmin", "geoMean", "mean", "median", "random"),
###			    method=c("maxORmin", "maxVar", "geoMean", "mean", "median", "random"),
			    idCol=1, byCol=2, absolute=TRUE, decreasing=TRUE,
			    trim=0, ...) {

	##check object
	if (! is.data.frame(object)) stop("filterRedundant works on data.frames")

	##stop if NA
	if ("na.rm" %in% names(match.call()) ) {
		stop("This function does not handle NA, please replace NA values with numbers")
	}

	##evaluate if method
	method <- match.arg(method)

	##check arguments: idCol
	if (is.character(idCol)) idCol <- which(colnames(object) %in% idCol)

	##check arguments: byCol
	if (is.character(byCol)) byCol <- which(colnames(object) %in% byCol)

	##check arguments: byCol and idCol
	if (idCol == byCol || idCol > ncol(object) || byCol > ncol(object)
	    || length(idCol)==0 || length(byCol)==0) {
		stop(paste("Provide valid and distinct 'idCol' and 'byCol' paramenters!",
			   "\n Valid 'idCol' and 'byCol' are of mode 'character' or 'numeric'."))
	}

	##check if there is redundancy
	if (length(object[,idCol]) != length(unique(object[,idCol]))) {

		## select maxOrMin for each feature/identifier
		if (method == "maxORmin") {
			mOut <- maxORminM(object, idCol, byCol, absolute, decreasing)
		## select most Variant for each feature/identifier
		}
### 		else if (method == "maxVar") {
### 			mOut <- maxVar(object, idCol, byCol, ...)
### 		}
### 		## use other methods
		else {
			## identifiers as indexes
 			ids <- object[, idCol]
			whichNumeric <- sapply(object, mode) == "numeric"
			mOut <- sapply(object , function(x, y) {
				##NUMERIC VALUES
				if (is.numeric(x)) {
					tapply(X=x, INDEX=y, FUN=function(num) {
						if (length(num) > 1) {
							## select method to remove redundancy
							num <- switch(method,
								      geoMean = geoMeanM(num, trim, ...),
								      mean = meanM(num, trim, ...),
								      random = randomM(num, ...),
								      median = medianM(num, ...),
								      )
						} else {
							num
						}
						return(num)
					}, ...)
					##NON-NUMERIC VALUES
				} else {
					tapply(X=x, INDEX=y, FUN=function(x) x[1] )
				}
			} , y=ids)
			##process numeric
			mOut <- data.frame(mOut, stringsAsFactors=FALSE)
			mOut[,whichNumeric] <- mOut[, whichNumeric] <- apply(mOut[, whichNumeric], 2, as.numeric)
		}
	}  else {
		warning("No redundant feature identifiers were found")
		mOut <- object
	}
	return(mOut)
}


### ###testing object
### mat <- data.frame(A=1:10, B=2:11, C=100:99,
### 		  D=200:209, neg1=-1*10:1, neg2=100:99*-2,
### 		  oo=letters[1:10], cc=sample(letters,10),
### 		  ID=c(LETTERS[1:5], LETTERS[1:5]), stringsAsFactors=FALSE)



###############################################################
###############################################################
### A new function to merge list of data.frames
### If all dataframes contain the same genes, than a quick do.call after reordering will do
### listOfDataFrames=list contaning matrices with common rownames
### idCol=list contaning matrices with common rownames

mergeData <- function(listOfDataFrames, idCol=1, byCol=2) {

	##is a list
	if (!is.list(listOfDataFrames)) {
		stop("Use a list containing matrices or data.frames to be merged")
	}

	##is the list of the correct length
	if (length(listOfDataFrames)<2) {
		stop("A minimun of 2 data.frames are needed")
	} else  {
		n <- length(listOfDataFrames)
	}

	##set the rownames using selected identifiers
	if (length(idCol)==1) {
		##message
		msg1 <- "Please provide a valid column name or index to select the desired 'dCol' identifiers column"
		##using column names as column identifers
		if (is.character(idCol)) {
			idCol <- sapply(listOfDataFrames,
					function(x,y) which(colnames(x)%in%y), y=idCol)
		##using column indexes as column identifers
		} else if (is.numeric(idCol)) {
			idCol <- rep(idCol, length(listOfDataFrames))
		}		

		## check if there are NA
		if (any(is.na(idCol))) {
			stop(msg1)
		}

		##set identifiers as rownames identifiers if non-redundant
		tmp <- mapply(x=listOfDataFrames, y=idCol, MoreArgs=list(msg=msg1),
			      FUN=function(x, y, msg){
				      if ( all(!duplicated(x[,y])) ) {
					      rownames(x) <- x[,y]
					      return(x)
				      } else {
					      stop(msg)
				      }
				      }, SIMPLIFY=FALSE)
		} else {
			stop(msg1)
		}

	##get ranking columns index and subset with it
	if (length(byCol)==1) {
		##message
		msg2 <- "Please provide a valid column name or index to select the desired 'byCol' value column"
		##using column names as column identifers
		if (is.character(byCol)) {
			byCol <- sapply(tmp,
					function(x, y) which(colnames(x)%in%y), y=byCol)
		##using column indexes as column identifers
		} else if (is.numeric(byCol)) {
			byCol <- rep(byCol, length(tmp))
		}

		## check if there are NA
		if (any(is.na(byCol))) {
			stop(msg2)
		}

		##keep selected columns for each data.frame in the list
		tmp <- mapply(x=tmp, y=byCol, MoreArgs=list(msg=msg2),
			      FUN=function(x, y, msg) {
				      if (y <= ncol(x)) {
					      x <- x[ , y, drop=FALSE]
					      return(x)
				      } else {
					      stop(msg)
				      }
			      },  SIMPLIFY=FALSE)
	} else {
		stop(msg2)
	}

	##to merge get all the unique rownames
	x <- tmp
	xID <- lapply(x, rownames)  #do not use sapply()
	allID <- unique(unlist(xID))
	sel <- apply(sapply(1:n, function(x,y,z) y%in%z[[x]], z=xID, y=allID), 1, sum) == n

	##ids intersection: common ids
	commonID <- allID[sel]

	## get all the unique of the matrices to merge
	xNms <- names(x)

	## subset with intesection genes
	x <- lapply(x, function(x,y) x[rownames(x)%in%y,,FALSE],y=commonID)

	## reorder by rownames
	x <- lapply(x, function(x) x <- x[order(rownames(x)),,FALSE])

	## make data.frame using cbind
	x <- do.call("cbind",x)

	##make new column names
	colnames(x) <- paste(xNms,colnames(x),sep=".")

	##add common ids as first column
	x <- data.frame(commonID=rownames(x),x,stringsAsFactors=FALSE)

	##return
	return(x)
}



###############################################################
###############################################################
####test for even numbers
is.even <- function(x) {
	x%%2 == 0
}

###############################################################
###############################################################
####test for odd numbers
is.odd <- function(x) {
	x%%2 == 1
}
