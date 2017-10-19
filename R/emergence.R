#==================================================================================
#==================================================================================
#==================================================================================
#' @title find emergence times for a species
#'
#' @description
#' @details
#' See Examples.
#'
#' @param futTs The future temperatures in the cell
#' @param focmaxT The T max of each species in the cell
#' @param runLength The number of consecutive years that temperatures must exceed focmaxT to qualify as emergence
#' @param prob The percent of subsequent years that must be emergent to qualify as emergence


# @keywords
#' @export
#'
# @examples
#'
#'
#' @return three emergence estimates.
#' \itemize{
#'  \item{ToEFirst }{year of first emergence - defined as when the temperature first exceeds focmaxT for a run of 'runLength' years }
#' \item{ToEPerm}{year of permanent emergence - if the temperature of the final year if emergent then it is when this final run of emergent years started}
#'  \item{ToEProb}{year when the probability of subsequent years being emergent equals or exceeds 'prob'}
# }

#' @author Alex Pigot <alex.pigot1@@gmail.com>,  Cory Merow
#' @note

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

findMoraToE<-function(focmaxT,futTs,yearsFut,runLength,prob){

  ToEFirst<-ToEPerm<-ToEProb<-NA

  emYears<-which(futTs>focmaxT)

  if(length(emYears)>0){

    x<-rep(0,length(yearsFut)) # create binary vector showing emergent (1) or non-emergent (0) years
    x[emYears]<-1

    #calculate the year when the temperature emerges for a period of at least 'runLength' years
    diffs <- x[-1L] != x[-length(x)] # find where the values change
    idx <- c(which(diffs), length(x)) # get the indexes, and then get the difference in subsequent indexes
    runs<-diff(c(0, idx))# calculate the length of the runs
    x2<-x[-which(c(NA,diff(x))==0)] # determine whether each run is a run of 1's or 0's
    yearsRunStart<-yearsFut[-which(c(NA,diff(x))==0)] # find the start year of each run
    runsTab<-data.frame(year=yearsRunStart,val=x2,run=runs) # store all this in a table
    runsTab<-runsTab[runsTab$val==1,] # remove runs of 0's
    runStart<-which(runsTab$run>=runLength)
    if(length(runStart)>0){
      ToEFirst<-min(runsTab$year[runStart]) #find the first year where there is a run of 1's equal to or greater than runLength
    }

    #calculate the time when the temperature emerges permamently
    if(x[length(x)]==1){
      nonemyears<-which(x==0)
      if(length(nonemyears)>0){
        ToEPerm<-1+yearsFut[max(nonemyears)]
      }else{
        ToEPerm<-yearsFut[1] # we need this in there in case all years in a cell are emergent - this can happen when we use a quantile to calculate focmaxT
      }
    }

    #calculate the time when the temperature emerges for a % of years given by prob
    #i.e. if prob = 0.8, when is the first year where 80% of subsequent years are emergent
    #n.b. subsequent years may drop below this threshold

    xProb<-rev(cumsum(rev(x)))/(length(yearsFut):1)
    xEm<-which(xProb>=prob)
    if(length(xEm)>0){
      ToEProb<-min(yearsFut[xEm])
    }
  }
  out<-list()
  out[[1]]<-ToEFirst
  out[[2]]<-ToEPerm
  out[[3]]<-ToEProb
  return(out)
}

#==================================================================================
#==================================================================================
#==================================================================================
#' @title find emergence times for the species in a cell
#'
#' @description
#' @details
#' See Examples.
#'
#' @param futTs The future temperatures in the cell
#' @param yearsFut The years for which we have future temperature data
#' @param focmaxT The T max of each species in the cell
#' @param permYearCut The number of years before the end of the projection for defining permanent emergence i.e if the projection runs to 2300 and permYearCut = 10, then all years past 2290 mus be in an emergent state
#' @param method Can be Mora or Henson

# @keywords
#' @export
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
#' @note

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

communityEmergenceTimes<-function(focmaxT,
                                  focSpec,
                                  futTs,
                                  yearsFut,
                                  permYearCut,
                                  method="Mora",
                                  runLength,
                                  prob){

  YearCut<-max(yearsFut)-permYearCut
  nSpec<-length(focmaxT)
  probEmyear<-permEmyear<-firstEmyear<-rep(NA,nSpec)

  for(it in 1:nSpec){
    if(method=="Mora"){ToE<-findMoraToE(focmaxT[it],futTs,yearsFut,runLength,prob)} #find timing of first and permanent emergence for each species using the Mora method
    #if(method=="Henson"){ToE<-findHensonToE(histTs,futTs,yearsHist,yearsFut)} #find timing of first and permanent emergence for each species using the Henson method

    firstEmyear[it]<-ToE[[1]]
    permEmyear[it]<-ToE[[2]]
    probEmyear[it]<-ToE[[3]]

    if(is.na(permEmyear[it])==FALSE){
      if(permEmyear[it]>=YearCut){permEmyear[it]<-NA} # only count permanent emergence if the last permYearCut years are all emergent
    }
    if(is.na(probEmyear[it])==FALSE){
      if(probEmyear[it]>=YearCut){probEmyear[it]<-NA} # only count prob emergence if it occurs prior to the last permYearCut years
    }
  }
  names(probEmyear)<-names(permEmyear)<-names(firstEmyear)<-focSpec
  Emyear<-list()
  Emyear[[1]]<-firstEmyear
  Emyear[[2]]<-permEmyear
  Emyear[[3]]<-probEmyear
  return(Emyear)
}

#==================================================================================
#==================================================================================
#==================================================================================
#' @title calculateDomainHorizonTimes
#'
#' @description this is a wrapper function for calculating emergence times across the domain for a set of different climate runs and speceis Tmax values
#' @details
#' See Examples.
#'
#' @param groupName
#' @param outFolder
#' @param TimeSeriesFolder
#' @param TempTimeSeries
#' @param rangeData
#' @param maxT
#' @param domain
#' @param printProgress
#' @param permYearCut
#' @param method
#' @param runLength
#' @param prob
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
#' @note

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

calculateDomainHorizonTimes<-function(groupName,
                                      outFolder,
                                      TimeSeriesFolder,
                                      TempTimeSeries,
                                      rangeData,
                                      maxT,
                                      domain,
                                      printProgress,
                                      permYearCut,
                                      method,
                                      runLength,
                                      prob){

  nCols<-length(TempTimeSeries)

  FirstEmergence<-PermEmergence<-ProbEmergence<-matrix(nrow=dim(domain)[1],ncol=nCols)

  for(i in 1:nCols){
    load(paste(TimeSeriesFolder,TempTimeSeries[i],sep=""))

    specT<-data.frame(Tip_Label=rownames(maxT),maxT=maxT[,i])

    output<-domainEmergence(domain,rangeData,specT,tempPostMat,printProgress,permYearCut,method,runLength,prob)

    outfile<-paste(outFolder,groupName,substr(TempTimeSeries[i], 12,nchar(TempTimeSeries[i])),sep="")

    save(output,file=outfile)
    print(i)
  }
}

#==================================================================================
#==================================================================================
#==================================================================================
#' @title find emergence times for a set of communities
#'
#' @description
#' @details
#' See Examples.
#'
#' @param domain Contains the cell ids (WorldID) we want to explore - we could just use a vector these rather than a big dataframe.......
#' @param rangeData Dataframe containing species names (Tip_Label) and cell ids (WorldID)
#' @param specT Dataframe containing species names (Tip_Label) and max T values
#' @param futTs The future temperatures in the cell
#' @param yearsFut The years for which we have future temperature data
#' @param permYearCut The number of years before the end of the projection for defining permanent emergence i.e if the projection runs to 2300 and permYearCut = 10, then all years past 2290 mus be in an emergent state
#' @param tempFutMat Matrix containing future temperature values for each cell in domain

# @keywords
#' @export
#'
# @examples
#'
#'
#' @return
#' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
#' @note

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


domainEmergence<-function(domain,
                          rangeData,
                          specT,
                          tempPostMat,
                          printProgress=TRUE,
                          permYearCut,
                          method="Mora",
                          runLength,prob){

  nCells<-dim(domain)[1]
  years<-as.numeric(colnames(tempPostMat))

  probEmergenceList<-permEmergenceList<-firstEmergenceList<-vector("list", nCells)
  probEmergenceSumList<-permEmergenceSumList<-firstEmergenceSumList<-vector("list", nCells)

  domain$probEmergence<-domain$permEmergence<-domain$firstEmergence<-domain$n<-domain$meanfocmaxT<-rep(NA,nCells)
  for(i in 1:nCells){

    focCell<-domain$WorldID[i] #pick a cell
    focSpec<-rangeData$Tip_Label[which(rangeData$WorldID==focCell)]#find species in that cell
    nSpec<-length(focSpec)
    domain$n[i]<-nSpec

    if(domain$n[i]>0){
      focmaxT<-specT$maxT[match(focSpec,specT$Tip_Label)]#find maxTs for species in the cell
      domain$meanfocmaxT[i]<-mean(focmaxT)

      futTs<-tempPostMat[i,] #extract future temperatures for cell
      yearsFut<-as.numeric(colnames(tempPostMat))
      Emyear<-communityEmergenceTimes(focmaxT,focSpec,futTs,yearsFut,permYearCut,method,runLength,prob)

      firstEmergenceList[[i]]<-Emyear[[1]]
      permEmergenceList[[i]]<-Emyear[[2]]
      probEmergenceList[[i]]<-Emyear[[3]]

      domain$firstEmergence[i]<-length(na.omit(Emyear[[1]]))
      domain$permEmergence[i]<-length(na.omit(Emyear[[2]]))
      domain$probEmergence[i]<-length(na.omit(Emyear[[3]]))
    }
    if(printProgress==TRUE){print(i)}
  }

  resultsList<-list()
  resultsList[[1]]<-domain
  resultsList[[2]]<-firstEmergenceList
  resultsList[[3]]<-permEmergenceList
  resultsList[[4]]<-probEmergenceList

  return(resultsList)
}


#' #==================================================================================
#' #==================================================================================
#' #==================================================================================
#' #' @title
#' #'
#' #' @description
#' #' @details
#' #' See Examples.
#' #'
#' #' @param
#' # @keywords
#' #' @export
#' #'
#' # @examples
#' #'
#' #'
#' #' @return a data.frame
#' #' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
#' #' @note
#'
#' # @seealso
#' # @references
#' # @aliases - a list of additional topic names that will be mapped to
#' # this documentation when the user looks them up from the command
#' # line.
#' # @family - a family name. All functions that have the same family tag will be linked in the documentation.

