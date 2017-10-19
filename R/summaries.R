#==================================================================================
#==================================================================================
#==================================================================================
#' @title extract emergence metrics for a set of runs
#'
#' @description
#' @details
#' See Examples.
#'
#' @param groupName
#' @param ET
#' @param TempTimeSeries
#' @param year
#' @param yearTr
#' @param EntropyRes
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


summariseHorizonTimes<-function(groupName,
                                ET,
                                TempTimeSeries,
                                year,
                                yearTr,
                                EntropyRes){

  nCols<-length(ET)

  firstEmerge<-matrix(ncol=nCols,nrow=dim(domain)[1])
  colnames(firstEmerge)<-ET
  probEmergeMeanYr<-permEmergeMeanYr<-firstEmergeMeanYr<-probEmergeProp<-permEmergeProp<-firstEmergeProp<-probEmerge<-permEmerge<-firstEmerge
  probEntropy<-permEntropy<-firstEntropy<-probAUC<-permAUC<-firstAUC<-probEmergeMeanYrTr<-permEmergeMeanYrTr<-firstEmergeMeanYrTr<-firstEmerge

  for(i in 1:nCols){

    load(paste(TimeSeriesFolder,TempTimeSeries[i],sep=""))

    years<-as.numeric(colnames(tempPostMat))

    load(paste(getwd(),"/HorizonTimes/",ET[i],sep=""))

    domain1<-output[[1]]
    firstEmergenceList<-output[[2]]
    permEmergenceList<-output[[3]]
    probEmergenceList<-output[[4]]

    #How many species emerge by a given year?
    firstEmerge[,i]<-NSpeciesEmerging(firstEmergenceList,year)
    permEmerge[,i]<-NSpeciesEmerging(permEmergenceList,year)
    probEmerge[,i]<-NSpeciesEmerging(probEmergenceList,year)

    firstEmergeProp[,i]<-firstEmerge[,i]/domain1$n
    permEmergeProp[,i]<-permEmerge[,i]/domain1$n
    probEmergeProp[,i]<-probEmerge[,i]/domain1$n

    #What is the mean emergence time within each cell?
    firstEmergenceListNoNA<-lapply(firstEmergenceList, function(x) x[!is.na(x)])
    firstEmergeMeanYr[,i]<-sapply(firstEmergenceListNoNA,mean)
    firstEmergenceListTr<-lapply(firstEmergenceList,function(x) replace(x,is.na(x),yearTr))
    firstEmergeMeanYrTr[,i]<-sapply(firstEmergenceListTr,mean)

    permEmergenceListNoNA<-lapply(permEmergenceList, function(x) x[!is.na(x)])
    permEmergeMeanYr[,i]<-sapply(permEmergenceListNoNA,mean)
    permEmergenceListTr<-lapply(permEmergenceList,function(x) replace(x,is.na(x),yearTr))
    permEmergeMeanYrTr[,i]<-sapply(permEmergenceListTr,mean)

    probEmergenceListNoNA<-lapply(probEmergenceList, function(x) x[!is.na(x)])
    probEmergeMeanYr[,i]<-sapply(probEmergenceListNoNA,mean)
    probEmergenceListTr<-lapply(probEmergenceList,function(x) replace(x,is.na(x),yearTr))
    probEmergeMeanYrTr[,i]<-sapply(probEmergenceListTr,mean)

    EmStats<-calcEmStats(firstEmergenceList,years,EntropyRes)
    firstAUC[,i]<-EmStats[[1]]
    firstEntropy[,i]<-EmStats[[2]]

    EmStats<-calcEmStats(permEmergenceList,years,EntropyRes)
    permAUC[,i]<-EmStats[[1]]
    permEntropy[,i]<-EmStats[[2]]

    EmStats<-calcEmStats(probEmergenceList,years,EntropyRes)
    probAUC[,i]<-EmStats[[1]]
    probEntropy[,i]<-EmStats[[2]]

    print(i)

    save(probEmerge,permEmerge,firstEmerge,probEmergeProp,permEmergeProp,firstEmergeProp,
         probEmergeMeanYr,permEmergeMeanYr,firstEmergeMeanYr,probEmergeMeanYrTr,permEmergeMeanYrTr,firstEmergeMeanYrTr,
         probAUC,permAUC,firstAUC,probEntropy,permEntropy,firstEntropy,
         file=paste(getwd(),"/HorizonResults/",groupName,"_",year,"_",yearTr,"_",EntropyRes,".rda",sep=""))
  }
}


#==================================================================================
#==================================================================================
#==================================================================================
#' @title calculate entropy and AUC fo emergence times
#'
#' @description
#' @details
#' See Examples.
#'
#' @param EmTimes
#' @param years
#' @param EntropyRes
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


calcEmStats<-function(EmTimes,
                      years,
                      EntropyRes){

  AUC<-Entropy<-rep(NA,length(EmTimes))
  for(i in 1:length(EmTimes)){

    ToE<-EmTimes[[i]]

    EmEvents<-table(match(ToE,years))
    yrEmEvents<-numeric(length(years))
    yrEmEvents[as.numeric(names(EmEvents))]<-as.numeric(EmEvents)

    EmCum<-cumsum(yrEmEvents)
    AUC[i]<-sum(EmCum)/length(ToE)

    nTP<-ceiling(length(years)/EntropyRes)
    TP<-rep(1:nTP,each=EntropyRes)
    TP<-TP[1:length(years)]

    EmEventsByPeriod<-sapply(split(yrEmEvents,TP),sum)

    Entropy[i]<-entropy(EmEventsByPeriod)
  }

  outList<-list()
  outList[[1]]<-AUC
  outList[[2]]<-Entropy
  return(outList)
}


#==================================================================================
#==================================================================================
#==================================================================================
#' @title NSpeciesEmerging
#'
#' @description using the species emergence times for each community calculate the number of species emerging by a particular year
#' @details
#' See Examples.
#'
#' @param emergenceTimes
#' @param year
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


NSpeciesEmerging<-function(emergenceTimes,
                           year){

  emergenceTimesNoNA<-lapply(emergenceTimes, function(x) x[!is.na(x)])
  emergenceBin<-lapply(emergenceTimesNoNA,function(x) replace(x,x<=year,1))
  emergenceBin<-lapply(emergenceBin,function(x) replace(x,x>year,0))

  nEmerged<-sapply(emergenceBin,sum)

  return(nEmerged)
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
#' #' @return
#' #' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
#' #' @note
#'
#' # @seealso
#' # @references
#' # @aliases - a list of additional topic names that will be mapped to
#' # this documentation when the user looks them up from the command
#' # line.
#' # @family - a family name. All functions that have the same family tag will be linked in the documentation.


