#==================================================================================
#==================================================================================
#==================================================================================
#' @title splitTimeSeries
#'
#' @description Takes a matrix with the time series of temperature for each cell and generates two matrices containing the pre and post inflection temperature values
#' @details
#' See Examples.
#'
#' @param InflectionTimes
#' @param tempMat
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
# @note

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

splitTimeSeries<-function(InflectionTimes,
                          tempMat){

  tempPreMat<-tempPostMat<-tempMat
  UIT<-na.omit(unique(InflectionTimes))

  for(i in 1:length(UIT)){
    focRows<-which(InflectionTimes==UIT[i])
    InflectionCol<-which(colnames(tempMat)==UIT[i])

    tempPreMat[focRows,InflectionCol:dim(tempPreMat)[2]]<-NA
    tempPostMat[focRows,1:(InflectionCol-1)]<-NA
  }
  output<-list()
  output[[1]]<-tempPreMat
  output[[2]]<-tempPostMat
  return(output)
}


#==================================================================================
#==================================================================================
#==================================================================================
#' @title calculateMaxT
#'
#' @description calculate the maximum temperature across species geographic ranges
#' @details
#' See Examples.
#'
#' @param groupName
#' @param outFolder
#' @param TimeSeriesFolder
#' @param TempTimeSeries
#' @param rangeCells
#' @param domain
#' @param quantWithin
#' @param quantAcross
#' @param StDevWithin
#' @param StDevAcross
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
#' @note wrapper function for 'extractMaxT' function

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


calculateMaxT<-function(groupName,
                        outFolder,
                        TimeSeriesFolder,
                        TempTimeSeries,
                        rangeCells,
                        domain,
                        quantWithin=NA,
                        quantAcross=NA,
                        StDevWithin=3,
                        StDevAcross=3){

  species<-names(rangeCells)
  nSpec<-length(species)
  nCols<-length(TempTimeSeries)
  maxT<-matrix(nrow=nSpec,ncol=nCols)
  colnames(maxT)<-TempTimeSeries
  rownames(maxT)<-species

  for(i in 1:nCols){
    load(paste(TimeSeriesFolder,TempTimeSeries[i],sep=""))
    maxT[,i]<-extractMaxT(rangeCells,domain,tempPreMat,quantWithin,quantAcross,StDevWithin,StDevAcross)
    print(i)
  }

  save(maxT,file=paste(outFolder,groupName,".rda",sep=""))

}

#==================================================================================
#==================================================================================
#==================================================================================
#' @title extractMaxT
#'
#' @description
#' @details
#' See Examples.
#'
#' @param ranges A list with each element containing the cell ids where a species occurs
#' @param domain A dataframe containing the cell id's
#' @param histT  A matrix containing the temperatures for each year (columns) for each cell (row) in 'domain'
#' @param quantWithin The quantile to use when calculating the maximum temperature within a cell
#' @param quantAcross The quantile to use when calculating the maximum temperature across cells
#' @param StDevWithin The number of standard deviations used to exclude extreme temperature values within a cell
#' @param StDevAcross The number of standard deviations used to exclude extreme temperature values across the range

# @keywords
#' @export
#'
# @examples
#'
#'
#' @return a data.frame
#' @author Alex Pigot <alex.pigot1@@gmail.com>, Cory Merow
# @note

# @seealso
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


extractMaxT<-function(ranges,
                      domain,
                      histT,
                      quantWithin,
                      quantAcross,
                      StDevWithin,
                      StDevAcross){

  nSpec<-length(ranges) # number of species

  maxTemp<-rep(NA,nSpec)

  if(is.na(StDevWithin)==FALSE & is.na(quantWithin)==FALSE | is.na(StDevAcross)==FALSE & is.na(quantAcross)==FALSE){
    print("Choose either quantile or standard deviations to remove extreme values, not both")
  }else{

    for(i in 1:nSpec){ # loop through each species

      focRows<-match(ranges[[i]],domain$WorldID) # match cells in range to rows in the domain table

      if(length(focRows)>0){ # if the range is represented on the domain

        histTFoc<-histT[focRows,] # extract the rows in the temperature matrix corresponding to these cells
        if(length(focRows)>1){ # remove any rows missing temperature data - these are typically coastal cells
          focRows<-which(is.na(rowQuantiles(histTFoc,na.rm=TRUE,1))==FALSE)
          histTFoc<-histTFoc[focRows,]
        }

        if(length(focRows)==1){ # if the species range is just a single cell
          histTFoc<-na.omit(histTFoc)
          if(length(histTFoc)>0){
            if(is.na(StDevWithin)==FALSE){
              histTFocMean<-mean(histTFoc)
              if(length(histTFoc)>1){
                histTFocSD<-sd(histTFoc)
                if(histTFocSD>0){
                  histTFocSD_upper<-histTFocMean+(histTFocSD*StDevWithin)
                  histTFoc<-histTFoc[which(histTFoc<histTFocSD_upper)]
                }
              }
              maxTemp[i]<-max(histTFoc)
            }
            if(is.na(quantWithin)==FALSE){
              maxTemp[i]<-as.numeric(quantile(histTFoc,quantWithin,na.rm=TRUE))
            }
          }
        }
        if(length(focRows)>1){	# if the species range has multiple cells then we have a matrix

          #find the maximum temperature within each cell
          if(is.na(StDevWithin)==FALSE){
            histTFocMean<-rowMeans(histTFoc,na.rm=TRUE)
            histTFocSD<-rowSds(histTFoc,na.rm=TRUE)
            histTFocSD_upper<-histTFocMean+(histTFocSD*StDevWithin)
            FindOutliers<-sweep(histTFoc,1,histTFocSD_upper)
            histTFoc[FindOutliers>0]<-NA
            histTFoc<-rowMaxs(histTFoc,na.rm=TRUE)
          }
          if(is.na(quantWithin)==FALSE){
            histTFoc<-rowQuantiles(histTFoc,quantWithin,na.rm=TRUE)
          }

          #find the maximum temperature across the range
          if(is.na(StDevAcross)==FALSE){
            histTFocMean<-mean(histTFoc)
            histTFocSD<-sd(histTFoc)
            if(histTFocSD>0){
              histTFocSD_upper<-histTFocMean+(histTFocSD*StDevAcross)
              histTFoc<-histTFoc[which(histTFoc<histTFocSD_upper)]
            }
            maxTemp[i]<-max(histTFoc)
          }
          if(is.na(quantAcross)==FALSE){
            maxTemp[i]<-as.numeric(quantile(histTFoc,quantAcross))
          }
        }
      }
      #print(i)
    }
  }
  return(maxTemp)
}

