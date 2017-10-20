#==================================================================================
#==================================================================================
#==================================================================================




#==================================================================================
#==================================================================================
#==================================================================================
#' @title findInflectionTimes
#'
#' @description loads in time series data and produces a matrix of inflection times across the domain
#'
#' @details
#' See Examples.
#'
#' @param
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


findInflectionTimes<-function(domain,
                              inputFolder,
                              TempTimeSeries,
                              method,
                              printProgress=FALSE,
                              parallelize=FALSE,
                              ...){

	# for testing
	# method="Henson"; parallelize = TRUE; sdDate=1960; nSD=3; sdThreshold=1; lowerBound=1950; upperBound=2010
  nrows<-dim(domain)[1]
  ncols<-length(TempTimeSeries)
  output<-matrix(nrow=nrows,ncol=ncols)
  colnames(output)<-TempTimeSeries
  if(!parallelize){
    for(i in 1:ncols){
      load(paste(inputFolder,TempTimeSeries[i],sep=""))
      output[,i]<-calcInflectionTimes(tempMat,method,printProgress,...)
      print(i)
    }
  } else {
    output=foreach(i = 1:ncols,.combine=cbind,.packages=c('pracma')) %dopar% {
      load(paste0(inputFolder,TempTimeSeries[i]))
      print(i)
      calcInflectionTimes(tempMat,method,printProgress,...)
    }
    #output=do.call('cbind',output)  
  }
  return(output)
}



#==================================================================================
#==================================================================================
#==================================================================================
#' @title calcInflectionTimes
#'
#' @description
#' @details
#' See Examples.
#'
#' @param
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

#CalcInflectionTimes
calcInflectionTimes<-function(tempMat,
															method,
															printProgress,
															...){

  nRows<-dim(tempMat)[1]
  InflectionYear<-rep(NA,nRows)

  for(i in 1:nRows){
    if(printProgress){print(i)}
    TsData<-data.frame(years=as.numeric(colnames(tempMat)),Ts=tempMat[i,])
    if(length(na.omit(TsData$Ts))>0){
      if(method=="Henson"){ InflectionYear[i]<-findHensonInflection(TsData,verbose=F,...)
      }
      if(method=="OPT"){InflectionYear[i]<-findOPTInflection(TsData)}
    }
  }
  return(InflectionYear)
}

#==================================================================================
#==================================================================================
#==================================================================================

#' @title Find inflections in temperature time series
#'
#' @description Inflection based on change exceding noise
#'
#' @details
#' See Examples.
#'
#' @param TsData dataframe?
#' @param option 1 or 2 for methods to calculate the cumulative gradient
#' @param sdDate end data for calculating the SD
#' @param nSD nuber of SD that cumulative gradient must excede the mean
#' @param sdThreshold distance, in SD, that cumulative gradient must excede the mean
#' @param lowerBound year to begin looking for inflection point
#' @param upperBound max year to use at inflection point
#' @param verbose logical
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

findHensonInflection<-function(TsData,
                               option=2,
                               sdDate=1980,
                               nSD=2,
                               sdThreshold=1,
                               lowerBound=1950,
                               upperBound=2010,
                               verbose=T){

  if(option==1){csTs<-as.numeric(cumsum(diff(TsData$Ts)))} #calculate the cummulative gradient in temperature using just the difference in temperature between years
  if(option==2){csTs<-as.numeric(cumsum(gradient(TsData$Ts,TsData$years)))} #calculate the cummulative gradient in temperature following Henson. Not exactly sure what this 'gradient' function is doing....
  # CM: its doing the diff as you do above except using a time step of 2 rather than 1.
  #if(option==3){csTs1<-as.numeric(cumsum(diff(TsData$Ts,lag=2)))}

  earlyT=sort(subset(TsData$Ts,TsData$years<sdDate))
  earlyTSD=sd(earlyT)
  earlyTMean=mean(earlyT)
  keep=abs(earlyT-earlyTMean)<(nSD*earlyTSD)
  if(verbose){
    if(any(!keep)) print(paste0(sum(!keep),' outliers tossed'))
  }
  earlyT=earlyT[keep]
  sd=sd(subset(TsData$Ts,TsData$years<sdDate))
  negVals<-which(csTs<(sd*sdThreshold)) # find which years had a negative cummulative gradient

  if(length(negVals)>0){
    inflection<-TsData$years[max(negVals)+1] # the last year where the cummulative gradient was -ve marks the inflection point
  } else {
    inflection<-TsData$years[1] # if no years were negative then the trend must have been continuously increasing (according to Henson anyway....)
  }
  if(is.na(inflection)) inflection=upperBound
  if(inflection<lowerBound) inflection=lowerBound
  if(inflection>upperBound) inflection=upperBound

  return(inflection)
}

#==================================================================================
#==================================================================================
#==================================================================================
# https://stats.stackexchange.com/questions/149627/piecewise-regression-with-constraints
#' @title Find inflection time
#'
#' @description
#'
#' @details
#' See Examples.
#'
#' @param
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


findOPTInflection<-function(TsData){

  y<-TsData$Ts
  x<-TsData$years
  #minT<-min(y)
  #maxT<-max(y)

  minT<-mean(y[1:50])
  maxT<-max(y[(length(y)-50):length(y)])


  fun <- function(par, x) {
    #set all y values to starting intercept
    y1 <- x^0 * par["i1"]
    #set values after second breakpoint to ending intercept
    y1[x >= par["x2"]] <- par["i2"]
    #which values are between breakpoints?
    r <- x > par["x1"] & x < par["x2"]
    #interpolate between breakpoints
    y1[r] <- par["i1"] + (par["i2"] - par["i1"]) / (par["x2"] - par["x1"]) * (x[r] - par["x1"])
    y1
  }

  #sum of squared residuals
  SSR <- function(par) {
    sum((y - fun(par, x))^2)
  }

  inflection<-optimx(par = c(x1 = 2000, x2 = 2100, i1 = minT, i2 = maxT), fn = SSR, method = "Nelder-Mead")[[1]]
  return(inflection)
}


#' @title findHensonToE
#'
#' @description function to find the emergence time and also inflection time (i.e. when the temperature time series starts to trend upwards) following the method of Henson (or as close as I can follow it anyway)
#'
#' @details
#' See Examples.
#'
#' @param
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


findHensonToE<-function(TsData,option){

  if(option==1){csTs<-as.numeric(cumsum(diff(TsData$Ts)))} #calculate the cummulative gradient in temperature using just the difference in temperature between years
  if(option==2){csTs<-as.numeric(cumsum(pracma::gradient(TsData$Ts,TsData$years)))} #calculate the cummulative gradient in temperature following Henson. Not exactly sure what this 'gradient' function is doing....
  #details of the gradient function can be found here
  # https://www.rdocumentation.org/packages/pracma/versions/1.9.9/topics/gradient

  negVals<-which(csTs<0) # find which years had a negative cummulative gradient

  if(length(negVals)>0){
    inflection<-TsData$years[max(negVals)+1] # the last year where the cummulative gradient was -ve marks the inflection point
  }else{
    inflection<-TsData$years[1] # if no years were negative then the trend must have been continuously increasing (according to Henson anyway....)
  }

  #use the years post the inflection to calculate the slope in temperature change
  toFit<-which(TsData$years>=inflection)
  yearsToFit<-TsData$years[toFit]
  TsToFit<-TsData$Ts[toFit]

  mod.gls <- gls(TsToFit ~ yearsToFit, correlation=corARMA(p=1), method="ML")

  intercept<-as.numeric(summary(mod.gls)$coef[1])
  slope<-as.numeric(summary(mod.gls)$coef[2])

  #use the years pre inflection to calculate the background noise in interannual temperatures
  preInflectionTemp<-TsData$Ts[which(TsData$years<inflection)]

  noise<-sd(preInflectionTemp)

  #use the estimates of the slope and noise to calculate the time of emergence
  ToE<-inflection+(2*noise/slope)

  out<-list()
  out[[1]]<-ToE
  out[[2]]<-NA
  out[[3]]<-inflection
  out[[4]]<-intercept
  out[[5]]<-slope
  out[[6]]<-noise
  out[[7]]<-preInflectionTemp

  names(out)<-c("ToE","","InflectionYear","Intercept","Slope","Noise","PreInflectionTemp")

  return(out)
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


