#==================================================================================
#==================================================================================
#==================================================================================
#' @title Community Horizon Plot
#'
#' @description
#' @details
#' See Examples.
#'
#' @param histTs
#' @param futTs
#' @param yearsHist
#' @param yearsFut
#' @param focmaxT
#' @param focSpec
#' @param permYearCut
#' @param method
#' @param runLength
#' @param prob
#' @param First logical;
#' @param Perm logical;
#' @param Prob logical;
#' @param addJitter logical;
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


CommunityHorizonPlot<-function(histTs,
                               futTs,
                               yearsHist,
                               yearsFut,
                               focmaxT,
                               focSpec,
                               permYearCut,
                               method,
                               runLength,
                               prob,
                               First=FALSE,
                               Perm=FALSE,
                               Prob=TRUE,
                               addJitter=TRUE){

  nSpec<-length(focmaxT)

  futTs<-futTs[which(yearsFut>max(yearsHist))] # only use years from future run that exceed all years in historic run
  yearsFut<-yearsFut[which(yearsFut>max(yearsHist))]
  Ts<-c(histTs,futTs)
  years<-c(yearsHist,yearsFut)

  ToE<-communityEmergenceTimes(focmaxT,focSpec,Ts,years,permYearCut,method="Mora",runLength,prob)

  firstEmyear<-ToE[[1]]
  permEmyear<-ToE[[2]]
  probEmyear<-ToE[[3]]

  FirstEmEvents<-table(match(firstEmyear,years))
  yearsFirstEmEvents<-rep(0,length(years))
  yearsFirstEmEvents[as.numeric(names(FirstEmEvents))]<-as.numeric(FirstEmEvents)
  FirstEmergeCum<-cumsum(yearsFirstEmEvents)

  PermEmEvents<-table(match(permEmyear,years))
  yearsPermEmEvents<-rep(0,length(years))
  yearsPermEmEvents[as.numeric(names(PermEmEvents))]<-as.numeric(PermEmEvents)
  PermEmergeCum<-cumsum(yearsPermEmEvents)

  ProbEmEvents<-table(match(probEmyear,years))
  yearsProbEmEvents<-rep(0,length(years))
  yearsProbEmEvents[as.numeric(names(ProbEmEvents))]<-as.numeric(ProbEmEvents)
  ProbEmergeCum<-cumsum(yearsProbEmEvents)

  mat <- matrix(ncol=2,nrow=3)
  mat[1,]<-rep(1,2)
  mat[2,]<-rep(2,2)
  mat[3,]<-rep(2,2)
  layout(mat)

  plot(0,0,type="n",ylim=c(0,nSpec),xlim=c(1850,2300),xlab="Year",ylab="N emerged")
  polygon(c(years,rev(years)),c(rep(0,length(years)),rev(ProbEmergeCum)),col="grey80",border=FALSE)
  lines(c(1850,2300),c(nSpec/2,nSpec/2),col="black")

  ymin<-min(focmaxT,Ts)
  ymax<-max(focmaxT,Ts)
  xmin<-min(years)
  xmax<-max(years)

  plot(years,Ts,ylim=c(ymin,ymax),type="n",xlab="Year",ylab="T")
  polygon(c(xmin-100,max(yearsHist),max(yearsHist),xmin-100),c(ymin-10,ymin-10,ymax+10,ymax+10),col="grey90",border=FALSE)
  lines(years,Ts,col="grey60")

  if(First==TRUE){toplot<-firstEmyear}
  if(Perm==TRUE){toplot<-permEmyear}
  if(Prob==TRUE){toplot<-probEmyear}

  focmaxTPlot<-focmaxT
  if(addJitter==TRUE){
    toplot<-jitter(toplot,amount=0.05)
    focmaxTPlot<-jitter(focmaxT,amount=0.05)
  }

  for(it in 1:nSpec){

    if(is.na(toplot[it])==FALSE){
      lines(c(min(years),toplot[it]),c(focmaxTPlot[it],focmaxTPlot[it]),col="red")
      lines(c(toplot[it],toplot[it]),c(0,focmaxTPlot[it]),col="black")
    }else{
      lines(c(xmin,xmax),c(focmaxTPlot[it],focmaxTPlot[it]),col="blue")
    }
  }

  #for(it in 1:nSpec){
  #	if(is.na(firstEmyear[it])==FALSE){
  #		lines(c(min(years),firstEmyear[it]),c(focmaxT[it],focmaxT[it]),col="red")
  #		lines(c(firstEmyear[it],firstEmyear[it]),c(0,focmaxT[it]),col="red")
  #		if(is.na(permEmyear[it])==FALSE){
  #			if(permEmyear[it]<YearCut){
  #				lines(c(firstEmyear[it],permEmyear[it]),c(focmaxT[it],focmaxT[it]),col="black")
  #				lines(c(permEmyear[it],permEmyear[it]),c(0,focmaxT[it]),col="black")
  #			}
  #		}
  #	}else{
  #		lines(c(xmin,xmax),c(focmaxT[it],focmaxT[it]),col="blue")
  #	}
  #}
}



#==================================================================================
#==================================================================================
#==================================================================================
#' @title map patterns of emergence
#'
#' @description
#' @details
#' See Examples.
#'
#' @param gridB
#' @param toplot
#' @param colors vector of colors, e.g. from `colorRampPalette`
#' @param brks2
#' @param legendPlot
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


# library(classInt)
mapIt<-function(gridB,
                toplot,
                colors=colorRampPalette(c("steelblue4","cadetblue2","yellow","red"))(10),
                brks2,
                legendPlot=TRUE){

  is.even<-function(x) x%%2==0
  is.odd<-function(x) x%%2!=0
  digits<-2
  colors_continent<-c("grey80")
  gridB@data$displaycont<-NA
  gridB@data$displaycont<-1
  gridB@data$display1<-toplot

  brks <- brks2$brks
  gridB@data$col1 <- findInterval(gridB@data$display1, brks, all.inside = TRUE)
  my.col.fr<-findColours(brks2,colors) # ramp colors based on classInts
  legtext <- names(attr(my.col.fr,"table"))  # declare label
  legtext <- substr(legtext,2, nchar(legtext)-1) # delete silly brackets
  tst <-  as.numeric(unlist(strsplit(legtext,",")))  # from here on  this all to standardize digits in legend text
  lower <- tst[is.odd(rep(1:length(tst)))]
  upper <- tst[is.even(rep(1:length(tst)))]
  t <- data.frame(lower,upper)
  t<-round(t,digits)
  legtext <- paste(format(t$lower,nsmall = digits),format(t$upper,nsmall = digits),sep="-")
  legcols <- attr(my.col.fr,"palette")
  plot(gridB, col = colors_continent[gridB$displaycont],border = NA)
  plot(gridB, col = colors[gridB$col1],border = NA,add=TRUE)
  if(legendPlot==TRUE){legend(-16000000,2500000,legend=legtext,fill=legcols,border= "black", cex=0.5, x.intersp = 0.6, y.intersp = 0.9,bg="white")}
}
#'
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

