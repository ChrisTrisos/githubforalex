#' @export
.calc50Quant<-function(x){
  quant<-quantile(na.omit(x),0.5)
  return(quant)
}

#' @export
#'
.calc025Quant<-function(x){
  quant<-quantile(na.omit(x),0.025)
  return(quant)
}

#' @export
.calc975Quant<-function(x){
  quant<-quantile(na.omit(x),0.975)
  return(quant)
}
