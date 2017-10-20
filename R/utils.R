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


#' @title Nice color palettes for ranges
#'
#' @description Nice color palettes for ranges
#'
#' @details
#' See Examples.
#'
#' @param x number of colors
#' @param bias as in colorRampPalette
# @keywords
#' @export
#'
# @examples
#'
#'
#' @return Returns a vector of colors
#' @author Cory Merow <cory.merow@@gmail.com>
# @note
#' @seealso colorRampPalette
# @references
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


cm.cols1=function(x,bias=1) {
  colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)
}
