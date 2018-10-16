#' @title Estimativa de maxima verossimilhanca dos parametros da distribuicao 
#' Exponencial (Exp) para dados censurados do tipo II.
#' 
#' @name estimaexp
#'
#' @description A funcao calcula a estimativa de maxima verossimilhanca dos parametros
#' da distribuicao Exponencial.
#'
#' @param x Um vetor numerico.
#' @param delta Um vetor numerico que e 1 se a observacao e nao censurada
#' e 0 se for censurada.
#'
#' @details Se o vetor \code{delta} nao for especificado, este assume por \emph{default}
#' um vetor apenas de 1 com comprimento igual ao do vetor \code{x}.
#'
#' @return A estimativa do parametro de taxa \eqn{\alpha}.
#'
#' @author Alan da Silva
#'
#' @seealso \code{\link{maxLik}}.
#'
#' @import VGAM
#' 
#' @importFrom maxLik maxLik
#' 
#' @export
estimaexp=function(x,delta=rep(1,times=length(x))){
  logvero <- function(parametro){
    a <- parametro[1]
    fx<- dexp(x, rate = a)
    log.vero <- sum(log(fx)*delta+(1-delta)*log(1-pexp(x,rate = a)))
    return(log.vero)  
  }
  
  # Chutes inicial: m?todo dos momentos
  alpha.inicial <- 1/mean(x)
  # Fun??o log-verossimilhan?a
  # Maximizando a fun??o de verossimilhan?a
  maxvero <- maxLik(logvero, start = c(alpha.inicial), method = "nr")
  # Par?metros estimados da BS
  alpha <- maxvero$estimate[1]
  return(c(alpha))
}
