#' @title Estimativa de maxima verossimilhanca dos parametros da distribuicao 
#' Birnbaum-Saunders (BS) para dados censurados do tipo II.
#' 
#' @name estimabs
#'
#' @description A funcao calcula a estimativa de maxima verossimilhanca dos parametros
#' da distribuicao Birnbaum-Saunders.
#'
#' @param x Um vetor numerico.
#' @param delta Um vetor numerico que e 1 se a observacao e nao censurada
#' e 0 se for censurada.
#'
#' @details Se o vetor \code{delta} nao for especificado, este assume por \emph{default}
#' um vetor apenas de 1 com comprimento igual ao do vetor \code{x}.
#'
#' @return A estimativa dos parametros de forma \eqn{\alpha} e o parametro de escala
#' \eqn{\beta}.
#'
#' @author Alan da Silva
#'
#' @seealso \code{\link{maxLik}}.
#' 
#' @references BIRNBAUM, Zygmund William; SAUNDERS, Sam C. 
#' A new family of life distributions. Journal of applied probability, 
#' v. 6, n. 2, p. 319-327, 1969.
#' @references BIRNBAUM, Zygmunt W.; SAUNDERS, Sam C.
#'  Estimation for a family of life distributions with 
#'  applications to fatigue. Journal of Applied Probability, v. 6, n. 2, p. 328-347, 1969.
#'
#' @import VGAM
#' 
#' @importFrom maxLik maxLik
#' 
#' @export
estimabs=function(x,delta=rep(1,times=length(x))){
  if((any(delta!=1)&&any(delta!=0))){
    stop("O vetor delta deve conter apenas 0 e 1")
  }
  if(length(x)!=length(delta)){
    stop("Os vetores x e delta devem apresentar comprimentos iguais")
  }
  logvero <- function(parametros){
    a <- parametros[1]
    b <- parametros[2]
    fx<- dbisa(x, scale = b, shape = a)
    log.vero <- sum(log(fx)*delta+(1-delta)*log(1-pbisa(x,b,a)))
    return(log.vero)  
  }
  n<-length(x)
  s <- sum(x)/n
  rr <- 1/((1/n)*sum(1/x))
  alpha.inicial <- (2*(sqrt(s/rr)-1))^(1/2)
  beta.inicial <- sqrt(s*rr)
  # Fun??o log-verossimilhan?a
  # Maximizando a fun??o de verossimilhan?a
  maxvero <- maxLik(logvero, start = c(alpha.inicial, beta.inicial), method = "nr")
  # Par?metros estimados da BS
  alpha <- maxvero$estimate[1]
  beta <- maxvero$estimate[2]
  return(c(alpha,beta))
}
