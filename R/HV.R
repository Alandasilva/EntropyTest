#' @title Estimativa da entropia populacional pelo estimador de Vasicek(1976)
#' @name HVmn
#'
#' @description A funcao calcula a estimativa da entropia populacional
#'     com base na amostra.
#'
#' @param x Um vetor numerico.
#'
#' @details O vetor \code{x} deve ter comprimento maior ou igual que 3.
#'
#' @return A estimativa da entropia com base no vetor \code{x}. \deqn{HV_{mn}=\sum _{i=1}^{n}log (n/2m)[x_{i+m} - x_{i-m}]}
#'
#' @author Alan da Silva
#'
#' @references VASICEK, Oldrich. A test for normality based on sample entropy. Journal of the Royal
#' Statistical Society. Series B (Methodological), p. 54-59, 1976
#'
#' @export
HVmn <- function(x){
  if(length(x)<3||is.character(x)){
    stop("x deve ser um vetor numerico de comprimento maior ou igual que 3")
  }
  x <- sort(x)
  n <- length(x)
  m <- round(sqrt(n)+0.5)
  vetor <- NULL
  for(i in 1:n){
    if((i-m) < 1){vetor[i] <- x[i+m]-min(x)}
    if((i+m) > n){vetor[i] <- max(x)-x[i-m]}
    if((i-m) >= 1 && (i+m) <= n){vetor[i] <- x[i+m]-x[i-m]}
  }
  estimate <- sum(log((n/(2*m))*vetor))/n
  return(estimate)
}

