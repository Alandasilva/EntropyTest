#' @title Estimativa da entropia populacional pelo estimador de Van Es(1992)
#' @name HVEmn
#'
#' @description A funcao calcula a estimativa da entropia populacional
#'     com base na amostra.
#'
#' @param x Um vetor numerico.
#'
#' @details O vetor \code{x} deve ter comprimento maior ou igual que 3.
#'
#' @return A estimativa da entropia com base no vetor \code{x}. 
#' \deqn{HVE_{mn}=1/(n-m) \sum _{i=1}^{n-m}log((n+1/m)[x_{i+m} - x_{i-m}]) + \sum _{k=m}^{n}(1/k) + log(m) - log(n+1)}
#'
#' @author Alan da Silva
#'
#' @references VAN ES, Bert. Estimating functionals related to a density by a class of statistics based on spacings. 
#' Scandinavian Journal of Statistics, p. 61-72, 1992.
#'
#' @export
  HVEmn <- function(x){
  if(length(x)<3||is.character(x)){
    stop("x deve ser um vetor numerico de comprimento maior ou igual que 3")
  }
  x <- sort(x)
  n <- length(x)
  m <- round(sqrt(n)+0.5)
  vetor <- NULL
  for(i in 1:(n-m)){
    vetor[i] <- x[i+m]-x[i]
  }
  vetor <- ((n+1)/m)*vetor
  vetor1 <- 1/(m:n)
  estimate <- (1/(n-m))*sum(log(vetor)) + sum(vetor1) + log(m) - log(n+1)
  return(estimate)
}

