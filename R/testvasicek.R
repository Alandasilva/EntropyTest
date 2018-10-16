#' @title Teste de normalidade de Vasicek(1976).
#' 
#' @name testvasicek
#'
#' @description A funcao calcula a estatistica de teste, o valor critico e o p-valor
#' da estatistica de teste de Vasicek.
#'
#' @param x e o vetor contendo a amostra.
#' @param nrep numero de repeticoes usadas na simulacao de Monte Carlo para a 
#' geracao da distribuicao amostral.
#' @param nivel nivel de significancia nominal para testar a hipotese de normalidade.
#'
#' @details Se  \code{nrep} nao for especificado, este assume por \emph{default} 10.000.
#' De maneira analoga, se \code{nivel} nao for especificado, entao nivel=0,05.
#' @details A estatistica de Vasicek(1976) testa a hipotese de normalidade dos dados.
#'
#' @return A  estatistica de Vasicek, o valor critico ao nivel de
#' significancia estabelecido na funcao, bem como o p-valor obtido.
#'
#' @author Alan da Silva
#' 
#' @references VASICEK, Oldrich. A test for normality based on sample entropy. 
#' Journal of the Royal Statistical Society. Series B (Methodological), p. 54-59, 1976.
#' 
#' @importFrom stats rnorm quantile
#' 
#' @export
testvasicek <- function(x,nivel=0.05,nrep=10000){
  n <- length(x)
  if((n<3)||is.character(n)){
    stop("n deve ser um inteiro maior ou igual a 3")
  }
  Emn<-NULL
  if (n<=8){m <- 1}
  if (9<=n&&n<=15){m <- 2}
  if (16<=n&& n<=35){m <- 3}
  if (36<=n&& n<=60){m <- 4}
  if (61<=n&& n<=80){m <- 5}
  if (81<=n&& n<=100){m <- 6}
  if (n>100){m <- round(sqrt(n)+0.5)}
  x<-sort(x)
  xbarra <- mean(x)
  s2 <- (1/(n))*sum((x-xbarra)^2)
  s <- sqrt(s2)
  vetor <- NULL
  for(i in 1:n){
    if((i-m) < 1){vetor[i] <- x[i+m]-min(x)}
    if((i+m) > n){vetor[i] <- max(x)-x[i-m]}
    if((i-m) >= 1 && (i+m) <= n){vetor[i] <- x[i+m]-x[i-m]}
  }
  Est.vasicek <- (n/(2*m*s))*(prod(vetor))^(1/n)
  Est.simulada <- vasicekcrit(n=n,nivel=nivel,nrep=nrep)
  p.valor <- mean(Est.simulada$estatistica<Est.vasicek)
  est.data <- list(Estatistica_Vasicek=Est.vasicek,p_valor=p.valor,valor.critico=Est.simulada$valor.critico)
  cat("Teste de normalidade de Vasicek", "\n")
  return(est.data)
}
