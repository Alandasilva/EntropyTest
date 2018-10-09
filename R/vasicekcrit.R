#' @title Valor critico da estatistica de Vasicek(1976).
#' 
#' @name vasicekcrit
#'
#' @description A funcao calcula a distribuicao amostral da estatistica de
#' Vasicek, bem como o seu valor critico.
#'
#' @param n tamanho da amostra.
#' @param nrep numero de repeticoes usadas na simulacao de Monte Carlo para a 
#' geracao da distribuicao amostral.
#' @param nivel nivel de significancia nominal para testar a hipotese de normalidade
#' da amostra.
#'
#' @details Se  \code{nrep} nao for especificado, este assume por \emph{default} 10.000.
#' De maneira analoga, se \code{nivel} nao for especificado, entao nivel=0,05.
#' @details A estatistica de Vasicek(1976) testa a hipotese de normalidade dos dados.
#'
#' @return A distribuicao amostral da estatistica de Vasicek e o valor critico ao nivel de
#' significancia estabelecido na funcao.
#'
#' @author Alan da Silva
#' 
#' @references VASICEK, Oldrich. A test for normality based on sample entropy. 
#' Journal of the Royal Statistical Society. Series B (Methodological), p. 54-59, 1976.
#' 
#' @importFrom stats rnorm quantile
#' 
#' @export
vasicekcrit=function(n=length(x),nrep=10000,nivel=0.05){
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
  for(j in 1:nrep){
    x<-rnorm(n)
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
    Emn[j] <- (n/(2*m*s))*(prod(vetor))^(1/n)
  }
  estatistica <- Emn
  valor.critico<- round(quantile(Emn, nivel), 4)
  mylist <- list(estatistica=estatistica,valor.critico=valor.critico)
  return(mylist)
}
