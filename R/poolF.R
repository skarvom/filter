#'Filtre simulaçoes do Pool
#'
#'Essa função retira as médias da RSA e o Wi relativo para o pool de espécies das simulações.
#'
#'@param df um data frame contendo resultado das simulações
#'@param meanRSA se definida com TRUE, A RSA das simualçõe será usada para a média, se for FALSE, a RSA será tirar a partir da média da abundância de todas as simualações, em geral, o primeiro método pode resultar em valores diferente de 1 quando somadas a RSA, já no segundo isso não ocorre.
#'
#'@export
poolF <- function(df) {
  library(dplyr)
  simPool <- df %>% filter(source=="pool")
  simPool <- select(simPool, simulation, abundance, wi)

  #agora ordenaremos na seguinte ordem: simulação (1 a x), a abundancia da simulação x será ordenada em decrescente e o wi acompanhará a posição especifica da sua abundância.
  simPool <- simPool[with(simPool, order(simulation, -abundance)),]

  #criaremos a matriz para receber as simulações
  matrizRSA <- matrix(0, nrow = max(table(simPool$simulation)),
                      ncol = max(simPool$simulation))

  #for loop para inserir as abundâncias de cada simulação do pool em uma coluna
  for (i in 1:ncol(matrizRSA)) {
    matrizRSA[1:length(simPool$abundance[simPool$simulation == i]), i] <- simPool$abundance[simPool$simulation == i]
  }

  #criaremos a mtriz para receber os valores de wi
  matrizWi <- matrix(0, nrow = max(table(simPool$simulation)),
                     ncol = max(simPool$simulation))

  #for loop para colocar os wi's de cada simulação em uma coluna
  for (z in 1:ncol(matrizWi)) {
    matrizWi[1:length(simPool$wi[simPool$simulation == z]), z] <- simPool$wi[simPool$simulation == z]

  }

  #colocando NA nos zeros
  matrizRSA[matrizRSA == 0] <- NA
  matrizWi[matrizWi == 0] <- NA
  #tirar e médias e desvios das abundâncias e do wi e o rank
  ABUNDANCIA <- rowMeans(matrizRSA, na.rm = T)
  SD_ABUNDANCIA <- apply(matrizRSA, 1, sd, na.rm=T)
  WI <- rowMeans(matrizWi, na.rm = T)
  SD_WI <- apply(matrizWi, 1, sd, na.rm=T)
  RANK <- 1:length(ABUNDANCIA)
  final <- cbind(RANK, ABUNDANCIA, SD_ABUNDANCIA, WI, SD_WI)
  final <- as.data.frame(final) #coloquei como dataframe pois ocorreu erro no próximo passo

  final <- final %>% filter(ABUNDANCIA != 0) #aqui estou retirando linhas cuja abundancia é 0
  #adicionar os valroes de sigma e migração
  final$SIGMA <- rep(df$sigma[1], length(final$ABUNDANCIA))
  final$MIG <- rep(df$migration[1], length(final$ABUNDANCIA))
  return(final)
}




