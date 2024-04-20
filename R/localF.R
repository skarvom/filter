#'Filtre comunidades locais
#'
#'Essa função filtra a média das comunidades locais produzidas em cada simulação.
#'
#'@param df um data frame contendo resultado das simulações
#'@param meanRSA se definida com TRUE, A RSA das simualçõe será usada para a média, se for FALSE, a RSA será tirar a partir da média da abundância de todas as simualações, em geral, o primeiro método pode resultar em valores diferente de 1 quando somadas a RSA, já no segundo isso não ocorre.
#'
#'@export
localF <- function(df, meanRSA = TRUE) {
  library(dplyr)
  #Filtrando o data frame apenas para local
  simLocais <- df %>% filter(source=="local.1")

  #criando matriz para abrigar as simulações realizadas
  matrizLocais <- matrix(0, nrow = max(table(simLocais$simulation)), #pega o número máximo de linhas da simulação com mais linhas
                         ncol = max(simLocais$simulation)) #pega o número máximo de simulações

  if (meanRSA == TRUE) {
    for (i in 1: max(simLocais$simulation)) {
      #Passa uma simulação para uma coluna, a abundância já irá em RSA
      matrizLocais[1:length(simLocais$abundance[simLocais$simulation == i]), i] <- simLocais$abundance[simLocais$simulation == i] / sum(simLocais$abundance[simLocais$simulation == i])
    }
    for (x in 1:ncol(matrizLocais)) {
      matrizLocais[,x] <- sort(matrizLocais[,x], decreasing = T)
    }

    #agora passar todos os zeros para NA
    matrizLocais[matrizLocais==0] <- NA

    #retirar a média e desvio da RSA de cada simulação
    RSA <- rowMeans(matrizLocais, na.rm = T)
    SD <- apply(matrizLocais, 1, sd, na.rm=T)
    final <- cbind(RSA, SD)
    final <- as.data.frame(final)
    #agora retirar as linhas que no qual não possua valor de RSA
    final <- final %>% filter(RSA != 0)

    #criando rank, sigma, migração,
    final$RANK <- 1:nrow(final)
    final$SIGMA <- rep(df$sigma[1], nrow(final))
    final$MIG <- rep(df$migration[1], nrow(final))
    print("Método de retirar RSA antes da média aplicado!")
    return(final)
  }

  if (meanRSA == FALSE) {
    for (i in 1: max(simLocais$simulation)) {
      matrizLocais[1:length(simLocais$abundance[simLocais$simulation == i]), i] <- simLocais$abundance[simLocais$simulation == i]
    }
    for (x in 1:ncol(matrizLocais)) {
      matrizLocais[,x] <- sort(matrizLocais[,x], decreasing = T)
    }

    #agora passar todos os zeros para NA
    matrizLocais[matrizLocais==0] <- NA

    #retirar a média e desvio da RSA de cada simulação
    media <- rowMeans(matrizLocais, na.rm = T)
    RSA <- media / sum(media, na.rm = T)
    desvio <- apply(matrizLocais, 1, sd, na.rm=T)
    SD <- desvio/ sum(desvio, na.rm = T)
    final <- cbind(RSA, SD)
    final <- as.data.frame(final)
    #agora retirar as linhas que no qual não possua valor de RSA
    final <- final %>% filter(RSA != 0)

    #criando rank, sigma, migração,
    final$RANK <- 1:nrow(final)
    final$SIGMA <- rep(df$sigma[1], nrow(final))
    final$MIG <- rep(df$migration[1], nrow(final))
    print("Método de retirar a média antes da RSA aplicado!")
    return(final)

  }

}


