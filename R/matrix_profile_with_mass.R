#' MASSpre
#'
#' Realiza a pré-computação do MASS, calculando a média móvel e o desvio padrão móvel (movstd,movmean)
#' 
#' @importFrom  zoo rollmean
#' @importFrom  zoo rollsumr
#' @param data.values se refere a um vetor numérico representando a série temporal
#' @param data.length se refere ao tamanho da série temporal passada
#' @param subsequence.length se refere ao tamanho da subsequência desejada, a mesma utilizada para 
#' buscar por discords ou motifs
#' @return Retorna uma lista contendo os dados no domínio na frequência, a média e o desvio padrão
#' pré-calculados das subsequências
#' @export
#' @examples 
#' tsLen <- 100
#' ts <- runif(n = 1000, min = 0, max = tsLen)
#' MASSPre(ts,1000,10)
MASSPre <- function(data.values,data.length, subsequence.length){
  data.freq <- stats::fft(c(data.values, numeric(length = subsequence.length)))
  data.roll_mean <- zoo::rollmean(data.values, subsequence.length)
  data.roll_sum <- zoo::rollsumr(data.values^2, subsequence.length)
  data.roll_sig <- sqrt((data.roll_sum / subsequence.length) - (data.roll_mean^ 2))
  return(list(freq = data.freq, rollmean = data.roll_mean, rollsig = data.roll_sig))
}


#' MASS
#'
#' Realiza o cálculo da distance profile para uma determina subsequência
#' 
#' @param query se a subsequência que se deseja construir a distance profile
#' @param data.length se refere ao tamanho da série temporal 
#' @param subs.length se refere ao tamanho da subsequência desejada, a mesma utilizada para 
#' buscar por discords ou motifs
#' @param ind refere-se a um inteiro que corresponde a trecho que essa subsequência foi
#' retirada
#' @param masspre refere-se ao resultado obtido atrabés da computação do \code{MASSPre}
#' @return Retorna um vetor numérico contendo a distance profile para aquela query
#' @examples 
#' tsLen <- 100
#' ts <- runif(n = 1000, min = 0, max = tsLen)
#' subLen <- 10
#' masspre <- MASSPre(ts,tsLen,subLen)
#' index <- 1
#' query <- ts[index:(index + subLen - 1)]
#' MASS(query,tsLen,subLen,index,masspre)
#' @export
MASS <- function(query, data.length, subs.length,ind,masspre){
  ifft <- function(x) stats::fft(x, inverse=TRUE) / length(x)
  query <- rev(query) #Reverse the query
  query <- c(query,numeric(data.length)) #Append zeros
  query.freq <- stats::fft(query) #Change to frequency domain
  product <- ifft(masspre$freq * query.freq)
  distance.profile <- 2 * (subs.length - (product[subs.length:data.length] -
                                            subs.length * masspre$rollmean * masspre$rollmean[ind])/
                             (masspre$rollsig * masspre$rollsig[ind]))
  return(distance.profile)
}


#' ComputeMatrixProfile
#'
#' Realiza o cálculo da matrix profile para uma determinada série temporal
#'
#' @param data Refere-se a um vetor numérico representando a série temporal
#' @param n se refere ao tamanho da subsequência desejada, a mesma utilizada para 
#' buscar por discords ou motifs
#' @return Retorna um data.frame contendo os valores calculados para a matrix profile
#' e o índice correspondente a cada distância
#' @examples 
#' ts <- runif(n = 1000, min = 0, max = 100)
#' subLen <- 10
#' ComputeMatrixProfile(ts,subLen)
#' @export
ComputeMatrixProfile <- function(data, n){

  IsInfIsNaN <- function(x) return(any(is.na(x) | is.infinite(x)))
  
  # Initialization
  exclusion.zone <- round(n * 0.5)
  serie.length <- length(data)
  mp.length <- serie.length - n + 1
  matrixProfile <- data.frame(index = numeric(length = mp.length),values = numeric(length = mp.length) + Inf)
  
  # Search for Inf e Nans
  isSkip <- logical(mp.length)
  for(i in 1:mp.length){
    isSkip[i] <- IsInfIsNaN(data[i:(i+n-1)])
  }
  data[is.na(data) | is.infinite(data)] <- 0
  
  massPre <- MASSPre(data,serie.length, n)
  idxOrder <- sample(mp.length) 
  
  for(idx in idxOrder){
    
    if(!isSkip[idx]){
      #Compute the distance profile
      query <- data[idx:(idx + n - 1)]
      distance.profile <- MASS(query, serie.length, n, idx, massPre)
      distance.profile <- sqrt(abs(distance.profile))
      distance.profile[isSkip] <- Inf
      
      #Apply exclusion zone
      excZoneStart <- max(1, idx - exclusion.zone)
      excZoneEnd <- min(mp.length, idx + exclusion.zone)
      distance.profile[excZoneStart:excZoneEnd] <- Inf
      
      #Update matrix profile
      updatePos <- distance.profile < matrixProfile$values
      matrixProfile[updatePos, ] <- list(idx, distance.profile[updatePos])
      idx.min <- which.min(distance.profile)
      matrixProfile[idx, ] <- list(idx.min, distance.profile[idx.min])
    }
  }
  matrixProfile$values[is.infinite(matrixProfile$values)] <- NA
  return(matrixProfile)
}


#' FindMotifs
#'
#' Extrai os motifs a partir da matrix profile
#'
#' @param data Refere-se a um vetor numérico representando a série temporal
#' @param subs.length se refere ao tamanho do motif que se deseja buscar (mesmo tamanho de
#' subsequência usado em \code{ComputeMatrixProfile})
#' @param num_motifs quantidade de motifs que se deseja encontrar
#' @param matrix.profile refere-se a matrix profile previamente calculada pelo método
#' \code{ComputeMatrixProfile}
#' @param num_neighbors refere-se a quantidade máxima de vizinhos a se
#' considerar por motifs
#' @param radius refere-se ao raio para considerar se uma subsequência caso com 
#' determinado motif
#' @return Retorna um lista contendo a matrix profile com as zonas de exclusão e 
#' a lista contendo os índices dos motifs
#' @examples 
#' ts <- runif(n = 1000, min = 0, max = 100)
#' subLen <- 10
#' mp <- ComputeMatrixProfile(ts,subLen)
#' motifs <- FindMotifs(ts,subLen,mp)
#' motifs$motifs
#' @export
FindMotifs <- function(data, subs.length, matrix.profile, num_motifs = 3, num_neighbors = 10, radius = 2){
  #Aux Function
  IsInfIsNaN <- function(x) return(any(is.na(x) | is.infinite(x)))
  
  data.length <- length(data)
  mp.length <- nrow(matrix.profile)
  exclusion.zone <- round(subs.length * 0.5)
  
  SetExclusionZone <- function(ind) list(start = max( 1, ind - exclusion.zone), 
                                         end   = min( mp.length, ind + exclusion.zone))
  
  # Search for Inf e Nans
  isSkip <- sapply(1:mp.length, function(i) data[i:( i + subs.length - 1 )])
  
  massPre <- MASSPre(data, data.length, subs.length)
  
  motifs <- list()
  
  mp.current <- matrix.profile$values
  #Encontrar 3 motifs mais significantes
  for( i in 1:num_motifs ){
    
    motif.ind <- which.min(mp.current)
    motif.pair <- sort(x = c(motif.ind, matrix.profile$index[motif.ind]))
    
    #Obter K-Motif e obter a distance profile para ele
    motif.ind <- motif.pair[1]
    query <- data[motif.ind:(motif.ind + subs.length - 1)]
    distance.profile <- MASS(query, data.length, subs.length, motif.ind, massPre)
    distance.profile <- abs(distance.profile)
    
    #Aplicar zona de exclusao
    motif.dist <- mp.current[motif.ind]^2
    distance.profile[distance.profile > motif.dist * radius] <- Inf
    
    
    motif.zone <- SetExclusionZone(motif.pair[1])
    distance.profile[motif.zone$start:motif.zone$end] <- Inf
    
    motif.zone <- SetExclusionZone(motif.pair[2])
    distance.profile[motif.zone$start:motif.zone$end] <- Inf
    
    #Ignorar trechos da serie com Inf e NaNs
    distance.profile[isSkip] <- Inf
    
    #Encontrar índices válidos de vizinhos para o motif atual
    distance <- sort(distance.profile, index.return=TRUE)
    distance.idx <- which(is.finite(distance$x))
    
    
    distance.idx <- distance$ix[distance.idx]
    distance.len <- length(distance.idx)
    
    motif.neighbor <- integer()
    motif.neighbor.ind <- 1
    
    #Buscar motifs
    for( j in 1:distance.len ){
      if(distance.len == 0) break
      add_to_list = !any(abs(distance.idx[j] - motif.neighbor) < exclusion.zone)
      if( add_to_list ){
        print(distance.idx[j])
        motif.neighbor[motif.neighbor.ind] <- distance.idx[j]
        motif.neighbor.ind <- motif.neighbor.ind + 1
      }
      if ( motif.neighbor.ind  > num_neighbors ){
        break
      }
    }
    
    #Aplicar zona de exclusao sobre os vizinhos dos motifs
    motifs[[i]] <- c(motif.pair,motif.neighbor)
    for(idx in motifs[[i]]){
      motif.zone <- SetExclusionZone(idx)
      mp.current[motif.zone$start:motif.zone$end] <- Inf
    }
    
  }
  
  return(list(motifs = motifs, matrix.profile = data.frame(index = matrix.profile$index,values = mp.current)))
  
}


#' FindDiscords
#'
#' Extrai os discords a partir da matrix profile
#'
#' @param mp.current Refere-se a matrix profile atual, a qual se deseja buscar pelos discords
#' @param subs.length se refere ao tamanho do discord que se deseja buscar (mesmo tamanho de
#' subsequência usado em \code{ComputeMatrixProfile})
#' @param num_discords quantidade de discords que se deseja encontrar
#' @return Retorna vetor indicando os índices das subsequências onde os discords foram encontrados
#' @examples
#' ts <- runif(n = 1000, min = 0, max = 100)
#' subLen <- 10
#' mp <- ComputeMatrixProfile(ts,subLen)
#' motifs <- FindMotifs(ts,subLen,mp)
#' mp.current <- motifs$matrix.profile
#' FindDiscords(mp.current,subLen)
#' @export
FindDiscords <- function(mp.current,subs.length, num_discords = 3){
  mp.current <- mp.current$values
  exclusion.zone <- round(subs.length*0.5)
  
  #Ordenar os valores dos maiores para os menores para obter os discords
  distance <- sort(mp.current, decreasing = TRUE,index.return = TRUE)
  distance.idx <- which(is.finite(distance$x))
  distance.idx <- distance$ix[distance.idx]
  distance.len <- length(distance.idx)
  
  #Buscar por discords  
  discords <- integer()
  discords.ind <- 1
  
  for( i in 1:distance.len ){
    add_to_list = !any(abs(distance.idx[i] - discords) < exclusion.zone)
    if( add_to_list ){
      discords[discords.ind] <- distance.idx[i]
      discords.ind <- discords.ind + 1
    }
    if ( discords.ind  > num_discords ){
      break
    }
  }
  
  return(discords)
}