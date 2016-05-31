#'makeGamete
#'
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
makeGamete <- function(geno, pos){
  diff <- diff(pos)
  rec <- (1 - exp(-2 * diff / 100)) / 2
  rec[rec < 0] <- 0.5
  rec <- c(0.5, rec)
  sample <- runif(length(rec))
  crossOver <- ((rec - sample) >= 0)
  selectHaplo <- cumsum(crossOver) %% 2
  selectHaplo <- selectHaplo + 1
  gamete <- geno[1, ]
  gamete[selectHaplo == 2] <- geno[2, selectHaplo == 2]
  return(gamete)
}

#'makeProgeny
#'
#'@param genoPat matrix of paternal haplotype
#'@param genoMat matrix of maternal haplotype
#'@param pos position of markers/QTLs
#'
#'@export
makeProgeny <- function(genoMat, genoPat, pos){
  progeny <- rbind(makeGamete(genoMat, pos), makeGamete(genoPat, pos))
  return(progeny)
}

#'makeProgenies
#'
#'@param parents ID of haplotypes that the parents harbor
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
makeProgenies <- function(parents, geno, pos){
  progenies <- matrix(NA, 2 * nrow(parents), ncol(geno))
  for(par in 1:nrow(parents)){
    genoMat <- geno[c(parents[par, 1] * 2 - 1, parents[par, 1] * 2), ]
    genoPat <- geno[c(parents[par, 2] * 2 - 1, parents[par, 2] * 2), ]
    progenies[c(par * 2 - 1, par * 2), ] <- makeProgeny(genoMat, genoPat, pos)
  }
  return(progenies)
}

#'DH
#'
#'@param genoParent matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
DH <- function(genoParent, pos){
  gamete <- makeGamete(genoParent, pos)
  progeny <- rbind(gamete, gamete)
  return(progeny)
}

#'DHs
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
DHs <- function(popSize, geno, pos){
  parentsSize <- nrow(geno) / 2
  if(popSize %% parentsSize == 0){
    rep <- popSize %/% parentsSize
    parent <- rep(1:parentsSize, rep)
    parents <- cbind(parent, parent)
    progenies <- matrix(NA, 2 * popSize, ncol(geno))
    for(par in 1:length(parent)){
      progenies[c(par * 2 - 1, par * 2), ] <- DH(geno[c(parent[par] * 2 - 1, parent[par] * 2), ], pos)
    }
  }else{
    if(popSize > parentsSize){
      rep <- popSize %/% parentsSize
      rem <- popSize %% parentsSize
      parent <- rep(1:parentsSize, rep)
      parent <- c(parent, sample(1:parentsSize, rem))
      parents <- cbind(parent, parent)
      progenies <- matrix(NA, 2 * popSize, ncol(geno))
      for(par in 1:length(parent)){
        progenies[c(par * 2 - 1, par * 2), ] <- DH(geno[c(parent[par] * 2 - 1, parent[par] * 2), ], pos)
      }
    }else{
      parent <- sample(1:parentsSize, popSize)
      parents <- cbind(parent, parent)
      progenies <- matrix(NA, 2 * popSize, ncol(geno))
      for(par in 1:length(parent)){
        progenies[c(par * 2 - 1, par * 2), ] <- DH(geno[c(parent[par] * 2 - 1, parent[par] * 2), ], pos)
      }
    }
  }
  return(list(progenies = progenies, pedigree = parents))
}

#'randomMate
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
randomMate <- function(popSize, geno, pos){
  parents <- t(sapply(rep(nrow(geno) / 2, popSize), sample, size=2))
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

#'randomMateAll
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
randomMateAll <- function(popSize, geno, pos){
  parent1 <- rep(1:(nrow(geno) / 2), popSize %/% (nrow(geno) / 2))
  if(popSize %% (nrow(geno) / 2) != 0) parent1 <- c(parent1, 1:(popSize %% (nrow(geno) / 2)))
  parent2 <- rep(NA, popSize)
  for(pollen in 1:popSize){
    parent2[pollen] <- sample((1:(nrow(geno) / 2))[-parent1[pollen]], size = 1)
  }
  parents <- cbind(parent1, parent2)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}

#'selfing
#'
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
#'@export
selfing <- function(geno, pos){
  parent1 <- parent2 <- 1:(nrow(geno) / 2)
  parents <- cbind(parent1, parent2)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}
