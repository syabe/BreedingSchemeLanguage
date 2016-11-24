#'makeGamete
#'
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
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
DH <- function(genoParent, pos){
  gamete <- makeGamete(genoParent, pos)
  progeny <- rbind(gamete, gamete)
  return(progeny)
}

#'makeDHs
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
makeDHs <- function(popSize, geno, pos){
  nPar <- nrow(geno) / 2
  nRep <- popSize %/% nPar
  rem <- popSize %% nPar
  parent <- rep(1:nPar, nRep)
  parent <- c(parent, sample(1:nPar, rem))
  progenies <- t(sapply(parent, function(par) DH(geno[c(parent[par] * 2 - 1, parent[par] * 2), ], pos)))
  return(list(progenies = progenies, pedigree = cbind(parent, parent)))
}

#'makeSelfs
#'
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
makeSelfs <- function(popSize, geno, pos){
  nPar <- nrow(geno) / 2
  nRep <- popSize %/% nPar
  rem <- popSize %% nPar
  parent <- rep(1:nPar, nRep)
  parent <- c(parent, sample(1:nPar, rem))
  progenies <- makeProgenies(cbind(parent, parent), geno, pos)
  return(list(progenies = progenies, pedigree = cbind(parent, parent)))
}

#'randomMate
#'
#'@param popSize population size
#'@param geno matrix of haplotypes
#'@param pos position of markers/QTLs
#'
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
randomMateAll <- function(popSize, geno, pos){
  nInd <- nrow(geno) / 2
  parent1 <- rep(1:nInd, length.out=popSize)
  parent2 <- sapply(parent1, function(par) sample((1:nInd)[-parent1[par]], size = 1))
  parents <- cbind(parent1, parent2)
  progenies <- makeProgenies(parents, geno, pos)
  return(list(progenies = progenies, pedigree = parents))
}
