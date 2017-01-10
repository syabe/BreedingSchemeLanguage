#'calcGenotypicValue
#'
#'@param geno matrix of haplotypes
#'@param mapData map data
#'
calcGenotypicValue <- function(geno, mapData){
  nPop <- nrow(geno) / 2
  gv1pos <- function(geno1pos, actType, effect){
    geno1pos <- as.matrix(geno1pos, nrow = 2)
    if (!is.null(mapData$domModel)){ # 1 dominant over -1; degree in actType
      mnMxGeno <- apply(geno1pos, 2, range)
      sapply(1:length(actType), function(i) c(1-actType[i], actType[i]) %*% mnMxGeno[,i])
    } else{ # Standard model
      coef <- ifelse(actType == 0, (geno1pos[1, ] + geno1pos[2, ])/2, -(geno1pos[1, ] * geno1pos[2, ]))
    }
    return(effect * prod(coef))
  }
  gv1ind <- function(genoVec, mapData){
    nQTL <- max(mapData$effectID)
    geno1ind <- rbind(genoVec[1:(length(genoVec) / 2)], genoVec[(length(genoVec) / 2 + 1):length(genoVec)])
    power <- 0
    for(i in 1:nQTL){
      power <- power + gv1pos(geno = geno1ind[, mapData$effectivePos[mapData$effectID == i]], actType = mapData$actionType[mapData$effectID == i], effect = mapData$effects[i, 1])
    }
    return(power)
  }
  genoVec <- cbind(geno[1:nPop * 2 - 1, ], geno[1:nPop * 2, ])
  gv <- apply(genoVec, 1, gv1ind, mapData = mapData)
  return(gv)
}

#'calcPhenotypicValue
#'
#'@param gv genotypic values
#'@param errorVar error variance
#'@param H2 a broad heritability
#'
calcPhenotypicValue <- function(gv, errorVar = NULL, H2 = NULL){
  if((!is.null(errorVar) & !is.null(H2)) | (is.null(errorVar) & is.null(H2))){
    stop("I cannot make phenotypic value!")
  }else{
    if(!is.null(errorVar)){
      pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
    }else{
      varG <- var(gv)
      errorVar <- varG * (1 - H2) / H2
      pv <- gv + rnorm(length(gv), 0, sqrt(errorVar))
    }
  }
  return(pv)
}
