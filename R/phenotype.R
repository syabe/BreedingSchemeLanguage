#'Evaluate the phenotypic value
#'
#'@param errorVar error variance
#'@param popID population ID to be evaluated (default: the latest population)
#'
#'@return phenotypic values and the all information created before (list)
#'@export
phenotype <- function(simEnv, errorVar = 1, popID = NULL){
  parent.env(simEnv) <- environment()
  phenotype.func <- function(data, errorVar, popID){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- rep(F, length(breedingData$GID))
    for(i in popID){
      tf[breedingData$popID == i] <- T
    }
    GID.now <- breedingData$GID[tf]
    gValue.now <- breedingData$gValue[tf]
    pValue <- calcPhenotypicValue(gv = gValue.now, errorVar = errorVar)
    if(is.null(breedingData$pValue)){
      breedingData$pValue <- pValue
      breedingData$error <- rep(errorVar, length(pValue))
      breedingData$phenoGID <- sort(GID.now)
    }else{
      breedingData$pValue <- c(breedingData$pValue, pValue)
      breedingData$error <- c(breedingData$error, rep(errorVar, length(pValue)))
      breedingData$phenoGID <- c(breedingData$phenoGID, sort(GID.now))
    }
    selCriterion <- list(popID = popID, criterion = "pheno")
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, phenotype.func, errorVar = errorVar, popID = popID)
      sfStop()
    }else{
      sims <- lapply(sims, phenotype.func, errorVar = errorVar, popID = popID)
    }
  })
}
