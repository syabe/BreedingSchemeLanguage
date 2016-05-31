#'Genomic prediction
#'
#'@param popID population ID to be predicted (default: the latest population)
#'@param trainingPopID population ID to be used for training a prediction model (default: all populations with phenotype data)
#'
#'@return predicted values and the all information created before (list)
#'
#'@export
predictBreedVal <- function(popID = NULL, trainingPopID = NULL){
  predict.func <- function(data, popID, trainingPopID){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    if(is.null(popID)){
      predictedPopID <- max(breedingData$popID)
    }else{
      predictedPopID <- popID
    }
    if(is.null(trainingPopID)){
      GID.train <- sort(unique(breedingData$phenoGID))
    }else{
      tf <- rep(F, length(breedingData$GID))
      for(i in trainingPopID){
        tf[breedingData$popID == i] <- T
      } # i
      GID.train <- breedingData$GID[tf]
    }
    GID.not.train <- breedingData$GID[-GID.train]
    GID <- breedingData$phenoGID
    y <- breedingData$pValue
    R <- breedingData$error
    if(max(breedingData$GID) > max(GID)){
      no.pheno <- max(breedingData$GID) - max(GID)
      GID <- c(GID, max(GID) + 1:no.pheno)
      y <- c(y, rep(NA, no.pheno))
      R <- c(R, rep(NA, no.pheno))
    }
    for(i in GID.not.train){
      y[GID == i] <- NA
    }
    data <- data.frame(GID = GID, y = y)
    K <- A.mat(score)
    if(length(GID) > length(unique(GID))){
      reduce <- T
    }else{
      reduce <- F
    }
    model <- kin.blup(data = data, geno = "GID", pheno = "y", K = K, R = R, reduce = reduce)
    predict <- as.numeric(model$g)
    predGID <- rownames(K)
    if(is.null(breedingData$predict)){
      breedingData$predict <- predict
      breedingData$predGID <- predGID
      breedingData$predNo <- rep(1, length(predGID))
    }else{
      breedingData$predict <- c(breedingData$predict, predict)
      breedingData$predGID <- c(breedingData$predGID, predGID)
      breedingData$predNo <- c(breedingData$predNo, rep(max(breedingData$predNo) + 1, length(predGID)))
    }
    selCriterion <- list(popID = predictedPopID, criterion = "pred")
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  if(nCore > 1){
    sfInit(parallel=T, cpus=nCore)
    lists <<- sfLapply(lists, predict.func, popID = popID, trainingPopID = trainingPopID)
    sfStop()
  }else{
    lists <<- lapply(lists, predict.func, popID = popID, trainingPopID = trainingPopID)
  }
}
