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
    if (is.null(popID)) popID <- max(breedingData$popID)
    if(is.null(trainingPopID)){
      GID.train <- sort(unique(breedingData$phenoGID))
    }else{
      tf <- breedingData$popID %in% trainingPopID
      GID.train <- breedingData$GID[tf]
    }
    GID.not.train <- breedingData$GID[-GID.train]
    GID <- breedingData$phenoGID
    y <- breedingData$pValue
    R <- breedingData$error
    # Put a cell in y for GID that do not have a phenotype
    GIDnoPheno <- setdiff(breedingData$GID, GID)
    nNoPheno <- length(GIDnoPheno)
    GID <- c(GID, GIDnoPheno)
    y <- c(y, rep(NA, nNoPheno))
    R <- c(R, rep(NA, nNoPheno))
    y[GID %in% GID.not.train] <- NA
    data <- data.frame(GID = GID, y = y)
    K <- A.mat(score)
    reduce <- (length(GID) > length(unique(GID)))
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
