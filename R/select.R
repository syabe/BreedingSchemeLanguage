#'Select individuals
#'
#'@param nSelect the number of selected individuals
#'@param popID population ID to be selected (default: When random=T, the last population. When random=F, it is the last evaluated population)
#'@param random assuming random selection or selection according to their features (T: random selection, F: selection of good individuals)
#'
#'@return information of the selected individuals and the all information created before (list)
#'
#'@export
select <- function(nSelect = 40, popID = NULL, random = F){
  select.func.random <- function(data, nSelect, popID){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    selCriterion <- data$selCriterion
    if(is.null(popID)){
      popID <- max(breedingData$popID)
    }
    tf <- rep(F, length(breedingData$GID))
    for(i in popID){
      tf[breedingData$popID == i] <- T
    }
    GID.now <- breedingData$GID[tf]
    selectedGID <- sample(GID.now, nSelect)
    popID.new <- max(breedingData$popID) + 1
    for(i in selectedGID){
      breedingData$popID[breedingData$GID == i] <- popID.new
    }
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  select.func <- function(data, nSelect, popID){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    selCriterion <- data$selCriterion
    criterion <- selCriterion$criterion
    if(is.null(popID)){
      popID <- selCriterion$popID
    }
    tf <- rep(F, length(breedingData$GID))
    for(i in popID){
      tf[breedingData$popID == i] <- T
    } # i
    GID.now <- breedingData$GID[tf]
    candValue <- NULL
    if(substr(criterion, 1, 5) == "pheno"){
      for(i in GID.now){
        error.now <- min(breedingData$error[breedingData$phenoGID == i])
        candValue <- c(candValue, mean(breedingData$pValue[(breedingData$phenoGID == i) & (breedingData$error == error.now)]))
      }
    }else{
      if(substr(criterion, 1, 4) == "pred"){
        for(i in GID.now){
          candValue <- c(candValue, breedingData$predict[(breedingData$predGID == i) & (breedingData$predNo == max(breedingData$predNo))])
        }
      }else{
        stop("Please define selection criterion in correct way!")
      }
    }
    order <- order(candValue, decreasing = T)
    selectedGID <- GID.now[order[1:nSelect]]
    popID.new <- max(breedingData$popID) + 1
    for(i in selectedGID){
      breedingData$popID[breedingData$GID == i] <- popID.new
    }
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  if(nCore > 1){
    if(random){
      sfInit(parallel=T, cpus=nCore)
      lists <<- sfLapply(lists, select.func.random, nSelect = nSelect, popID = popID)
      sfStop()
    }else{
      sfInit(parallel=T, cpus=nCore)
      lists <<- sfLapply(lists, select.func, nSelect = nSelect, popID = popID)
      sfStop()
    }
  }else{
    if(random){
      lists <<- lapply(lists, select.func.random, nSelect = nSelect, popID = popID)
    }else{
      lists <<- lapply(lists, select.func, nSelect = nSelect, popID = popID)
    }
  }
}
