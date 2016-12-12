#'Select individuals
#'
#'@param nSelect the number of selected individuals
#'@param popID population ID to be selected (default: When random=T, the last population. When random=F, it is the last evaluated population)
#'@param random assuming random selection or selection according to their features (T: random selection, F: selection of good individuals)
#'
#'@return information of the selected individuals and the all information created before (list)
#'
#'@export
select <- function(simEnv, nSelect = 40, popID = NULL, random = F){
  parent.env(simEnv) <- environment()
  select.func <- function(data, nSelect, popID, random=FALSE){
    mapData <- data$mapData
    breedingData <- data$breedingData
    score <- data$score
    selCriterion <- data$selCriterion
    criterion <- selCriterion$criterion
    if(is.null(popID)){
      popID <- selCriterion$popID
    }
    tf <- breedingData$popID %in% popID
    GID.now <- breedingData$GID[tf]
    if (random){
      selectedGID <- sample(GID.now, nSelect)
    } else{
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
    }#END not random selection
    popID.new <- max(breedingData$popID) + 1
    for(i in selectedGID){
      breedingData$popID[breedingData$GID == i] <- popID.new
    }
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  } #END select.func
  with(simEnv, {
    if(nCore > 1){
      sfInit(parallel=T, cpus=nCore)
      sims <- sfLapply(sims, select.func, nSelect=nSelect, popID=popID, random=random)
      sfStop()
    }else{
      sims <- lapply(sims, select.func, nSelect = nSelect, popID = popID, random=random)
    }
  })
}
