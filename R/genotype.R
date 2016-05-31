#'Genotype markers
#'
#'@return marker genotype and the all information created before (list)
#'
#'@export
genotype <- function(){
  genotype.func <- function(data){
    mapData <- data$mapData
    breedingData <- data$breedingData
    selCriterion <- data$selCriterion
    geno <- breedingData$geno
    nPop <- nrow(geno) / 2
    score <- geno[1:nPop * 2 - 1, ] + geno[1:nPop * 2,]
    score <- score / 2
    score <- score[, mapData$markerPos]
    rownames(score) <- breedingData$GID
    return(list(mapData = mapData, breedingData = breedingData, score = score, selCriterion = selCriterion))
  }
  if(nCore > 1){
    sfInit(parallel=T, cpus=nCore)
    lists <<- sfLapply(lists, genotype.func)
    sfStop()
  }else{
    lists <<- lapply(lists, genotype.func)
  }
}
