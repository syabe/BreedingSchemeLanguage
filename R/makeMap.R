makeMap <- function(map, nLoci, nMarkers, nQTL, propDomi, interactionMean, varEffects = 1){
  nEffectiveLoci <- 1 + rpois(n = nQTL, lambda = interactionMean)
  posEffectiveLoci <- sample(1:nLoci, sum(nEffectiveLoci))
  actionType <- rbinom(sum(nEffectiveLoci), 1, propDomi)
  effectID <- NULL
  for(i in 1:nQTL){
    effectID <- c(effectID, rep(i, nEffectiveLoci[i]))
  } # i
  if(length(varEffects) == 1){
    vEffect <- varEffects
  }else{
    vEffect <- varEffects[1]
  }
  effect <- rnorm(nQTL, 0, sqrt(varEffects))
  effects <- as.matrix(effect, nQTL, 1)
  rownames(effects) <- 1:nQTL
  colnames(effects) <- "trait1"
  if(nLoci - length(posEffectiveLoci) < nMarkers) print("Warning: Number of markers is not enough!")
  mrkPos <- sort(sample((1:nLoci)[-posEffectiveLoci], nMarkers, replace = F))
  return(list(map = map, markerPos = mrkPos, effectID = effectID, effectivePos = posEffectiveLoci, actionType = actionType, effects = effects))
}
