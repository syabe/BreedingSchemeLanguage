#' Save the results
#'
#' @param summarize if T a result averaged over all the replications is saved, if F each replication's result is saved
#' @param directory the directory to which the output will be saved (Enclose the phrase in double quotation!) (default: the current directory)
#' @param saveDataFileName the file name to save the simulated data with double-quotation, like "result1_1". (default: "BSLoutput")
#'
#' @return The simulation results (The output data was saved as BSLoutput.RData. After you load the data in R, you can find the data named as BSLoutput.)
#'
#' @export
outputResults <- function(summarize = T, directory = NULL, saveDataFileName = "BSLoutput"){
  if(summarize){
    breedingData <- lists[[1]]$breedingData
    popID <- sort(unique(breedingData$popIDsel))
    phenotype(popID = 0:max(breedingData$popID))
    muSim <- matrix(NA, length(popID), nSim)
    varSim <- matrix(NA, length(popID), nSim)
    for(sim in 1:nSim){
      breedingData <- lists[[sim]]$breedingData
      mu <- rep(NA, length(popID))
      var <- rep(NA, length(popID))
      for(i in 1:length(popID)){
        GID.now <- breedingData$GID[breedingData$popIDsel == popID[i]]
        g.now <- NULL
        for(j in GID.now){
          g.now <- c(g.now, breedingData$gValue[breedingData$phenoGID == j][1])
        }
        mu[i] <- mean(g.now)
        var[i] <- var(g.now)
      }
      muSim[, sim] <- mu
      varSim[, sim] <- var
    }
    dataSim <- cbind(muSim, varSim)
    rownames(dataSim) <- popID
    colnames(dataSim) <- c(paste("mu", 1:nSim, sep = ""), paste("var", 1:nSim, sep = ""))
  }else{
    dataSim <- lists
  }
  BSLoutput <- dataSim
  if(is.null(directory)){
    save(BSLoutput, file = "BSLoutput.RData")
  }else{
    save(BSLoutput, file = paste(directory, "/", saveDataFileName, ".RData", sep = ""))
  }
}
