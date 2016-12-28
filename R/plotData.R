#'Plot the results
#'
#'@param ymax the maximum value of y-axis (default: the maximun value in the data)
#'@param add if T a new result will be added to an existing plot, if F a new plot will be drawn (default)
#'@param addDataFileName the name to save the summarized data for the next simulation with double-quotation, like "plot1_1". (default: "plotData")
#'@param popID vector of the population IDs you want plotted
#'
#'@return A ggplot object of the simulation results
#'
#'@export
plotData <- function(simEnv, ymax = NULL, add = F, addDataFileName = "plotData", popIDplot = NULL){
  plotBase <- is.null(popIDplot)
  if (plotBase) popIDplot <- sort(unique(simEnv$sims[[1]]$breedingData$popIDsel))
  
  getMeans <- function(sim){
    breedingData <- sim$breedingData
    if (plotBase) popID <- breedingData$popIDsel
    else popID <- breedingData$popID
    mu <- tapply(breedingData$gValue, as.factor(popID), mean)
    mu <- mu[as.character(popIDplot)]
  }
  muSim <- t(sapply(simEnv$sims, getMeans))
  
  muSim <- muSim - muSim[, 1]
  g <- NULL
  group <- NULL
  col <- NULL
  size <- NULL
  nGenPlot <- length(popIDplot)
  for(sim in 1:simEnv$nSim){
    g <- c(g, muSim[sim, ])
    group <- c(group, rep(sim, nGenPlot))
    size <- c(size, rep(1, nGenPlot))
  }
  g <- c(g, apply(muSim, 2, mean))
  group <- c(group, rep(simEnv$nSim + 1, nGenPlot))
  size <- c(size, rep(2, nGenPlot))
  data <- data.frame(g = g, popID = rep(0:(nGenPlot - 1), simEnv$nSim + 1), size = size, group = group, scheme = rep(1, length(g)))
  if (add){
    load(file = paste(addDataFileName, ".RData", sep = ""))
    data$scheme <- data$scheme + max(data.previous$scheme)
    data$group <- data$group + max(data.previous$group)
    data <- rbind(data, data.previous)
  }
  data.previous <- data
  save(data.previous, file = paste(addDataFileName, ".RData", sep = ""))
  p <- ggplot(data = data, aes(x = popID, y = g))
  p <- p + geom_line(aes(size = factor(size), linetype = factor(scheme), group = factor(group)))
  if (is.null(ymax)) {
    p <- p + ylim(min(data$g), max(data$g))
  }
  else {
    p <- p + ylim(min(data$g), ymax)
  }
  p <- p + scale_size_manual(name = "", values = c(0.2, 3), labels = c("Each", "Mean"))
  p <- p + labs(title = "", x = "Generation", y = "Genetic improvement")
  p <- p + guides(size = guide_legend("Lines"))
  p <- p + guides(linetype = guide_legend("Scheme"))
  print(p)
}
