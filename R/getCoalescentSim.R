#'getCalescentSim
#'
#'@useDynLib BreedingSchemeLanguage mainR
#'
#'@param nPopsSamples -pop[2] (GENOME)
#'@param effPopSize -N (GENOME)
#'@param nChr -c (GENOME)
#'@param nPiecesPerChr -pieces (GENOME)
#'@param recBTpieces -rec (GENOME)
#'@param nMrkOrMut parameter tocalculate -s (GENOME)
#'@param minMAF -maf (GENOME)
#'@param seed -seed (GENOME)
#'@param tree -tree (GENOME)
#'
#'@export
getCoalescentSim <- function(nPopsSamples=NULL, effPopSize=100, nChr=7, nPiecesPerChr=15000, recBTpieces=0.0001, nMrkOrMut=100, minMAF=0.01, seed=as.integer(Sys.time()), tree=0){

systemCall <- rep(NA, 11)
systemCall <- rep(NA, 10)
nPopsSamples <- effPopSize
systemCall[1] <- 1     # pop1
systemCall[2] <- nPopsSamples     #pop2
systemCall[3] <- effPopSize     # N
systemCall[4] <- nChr     #c
systemCall[5] <- nPiecesPerChr     # pieces
systemCall[6] <- recBTpieces     # rec
systemCall[7] <- minMAF     #maf
systemCall[8] <- seed     # seed
systemCall[9] <- tree     # tree
systemCall[10] <- round(8 * nMrkOrMut / nChr)     # s

doGenome <- function(call){
  .Call('mainR',
        par1 = call[1],
        par2 = call[2],
        par3 = as.character(call[3]),
        par4 = call[4],
        par5 = call[5],
        par6 = as.character(call[6]),
        par7 = call[7],
        par8 = call[8],
        par9 = call[9],
        par10 = call[10])
}

if(file.exists("genomeOUT.txt")){
  file.remove("genomeOUT.txt")
}
sink("genomeOUT.txt")
doGenome(systemCall)
sink()
genOut <- readLines("genomeOUT.txt", n=-1)
file.remove("genomeOUT.txt")
strtPos <- grep("SNP genetic position", genOut)
map <- NULL
for (chr in 1:nChr){
  posVec <- round(as.numeric(strsplit(genOut[strtPos[chr]+1], split=" ")[[1]]) * 100 * (nPiecesPerChr - 1) * recBTpieces, 6)
  if (length(posVec > 0)) map <- rbind(map, cbind(chr, posVec))
}
colnames(map) <- c("Chr", "Pos")
strtInd <- NULL
for (pop in 1:length(nPopsSamples)){
  strtInd <- c(strtInd, grep(paste("POP", pop, ":", sep=""), genOut))
}
genOut <- substring(genOut[strtInd], 7)
nSNP <- nchar(genOut[1])
markers <- t(array(as.numeric(unlist(strsplit(genOut, split=""))), c(nSNP, sum(nPopsSamples))))
freq <- colMeans(markers)
maf <- abs((freq > 0.5) - freq)
markers <- markers[, maf >= minMAF]
map <- map[maf >= minMAF,]
maf <- maf[maf >= minMAF]
if (length(nMrkOrMut) == 1){
  if (ncol(markers) < nMrkOrMut){
    print("Warning! Settings such that fewer markers simulated than demanded")
    print(paste("There were", ncol(markers), "markers"))
  } else{
    keepMrk <- sort(sample(ncol(markers), nMrkOrMut))
    markers <- markers[, keepMrk]
    map <- map[keepMrk, ]
    maf <- maf[keepMrk]
  }
}
nMrk <- ncol(markers)
for (chr in 1:nChr){
  mrkThisChr <- map[,"Chr"] == chr
  uniqPos <- unique(map[mrkThisChr,"Pos"])
  for (pos in uniqPos){
    mrkAtPos <- which(mrkThisChr & map[,"Pos"] == pos)
    if (length(mrkAtPos) > 1) map[mrkAtPos,"Pos"] <- map[mrkAtPos,"Pos"] + sort(round(runif(length(mrkAtPos), 0, 100 * recBTpieces), 6))
  }
}
map <- as.data.frame(map)
return(list(markers=markers, map=map))
}
