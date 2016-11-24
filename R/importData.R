#'Import phased data from hapmap format file defining founder haplotypes
#'
#'@param fileName if null read a file called "founderHaplotypes.hmp" (default) else read
#'the specified file.
#'
#'@return NULL
#'
#'@export
#'
importData <- function(fileName = "founderHaplotypes.hmp"){
  hm <- read.table(fileName)
  BreedingSchemeLanguage::importedData <- phasedHapMap2mat(hm)
}

#' Transform a data.frame with a hapmap data in it into a marker dosage and map list
#' 
#' @param hm The data.frame
#' @param firstCol Number of the first column you want in the marker dosage file
#' @param lastCol Number of the last column you want in the marker dosage file
#' 
#' @return list with markers and map objects
#' 
# This function will not allow missing data so anything missing will be
# set to the major allele
# Assume marker names in column 1 and genotypes (e.g. "A/G") in column 2
# Assume chromosome and position are in columns 3 and 4 of the hapmap
# and that position is in cM: WARNING cM is probably not standard hapmap
phasedHapMap2mat <- function(hm, firstCol=12, lastCol=ncol(hm)){
  nInd <- lastCol - 11
  nMrk <- nrow(hm)
  mrkNames <- as.character(hm[[1]])
  map <- data.frame(Chr=hm[[3]], Pos=hm[[4]])
  hm <- apply(hm[c(2, firstCol:lastCol)], 2, as.character)
  if (nchar(hm[1, 2]) == 1){
    hapMap2num <- function(vec){ # Convert to numeric if one character codes
      hets <- c("W","R","Y","K","M","S") # NOTE: hets become missing here
      missings <- c("-", "0","N")
      alleles <- unlist(strsplit(vec[1], "/"))
      vec <- vec[-1]
      vec[vec %in% c(hets, missings)] <- NA
      vec[alleles[1]==vec] <- 1
      vec[alleles[2]==vec] <- 0
      return(suppressWarnings(as.numeric(vec)))
    }
  } else{ # Convert to numeric if two character codes
    hapMap2num <- function(vec){
      missings <- c("-", "0","N")
      alleles <- unlist(strsplit(vec[1], "/"))
      vec <- vec[-1]
      codes <- c(alleles, missings)
      numCodes <- c(1, 0, rep(NA, length(missings)))
      gam1 <- sapply(substr(vec, 1, 1), function(code) numCodes[code == codes])
      gam2 <- sapply(substr(vec, 2, 2), function(code) numCodes[code == codes])
      return(cbind(gam1, gam2))
    }
  }
  res <- apply(hm, 1, hapMap2num) 
  return(list(markers=res, map=map))
}

#' Generate a data.frame with a hapmap data in it to test phasedHapMap2mat
#' 
#' @param nInd The number of individuals with marker data
#' @param nMrk The number of markers
#' @param nChr The number of chromosomes the species has
#' @param lenChr The length of the chromosomes (all equal) in cM
#' @param maf The desired minor allele frequency of each marker
#' @param nCharCode Whether the genotype codes should be one or two characters
#' 
#' @return A data.frame with hapmap data in it
#' 
# Quick function to create a HapMap data.frame for testing
#
simHapMap <- function(nInd=20, nMrk=50, nChr=7, lenChr=150, maf=runif(nMrk), nCharCode=2){
  nucl <- c("A", "C", "G", "T")
  gt <- replicate(nMrk, sample(4, 2))
  alleles <- paste(nucl[gt[1,]], nucl[gt[2,]], sep="/")
  chr <- sort(sample(nChr, nMrk, replace=T))
  pos <- round(lenChr * runif(nMrk), 2)
  pos <- pos[order(chr, pos)]
  hapmap <- data.frame(locus=paste("locus", 1:nMrk, sep=""), alleles=alleles, chrom=chr, pos=pos)
  if (!(nCharCode %in% 1:2)) nCharCode <- 2
  nInd <- nInd * 2 / nCharCode
  gam <- matrix(rbinom(nInd*nMrk, 1, 1 - maf) + 1, nMrk, nInd)
  gam1 <- apply(gam, 2, function(gamVec) nucl[gt[cbind(gamVec, 1:nMrk)]])
  if (nCharCode == 2){
    gam <- matrix(rbinom(nInd*nMrk, 1, 1 - maf) + 1, nMrk, nInd)
    gam2 <- apply(gam, 2, function(gamVec) nucl[gt[cbind(gamVec, 1:nMrk)]])
  } else gam2 <- NULL
  hm <- matrix(paste(gam1, gam2, sep=""), nMrk, nInd)
  hapmap <- cbind(hapmap, matrix(NA, nMrk, 7), hm)
  colnames(hapmap)[12:(nInd+11)] <- paste("ind", 1:nInd, sep="")
  return(hapmap)
}
