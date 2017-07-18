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

########################################
kin.blup4.4 <- function(data,geno,pheno,GAUSS=FALSE,K=NULL,fixed=NULL,covariate=NULL,PEV=FALSE,n.core=1,theta.seq=NULL,reduce=FALSE,R=NULL)  {

	make.full <- function(X) {
		svd.X <- svd(X)
		r <- max(which(svd.X$d>1e-8))
		return(as.matrix(svd.X$u[,1:r]))
	}
	
	names <- colnames(data)	
	ypos <- match(pheno,names)
	if (is.na(ypos)) {
		stop("Phenotype name does not appear in data.")
	} else {
		y <- data[,ypos]
	}

	if (!is.null(R)&(length(R)!=length(y))) {
		stop("Length of R does not equal length of y")
	}
	
	not.miss <- which(!is.na(y))
	if (length(not.miss)<length(y)) {
		data <- data[not.miss,]
		y <- y[not.miss]
		if (!is.null(R)) {R <- R[not.miss]}
	} 
	n <- length(y)
   
	X <- matrix(1,n,1)		
	if (!is.null(fixed)) {
		p <- length(fixed)
		for (i in 1:p) {
			xpos <- match(fixed[i],names)
			xx <- factor(data[,xpos])	
			if (length(unique(xx)) > 1) {X <- cbind(X,model.matrix(~x-1,data.frame(x=xx)))}
		}
	}
	if (!is.null(covariate)) {
		p <- length(covariate)
		for (i in 1:p) {
			xpos <- match(covariate[i],names)
			X <- cbind(X,data[,xpos])
		}
	}	

	gid.pos <- match(geno,names)
	if (is.na(gid.pos)) {stop("Genotype name does not appear in data.")}

	not.miss.gid <- as.character(unique(data[,gid.pos]))
	if (is.null(K)) {
		if (reduce) {print("reduce=TRUE is not valid for independent genotypes. Proceeding without reduction.")}
		gid <- not.miss.gid
		v <- length(gid)
		Z <- matrix(0,n,v)
		colnames(Z) <- gid
		Z[cbind(1:n,match(data[,gid.pos],gid))] <- 1
		X2 <- make.full(X)
		ans <- mixed.solve(y=y,X=X2,Z=Z,SE=PEV)
		resid <- y-X2%*%ans$beta-Z%*%ans$u
		if (PEV) {
			return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,PEV=ans$u.SE^2,resid=resid))
		} else {
			return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u,resid=resid))
		}
	} else {

		if (class(K)=="dist") {K <- as.matrix(K)}
		gid <- rownames(K)
		ix.pheno <- match(not.miss.gid,gid)
		miss.pheno.gid <- which(is.na(ix.pheno))
		if (length(miss.pheno.gid)>0) {
			stop(paste("The following lines have phenotypes but no genotypes:",paste(not.miss.gid[miss.pheno.gid],collapse=" ")))
		}
		miss.gid <- setdiff(gid,not.miss.gid)
		ix <- c(ix.pheno,match(miss.gid,gid))
		K <- K[ix,ix]
		v <- length(not.miss.gid)
		Z <- matrix(0,n,v)
		Z[cbind(1:n,match(data[,gid.pos],not.miss.gid))] <- 1
	
		if (!is.null(R)) {
			sqrt.R <- sqrt(R)
			X2 <- X/sqrt.R
			y2 <- y/sqrt.R
			Z2 <- Z/sqrt.R
		} else {
			X2 <- X
			y2 <- y
			Z2 <- Z
		}

		if ((n > v)&(reduce)) {
			#transform
			w <- sqrt(diag(crossprod(Z2)))
			X2 <- make.full(crossprod(Z2,X2)/w)
			y2 <- crossprod(Z2,y2)/w
			Z2 <- cbind(diag(w),matrix(0,v,nrow(K)-v))
			reduced <- TRUE
		} else {
			X2 <- make.full(X2)
			Z2 <- cbind(Z2,matrix(0,n,nrow(K)-v))			
			reduced <- FALSE
		}

		rm(X,Z,y)
	
		if (!GAUSS) {
			ans <- mixed.solve(y=y2,X=X2,Z=Z2,K=K,SE=PEV)
			ix <- match(gid,rownames(ans$u))
			if (reduced) {
				resid <- NULL
			} else {			
				resid <- y2-X2%*%ans$beta-Z2%*%ans$u
			}
			if (PEV) {
				return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid))
			} else {
				return(list(Vg=ans$Vu,Ve=ans$Ve,g=ans$u[ix],resid=resid))
			}
		
		} else {

		if (is.null(theta.seq)) {
			theta <- setdiff(seq(0,max(K),length.out=11),0)
		} else {
			theta <- theta.seq
		}
		n.profile <- length(theta)		
		ms.fun <- function(theta) {
	    	soln <- list()
    		n.t <- length(theta)
	    	for (i in 1:n.t) {
   			   soln[[i]] <- mixed.solve(y=y2,X=X2,Z=Z2,K=exp(-(K/theta[i])^2),SE=PEV)
    		}
	    	return(soln)
  		}

		if (n.core > 1) {
    		library(parallel)
    		it <- split(theta,factor(cut(theta,n.core,labels=FALSE)))
	    	soln <- unlist(mclapply(it,ms.fun,mc.cores=n.core),recursive=FALSE)
  		} else {
    		soln <- ms.fun(theta)
	  	}      

 		LL <- rep(0,n.profile)
	  	for (i in 1:n.profile) {LL[i] <- soln[[i]]$LL}
		ans <- soln[[which.max(LL)]]	
		profile <- cbind(theta,LL)
		ix <- match(gid,rownames(ans$u))

		if (reduced) {
			resid <- NULL
		} else {			
			resid <- y2-X2%*%ans$beta-Z2%*%ans$u
		}

		if (PEV) {
			return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],PEV=ans$u.SE[ix]^2,resid=resid))
		} else {
			return(list(Vg=ans$Vu,Ve=ans$Ve,profile=profile,g=ans$u[ix],resid=resid))
		}	

		} #else GAUSS
	} #else is.null(K)
} #kin.blup
##########################


    model <- kin.blup4.4(data = data, geno = "GID", pheno = "y", K = K, R = R, reduce = reduce)
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
