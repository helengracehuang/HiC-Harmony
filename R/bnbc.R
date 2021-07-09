bnbc <- function(cg, batch, threshold=NULL, step=NULL,
                 downSample=FALSE, qn=TRUE, Harmony=FALSE, vars_use = c('Batch'), 
                 scalePCA=TRUE, theta=c(2), lambda=c(1), ComBat=FALSE, nbands=NULL, mod=NULL,
                 mean.only=FALSE, tol=5, bstart=1, verbose = TRUE){
    nppl <- ncol(cg)
    tacts <- contacts(cg)
    if(is.null(nbands)){
        nbands <- distanceIdx(cg, threshold = threshold, step = step)
    }
    stopifnot(bstart >= 1, bstart <= nbands)
    mat.list <- list()
    if(verbose) pb <- txtProgressBar(style = 3)
    for (ii in bstart:nbands){
        if (verbose && ii %% 50 == 0){ setTxtProgressBar(pb, (ii - bstart) / nbands) }
        mat <- getBandMatrix(cg, ii)
        mat.good <- seq_len(nrow(mat))
        
        if (downSample & min(colSums(mat))>0) {
            temp_mat <- c()
            for (kk in 1:length(mat[1,])) {
                sp <- sample(length(mat[,kk]), size=min(colSums(mat)), replace=TRUE, prob = mat[,kk])
                temp <- rep(0, length(mat[,kk]))
                for (jj in 1:length(sp)) {
                    temp[sp[jj]] = temp[sp[jj]] + 1
                }
                temp_mat <- c(temp_mat, temp)
            }
            mat <- matrix(temp_mat, ncol=length(mat[1,]), dimnames=list(c(1:length(mat[,1])),cg@colData@rownames))
        }
        
        if (qn) {
            mat[,c(1, 3, 5, 7, 9, 11, 13, 15)] <- normalize.quantiles(mat[,c(1, 3, 5, 7, 9, 11, 13, 15)], copy = FALSE)
            mat[,c(2, 4, 6, 8, 10, 12, 14, 16)] <- normalize.quantiles(mat[,c(2, 4, 6, 8, 10, 12, 14, 16)], copy = FALSE)
        }
        
        if (Harmony) {
            mat[,c(1, 3, 5, 7, 9, 11, 13, 15)] <- normalize.quantiles(mat[,c(1, 3, 5, 7, 9, 11, 13, 15)], copy = FALSE)
            mat[,c(2, 4, 6, 8, 10, 12, 14, 16)] <- normalize.quantiles(mat[,c(2, 4, 6, 8, 10, 12, 14, 16)], copy = FALSE)
            
            pca = prcomp(mat, scale. = scalePCA)
            scores = as.data.frame(pca$rotation)
            harmony_embeddings <- harmony::HarmonyMatrix(scores[,1:8], meta_data=cg@colData, 
                                                         vars_use=vars_use, theta=theta, lambda=lambda, do_pca=FALSE)
            total_embeddings <- cbind(as.data.frame(harmony_embeddings), scores[,9:16])
            total_embeddings <- data.matrix(total_embeddings)
            
            if (scalePCA) {
                reverse = t(t(pca$x %*% t(total_embeddings)) * pca$scale + pca$center)
            }
            else {
                reverse = t(t(pca$x %*% t(total_embeddings)) + pca$center)
            }
            
            mat[,c(1, 3, 5, 7, 9, 11, 13, 15)] <- normalize.quantiles.use.target(reverse[,c(1, 3, 5, 7, 9, 11, 13, 15)],mat[,1],subset=NULL)
            mat[,c(2, 4, 6, 8, 10, 12, 14, 16)] <- normalize.quantiles.use.target(reverse[,c(2, 4, 6, 8, 10, 12, 14, 16)],mat[,2],subset=NULL)
            mat <- round(mat, digits = 5)
            
        }
        
        if (ComBat) {
            if(!mean.only){
                batchvars <- bandLevelBatchVars(mat, batch)
                mat.good <- abs(rowMeans(mat)) > 0 &
                    round(rowMeans(batchvars), tol)  > 0
            }
            tryCatch({
                suppressMessages({
                    mat[mat.good,] <- ComBat(mat[mat.good,], batch, mod=mod,
                                             mean.only=mean.only)
                })
            }, error=warning)
        }
        
        mat[!mat.good,] <- 0
        tacts <- updateBand(tact_list=tacts,
                            idx=getBandIdx(nrow(tacts[[1]]), ii)-1,
                            band=mat)
    }
    if(verbose) close(pb)
    make.sym <- function(mat){
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
        mat
    }
    tacts <- lapply(tacts, make.sym)
    new.cg <- cg
    contacts(new.cg) <- tacts
    new.cg
}
