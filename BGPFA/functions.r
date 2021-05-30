
monteCarloP <- function(x, pn='p') {
    if(pn == 'n') {
        res <- (1 + sum(x >= 0, na.rm=TRUE)) / (1 + length(x))
    } else if(pn == 'p') {
        res <- (1 + sum(x <= 0, na.rm=TRUE)) / (1 + length(x))
    }
} # https://arxiv.org/pdf/1603.05766.pdf

orient_axes <- function(extracted, reference = NULL) {
    K <- dim(extracted)[[2]]
    if(is.null(reference)) reference <- extracted[,,sample(dim(extracted)[[3]],1)]
    oriented <- array(NA, dim=dim(extracted))
    which_match <- array(NA,dim(extracted)[c(2,3)])
    reflected <- array(NA,dim(extracted)[c(2,3)])
    order_axes <- order(sqrt(colSums(reference^2)))
    for(draw in 1:dim(extracted)[[3]]) {
        available <- rep(TRUE,K)
        for(axis in order_axes) {
            idx <- (1:K)[available]
            cors <- sapply(idx, function(k) cor(extracted[,k,draw], reference[,axis]))
            which_match[axis,draw] <- idx[which.max(abs(cors))]
            available[which_match[axis,draw]] <- FALSE
            reflected[axis,draw] <- sign(cors[which.max(abs(cors))])
            oriented[,axis,draw] <- reflected[axis,draw] * extracted[,which_match[axis,draw],draw]
        }
    }
    return(list(oriented=oriented, which_match=which_match, reflected=reflected))
} # swap axes such that they best match a reference (by default a random draw). If axes are degenerate during model fitting, this should orient them. Can this be made probabilistic rather than iterating max fit?

no_reflections <- function(extracted, reference = NULL) {
    K <- dim(extracted)[[2]]
    if(is.null(reference)) reference <- extracted[,,sample(dim(extracted)[[3]],1)]
    oriented <- array(NA, dim=dim(extracted))
    which_match <- array(NA,dim(extracted)[c(2,3)])
    reflected <- array(NA,dim(extracted)[c(2,3)])
    order_axes <- order(sqrt(colSums(reference^2)))
    for(draw in 1:dim(extracted)[[3]]) {
        available <- rep(TRUE,K)
        for(axis in order_axes) {
            idx <- (1:K)[available]
            cors <- sapply(idx, function(k) cor(extracted[,k,draw], reference[,axis]))
            which_match[axis,draw] <- idx[which.max(abs(cors))]
            available[which_match[axis,draw]] <- FALSE
            reflected[axis,draw] <- sign(cors[which.max(abs(cors))])
            oriented[,axis,draw] <- reflected[axis,draw] * extracted[,which_match[axis,draw],draw]
        }
    }
    return(list(oriented=oriented, which_match=which_match, reflected=reflected))
} # swap axes such that they best match a reference (by default a random draw). If axes are degenerate during model fitting, this should orient them. Can this be made probabilistic rather than iterating max fit?


prob_axes <- function(extracted, maxmult) {
    K <- dim(extracted)[[2]]
    reference <- array(NA, dim=c(dim(extracted)[[1]],0))
    which_match <- array(NA,dim(extracted)[c(2,3)])
    reflected <- array(NA,dim(extracted)[c(2,3)])
    oriented <- array(NA, dim=c(dim(extracted)[1],0,dim(extracted)[3]))
    refcount <- NULL
    for(draw in sample(1:dim(extracted)[[3]])) {
        K_ref <- dim(reference)[[2]]
        order_axes <- order(sqrt(colSums(extracted[,,draw]^2)))
        available <- rep(TRUE,K_ref)
        for(axis in order_axes) {
            sub <- extracted[,axis,draw]
            idx <- (1:K_ref)[available]
            cors <- cor(sub, reference[,idx])
            if(K_ref < (maxmult * K)) {
                cors <- unlist(c(cors,0.5))
            }
            wm <- which.max(abs(cors))
            if(wm == (length(idx) + 1)) {
                K_ref <- K_ref + 1
                reference <- cbind(reference,sub)
                refcount <- c(refcount,1)
                which_match[axis,draw] <- K_ref
                available <- c(available, FALSE)
                reflected[axis,draw] <- sign(cors[wm])
                oriented <- abind:::abind(oriented, array(0, dim=c(dim(oriented)[1],1,dim(oriented)[3])), along=2)
                oriented[,K_ref,draw] <- reflected[axis,draw] * sub
            } else {
                reference[,idx[wm]] <- (refcount[idx[wm]] * reference[,idx[wm]] + sign(cors[wm]) * sub) / (refcount[idx[wm]]+1)
                refcount[idx[wm]] <- refcount[idx[wm]] + 1
                which_match[axis,draw] <- idx[wm]
                available[which_match[axis,draw]] <- FALSE
                reflected[axis,draw] <- sign(cors[wm])
                oriented[,idx[wm],draw] <- reflected[axis,draw] * sub
            }
        }
    }
    return(list(oriented=oriented, which_match=which_match, reflected=reflected, refcount=refcount, reference=reference))
}# use correlation with reference to determine a probability of a match rather than deterministically giving it one. reference could be added to, such that there is always a probability that the axis doesn't exist at all during a draw


mytriplot <- function(object1, object2, object3, a1 = 1, a2 = 2, fact, labs = rep('',nrow(object2)), nv = nrow(object2), normalize = FALSE, vars = NULL, rdavars = NULL, pca = FALSE, center=FALSE, sigs = NULL, effsBoth = FALSE, cex1= 0.3, xlim = NULL, ylim = NULL) {

    if(pca) {
        pcob <- prcomp(object1, center=center)
        ob1 <- pcob$x[,c(a1,a2)]
        ob2 <- (object2 %*% pcob$rotation)[,c(a1,a2)]
        ob3 <- (object3 %*% pcob$rotation)[,c(a1,a2)]
    } else {
        ob1 <- apply(object1[,c(a1,a2)], 2, function(x) x - mean(x))
        ob2 <- object2[,c(a1,a2)]
        ob3 <- object3[,c(a1,a2)]
    }

    if(is.null(vars) & is.null(sigs)) {
        effs <- setNames(apply(ob2, 1, function(x) sqrt(sum(x^2))), labs)
        if(effsBoth) {
            keep <- tail(order(effs), nv)
            ob2 <- ob2[keep, ]
        } else {
            keep1 <- tail(order(abs(ob2[,1])),nv)
            keep2 <- tail(order(abs(ob2[,2])),nv)
            keep <- union(keep1, keep2)
            res <- list(rev(setNames(ob2[keep1,1], labs[keep1])), rev(setNames(ob2[keep2,2], labs[keep2])))
            ob2 <- ob2[keep, ]
        }
        newlabs = labs[keep]
    } else if(!is.null(vars)) {
        ob2 <- array(ob2[labs %in% vars,], dim=c(length(vars), 2))
        newlabs <- labs[labs %in% vars]
        res <- newlabs
    } else {
        mysigs <- apply(sigs[,c(a1,a2)], 1, function(x) as.logical(sum(x)))
        ob2 <- array(ob2[mysigs,], dim=c(sum(mysigs), 2))
        newlabs <- labs[mysigs]
        res <- newlabs
    }
    if(is.null(rdavars)) {
        effsRDA <- apply(ob3, 1, function(x) sqrt(sum(x^2)))
        if(effsBoth) {
            keepRDA <- tail(order(effsRDA), nv)
            ob3 <- ob3[keepRDA, ]
        } else {
            keep1 <- tail(order(abs(ob3[,1])),nv)
            keep2 <- tail(order(abs(ob3[,2])),nv)
            keep <- union(keep1, keep2)
            res <- list(res, list(rev(setNames(ob3[keep1,1], rownames(ob3)[keep1])), rev(setNames(ob3[keep2,2], rownames(ob3)[keep2]))))
            ob3 <- ob3[keep, ]
        }
        newnames <- rownames(ob3)
    } else {
        ob3 <- ob3[rdavars,]
        newnames = rownames(ob3)
        res <- list(res, newnames)
    }

    if(normalize) {
        ob1[,1] <- ob1[,1] / max(abs(ob1[,1]))
        ob1[,2] <- ob1[,2] / max(abs(ob1[,2]))
        ob2[,1] <- ob2[,1] / max(abs(ob2[,1]))
        ob2[,2] <- ob2[,2] / max(abs(ob2[,2]))
        ob3[,1] <- ob3[,1] / max(abs(ob3[,1]))
        ob3[,2] <- ob3[,2] / max(abs(ob3[,2]))
    }

    if(is.factor(fact)) {
        cols <- fact
    } else {
        cols <- colorRampPalette(brewer.pal(8,"RdYlGn"))(length(fact))[findInterval(fact, sort(fact))]
        cols[is.na(cols)] <- rgb(0, 0, 0, max = 255, alpha = 75)
    }

    if(!is.null(xlim)) {
        offplotxp <- sapply(ob2[,1],function(x) x > xlim[[2]])
        if(sum(offplotxp) > 1) {
            ob2[offplotxp,2] <- apply(ob2[offplotxp,],1,function(x) x[[2]] * xlim[[2]] / x[[1]])
            ob2[offplotxp,1] <- xlim[[2]]
        } else if(sum(offplotxp) == 1) {
            ob2[offplotxp,2] <- ob2[offplotxp,2] * xlim[[2]] / ob2[offplotxp,1]
            ob2[offplotxp,1] <- xlim[[2]]
        }
        offplotxn <- sapply(ob2[,1],function(x) x < xlim[[1]])
        if(sum(offplotxn) > 1) {
            ob2[offplotxn,2] <- apply(ob2[offplotxn,],1,function(x) x[[2]] * xlim[[1]] / x[[1]])
            ob2[offplotxn,1] <- xlim[[1]]
        } else if(sum(offplotxn) == 1) {
            ob2[offplotxn,2] <- ob2[offplotxn,2] * xlim[[1]] / ob2[offplotxn,1]
            ob2[offplotxn,1] <- xlim[[1]]
        }
        offplotx <- offplotxp | offplotxn
    } else {
        offplotx <- rep(FALSE,nrow(ob2))
    }
    if(!is.null(ylim)) {
        offplotyp <- sapply(ob2[,2],function(x) x > ylim[[2]])
        if(sum(offplotyp) > 1) {
            ob2[offplotyp,1] <- apply(ob2[offplotyp,],1,function(x) x[[1]] * ylim[[2]] / x[[2]])
            ob2[offplotyp,2] <- ylim[[2]]
        } else if(sum(offplotyp) == 1) {
            ob2[offplotyp,1] <- ob2[offplotyp,1] * ylim[[2]] / ob2[offplotyp,2]
            ob2[offplotyp,2] <- ylim[[2]]
        }
        offplotyn <- sapply(ob2[,1],function(x) x < xlim[[1]])
        if(sum(offplotyn) > 1) {
            ob2[offplotyn,1] <- apply(ob2[offplotyn,],1,function(x) x[[1]] * xlim[[1]] / x[[2]])
            ob2[offplotyn,2] <- ylim[[1]]
        } else if(sum(offplotyn) == 1) {
            ob2[offplotyn,1] <- ob2[offplotyn,1] * xlim[[1]] / ob2[offplotyn,2]
            ob2[offplotyn,2] <- ylim[[1]]
        }
        offploty <- offplotyp | offplotyn
    } else {
        offploty <- rep(FALSE,nrow(ob2))
    }
    offplot <- offplotx | offploty


    if(is.null(xlim)) xlim = range(c(ob1[,1],ob2[,1],ob3[,1]))
    if(is.null(ylim)) ylim = range(c(ob1[,2],ob2[,2],ob3[,2]))
    plot(ob1,
         col  = cols,
         pch  = 16,
         xlim = xlim,
         ylim = ylim,
         cex  = 0.5)
    if(is.factor(fact)) {
        legend("bottomright", legend = levels(fact), col = 1:nlevels(fact), pch = 16)
    }

    if(any(!offplot)) text(ob2[!offplot,], labels = newlabs[!offplot], cex=cex1)
    if(any(offplot)) text(ob2[offplot,], labels = newlabs[offplot], cex=cex1, col='red')
    text(ob3, labels = newnames, cex=cex1, col = 'orange')
    return(res)
}

findIntRuns <- function(run){
    rundiff <- c(1, diff(run))
    difflist <- split(run, cumsum(rundiff!=1))
    unlist(lapply(difflist, function(x){
        if(length(x) %in% 1:2) as.character(x) else paste0(x[1], "-", x[length(x)])
    }), use.names=FALSE)
} ## https://stackoverflow.com/questions/16911773/collapse-runs-of-consecutive-numbers-to-ranges

filter_large_stan_csv <- function(infile, outfile, params) {
    hi <- data.table:::fread(infile,nrows=1,skip='lp__',sep=',')
    s <- paste(c('lp__',params), collapse='|')
    cols <- findIntRuns(grep(s,colnames(hi)))
    cols <- paste(cols, collapse=',')
    system(paste0('cut -f ',cols,' -d \',\' ',infile,' > ',outfile))
}
