mytriplot <- function(object1, object2, object3, a1 = 1, a2 = 2, fact, labs = rep('',nrow(object2)), nv = nrow(object2), normalize = FALSE, vars = NULL, rdavars = NULL, pca = FALSE, center=FALSE, effsBoth = FALSE, cex1= 0.3, xlim = NULL, ylim = NULL) {

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

    if(normalize) {
        ob1[,1] <- ob1[,1] / max(abs(ob1[,1]))
        ob1[,2] <- ob1[,2] / max(abs(ob1[,2]))
        ob2[,1] <- ob2[,1] / max(abs(ob2[,1]))
        ob2[,2] <- ob2[,2] / max(abs(ob2[,2]))
        ob3[,1] <- ob3[,1] / max(abs(ob3[,1]))
        ob3[,2] <- ob3[,2] / max(abs(ob3[,2]))
    }

    if(is.null(vars)) {
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
    } else {
        ob2 <- array(ob2[labs %in% vars,], dim=c(length(vars), 2))
        newlabs <- labs[labs %in% vars]
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
