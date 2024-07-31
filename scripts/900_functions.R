do_MEM <- function(live.mats, markers) {
    require(cluster)
    .mats.med <- lapply(live.mats,function(mt,ch) {
        mm <- apply(mt[,ch],2,median)
        return(mm)
    },with(markers,name[type %in% c("phenotypic","functional")]))
    
    .mats.iqr <- lapply(live.mats,function(mt,ch) {
        mi <- apply(mt[,ch],2,IQR)
        if(min(mi)==0) {
            thresh <- min(mi[mi>0])
            thresh <- 10^floor(log10(thresh))
            mi <- replace(mi,mi==0,thresh)
        }
        return(mi)
    },with(markers,name[type %in% c("phenotypic","functional")]))
    
    .ref <- pam(bind_rows(.mats.med),1)
    .ref <- .ref$id.med
    
    .med.diffs <- sweep(bind_rows(.mats.med),2,.mats.med[[.ref]],`-`)
    .iqr.ratio <- sweep(bind_rows(.mats.iqr),2,.mats.iqr[[.ref]],`/`)
    .iqr.ratio <- 1/.iqr.ratio
    
    .mem <- abs(.med.diffs)+.iqr.ratio-1
    .mem[.med.diffs<0] <- -.mem[.med.diffs<0]
    .mem <- asinh(.mem)
    rownames(.mem) <- names(.mats.med)
    .mem <- .mem[-.ref,]
    
    return(list(MED = .mats.med,
                IQR = .mats.iqr,
                REF = .ref,
                MEM = .mem))
}

do_optim_wrap <- function(live.sim, markers) {
    require(Radviz)
    .S <- make.S(markers)
    .optim <- lapply(seq(1,10),function(i) {
        .optim.cur.cells <- do.optimRadviz(.S,live.sim,top=75)
        return(list(best=in.da(make.S(get.optim(.optim.cur.cells)),live.sim),
                    springs=get.optim(.optim.cur.cells)))
    })
    .optim.best <- lapply(.optim,function(x) x$best)
    make.S(.optim[[which.max(.optim.best)]]$springs)
}

do_shave <- function(v, fun = range, na.rm = T) {
    if (na.rm) {
        clean_v <- v[!is.na(v)]
    }
    else {
        clean_v <- v
    }
    range_v <- match.fun(fun)(clean_v)
    min_v <- min(range_v)
    max_v <- max(range_v)
    clean_v[clean_v<min_v] <- min_v
    clean_v[clean_v>max_v] <- max_v
    return(clean_v)
}
