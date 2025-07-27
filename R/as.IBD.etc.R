as.IBD <-
function (x='sibs', targets=NULL, ped=FALSE) 
{
#if(is.null(targets) && is.list(x)) {targets<-x[[2]]; x<-x[[1]]}
if(is.null(targets)) targets<-attr(x,'targets')
# process character argument to numeric or list
if (is.character(x)) 
        {
if((!is.null(dim(x)))&&dim(x)[2]==3)
{
id<-unique(as.vector(x))
nid<-length(id)
fid<-mid<-sex<-rep(0,nid)
for(i in 1:nrow(x))
{
j<-match(x[i,1],id)
mid[j]<-x[i,2]
fid[j]<-x[i,3]
}
for(j in 1:nid)
{
if(id[j]%in%mid) sex[j]<-2 else if(id[j]%in%fid) sex[j]<-1
}
x<-ped(id,fid,mid,sex)
} else {
	  x <- switch(pmatch(x, c("sibs", "parent-child", "half-sibs", 
            "cousins", "half-cousins", "second-cousins", "double-first-cousins", 
            "quadruple-half-first-cousins", "3cousins-cyclic", 
            "3cousins-star", "trio")), c(0.25, 0.5, 0.25), c(0, 
            1, 0), c(0.5, 0.5, 0), c(0.75, 0.25, 0), c(0.875, 
            0.125, 0), c(0.9375, 0.0625, 0), c(0.5625, 0.375, 
            0.0625), c(17, 14, 1)/32, list(pr = c(0.015625, 0.046875, 
            0.046875, 0.140625, 0.046875, 0.140625, 0.140625, 
            0.421875), patt = matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 
            2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 3, 3, 3, 3, 3, 
            3, 3, 3, 4, 4, 4, 4, 2, 2, 3, 4, 1, 1, 3, 5, 3, 4, 
            4, 5, 3, 5, 5, 6), 8, 6)), list(pr = c(0.0625, 0.1875, 
            0.1875, 0.1875, 0.375), patt = matrix(c(1, 1, 1, 
            1, 1, 2, 2, 2, 2, 2, 1, 1, 3, 3, 3, 3, 3, 4, 4, 4, 
            1, 4, 1, 3, 5, 4, 5, 5, 5, 6), 5, 6)), list(patt = c(1, 
            2, 3, 4, 1, 3)))
	  if(is.null(x)) stop('invalid character argument to as.IBD')
}
	  }

# process numeric to list
if(is.numeric(x))
	{
	if(is.null(dim(x)))
	{
	if(length(x)==3) {
		kappa <- x
            patt <- matrix(c(1, 2, 3, 4, 1, 2, 3, 1, 1, 2, 2, 
                1), 3, 4, byrow = TRUE)
            wk <- which(kappa > 0)
            x <- list(pr = kappa[wk], patt = patt[wk, , drop = FALSE])
	} else if(length(x)==9) {
		Delta <- x
            patt <- matrix(c(1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 
                2, 1, 1, 2, 3, 1, 2, 1, 1, 1, 2, 3, 3, 1, 2, 
                1, 2, 1, 2, 1, 3, 1, 2, 3, 4), 9, 4, byrow = TRUE)
            wk <- which(Delta > 0)
            x <- list(pr = Delta[wk], patt = patt[wk, , drop = FALSE])
	} else if(length(x)%%2==0) {
	x<-list(pr=1,patt=matrix(x,1,length(x)))
	} else stop('invalid numeric vector argument to as.IBD')
} else if(length(dim(x))!=2) stop('invalid numeric array argument to as.IBD') else
x<-list(pr=rep(1,nrow(x))/nrow(x),patt=x)
}

# x should now be a list - but this may be ped or pedList, or a (ped/pedList, targets) list
if(is.ped(x)||is.pedList(x))
{
	return(pedigreeIBD(x,targets,ped=ped))
} else {
	if(is.ped(x[[1]])||is.pedList(x[[1]]))
	{
		if(is.null(targets)) targets<-x[[2]]
		return(pedigreeIBD(x[[1]],targets))
	} 
}
# fix a special case
if(is.null(x$patt)) stop('invalid x')
if(is.null(dim(x$patt))) x$patt<-matrix(x$patt,nrow=1)
if(is.null(x$pr)) x$pr<-rep(1,nrow(x$patt))/nrow(x$patt)
# now has to be (pr,patt) list: just return suitably annotated
if(is.null(targets)) targets<-letters[1:(ncol(x$patt)/2)]
structure(x,class='IBD',targets=targets)
} # end of function as.IBD

convertIBD <- function (x='sibs', targets=NULL, ped=FALSE) as.IBD(x,targets,ped)

pedigreeIBD <-
function (x, targets, cond = TRUE, ped=FALSE, quiet = TRUE, verbose = FALSE) 
{
    if (is.ped(x)) {
	  if(missing(targets)||is.null(targets)) targets<-x$ID
        pdf <- as.data.frame(x)
        id <- pdf$id
        trim <- (id %in% targets) | (id %in% ancestors(x, targets))
        id <- id[trim]
        fid <- pdf$fid[trim]
        mid <- pdf$mid[trim]
        sex <- pdf$sex[trim]
    }
    else {
	  if(missing(targets)||is.null(targets)) targets<-unlist(lapply(x,function(x) x$ID))
        id <- fid <- mid <- sex <- NULL
        for (ic in 1:length(x)) {
            pdf <- as.data.frame(x[[ic]])
            idi <- pdf$id
            fidi <- pdf$fid
            midi <- pdf$mid
            sexi <- pdf$sex
            trim <- (idi %in% targets) | (idi %in% ancestors(x[[ic]], 
                targets[targets %in% idi]))
            id <- c(id, idi[trim])
            fid <- c(fid, fidi[trim])
            mid <- c(mid, midi[trim])
            sex <- c(sex, sexi[trim])
        }
    }
    if(!all(targets%in%id)) stop('invalid targets')
    x <- list(ID = id, FIDX = match(fid, id, nomatch = 0), MIDX = match(mid, 
        id, nomatch = 0), SEX = sex)
    ida <- id[1]
    pattern <- matrix(1:2, 1, 2)
    mlp <- 2
    count <- 1
    for (i in 2:length(id)) {
        idnext <- id[i]
        iid <- match(idnext, id)
        ida <- c(ida, idnext)
        wf <- match(fid[iid], ida)
        if (is.na(wf)) {
            pattern <- cbind(pattern, max(pattern) + 1)
        }
        else {
            pattern <- rbind(cbind(pattern, pattern[, 2 * wf - 
                1]), cbind(pattern, pattern[, 2 * wf]))
            count <- c(count, count)
        }
        wm <- match(mid[iid], ida)
        if (is.na(wm)) {
            pattern <- cbind(pattern, max(pattern) + 1)
        }
        else {
            pattern <- rbind(cbind(pattern, pattern[, 2 * wm - 
                1]), cbind(pattern, pattern[, 2 * wm]))
            count <- c(count, count)
        }
        keep <- prune(x, targets, ida, verbose) | (ida %in% targets)
        ida <- ida[keep]
        pattern <- pattern[, rep(keep, rep(2, length(keep))), 
            drop = FALSE]
        if (cond) 
            pattern <- t(apply(pattern, 1, function(r) minimalPattern(r)))
        else pattern <- t(apply(pattern, 1, function(r) match(r, 
            unique.default(r))))
        o <- statnet.common::order(pattern)
        pattern <- pattern[o, , drop = FALSE]
        nrp <- nrow(pattern)
        w <- rep(0, nrp)
        if (nrp > 1) 
            for (j in 1:(nrp - 1)) w[j] <- as.integer(any(pattern[j, 
                ] != pattern[j + 1, ]))
        w[nrp] <- 1
        cumct <- cumsum(count[o])
        count <- diff(c(0, cumct[which(w == 1)]))
        pattern <- pattern[w == 1, , drop = FALSE]
        mlp <- max(mlp, length(pattern))
        if (verbose) {
            message(paste(idnext, dim(pattern), "\n"))
            dimnames(pattern) <- list(rep(" ", nrow(pattern)), 
                as.vector(rbind(ida, " ")))
            message(pattern)
        }
    }
    dimnames(pattern) <- list(rep(" ", nrow(pattern)), 
        as.vector(rbind(ida, " ")))
    mti <- match(targets, ida)
    pattern <- pattern[, as.vector(rbind(2 * mti - 1, 2 * mti))]
    if (!quiet) 
        print(cbind(pr = count/sum(count), pattern = pattern))
if(ped) IBD <- structure(list(pr = count/sum(count), patt = pattern), 
        ped = ped(id, fid, mid, sex), mlp = mlp, targets=targets, class='IBD') else
	  IBD <- structure(list(pr = count/sum(count), patt = pattern), 
        mlp = mlp, targets=targets, class='IBD')
    IBD
}
