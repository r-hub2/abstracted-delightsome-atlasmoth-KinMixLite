aca2gt <-
function (mixture, aca) 
{
    nm <- length(names(aca))
    allele1 <- allele2 <- rep(0, nm)
    for (j in 1:nm) {
        m <- names(aca)[j]
        alleles <- mixture$data[[m]]$allele
        w <- which(aca[[m]] > 0)
        if (length(w) == 1) 
            allele1[j] <- allele2[j] <- alleles[w]
        else {
            allele1[j] <- alleles[w[1]]
            allele2[j] <- alleles[w[2]]
        }
    }
    data.frame(marker = names(aca), allele1, allele2)
}
addU <-
function (IBD, nU) 
{
    extra <- matrix(0, nrow(IBD$patt), 2 * nU)
    m <- apply(IBD$patt, 1, max)
    for (r in 1:nrow(IBD$patt)) extra[r, ] <- m[r] + (1:(2 * 
        nU))
    IBD$patt <- cbind(IBD$patt, extra)
    IBD
}
as.gt <-
function (res, ind) 
{
#    z <- subset(res, Sim == ind)
    z<-res[res['Sim']==ind,]
    allele <- z$Allele
    allele[allele == "X"] <- 0
    allele[allele == "Y"] <- 1
    allele <- matrix(as.numeric(allele), ncol = 2, byrow = TRUE)
    data.frame(marker = unique(z$Marker), allele1 = allele[, 
        1], allele2 = allele[, 2])
}
binary <-
function (n, d = NA) 
{
    res <- NULL
    while (ifelse(is.na(d), n > 0, d != 0)) {
        d <- abs(d) - 1
        res <- c(n%%2, res)
        n <- n%/%2
    }
    res
}
checkpeaks <-
function (x, db, fix = 0.003) 
{
    first <- TRUE
    if ("height" %in% names(x)) {
        for (m in unique(x$marker)) {
            a <- x[x$marker == m, "allele"]
            adb <- db[db$marker == m, "allele"]
            w <- !a %in% adb
            if (any(w)) {
                if (first) {
                  message("peaks at alleles not in database\n")
                  first <- FALSE
                }
                message(paste(m, ":", a[w], ":", x[x$marker == m, "height"][w], 
                  "\n"))
                if (fix > 0) {
                  db <- merge(db, data.frame(marker = m, allele = a[w], 
                    frequency = fix), all = T)
                  f <- db[db$marker == m, ]$frequency
                  db[db$marker == m, ]$frequency <- f/sum(f)
                }
            }
        }
    }
    else {
        for (m in unique(x$marker)) {
            a1 <- x[x$marker == m, "allele1"]
            a2 <- x[x$marker == m, "allele2"]
            adb <- db[db$marker == m, "allele"]
            for (a in list(a1, a2)) {
                w <- !a %in% adb
                if (any(w)) {
                  if (first) {
                    message("observed alleles not in database\n")
                    first <- FALSE
                  }
                  message(paste(m, ":", a[w], "\n"))
                  if (fix > 0) 
                    db <- merge(db, data.frame(marker = m, allele = a[w], 
                      frequency = fix), all = T)
                }
            }
            f <- db[db$marker == m, ]$frequency
            db[db$marker == m, ]$frequency <- f/sum(f)
        }
    }
    invisible(db)
}
dbetabinom <-
function (x, n, alpha, beta) 
{
    lx <- max(length(x), length(n), length(alpha), length(beta))
    x <- rep_len(x, length.out = lx)
    n <- rep_len(n, length.out = lx)
    alpha <- rep_len(alpha, length.out = lx)
    beta <- rep_len(beta, length.out = lx)
    suppressWarnings(res <- choose(n, x) * exp(lgamma(alpha + x) + lgamma(beta + 
        n - x) + lgamma(alpha + beta) - (lgamma(alpha + beta + 
        n) + lgamma(alpha) + lgamma(beta))))
    res[is.na(res)] <- (n == x)[is.na(res)]
    res
}
delete.DQnodes <-
function (mixture, which = "DQ") 
{
    for (m in mixture$markers) {
        d <- mixture$domains[[m]]
        z <- get.nodes(d)
        if (!is.na(pmatch("D", unlist(strsplit(which, ""))))) 
            for (n in z[substring(z, 1, 1) == "D"]) delete.node(d, 
                n)
        if (!is.na(pmatch("Q", unlist(strsplit(which, ""))))) 
            for (n in z[substring(z, 1, 1) == "Q"]) delete.node(d, 
                n)
    }
}
extract <-
function (db, epg, gts, C, markers = sort(unique(db$marker))) 
{
    z <- list()
    for (m in markers) if (m != "AMEL") {
#        e <- subset(epg, marker == m)
        e <- epg[epg['marker']==m,]
        z[[m]]$R <- sort(e$allele[e$height > C])
        for (gt in gts) {
            q <- get(gt)
            q <- sort(unlist(q[q$marker == m, 2:3]))
            names(q) <- NULL
            z[[m]][[gt]] <- q
        }
 #       dbm <- subset(db, marker == m)
        dbm <- db[db['marker']==m,]
        z[[m]]$alleles <- dbm$allele
        z[[m]]$frequency <- dbm$freq
    }
    z
}
fuzz <-
function (mixture, eps = 1e-08) 
{
    for (m in mixture$markers) {
        d <- mixture$domains[[m]]
        z <- get.nodes(d)
        for (n in z[substring(z, 1, 1) == "n"]) {
            tab <- get.table(d, n)
            tab$Freq <- pmax(eps, tab$Freq)
            set.table(d, n, tab, type = "cpt")
        }
    }
}
gt2aca <-
function (mixture, gt, eps = 0) 
{
    logGt <- 0
    aca <- list()
    namesaca <- NULL
    for (m in 1:length(mixture$markers)) {
        md <- match(mixture$markers[m], as.character(gt$marker))
        if (!is.na(md)) {
            a1 <- match(gt$allele1[md], mixture$data[[m]]$allele)
            a2 <- match(gt$allele2[md], mixture$data[[m]]$allele)
            if (is.na(a1) || is.na(a2)) 
                warning("invalid allele in genotype")
            else {
                ac <- rep(0, length(mixture$data[[m]]$allele))
                ac[a1] <- ac[a1] + 1
                ac[a2] <- ac[a2] + 1
                if (eps > 0) {
                  AC <- matrix(eps, length(ac), 3)
                  AC[cbind(1:length(ac), ac + 1)] <- 1
                  aca <- c(aca, list(AC))
                }
                else {
                  aca <- c(aca, list(ac))
                }
                namesaca <- c(namesaca, mixture$markers[m])
                fq <- prod(mixture$data[[m]]$freq[c(a1, a2)])
                if (a1 != a2) 
                  fq <- 2 * fq
                logGt <- logGt + log(fq)
            }
        }
    }
    names(aca) <- namesaca
    attr(aca, "logGt") <- logGt
    aca
}
intoMix <-
function (res) 
{
#    res0 <- subset(res00, !is.na(DNA))
    res0 <- res[!is.na(res['DNA'])]
    res0$Sim <- match(res0$Sim, unique(res0$Sim))
    res0
}
loop.rpt.IBD <-
function (listdata, pars, IBD, typed.gts = NULL, inds = 1, jtyped = ncol(IBD$patt)/2 - 
    length(typed.gts) + seq_along(typed.gts), jcontr = seq_along(inds), 
    targets = NULL, contribs=NULL, quiet = FALSE, verbose = FALSE, 
    presence.only = FALSE, ...) 
{
#------------#
IBD <- convertIBD(IBD)

if(any(!(contribs%in%targets))) { 
	extra<-setdiff(contribs, targets)
	typed<-names(typed.gts)
	keep<-intersect(extra,typed)
	drop<-setdiff(extra,typed)
if(length(keep)>0) 
	{
	cat('warning, in future put',keep,'in reference.profiles\n',
	'here, we extend IBD and targets to include it/them\n')
	targets<-c(targets,keep)
	pattkeep<-matrix(max(IBD$patt)+(1:(2*length(keep))),nrow(IBD$patt),2*length(keep),byrow=TRUE)
	IBD$patt<-cbind(IBD$patt,pattkeep)
	attr(IBD,'targets')<-targets
	}
if(length(drop)>0) cat('warning, in future omit',drop,':',
	'here, just ignored\n')
contribs<-setdiff(contribs,drop) 
}

    if (!(is.null(targets)&&is.null(contribs))) {
        jcontr <- match(contribs, targets)
	  inds<-which(!is.na(jcontr))
        jcontr <- jcontr[!is.na(jcontr)]
        jtyped <- match(names(typed.gts), targets)
    }
    if (length(jtyped) != length(typed.gts)) 
        stop("jtyped and typed.gts incompatible\n")
    if (length(jcontr) != length(inds)) 
        stop("jcontr and inds incompatible\n")
    if (!all(jtyped %in% (1:ncol(IBD$patt)/2)) || !all(jcontr %in% 
        (1:ncol(IBD$patt)/2))) 
        stop("jtyped, jcontr and IBD incompatible;\n",jtyped,"\n",jcontr,"\n",ncol(IBD$patt)/2,"\n")
#------------#
    sumlogLR <- 0
    logLR<-NULL
    for (t in 1:length(listdata)) {
        data <- listdata[[t]]
        for (m in setdiff(unique(data$marker), "AMEL")) if (all(unlist(lapply(typed.gts, 
            function(x) {
                m %in% x$marker
            })))) {
            listdatam <- lapply(listdata, function(d) subset(d, 
                d$marker == m))
            mixDm <- DNAmixture(listdatam, ...)
            baseline <- protected(logL(mixDm, presence.only)(pars))
            LR <- 0
            wtsum <- 0
            for (r in 1:nrow(IBD$patt)) {
                mixDmr <- DNAmixture(listdatam, triangulate = F, 
                  compile = F, ...)
                wt <- IBD$pr[r] * rpt.IBD(mixDmr, list(patt = IBD$patt[r, 
                  ]), typed.gts, inds, jtyped, jcontr)
                if (length(wt) > 0 && wt > 0) {
                  log10LR <- (protected(logL(mixDmr, presence.only)(pars)) - 
                    baseline)/log(10)
cat(r,m,protected(logL(mixDmr, presence.only)(pars)),'\n')
                  LR <- LR + wt * 10^log10LR
                  wtsum <- wtsum + wt
                }
            }
            LR <- LR/wtsum
            if (verbose) 
                message(paste(m, "log10 LR", log10(LR), "; LR", 
                  LR, "\n"))
		logLR<-c(logLR,log10(LR))
            sumlogLR <- sumlogLR + log10(LR)
        }
    }
    if (verbose) 
        message(paste("all trace all marker overall log10 LR", sumlogLR, 
            "LR", 10^sumlogLR, "\n"))
    names(sumlogLR)<-'sumlog10LR'
    invisible(structure(sumlogLR,log10LR=logLR))
}
make.profile <-
function (gt, name = "K") 
{
    marker <- NULL
    allele <- NULL
    count <- NULL
    for (i in 1:nrow(gt)) {
        a1 <- gt$allele1[i]
        a2 <- gt$allele2[i]
        if (a1 == a2) {
            marker <- c(marker, gt$marker[i])
            allele <- c(allele, a1)
            count <- c(count, 2)
        }
        else {
            marker <- c(marker, gt$marker[i], gt$marker[i])
            allele <- c(allele, a1, a2)
            count <- c(count, 1, 1)
        }
    }
    out <- data.frame(marker = marker, allele = allele, count = count)
    names(out)[3] <- name
    out
}
new.node<-function(domain,patt)
{
# returns string of form paste0(patt,i) that is not already the name of a node in domain
z<-grep(patt,get.nodes(domain),value=TRUE)
for(i in c(1:9,sample(10:99,10))) if(!(i%in%substring(z,nchar(patt)+1))) return(paste0(patt,i))
stop('too many attempts in new.node')
}
plot.IBD<-function (x, labels = NULL, probs = NULL, order = NULL, 
    colrs = c("black", "red", "blue"), digits = 3, nr = ceiling(sqrt(np)), ...) 
{
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
    x<-as.IBD(x)
    if(is.null(probs)) probs<-x$pr
    x<-x$patt
    np <- nrow(x)
    if (!is.null(labels) && is.na(labels)) 
        labels <- apply(x, 1, function(z) paste(z, collapse = " "))
    chcol <- -1
    if (!is.null(order)) {
        if (length(order) > 1) {
            o <- order(order)
            chcol <- which(labels[o][-1] != labels[o][-np])
        }
        else {
            switch(pmatch(order, c("pattern", "probs", "labels")), 
                {
                  nout <- apply(x, 1, function(p) sum(0 != 
                    apply(matrix(p, nrow = 2), 2, diff)))
                  ndist <- apply(x, 1, max)
                  maxf <- apply(x, 1, function(p) max(tabulate(p, 
                    6)))
                  z <- apply(x, 1, function(p) sum(apply(matrix(apply((outer(p, 
                    p, "-") + diag(6)) == 0, 1, sum), 2, 3), 
                    2, min)))
                  o <- statnet.common::order(cbind(nout, ndist, 
                    maxf, z))
                  zz <- cbind(nout, ndist, maxf, z)[o, ]
                  chcol <- which(apply(apply(zz, 2, diff) != 
                    0, 1, any))
                }, {
                  o <- order(probs, decreasing = TRUE)
                  chcol <- which(diff(probs[o]) != 0)
                }, {
                  if (is.null(labels)) o <- 1:np else {
                    o <- order(labels)
                    chcol <- which(labels[o][-1] != labels[o][-np])
                  }
                })
            x <- x[o, ]
            if (!is.null(labels)) 
                labels <- labels[o]
            if (!is.null(probs)) 
                probs <- probs[o]
        }
    }
    nc <- ceiling(np/nr)
    opars <- par(mfrow = c(nr, nc), mar = c(0.5, 0.5, 0.5, 0.5))
    nind <- ncol(x)/2
    y <- rep((nind - 1):0, each = 2)
    xp <- rep(c(0, 1), nind)
    curve <- spline(seq(0, 1, len = 5), c(0, 2/3, 1, 2/3, 0))
    sc <- c(0, 0.15, 0.4)
    icol <- length(colrs)
    for (f in 1:np) {
        patt <- x[f, ]
        plot(xp, y, xlim = c(-0.5, 1.5), ylim = c(-0.9, nind - 
            0.5), type = "n", xaxt = "n", yaxt = "n", bty = "n", 
            asp = 1)
        xo <- yo <- 0
        for (i in 1:(2 * nind - 1)) for (j in (i + 1):(2 * nind)) if (patt[i] == 
            patt[j]) {
            if (i%%2 == 1 && j%%2 == 1 & abs(i - j) > 2) 
                lines(xo + xp[j] - sc[abs(i - j)/2] * curve$y, 
                  yo + y[j] + (y[i] - y[j]) * curve$x, col = colrs[icol])
            else if (i%%2 == 0 && j%%2 == 0 & abs(i - j) > 2) 
                lines(xo + xp[j] + sc[abs(i - j)/2] * curve$y, 
                  yo + y[j] + (y[i] - y[j]) * curve$x, col = colrs[icol])
            else lines(xo + xp[c(i, j)], yo + y[c(i, j)], col = colrs[icol])
        }
        points(xo + xp, yo + y, pch = 16)
        yt <- 0.2
        if (!is.null(labels)) {
            yt <- yt - 0.5
            text(0.5, yt, labels[f], cex = 1)
        }
        if (!is.null(probs)) {
            yt <- yt - 0.5
            text(0.5, yt, signif(probs[f], digits), cex = 1)
        }
        if (f %in% chcol) 
            icol <- (icol%%length(colrs)) + 1
    }
    par(opars)
}
print.IBD <-
function(x,...)
{
nr<-length(x$pr)
nc<-length(x$patt)/nr
res<-cbind(x$pr,matrix(x$patt,nr,nc))
if(!is.null(attr(x,'targets'))) targets<-attr(x,'targets')
prmatrix(res,rowlab=rep('',nr),collab=c('pr',as.vector(rbind(targets,rep('',nc/2)))))
}
print.IBD <-
function(x,targets=letters[1:(nc/2)],...)
{
IBD<-convertIBD(x)
nr<-length(IBD$pr)
nc<-length(IBD$patt)/nr
res<-cbind(IBD$pr,matrix(IBD$patt,nr,nc))
if(!is.null(attr(IBD,'targets'))) targets<-attr(IBD,'targets')
prmatrix(cbind(IBD$pr,IBD$patt),rowlab=rep('',nr),collab=c('pr',as.vector(rbind(targets,rep('',nc/2)))))
}
print.tablesize <-
function (x, ...) 
{
    cat("total table size", x, "\n")
    invisible(x)
}
process.patterns <-
function (IBD, igt, jtyped, jcontr, q, keeplabels) 
{
    ntyped <- length(jtyped)
    nperms <- 2^ntyped
    apa <- NULL
    ape <- NULL
    pdenom <- NULL
    type <- NULL
    val <- NULL
    npp <- 0
    rlen <- 2 * length(jcontr)
    for (ipatt in 1:length(IBD$pr)) {
        kt <- IBD$patt[ipatt, as.vector(outer(c(-1, 0), 2 * 
            jtyped, "+"))]
        ukt <- unique(kt)
        ap <- rep(0, max(IBD$patt))
        for (iperm in 1:nperms) {
            ipp <- nperms * (ipatt - 1) + iperm
            possible <- TRUE
            if (ntyped != 0) {
                w <- binary(iperm - 1, ntyped)
                pkt <- kt[as.vector(rbind(w, 1 - w)) + rep(2 * 
                  (1:ntyped) - 1, rep(2, ntyped))]
                ap[pkt] <- igt
                possible <- all(ap[pkt] == igt)
            }
            if (possible) {
                apa <- c(apa, ipatt)
                ape <- c(ape, iperm)
                type <- rbind(type, rep(" ", rlen))
                val <- rbind(val, rep(0, rlen))
                if (is.null(dim(q))) 
                  pqa <- prod(q[ap[ukt]])
                else pqa <- prod(q[cbind(ukt, ap[ukt])])
                pdenom <- c(pdenom, IBD$pr[ipatt] * pqa)
                npp <- npp + 1
                for (ij in 1:(2 * length(jcontr))) {
                  j <- as.vector(outer(c(-1, 0), 2 * jcontr, 
                    "+"))[ij]
                  ku <- IBD$patt[ipatt, j]
                  for (ik in 1:length(ku)) {
                    ia <- match(ku[ik], ukt)
                    if (is.na(ia)) {
                      type[npp, ij] <- "d"
                      val[npp, ij] <- ku[ik]
                    }
                    else {
                      type[npp, ij] <- "c"
                      val[npp, ij] <- ap[ukt[ia]]
                    }
                  }
                }
            }
        }
    }
    res <- NULL
    if (!is.null(val)) {
        oval <- val
        if (keeplabels) {
            todraw <- NULL
            for (ir in 1:nrow(val)) {
                w <- which(type[ir, ] == "d")
                vw <- val[ir, w]
                todraw <- c(todraw, vw)
            }
            todraw <- unique(todraw)
        }
        else {
            ntodraw <- 0
            for (ir in 1:nrow(val)) {
                w <- which(type[ir, ] == "d")
                vw <- val[ir, w]
                uvw <- unique(vw)
                oval[ir, w] <- match(vw, uvw)
                ntodraw <- max(ntodraw, match(vw, uvw))
            }
            if (ntodraw > 0) 
                todraw <- 1:ntodraw
            else todraw <- NULL
        }
        res <- data.frame(apa, ape, pdenom, type = type, oval = oval, 
            stringsAsFactors = FALSE)
        attr(res, "todraw") <- todraw
        attr(res, "npp") <- npp
    }
    res
}
protected <-
function (x, default = -Inf) 
{
    if (is.numeric(tryCatch.W.E(x)$val)) 
        x
    else default
}
protected.mixML <-
function (mixture, pars, constraints = NULL, phi.eq = FALSE, 
    val = NULL, trace = FALSE, order.unknowns = TRUE, default = -999999, 
    ...) 
{
    R <- mixture$ntraces
    k <- mixture$k
    U <- mixture$U
    K <- mixture$K
    contr <- c(U, K)
    n.unknown <- mixture$n.unknown
    x2phi <- function(x) {
        phi <- if (phi.eq) {
            rep(list(tail(x, -3 * R)), R)
        }
        else {
            split(tail(x, -3 * R), rep(1:R, each = k))
        }
        mapply(function(th, contributor) {
            names(th) <- contributor
            th
        }, phi, rep(list(contr), times = R), SIMPLIFY = FALSE)
    }
    x2phiU <- function(x) {
        lapply(x2phi(x), function(th) {
            head(th, n.unknown)
        })
    }
    x2arr <- function(x) {
        arr <- array(list(NULL), dim = c(R, 4), dimnames = list(NULL, 
            c("rho", "eta", "xi", "phi")))
        arr[1:(3 * R)] <- as.list(head(x, 3 * R))
        arr[, "phi"] <- x2phi(x)
        arr
    }
    parlist2x <- function(parlist) {
        rex <- unlist(parlist[, 1:3], use.names = FALSE)
        if (phi.eq) {
            phi <- parlist[[1, 4]][contr]
        }
        else {
            phi <- unlist(lapply(parlist[, 4], function(x) x[contr]), 
                use.names = FALSE)
        }
        c(rex, phi)
    }
    logl <- logL(mixture)
    funvals <- numeric(0)
    minus.loglikelihood <- function(x) {
        xs <- x2arr(x)
        if (trace) 
            message(paste("xs",xs))
        val <- -protected(logl(xs), default)
        if (trace) 
            message(paste("-val",-val))
        val
    }
    lb <- rep(0, times = 3 * R + ifelse(phi.eq, k, k * R))
    ub <- rep(c(Inf, 1), times = c(2 * R, R + ifelse(phi.eq, 
        k, k * R)))
    if (phi.eq) {
        phi.sum.constraint <- function(x) {
            sum(tail(x, -3 * R))
        }
        eqB <- 1
    }
    else {
        phi.sum.constraint <- function(x) {
            sapply(x2phi(x), sum)
        }
        eqB <- rep(1, R)
    }
    if (!missing(constraints)) {
        eqfun <- function(x) {
            c(phi.sum.constraint(x), do.call(constraints, list(x2arr(x))))
        }
        eqB <- c(eqB, val)
    }
    else {
        eqfun <- phi.sum.constraint
    }
    if (n.unknown > 1 & order.unknowns) {
        phi.symmetry <- function(x) {
            phiU <- x2phiU(x)[[1]]
            diff(phiU, lag = 1)
        }
        n.diffs <- (n.unknown - 1)
        ineqLB <- rep(-1, n.diffs)
        ineqUB <- rep(0, n.diffs)
    }
    else phi.symmetry <- ineqLB <- ineqUB <- NULL
    x0 <- parlist2x(pars)
    soln <- solnp(x0, fun = minus.loglikelihood, LB = lb, UB = ub, 
        eqfun = eqfun, eqB = eqB, ineqfun = phi.symmetry, ineqLB = ineqLB, 
        ineqUB = ineqUB, control = list(trace = 0, delta = 1e-07, 
            tol = 1e-08), ...)
    est <- x2arr(soln$pars)
    class(est) <- "mixpar"
    val <- -tail(soln$value, 1)
    out <- list(mle = est, lik = val, funvals = funvals, starting.point = x0, 
        minimization.output = soln)
}
prune <-
function (x, targets, S, verbose = FALSE) 
{
    id <- x$ID
    ifid <- x$FIDX
    imid <- x$MIDX
    T <- NULL
    stack <- targets
    repeat {
        now <- stack[1]
        stack <- stack[-1]
        if (!now %in% S) {
            inow <- match(now, id)
            f <- ifid[inow]
            if (f != 0) {
                cf <- id[f]
                if (cf %in% S) 
                  T <- unique(c(T, cf))
                else if (!cf %in% stack) 
                  stack <- c(stack, cf)
            }
            m <- imid[inow]
            if (m != 0) {
                cm <- id[m]
                if (cm %in% S) 
                  T <- unique(c(T, cm))
                else if (!cm %in% stack) 
                  stack <- c(stack, cm)
            }
        }
        if (length(stack) == 0) 
            break
    }
    if (verbose) 
        message(paste("can delete", S[!S %in% T], "\nmust keep ", S[S %in% 
            T], "\n"))
    S %in% T
}
read.db <-
function (file = "db.csv", sep = ifelse(Sys.getenv("USERNAME") == 
    "Julia", ";", ",")) 
{
    z <- read.csv(file, stringsAsFactors = FALSE, sep = sep)
    data.frame(marker = as.factor(z$marker), allele = as.numeric(z$allele), 
        height = z$frequency)
}
read.epg <-
function (file = "epg.csv", sep = ifelse(Sys.getenv("USERNAME") == 
    "Julia", ";", ",")) 
{
    z <- read.csv(file, stringsAsFactors = FALSE, sep = sep)
    data.frame(marker = as.factor(z$marker), allele = as.numeric(z$allele), 
        height = z$height)
}
replace.tables.for.UAF <- function (mixture, M, compile = TRUE) rpt.UAF(mixture, M, compile) 

replace.Ui.tables <-
function (mixture, aca, ind = 1) 
{
    if (is.data.frame(aca)) 
        aca <- gt2aca(mixture, aca)
    for (m in names(mixture$domains)) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            nC <- aca[[m]]
            na <- length(nC)
            if (max(aca[[m]]) == 2) {
                am <- which(aca[[m]] == max(aca[[m]]))
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a] - (am == 1), 
                      1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a] - (am == a), 
                      pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                        0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, tab[, S2am1]), 1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
            else {
                ab <- which(aca[[m]] == max(aca[[m]]))
                am <- ab[1]
                bm <- ab[2]
                add.node(d, "ibdyet", subtype = "boolean")
                tab <- get.table(d, "ibdyet")
                tab$Freq <- c(0.5, 0.5)
                set.table(d, "ibdyet", tab, type = "cpt")
                add.edge(d, paste("n_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("n_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == am) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                          0, s, 1))
                    }
                  }
                  else if (a == bm) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], 1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], pmax(0, 1 - tab[, S2am1]), 
                        q[a]/ifelse(s > 0, s, 1))
                    }
                  }
                  else if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a], 1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a], pmax(0, 1 - 
                      tab[, S2am1]), q[a]/ifelse(s > 0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                  if (a == bm) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
        }
        compile(d)
    }
}
require.compiled <-
function (mixture) 
{
    for (m in mixture$markers) {
        d <- mixture$dom[[m]]
        if (!unclass(summary(d))$domain$compiled) {
            compile(d)
            message(paste(m, "compiled\n"))
        }
        else message(paste(m, "already compiled\n"))
    }
}
require.edge <-
function (d, n, p) 
{
    if (!p %in% get.parents(d, n)) 
        add.edge(d, n, p)
}
rGTs <-
function (nreps, IBD, db, DNA, sex = rep(0, ncontr), nU = 0) 
{
    IBD <- convertIBD(IBD)
    if (nU > 0) 
        IBD <- addU(IBD, nU)
    ncontr <- ncol(IBD$patt)/2
    res <- list()
    markers <- as.character(unique(db$marker))
    nmark <- length(markers)
    z <- matrix(0, nmark, ncol(IBD$patt))
    dimnames(z) <- list(markers, NULL)
    for (ir in 1:nreps) {
        for (m in markers) {
            dbm <- db[db$marker == m, ]
            zi <- IBD$patt[sample(length(IBD$pr), 1, replace = TRUE, 
                prob = IBD$pr), ]
            g <- sample(dbm$allele, max(IBD$patt), replace = TRUE, 
                prob = dbm$frequency)
            z[m, ] <- g[zi]
        }
        Sim <- rep(1:ncontr, each = 2 * nmark)
        Marker <- rep(markers, each = 2, times = ncontr)
        Allele <- as.character(as.vector(aperm(array(z, c(nmark, 
            2, ncontr)), c(2, 1, 3))))
        w <- Marker == "AMEL"
        if (!is.null(w)) {
            amel <- rep("X", 2 * ncontr)
            for (ic in 1:ncontr) amel[2 * ic] <- switch(sex[ic] + 
                1, sample(c("X", "Y"), 1), "Y", "X")
            Allele[w] <- amel
        }
        sn <- rep("rGTs", 2 * ncontr * nmark)
        res[[ir]] <- data.frame(Sim, Sample.Name = sn, Marker, 
            Allele, DNA = DNA[Sim], stringsAsFactors = FALSE)
    }
    res
}
rni <-
function (seed = 0) 
{
    if (seed == 0) {
        tt <- Sys.time()
        st <- strftime(tt, "%s")
        t <- as.numeric(paste(substring(st, nchar(st) - 5), substring(as.character(tt), 
            21, 23), sep = ""))
        set.seed(t)
        seed <- sample(999999999, 1)
        message(paste("to rerun type rni(", seed, ")\n", sep = ""))
    }
    set.seed(seed)
    invisible(seed)
}
rpt.AMEL <-
function (mixture, sex, compile = TRUE) 
{
    if ("AMEL" %in% mixture$markers) {
        d <- mixture$domains[["AMEL"]]
        for (i in seq_along(sex)) {
            ni1 <- paste("n", i, 1, sep = "_")
            if (ni1 %in% get.nodes(d)) {
                tab <- get.table(d, ni1)
                freq <- switch(sex[i] + 1, c(0, 0.5, 0.5), c(0, 
                  1, 0), c(0, 0, 1))
                tab$Freq <- freq
                set.table(d, ni1, tab, type = "cpt")
            }
        }
    }
    if (compile) 
        for (d in mixture$domains) compile(d)
}
rpt.IBD <-
function (mixture, IBD='parent-child', typed.gts = NULL, inds = 1, jtyped = ncol(IBD$patt)/2 - 
    length(typed.gts) + seq_along(typed.gts), jcontr = seq_along(inds), 
    targets = NULL, contribs = NULL, quiet = FALSE, all.freq = NULL, 
    compile = TRUE) 
{
#------------#
IBD <- convertIBD(IBD)

if(any(!(contribs%in%targets))) { 
	extra<-setdiff(contribs, targets)
	typed<-names(typed.gts)
	keep<-intersect(extra,typed)
	drop<-setdiff(extra,typed)
if(length(keep)>0) 
	{
	cat('warning, in future put',keep,'in reference.profiles\n',
	'here, we extend IBD and targets to include it/them\n')
	targets<-c(targets,keep)
	pattkeep<-matrix(max(IBD$patt)+(1:(2*length(keep))),nrow(IBD$patt),2*length(keep),byrow=TRUE)
	IBD$patt<-cbind(IBD$patt,pattkeep)
	attr(IBD,'targets')<-targets
	}
if(length(drop)>0) cat('warning, in future omit',drop,':',
	'here, just ignored\n')
contribs<-setdiff(contribs,drop) 
}

    if (!(is.null(targets)&&is.null(contribs))) {
        jcontr <- match(contribs, targets)
	  inds<-which(!is.na(jcontr))
        jcontr <- jcontr[!is.na(jcontr)]
        jtyped <- match(names(typed.gts), targets)
    }
    if (length(jtyped) != length(typed.gts)) 
        stop("jtyped and typed.gts incompatible\n")
    if (length(jcontr) != length(inds)) 
        stop("jcontr and inds incompatible\n")
    if (!all(jtyped %in% (1:ncol(IBD$patt)/2)) || !all(jcontr %in% 
        (1:ncol(IBD$patt)/2))) 
        stop("jtyped, jcontr and IBD incompatible;\n",jtyped,"\n",jcontr,"\n",ncol(IBD$patt)/2,"\n")
#------------#

    ntyped <- length(jtyped)
    ap <- rep(0, max(IBD$patt))
    ptypedgts <- NULL
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (all(unlist(lapply(typed.gts, function(x) {
            m %in% x$marker
        })))) {
            keeplabels <- FALSE
            if (missing(all.freq) || is.null(all.freq)) {
                q <- mixture$data[[m]]$freq
                q <- q/sum(q)
            }
            else {
                if (is.data.frame(all.freq)) {
#                  dbm <- subset(all.freq, marker == m) 
                  dbm <- all.freq[all.freq['marker']==m,]
                  q <- dbm$frequency[match(mixture$data[[m]]$allele, 
                    dbm$allele)]
                  q <- q/sum(q)
                }
                else {
                  keeplabels <- TRUE
                  q <- matrix(0, max(IBD$patt), nrow(mixture$data[[m]]))
                  for (k in 1:max(IBD$patt)) {
 #                   dbm <- subset(all.freq[[k]], marker == m)
                    dbm <- all.freq[all.freq[[k]]['marker']==m,]
                    q[k, ] <- dbm$frequency[match(mixture$data[[m]]$allele, 
                      dbm$allele)]
                    q[k, ] <- q[k, ]/sum(q[k, ])
                  }
                }
            }
            igt <- NULL
            if (ntyped > 0) 
                for (l in 1:ntyped) {
                  gt <- typed.gts[[l]]
                  a1 <- match(gt$allele1[gt$marker == m], mixture$data[[m]]$allele)
                  a2 <- match(gt$allele2[gt$marker == m], mixture$data[[m]]$allele)
                  igt <- c(igt, a1, a2)
                }
            z <- process.patterns(IBD, igt, jtyped, jcontr, q, 
                keeplabels)
            if (is.null(z)) {
                ptypedgts <- c(ptypedgts, 0)
                names(ptypedgts)[length(ptypedgts)] <- m
            }
            else {
                nodes <- get.nodes(d)
                for (i in inds) {
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("n", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                  }
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("S", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                    delete.node(d, n)
                  }
                }
		    pattperm<-new.node(d,"pp_")
                add.node(d, pattperm, states = 1:attr(z, 
                  "npp"))
                tab <- get.table(d, pattperm)
                f <- as.vector(z$pdenom)
                ptypedgts <- c(ptypedgts, sum(f))
                names(ptypedgts)[length(ptypedgts)] <- m
                tab$Freq <- f/sum(f)
                set.table(d, pattperm, tab)
                na <- nrow(mixture$data[[m]])
# count how many nodes "np*_1" already exist
npn<-grep("np",get.nodes(d),value=TRUE)
koff<-length(npn)
if(koff>0) koff<-sum(sapply(strsplit(npn,'_'),function(x) {x[2]=="1"}))
                for (k in attr(z, "todraw")) {
                  if (is.null(dim(q))) {
                    qstar <- q/rev(cumsum(rev(q)))
                  }
                  else {
                    qstar <- q[k, ]/rev(cumsum(rev(q[k, ])))
                  }
                  for (a in 1:na) {
                    npka <- paste("np", koff+k, "_", a, 
                      sep = "")
                    add.node(d, npka, states = 0:1)
                    if (a == 1) {
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[1], qstar[1])
                    }
                    else {
                      Spka1 <- ifelse(a == 2, paste("np", 
                        koff+k, "_", 1, sep = ""), paste("Sp", 
                        koff+k, "_", a - 1, sep = ""))
                      add.edge(d, npka, Spka1)
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[a], qstar[a], 1, 
                        0)
                    }
                    set.table(d, npka, tab)
                    if (a != 1 & a != na) {
                      Spka <- paste("Sp", koff+k, "_", 
                        a, sep = "")
                      Spka1 <- ifelse(a == 2, paste("np", 
                        koff+k, "_", 1, sep = ""), paste("Sp", 
                        koff+k, "_", a - 1, sep = ""))
                      add.node(d, Spka, states = 0:1)
                      add.edge(d, Spka, npka)
                      add.edge(d, Spka, Spka1)
                      tab <- get.table(d, Spka)
                      tab$Freq <- as.numeric(tab[, Spka] == pmax(tab[, 
                        npka], tab[, Spka1]))
                      set.table(d, Spka, tab)
                    }
                  }
                }
                for (l in 1:length(inds)) {
                  i <- inds[l]
                  j <- l
                  t <- cbind(z[[paste("type.", 2 * j - 
                    1, sep = "")]], z[[paste("type.", 
                    2 * j, sep = "")]])
                  npj <- cbind(paste("np", koff+z[[paste("oval.", 
                    2 * j - 1, sep = "")]], sep = ""), 
                    paste("np", koff+z[[paste("oval.", 
                      2 * j, sep = "")]], sep = ""))
                  ov <- cbind(z[[paste("oval.", 2 * j - 
                    1, sep = "")]], z[[paste("oval.", 
                    2 * j, sep = "")]])
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    add.edge(d, nia, pattperm)
                    for (ipp in 1:attr(z, "npp")) {
                      if (t[ipp, 1] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 1], 
                          a, sep = "_"))
                      if (t[ipp, 2] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 2], 
                          a, sep = "_"))
                    }
                  }
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    tab <- get.table(d, nia)
                    for (ipp in 1:attr(z, "npp")) {
                      w <- tab[, pattperm] == ipp
                      value <- rep(0, sum(w))
                      if (t[ipp, 1] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          1], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 1])
                      if (t[ipp, 2] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          2], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 2])
                      tab$Freq[w] <- as.numeric(tab[w, nia] == 
                        value)
                    }
                    set.table(d, nia, tab)
                  }
                }
            }
        }
        if (compile) {
            compile(d)
        }
    }
    if ("AMEL" %in% mixture$markers) {
        ptypedgts <- c(ptypedgts, 1)
        names(ptypedgts)[length(ptypedgts)] <- "AMEL"
        if (compile) {
            compile(mixture$domains[["AMEL"]])
        }
    }
    invisible(ptypedgts)
}
rpt.typed.child <-
function (mixture, aca, ind = 1) 
{
    if (is.data.frame(aca)) 
        aca <- gt2aca(mixture, aca)
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            nC <- aca[[m]]
            na <- length(nC)
            if (max(aca[[m]]) == 2) {
                am <- which(aca[[m]] == max(aca[[m]]))
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a] - (am == 1), 
                      1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a] - (am == a), 
                      pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                        0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1, tab[, S2am1]), 1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
            else {
                ab <- which(aca[[m]] == max(aca[[m]]))
                am <- ab[1]
                bm <- ab[2]
                add.node(d, "ibdyet", subtype = "boolean")
                tab <- get.table(d, "ibdyet")
                tab$Freq <- c(0.5, 0.5)
                set.table(d, "ibdyet", tab, type = "cpt")
                add.edge(d, paste("n_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", am, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("n_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                add.edge(d, paste("S_", ind, "_", bm, sep = ""), 
                  "ibdyet")
                for (a in 1:na) {
                  n2a <- paste("n_", ind, "_", a, sep = "")
                  S2a <- paste("S_", ind, "_", a, sep = "")
                  S2am1 <- paste("S_", ind, "_", a - 1, sep = "")
                  tab <- get.table(d, n2a)
                  if (a == am) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - tab[, "ibdyet"], 
                        pmax(0, 1 - tab[, S2am1]), q[a]/ifelse(s > 
                          0, s, 1))
                    }
                  }
                  else if (a == bm) {
                    if (a == 1) {
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], 1, q[1])
                    }
                    else {
                      tab <- get.table(d, n2a)
                      s <- sum(tail(q, -(a - 1)))
                      tab$Freq <- dbinom(tab[, n2a] - 1 + tab[, 
                        "ibdyet"], pmax(0, 1 - tab[, S2am1]), 
                        q[a]/ifelse(s > 0, s, 1))
                    }
                  }
                  else if (a == 1) {
                    tab$Freq <- dbinom(tab[, n2a], 1, q[1])
                  }
                  else {
                    tab <- get.table(d, n2a)
                    s <- sum(tail(q, -(a - 1)))
                    tab$Freq <- dbinom(tab[, n2a], pmax(0, 1 - 
                      tab[, S2am1]), q[a]/ifelse(s > 0, s, 1))
                  }
                  set.table(d, n2a, tab, type = "cpt")
                  if (a == am) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                  if (a == bm) {
                    tab <- get.table(d, S2a)
                    if (a == 1) {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], 0), 1, 0)
                    }
                    else {
                      tab$Freq <- ifelse(tab[, S2a] == pmax(tab[, 
                        n2a] - 1 + tab[, "ibdyet"], tab[, S2am1]), 
                        1, 0)
                    }
                    set.table(d, S2a, tab, type = "cpt")
                  }
                }
            }
        }
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
rpt.typed.parents <-
function (mixture, Mgt, Fgt, ind = 1, compile = TRUE) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (m %in% Mgt$marker && m %in% Fgt$marker) {
            z <- get.nodes(d)
            for (n in z[substring(z, 1, 3 + floor(log10(ind))) == 
                paste("S", ind, sep = "_")]) delete.node(d, n)
            na <- nrow(mixture$data[[m]])
            ma1 <- match(Mgt$allele1[Mgt$marker == m], mixture$data[[m]]$allele)
            ma2 <- match(Mgt$allele2[Mgt$marker == m], mixture$data[[m]]$allele)
            fa1 <- match(Fgt$allele1[Fgt$marker == m], mixture$data[[m]]$allele)
            fa2 <- match(Fgt$allele2[Fgt$marker == m], mixture$data[[m]]$allele)
            add.node(d, "patt", states = 1:4)
            tab <- get.table(d, "patt")
            tab$Freq <- rep(0.25, 4)
            set.table(d, "patt", tab, type = "cpt")
            for (a in 1:na) {
                nia <- paste("n_", ind, "_", a, sep = "")
                add.edge(d, nia, "patt")
                tab <- get.table(d, nia)
                sum <- ifelse(tab[, "patt"] == 1, (a == ma1) + 
                  (a == fa1), 0) + ifelse(tab[, "patt"] == 2, 
                  (a == ma2) + (a == fa1), 0) + ifelse(tab[, 
                  "patt"] == 3, (a == ma1) + (a == fa2), 0) + 
                  ifelse(tab[, "patt"] == 4, (a == ma2) + (a == 
                    fa2), 0)
                tab$Freq <- ifelse((tab[, nia] == sum), 1, 0)
                set.table(d, nia, tab, type = "cpt")
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
rpt.typed.relative <-
function (mixture, Rgt, IBD = c(0.25, 0.5, 0.25), ind = 1, 
    compile = TRUE) 
{
    if (is.list(IBD)) {
        if (ncol(IBD$patt) != 4) 
            stop("can only specify 2-person pattern")
        if (any(apply(matrix(t(IBD$patt), nrow = 2), 2, diff) == 
            0)) 
            stop("this function does not handle autozygosity")
        kappas <- rep(0, 3)
        for (j in 1:nrow(IBD$patt)) {
            patt <- IBD$patt[j, ]
            k <- 1 + sum(outer(patt[1:2], patt[3:4], "=="))
            kappas[k] <- kappas[k] + IBD$pr[j]
        }
        IBD <- kappas
    }
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (m %in% Rgt$marker) {
            z <- get.nodes(d)
            for (n in z[substring(z, 1, 3 + floor(log10(ind))) == 
                paste("S", ind, sep = "_")]) delete.node(d, n)
            ra1 <- match(Rgt$allele1[Rgt$marker == m], mixture$data[[m]]$allele)
            ra2 <- match(Rgt$allele2[Rgt$marker == m], mixture$data[[m]]$allele)
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            qstar <- q/rev(cumsum(rev(q)))
            na <- length(q)
            add.node(d, "patt", states = 1:4)
            tab <- get.table(d, "patt")
            tab$Freq <- c(IBD[1], IBD[2]/2, IBD[2]/2, 
                IBD[3])
            set.table(d, "patt", tab, type = "cpt")
            for (a in 1:na) {
                nsia <- paste("ns_", ind, "_", a, sep = "")
                add.node(d, nsia, states = 0:2)
                add.edge(d, nsia, "patt")
                if (a != 1) {
                  Siam1 <- ifelse(a == 2, paste("ns_", ind, "_1", 
                    sep = ""), paste("S_", ind, "_", a - 1, sep = ""))
                  add.edge(d, nsia, Siam1)
                  tab <- get.table(d, nsia)
                  tsofar <- tab[, Siam1]
                }
                else {
                  tab <- get.table(d, nsia)
                  tsofar <- rep(0, nrow(tab))
                }
                trials <- pmax(0, c(2, 1, 1, 0)[tab[, "patt"]] - 
                  tsofar)
                tab$Freq <- dbinom(tab[, nsia], trials, qstar[a])
                set.table(d, nsia, tab, type = "cpt")
                if (a != 1 & a != na) {
                  Siam1 <- ifelse(a == 2, paste("ns_", ind, "_1", 
                    sep = ""), paste("S_", ind, "_", a - 1, sep = ""))
                  Sia <- paste("S_", ind, "_", a, sep = "")
                  add.node(d, Sia, states = 0:2)
                  add.edge(d, Sia, nsia)
                  add.edge(d, Sia, Siam1)
                  tab <- get.table(d, Sia)
                  sum <- tab[, nsia] + tab[, Siam1]
                  tab$Freq <- ifelse((tab[, Sia] == sum) | (sum > 
                    2), 1, 0)
                  set.table(d, Sia, tab, type = "cpt")
                }
                nia <- paste("n_", ind, "_", a, sep = "")
                add.edge(d, nia, nsia)
                add.edge(d, nia, "patt")
                tab <- get.table(d, nia)
                sum <- ifelse(tab[, "patt"] == 1, tab[, nsia], 
                  0) + ifelse(tab[, "patt"] == 2, tab[, nsia] + 
                  (a == ra1), 0) + ifelse(tab[, "patt"] == 3, 
                  tab[, nsia] + (a == ra2), 0) + ifelse(tab[, 
                  "patt"] == 4, tab[, nsia] + (a == ra1) + (a == 
                  ra2), 0)
                tab$Freq <- ifelse((tab[, nia] == sum) | (sum > 
                  2), 1, 0)
                set.table(d, nia, tab, type = "cpt")
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
rpt.typed.relatives <-
function (mixture, IBD="parent-child", typed.gts = NULL, inds = 1, jtyped = ncol(IBD$patt)/2 - 
    length(typed.gts) + seq_along(typed.gts), jcontr = seq_along(inds), 
    targets = NULL, contribs = NULL, quiet = FALSE, all.freq = NULL, 
    compile = TRUE) 
{
    if (!is.null(targets)) {
        jcontr <- match(contribs, targets)
	  inds<-which(!is.na(jcontr))
        jcontr <- jcontr[!is.na(jcontr)]
        jtyped <- match(names(typed.gts), targets)
    }
    IBD <- convertIBD(IBD)
    if (length(jtyped) != length(typed.gts)) 
        stop("jtyped and typed.gts incompatible")
    if (length(jcontr) != length(inds)) 
        stop("jcontr and inds incompatible")
    if (!all(jtyped %in% (1:ncol(IBD$patt)/2)) || !all(jcontr %in% 
        (1:ncol(IBD$patt)/2))) 
        stop("jtyped, jcontr and IBD incompatible")
    ntyped <- length(jtyped)
    ap <- rep(0, max(IBD$patt))
    ptypedgts <- NULL
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (all(unlist(lapply(typed.gts, function(x) {
            m %in% x$marker
        })))) {
            keeplabels <- FALSE
            if (missing(all.freq) || is.null(all.freq)) {
                q <- mixture$data[[m]]$freq
                q <- q/sum(q)
            }
            else {
                if (is.data.frame(all.freq)) {
#                  dbm <- subset(all.freq, marker == m)
                  dbm <- all.freq[all.freq['marker']==m,]
                  q <- dbm$frequency[match(mixture$data[[m]]$allele, 
                    dbm$allele)]
                  q <- q/sum(q)
                }
                else {
                  keeplabels <- TRUE
                  q <- matrix(0, max(IBD$patt), nrow(mixture$data[[m]]))
                  for (k in 1:max(IBD$patt)) {
#                    dbm <- subset(all.freq[[k]], marker == m)
                    dbm <- all.freq[all.freq[[k]]['marker']==m,]
                    q[k, ] <- dbm$frequency[match(mixture$data[[m]]$allele, 
                      dbm$allele)]
                    q[k, ] <- q[k, ]/sum(q[k, ])
                  }
                }
            }
            igt <- NULL
            if (ntyped > 0) 
                for (l in 1:ntyped) {
                  gt <- typed.gts[[l]]
                  a1 <- match(gt$allele1[gt$marker == m], mixture$data[[m]]$allele)
                  a2 <- match(gt$allele2[gt$marker == m], mixture$data[[m]]$allele)
                  igt <- c(igt, a1, a2)
                }
            z <- process.patterns(IBD, igt, jtyped, jcontr, q, 
                keeplabels)
            if (is.null(z)) {
                ptypedgts <- c(ptypedgts, 0)
                names(ptypedgts)[length(ptypedgts)] <- m
            }
            else {
                nodes <- get.nodes(d)
                for (i in inds) {
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("n", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                  }
                  for (n in nodes[substring(nodes, 1, 3 + floor(log10(i))) == 
                    paste("S", i, sep = "_")]) {
                    par <- get.parents(d, n)
                    for (p in par) delete.edge(d, n, p)
                    delete.node(d, n)
                  }
                }
                add.node(d, "pattperm", states = 1:attr(z, 
                  "npp"))
                tab <- get.table(d, "pattperm")
                f <- as.vector(z$pdenom)
                ptypedgts <- c(ptypedgts, sum(f))
                names(ptypedgts)[length(ptypedgts)] <- m
                tab$Freq <- f/sum(f)
                set.table(d, "pattperm", tab)
                na <- nrow(mixture$data[[m]])
                for (k in attr(z, "todraw")) {
                  if (is.null(dim(q))) {
                    qstar <- q/rev(cumsum(rev(q)))
                  }
                  else {
                    qstar <- q[k, ]/rev(cumsum(rev(q[k, ])))
                  }
                  for (a in 1:na) {
                    npka <- paste("np", k, "_", a, 
                      sep = "")
                    add.node(d, npka, states = 0:1)
                    if (a == 1) {
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[1], qstar[1])
                    }
                    else {
                      Spka1 <- ifelse(a == 2, paste("np", 
                        k, "_", 1, sep = ""), paste("Sp", 
                        k, "_", a - 1, sep = ""))
                      add.edge(d, npka, Spka1)
                      tab <- get.table(d, npka)
                      tab$Freq <- c(1 - qstar[a], qstar[a], 1, 
                        0)
                    }
                    set.table(d, npka, tab)
                    if (a != 1 & a != na) {
                      Spka <- paste("Sp", k, "_", 
                        a, sep = "")
                      Spka1 <- ifelse(a == 2, paste("np", 
                        k, "_", 1, sep = ""), paste("Sp", 
                        k, "_", a - 1, sep = ""))
                      add.node(d, Spka, states = 0:1)
                      add.edge(d, Spka, npka)
                      add.edge(d, Spka, Spka1)
                      tab <- get.table(d, Spka)
                      tab$Freq <- as.numeric(tab[, Spka] == pmax(tab[, 
                        npka], tab[, Spka1]))
                      set.table(d, Spka, tab)
                    }
                  }
                }
                for (l in 1:length(inds)) {
                  i <- inds[l]
                  j <- l
                  t <- cbind(z[[paste("type.", 2 * j - 
                    1, sep = "")]], z[[paste("type.", 
                    2 * j, sep = "")]])
                  npj <- cbind(paste("np", z[[paste("oval.", 
                    2 * j - 1, sep = "")]], sep = ""), 
                    paste("np", z[[paste("oval.", 
                      2 * j, sep = "")]], sep = ""))
                  ov <- cbind(z[[paste("oval.", 2 * j - 
                    1, sep = "")]], z[[paste("oval.", 
                    2 * j, sep = "")]])
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    add.edge(d, nia, "pattperm")
                    for (ipp in 1:attr(z, "npp")) {
                      if (t[ipp, 1] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 1], 
                          a, sep = "_"))
                      if (t[ipp, 2] == "d") 
                        require.edge(d, nia, paste(npj[ipp, 2], 
                          a, sep = "_"))
                    }
                  }
                  for (a in 1:na) {
                    nia <- paste("n", i, a, sep = "_")
                    tab <- get.table(d, nia)
                    for (ipp in 1:attr(z, "npp")) {
                      w <- tab[, "pattperm"] == ipp
                      value <- rep(0, sum(w))
                      if (t[ipp, 1] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          1], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 1])
                      if (t[ipp, 2] == "d") 
                        value <- value + tab[w, paste(npj[ipp, 
                          2], a, sep = "_")]
                      else value <- value + (a == ov[ipp, 2])
                      tab$Freq[w] <- as.numeric(tab[w, nia] == 
                        value)
                    }
                    set.table(d, nia, tab)
                  }
                }
            }
        }
        if (compile) {
            compile(d)
        }
    }
    if ("AMEL" %in% mixture$markers) {
        ptypedgts <- c(ptypedgts, 1)
        names(ptypedgts)[length(ptypedgts)] <- "AMEL"
        if (compile) {
            compile(mixture$domains[["AMEL"]])
        }
    }
    invisible(ptypedgts)
}
rpt.UAF <-
function (mixture, M, compile = TRUE) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        na <- nrow(mixture$data[[m]])
        q <- mixture$data[[m]]$freq
        q <- q/sum(q)
        alpha <- M * q
        beta <- c(rev(cumsum(rev(alpha)))[-1], 0)
        z <- get.nodes(d)
        n.unknown <- length(grep("^n_.*_1$", z))
        for (n in z[substring(z, 1, 1) == "n"]) for (par in get.parents(d, 
            n)) delete.edge(d, n, par)
        n <- outer(1:n.unknown, 1:na, function(i, a) paste("n", 
            i, a, sep = "_"))
        S <- outer(1:n.unknown, 1:na, function(i, a) paste("S", 
            i, a, sep = "_"))
        T <- outer(1:n.unknown, 1:na, function(i, a) paste("T", 
            i, a, sep = "_"))
        U <- outer(1:n.unknown, 1:na, function(i, a) paste("U", 
            i, a, sep = "_"))
        for (a in 1:na) {
            if (n.unknown > 1) 
                for (i in 1:(n.unknown - 1)) {
                  add.node(d, T[i, a], states = 0:(2 * i), subtype = "numbered")
                  add.node(d, U[i, a], states = 0:(2 * i), subtype = "numbered")
                }
            for (i in 1:n.unknown) {
                if (a == 1) {
                  if (i == 1) {
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], rep(2, 
                      lnval), rep(alpha[1], lnval), rep(beta[1], 
                      lnval))
                    set.table(d, n[i, a], tab)
                  }
                  else {
                    add.edge(d, n[i, a], T[i - 1, a])
                    add.edge(d, n[i, a], U[i - 1, a])
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], rep(2, 
                      lnval), alpha[a] + tab[, T[i - 1, a]], 
                      (beta[a] + pmax(0, tab[, U[i - 1, a]])))
                    set.table(d, n[i, a], tab)
                  }
                }
                else {
                  if (i == 1) {
                    add.edge(d, n[i, a], S[i, a - 1])
                    tab <- get.table(d, n[i, a])
                    lnval <- nrow(tab)
                    tab$Freq <- dbetabinom(tab[, n[i, a]], 2 - 
                      tab[, S[i, a - 1]], rep(alpha[a], lnval), 
                      rep(beta[a], lnval))
                    set.table(d, n[i, a], tab)
                  }
                  else {
                    add.edge(d, n[i, a], T[i - 1, a])
                    add.edge(d, n[i, a], U[i - 1, a])
                    add.edge(d, n[i, a], S[i, a - 1])
                    tab <- get.table(d, n[i, a])
                    tab$Freq <- dbetabinom(tab[, n[i, a]], 2 - 
                      tab[, S[i, a - 1]], alpha[a] + tab[, T[i - 
                      1, a]], (beta[a] + pmax(0, tab[, U[i - 
                      1, a]])))
                    set.table(d, n[i, a], tab)
                  }
                }
                if (i < n.unknown) {
                  if (i == 1) {
                    add.edge(d, T[i, a], n[i, a])
                    tab <- get.table(d, T[i, a])
                    tab$Freq <- ifelse(tab[, T[i, a]] == tab[, 
                      n[i, a]], 1, 0)
                    set.table(d, T[i, a], tab)
                  }
                  else {
                    add.edge(d, T[i, a], n[i, a])
                    add.edge(d, T[i, a], T[i - 1, a])
                    tab <- get.table(d, T[i, a])
                    sumvals <- tab[, n[i, a]] + tab[, T[i - 1, 
                      a]]
                    tab$Freq <- ifelse(sumvals <= 2 * i, ifelse(tab[, 
                      T[i, a]] == sumvals, 1, 0), 1/(1 + 2 * 
                      i))
                    set.table(d, T[i, a], tab)
                  }
                }
                if (i < n.unknown) {
                  if (i == 1) {
                    add.edge(d, U[i, a], S[i, a])
                    tab <- get.table(d, U[i, a])
                    tab$Freq <- ifelse(tab[, U[i, a]] == 2 - 
                      tab[, S[i, a]], 1, 0)
                    set.table(d, U[i, a], tab)
                  }
                  else {
                    add.edge(d, U[i, a], S[i, a])
                    add.edge(d, U[i, a], U[i - 1, a])
                    tab <- get.table(d, U[i, a])
                    sumvals <- 2 - tab[, S[i, a]] + tab[, U[i - 
                      1, a]]
                    tab$Freq <- ifelse(sumvals <= 2 * i, ifelse(tab[, 
                      U[i, a]] == sumvals, 1, 0), 1/(1 + 2 * 
                      i))
                    set.table(d, U[i, a], tab)
                  }
                }
            }
        }
        if (compile) 
            compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
size <-
function (mixture) 
{
    if((.Platform$OS.type=="windows")&&("package:RHugin"%in%search())) 
		res<-0 else res <- sum(unlist(lapply(mixture$dom, function(d) sizedomain(d)))) 
    class(res) <- "tablesize"
    res
}
sizedomain<-function(domain)
{
if(is(domain,'gRaven')) 
{
sum(sapply(summary(domain,jt=TRUE)$jt$cliques,
	function(x) prod(sapply(x,function(x) length(get.states(domain, x))))))
} else if(is(domain,'RHuginDomain')){
sum(unlist(lapply(summary(domain, jt = TRUE)$jt, function(j) j$size)))
} else NULL
}

tryCatch.W.E <-
function (expr) 
{
    W <- NULL
    w.handler <- function(w) {
        W <<- w
        invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e), 
        warning = w.handler), warning = W)
}
wlr <-
function (sep, Cgt, db, ind = 1, Mgt = NULL) 
{
    LR <- 1
    for (m in Cgt$marker) {
        sepm <- summary(sep)[[m]]
        pa1 <- sepm[[paste("U", ind, ".1", sep = "")]]
        pa2 <- sepm[[paste("U", ind, ".2", sep = "")]]
        Prob <- sepm$Prob
        if (is.null(Mgt)) {
            ca1 <- Cgt[Cgt$marker == m, ]$allele1
            ca2 <- Cgt[Cgt$marker == m, ]$allele2
            dbm <- db[db$marker == m, ]
            qca1 <- dbm$freq[dbm$allele == ca1]
            qca2 <- dbm$freq[dbm$allele == ca2]
            Den <- qca1 * qca2
            if (ca1 != ca2) 
                Den <- 2 * Den
            if (ca1 == ca2) {
                Num <- 0.5 * qca1 * ((pa1 == ca1) + (pa2 == ca1))
            }
            else {
                Num <- 0.5 * qca1 * ((pa1 == ca2) + (pa2 == ca2)) + 
                  0.5 * qca2 * ((pa1 == ca1) + (pa2 == ca1))
            }
            LRmfgt <- Num/Den
        }
        else {
            ca1 <- Cgt[Cgt$marker == m, ]$allele1
            ca2 <- Cgt[Cgt$marker == m, ]$allele2
            ma1 <- Mgt[Mgt$marker == m, ]$allele1
            ma2 <- Mgt[Mgt$marker == m, ]$allele2
            w <- match(ca1, c(ma1, ma2))
            if (is.na(w)) {
                w <- match(ca2, c(ma1, ma2))
                if (!is.na(w)) {
                  t <- ca1
                  ca1 <- ca2
                  ca2 <- t
                }
            }
            if (is.na(w)) {
                LRmfgt <- 0
            }
            else {
                if (w == 2) {
                  t <- ma1
                  ma1 <- ma2
                  ma2 <- t
                }
                dbm <- db[db$marker == m, ]
                qca1 <- dbm$freq[dbm$allele == ca1]
                qca2 <- dbm$freq[dbm$allele == ca2]
                Num <- 0.5 * ((pa1 == ca2) + (pa2 == ca2))
                Den <- qca2
                if (ca2 == ma2 && ca1 != ca2) {
                  Num <- Num + 0.5 * ((pa1 == ca1) + (pa2 == 
                    ca1))
                  Den <- Den + qca1
                }
                LRmfgt <- Num/Den
            }
        }
        LRm <- sum(LRmfgt * Prob)/sum(Prob)
        LR <- LR * LRm
    }
    LR
}
write.epg <-
function (res, file = "epg.csv", C = 0) 
{
#    epg <- subset(res[, c("Marker", "Allele", "Height")], Height > C)
    epg <- res[, c("Marker", "Allele", "Height")]
    epg<-epg[epg['Height']>C,]
    names(epg) <- tolower(names(epg))
    epg$allele[epg$allele == "X"] <- 0
    epg$allele[epg$allele == "Y"] <- 1
    epg$allele <- as.numeric(epg$allele)
    if (file != "") 
        write.csv(epg, file, row.names = FALSE)
    else epg
}
