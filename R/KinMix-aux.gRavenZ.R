# additional functions needed to implement ALN and MBN methods

add.child.meiosis.nodes <-
function (mixture, aca, ind = 1) 
{
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            nC <- aca[[m]]
            na <- length(nC)
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            for (a in 1:na) {
                Cpa <- paste("Cp", a, sep = "")
                Fna <- paste("n_", ind, "_", a, sep = "")
                ga <- paste("g", a, sep = "")
                gam1 <- paste("g", a - 1, sep = "")
                Cma <- paste("Cm", a, sep = "")
                CmSa <- paste("CmS", a, sep = "")
                CmSam1 <- ifelse(a == 2, "Cm1", paste("CmS", 
                  a - 1, sep = ""))
                Cna <- paste("Cn", a, sep = "")
                add.node(d, Cpa, states = 0:1, subtype = "numbered")
                add.edge(d, Cpa, Fna)
                if (a == 1) {
                  tab <- get.table(d, Cpa)
                  tab$Freq <- dbinom(tab[, Cpa], 1, tab[, Fna]/2)
                }
                else {
                  add.edge(d, Cpa, gam1)
                  tab <- get.table(d, Cpa)
                  probs <- (tab[, Fna] == 1) * tab[, gam1]/2 + 
                    (tab[, Fna] == 2)
                  tab$Freq <- dbinom(tab[, Cpa], 1, probs)
                }
                set.table(d, Cpa, tab, type = "cpt")
                if (a < na) {
                  add.node(d, ga, states = 0:2, subtype = "numbered")
                  if (a == 1) {
                    add.edge(d, ga, Cpa)
                    add.edge(d, ga, Fna)
                    tab <- get.table(d, ga)
                    tab$Freq <- ((tab[, ga] == 2) * (tab[, Cpa] == 
                      0) + (tab[, ga] == 0) * (tab[, Cpa] == 
                      1)) * (tab[, Fna] >= 1) + (tab[, ga] == 
                      1) * (tab[, Fna] == 0)
                  }
                  else {
                    add.edge(d, ga, gam1)
                    add.edge(d, ga, Cpa)
                    add.edge(d, ga, Fna)
                    tab <- get.table(d, ga)
                    tab$Freq <- ((tab[, ga] == 2) * (tab[, Cpa] == 
                      0) + (tab[, ga] == 0) * (tab[, Cpa] == 
                      1)) * (tab[, Fna] >= 1) * (tab[, gam1] == 
                      1) + (tab[, ga] == tab[, gam1]) * ((tab[, 
                      Fna] == 0) + (tab[, Fna] >= 1) * (tab[, 
                      gam1] != 1))
                  }
                  set.table(d, ga, tab, type = "cpt")
                }
                add.node(d, Cma, states = 0:1, subtype = "numbered")
                if (a == 1) {
                  tab <- get.table(d, Cma)
                  tab$Freq <- dbinom(tab[, Cma], 1, q[1])
                }
                else {
                  add.edge(d, Cma, CmSam1)
                  tab <- get.table(d, Cma)
                  s <- sum(tail(q, -(a - 1)))
                  tab$Freq <- dbinom(tab[, Cma], 1 - tab[, CmSam1], 
                    q[a]/ifelse(s > 0, s, 1))
                }
                set.table(d, Cma, tab, type = "cpt")
                if (a > 1 && a < na) {
                  add.node(d, CmSa, states = 0:1, subtype = "numbered")
                  add.edge(d, CmSa, Cma)
                  add.edge(d, CmSa, CmSam1)
                  tab <- get.table(d, CmSa)
                  tab$Freq <- ifelse(tab[, CmSa] == pmax(tab[, 
                    Cma], tab[, CmSam1]), 1, 0)
                  set.table(d, CmSa, tab, type = "cpt")
                }
                add.node(d, Cna, states = 0:2, subtype = "numbered")
                add.edge(d, Cna, Cpa)
                add.edge(d, Cna, Cma)
                tab <- get.table(d, Cna)
                tab$Freq <- ifelse(tab[, Cna] == tab[, Cma] + 
                  tab[, Cpa], 1, 0)
                set.table(d, Cna, tab, type = "cpt")
            }
        }
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
add.motherchild.likd.node <-
function (mixture, Cgt, Mgt, db, ind = 1) 
{
    revid <- list()
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (m %in% Cgt$marker) {
            add.node(d, "Rlikd", subtype = "boolean")
            datam <- mixture$data[[m]]
            q <- datam$freq
            q <- q/sum(q)
            ca1 <- match(Cgt[Cgt$marker == m, ]$allele1, datam$allele)
            ca2 <- match(Cgt[Cgt$marker == m, ]$allele2, datam$allele)
            ma1 <- match(Mgt[Mgt$marker == m, ]$allele1, datam$allele)
            ma2 <- match(Mgt[Mgt$marker == m, ]$allele2, datam$allele)
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
                tab <- get.table(d, "Rlikd")
                tab$Freq <- c(1, 0)
            }
            else {
                if (w == 2) {
                  t <- ma1
                  ma1 <- ma2
                  ma2 <- t
                }
                if (ca2 == ma2 && ca1 != ca2) {
                  Fna <- paste("n_", ind, "_", ca1, sep = "")
                  Fnb <- paste("n_", ind, "_", ca2, sep = "")
                  add.edge(d, "Rlikd", Fna)
                  add.edge(d, "Rlikd", Fnb)
                  tab <- get.table(d, "Rlikd")
                  lrs <- 0.5 * (tab[, Fna] + tab[, Fnb])/(q[ca1] + 
                    q[ca2])
                  mlr <- max(lrs)
                  lrs <- lrs/mlr
                  odd <- 1 == (1:length(lrs))%%2
                  lrs[odd] <- 1 - lrs[odd]
                  tab$Freq <- lrs
                  revid[[m]] <- mlr
                }
                else {
                  Fna <- paste("n_", ind, "_", ca2, sep = "")
                  add.edge(d, "Rlikd", Fna)
                  tab <- get.table(d, "Rlikd")
                  lrs <- 0.5 * tab[, Fna]/q[ca2]
                  mlr <- max(lrs)
                  lrs <- lrs/mlr
                  odd <- 1 == (1:length(lrs))%%2
                  lrs[odd] <- 1 - lrs[odd]
                  tab$Freq <- lrs
                  revid[[m]] <- mlr
                }
            }
            set.table(d, "Rlikd", tab, type = "cpt")
        }
        d$Revid<-revid[m]
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
}
add.relative.likd.node <-
function (mixture, aca, ind = 1) 
{
    revid <- list()
    for (m in setdiff(mixture$markers, "AMEL")) {
        d <- mixture$domains[[m]]
        if (!is.na(match(m, names(aca)))) {
            q <- mixture$data[[m]]$freq
            q <- q/sum(q)
            add.node(d, "Rlikd", subtype = "boolean")
            if (max(aca[[m]]) == 2) {
                a <- which(aca[[m]] == max(aca[[m]]))
                Fna <- paste("n_", ind, "_", a, sep = "")
                add.edge(d, "Rlikd", Fna)
                tab <- get.table(d, "Rlikd")
                lrs <- 0.5 * tab[, Fna]/q[a]
                mlr <- max(lrs)
                lrs <- lrs/mlr
                odd <- 1 == (1:length(lrs))%%2
                lrs[odd] <- 1 - lrs[odd]
                tab$Freq <- lrs
                revid[[m]] <- mlr
            }
            else {
                ab <- which(aca[[m]] == max(aca[[m]]))
                a <- ab[1]
                b <- ab[2]
                Fna <- paste("n_", ind, "_", a, sep = "")
                Fnb <- paste("n_", ind, "_", b, sep = "")
                add.edge(d, "Rlikd", Fna)
                add.edge(d, "Rlikd", Fnb)
                tab <- get.table(d, "Rlikd")
                lrs <- 0.25 * (tab[, Fna]/q[a] + tab[, Fnb]/q[b])
                mlr <- max(lrs)
                lrs <- lrs/mlr
                odd <- 1 == (1:length(lrs))%%2
                lrs[odd] <- 1 - lrs[odd]
                tab$Freq <- lrs
                revid[[m]] <- mlr
            }
            set.table(d, "Rlikd", tab, type = "cpt")
        }
        d$Revid<-revid[m]
        compile(d)
    }
    if ("AMEL" %in% mixture$markers) 
        compile(mixture$domains[["AMEL"]])
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
logL.UKX <-
function (mixture, expr.extra.findings, initialize = FALSE) 
{
    if (initialize) 
        lapply(mixture$domains, initialize.domain)
    C <- mixture$C
    n.unknown <- mixture$n.unknown
    U <- mixture$U
    K <- mixture$K
    n_K <- lapply(mixture$data, function(d) d[mixture$K])
    function(pararray) {
        logL.m <- function(m) {
            domain <- mixture$domains[[m]]
            d <- mixture$data[[m]]
            rho <- pararray[, "rho"]
            xi <- pararray[, "xi"]
            eta <- pararray[, "eta"]
            phi <- pararray[, "phi"]
            for (r in seq_len(mixture$ntraces)) {
                if (!all(is.na(d[, r + 1]))) {
                  evidence <- setCPT.O(domain, rho[[r]], xi[[r]], 
                    eta[[r]], phi[[r]][c(U, K)], d[, r + 1], 
                    C[[r]], n.unknown, n_K[[m]], attr(domain, 
                      "O")[[r]], d$gets_stutter, d$can_stutter, 
                    d$stutter.from)
                  lapply(seq_along(attr(domain, "O")[[r]]), function(i) set.finding(domain, 
                    attr(domain, "O")[[r]][i], evidence[, i]))
                }
            }
		name<-names(mixture$domains)[m]
		for(e in expr.extra.findings) eval(e)
            if("RHuginDomain"%in%class(domain)) propagate(domain)
            get.normalization.constant(domain, log = TRUE)
        }
        sum(sapply(seq_along(mixture$domains), logL.m))
    }
}
logLX <-
function (mixture, expr.extra.findings, presence.only = FALSE, initialize = FALSE) 
{
    if (mixture$n.unknown > 0) {
        if (presence.only) 
            logLpres.UK(mixture, initialize = initialize)
        else logL.UKX(mixture, expr.extra.findings, initialize = initialize)
    }
    else {
        if (presence.only) 
            logLpres.K(mixture)
        else logL.K(mixture)
    }
}

expr.make.findings<-function (z) 
{
# version of make.findings that returns vector of expressions instead of a function
    res<-NULL
    if (length(z) > 0) 
        for (j in 1:length(z)) {
            if (z[[j]][1] == "Male") {
                ind <- z[[j]]$ind
		    expr.as.text<-paste0("if(name==\"AMEL\"){\n", 
			"\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_1\",sep=\"\"), c(0,1,0))\n", sep = "", 
                  "\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_2\",sep=\"\"), c(0,1,0))\n", sep = "", 
                  "}")
            }
            else if (z[[j]][1] == "Female") {
                ind <- z[[j]]$ind
		    expr.as.text<-paste0("if(name==\"AMEL\"){\n", 
			"\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_1\",sep=\"\"), c(0,0,2))\n", sep = "", 
                  "\tset.finding(domain, paste(\"n_\",", ind, 
                  ",\"_2\",sep=\"\"), c(2,0,0))\n", sep = "", 
                  "}")
            }
            else if (z[[j]][1] == "Aca") {
                ind <- z[[j]]$ind
                aca <- z[[j]]$aca
		    expr.as.text<-paste0("w<-", aca, "[[name]]\n", "if(!is.null(w))\n", 
                  "\t{\n", "\tif(is.null(dim(w))) for(k in 1:length(w))\n", 
                  "\tset.finding(domain,paste(\"n_\",", ind, 
                  ",\"_\",k,sep=\"\"),w[k])\n", "\telse for(k in 1:nrow(w))\n", 
                  "\tset.finding(domain,paste(\"n_\",", ind, 
                  ",\"_\",k,sep=\"\"),w[k,])\n", "\t}")
            }
            else if (z[[j]][1] == "Caca") {
                aca <- z[[j]]$aca
		    expr.as.text<-paste0("if(!is.null(", aca, "[[name]]))", " for(k in 1:length(", 
                  aca, "[[name]]))\n", "\tset.finding(domain,paste(\"Cn\",k,sep=\"\"),", 
                  aca, "[[name]][k])")
            }
            else if (z[[j]][1] == "Denom") {
		    expr.as.text<-paste0("\tset.finding(domain, \"denom\", 1)")
		    }
            else if (z[[j]][1] == "Rlikd") {
                aca <- z[[j]]$aca
                cgt <- z[[j]]$cgt
                evid <- z[[j]]$evid
                expr.as.text<-paste0("if(!is.null(domain$", evid, "))\n", 
			"\tset.finding(domain,\"Rlikd\",c(0,domain$",evid, "[[name]]))")
            }
	  res<-c(res,expr.as.text)
        }
	str2expression(res)
}
mixMLX <-
function (mixture, expr.extra.findings, pars, constraints = NULL, phi.eq = FALSE, 
    val = NULL, trace = FALSE, order.unknowns = TRUE, initialize = FALSE, ...) 
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
    logl <- logLX(mixture,expr.extra.findings,initialize=initialize)
    funvals <- numeric(0)
    minus.loglikelihood <- function(x) {
        xs <- x2arr(x)
        if (trace) 
            print(xs)
        val <- -logl(xs)
        if (trace) 
            print(-val)
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
        ineqUB = ineqUB, control = list(trace = 0), ...)
    est <- x2arr(soln$pars)
    class(est) <- "mixpar"
    val <- -tail(soln$value, 1)
    out <- list(mle = est, lik = val, funvals = funvals, starting.point = x0, 
        minimization.output = soln)
}



