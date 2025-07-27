logLm<-function (mixture, presence.only = FALSE, initialize = TRUE) 
{
    if (mixture$n.unknown > 0) {
        if (presence.only) 
            logLpresm.UK(mixture, initialize = initialize)
        else logLm.UK(mixture, initialize = initialize)
    }
    else {
        if (presence.only) 
            logLpresm.K(mixture)
        else logLm.K(mixture)
    }
}

logLpresm.UK<-function (mixture, initialize = TRUE) 
{
    if (initialize) 
        lapply(mixture$domains, initialize.domain)
    C <- mixture$C
    n.unknown <- mixture$n.unknown
    U <- mixture$U
    K <- mixture$K
    n_K <- lapply(mixture$data, function(d) subset(d, select = K))
    function(pararray) {
        logLpres.m <- function(m) {
            domain <- mixture$domains[[m]]
            d <- mixture$data[[m]]
            rho <- pararray[, "rho"]
            xi <- pararray[, "xi"]
            eta <- pararray[, "eta"]
            phi <- pararray[, "phi"]
            for (r in seq_len(mixture$ntraces)) {
                if (!all(is.na(d[, r + 1]))) {
                  setCPT.D(domain, rho[[r]], xi[[r]], eta[[r]], 
                    phi[[r]][c(U, K)], C[[r]], n.unknown, n_K[[m]], 
                    attr(domain, "D")[[r]], d$gets_stutter, d$can_stutter, 
                    d$stutter.from)
                  mapply(function(x, finding) set.finding(domain, 
                    x, finding), attr(domain, "D")[[r]], ifelse(d[, 
                    r + 1] >= C[[r]], 0, 1))
                }
            }
            propagate(domain)
            get.normalization.constant(domain, log = TRUE)
        }
        sapply(seq_along(mixture$domains), logLpres.m)
    }
}

logLpresm.K<-function (mixture) 
{
    domains <- mixture$domains
    C <- mixture$C
    k <- mixture$k
    data <- mixture$data
    n_K <- lapply(data, function(d) subset(d, select = mixture$K))
    function(pars) {
        rhos <- pars[, "rho"]
        etas <- pars[, "eta"]
        xis <- pars[, "xi"]
        phis <- pars[, "phi"]
        logLpres.m <- function(m) {
            d <- data[[m]]
            heights <- d[, seq_len(mixture$ntraces) + 1, drop = FALSE]
            can_stutter <- d$can_stutter
            gets_stutter <- d$gets_stutter
            stutter.from <- d$stutter.from
            n_K <- as.matrix(n_K[[m]])
            one.allele <- function(a) {
                one.trace <- function(rho, eta, xi, phi, height, 
                  threshold) {
                  phi <- phi[mixture$K]
                  if (is.na(height)) {
                    return(0)
                  }
                  shape <- rho * (1 - xi * can_stutter[a]) * 
                    (n_K[a, ] %*% phi)
                  if (gets_stutter[a]) {
                    st <- stutter.from[a]
                    shape <- shape + rho * xi * (n_K[st, ] %*% 
                      phi)
                  }
                  shape <- as.numeric(shape)
                  if (height == 0) 
                    pgamma(threshold, shape = shape, scale = eta, 
                      log.p = TRUE, lower.tail = TRUE)
                  else pgamma(threshold, shape = shape, scale = eta, 
                    log.p = TRUE, lower.tail = FALSE)
                }
                sum(mapply(FUN = one.trace, rhos, etas, xis, 
                  phis, heights[a, ], C, SIMPLIFY = TRUE))
            }
            sum(sapply(seq_len(nrow(d)), one.allele))
        }
        sapply(seq_along(mixture$markers), logLpres.m)
    }
}

logLm.UK<-function (mixture, initialize = TRUE) 
{
    if (initialize) 
        lapply(mixture$domains, initialize.domain)
    C <- mixture$C
    n.unknown <- mixture$n.unknown
    U <- mixture$U
    K <- mixture$K
    n_K <- lapply(mixture$data, function(d) subset(d, select = K))
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
            propagate(domain)
            get.normalization.constant(domain, log = TRUE)
        }
        sapply(seq_along(mixture$domains), logL.m)
    }
}

logLm.K<-function (mixture) 
{
    C <- mixture$C
    k <- mixture$k
    data <- mixture$data
    n_K <- lapply(data, function(d) subset(d, select = mixture$K))
    function(pars) {
        rhos <- pars[, "rho"]
        etas <- pars[, "eta"]
        xis <- pars[, "xi"]
        phis <- pars[, "phi"]
        logL.m <- function(m) {
            d <- data[[m]]
            heights <- d[, seq_len(mixture$ntraces) + 1, drop = FALSE]
            can_stutter <- d$can_stutter
            gets_stutter <- d$gets_stutter
            stutter.from <- d$stutter.from
            n_K <- as.matrix(n_K[[m]])
            one.allele <- function(a) {
                one.trace <- function(rho, eta, xi, phi, height, 
                  threshold) {
                  phi <- phi[mixture$K]
                  if (is.na(height)) {
                    return(0)
                  }
                  shape <- rho * (1 - xi * can_stutter[a]) * 
                    (n_K[a, ] %*% phi)
                  if (gets_stutter[a]) {
                    st <- stutter.from[a]
                    shape <- shape + rho * xi * (n_K[st, ] %*% 
                      phi)
                  }
                  shape <- as.numeric(shape)
                  if (height == 0) 
                    pgamma(threshold, shape = shape, scale = eta, 
                      log.p = TRUE)
                  else dgamma(height, shape = shape, scale = eta, 
                    log = TRUE)
                }
                sum(mapply(FUN = one.trace, rhos, etas, xis, 
                  phis, heights[a, ], C, SIMPLIFY = TRUE))
            }
            sum(sapply(seq_len(nrow(d)), one.allele))
        }
        sapply(seq_along(mixture$markers), logL.m)
    }
}
