## ph_fitting.R
#.
#. initiate required packages and functions to fit PH-family distribution
#.  PH-family : PH, logPH, scaled logPH...
#.
#.  also some models from precedent studies 
#.

###################################################################################################

suppressPackageStartupMessages({
  
  library(actuar)
  library(data.table)
  library(evir)
  library(gamlss)
  library(gamlss.dist)
  library(GoFKernel)
  library(magrittr)
  library(dplyr)
  
})

###################################################################################################
### 1. Functions for scaled logPH dist'n and PH dist'n
###################################################################################################

#. back transform, functions for scaled logPH distribution
#. these functions are depend on package 'actuar' 

dPH <- function(x, fit) {
  
  prob <- fit$alpha ; rates <- fit$T
  dphtype(x, prob = prob, rates = rates)
  
}

qPH <- function(p, fit) {
  
  prob <- fit$alpha ; rates <- fit$T
  qphtype(p, prob = prob, rates = rates)
  
}

# density funtcion of scaled logPH distribution. 
dscaled.logph <- 
  
  function(x, prob, rates, c) {
    
    dphtype(log(x) + c, prob = prob, rates = rates) / x
    
  }

# distribution funtcion of scaled logPH distribution. 
pscaled.logph <- 
  
  function(q, prob, rates, c) {
    
    pphtype(log(q) + c, prob = prob, rates = rates)
    
  }

# quantile funtcion of scaled logPH distribution. 
qscaled.logph <- 
  
  function(p, prob, rates, c) {
    
    f <- function(x) pscaled.logph(x, prob = prob, rates = rates, c = c)
    inv.f <- inverse(f, lower = 0, upper = 1e+200)
    sapply(p, inv.f)
    
  }


# random number generator of scaled logPH distribution. 
rscaled.logph <- 
  
  function(n, prob, rates, c) {
    
    exp(rphtype(n = n, prob = prob, rates = rates) - c)
    
  }

# quantile function of phase-type distribution
qphtype <- 
  
  function(p, prob, rates) {
    
    f <- function (x) pphtype(x, prob = prob, rates = rates)
    inv.f <- inverse(f, lower = 0, upper = 1e+300)
    sapply(p, inv.f)
    
  }


# phase type, log likelihood function
lphtype <- 
  
  function(x, prob, rates) {
    
    dphtype(x = x, prob = prob, rates = rates) %>% log() %>% sum()
    
  }

# phase type, log likelihood function
lscaled.logph <- 
  
  function(x, prob, rates, c) {
    
    dscaled.logph(x = x, prob = prob, rates = rates, c = c) %>% log() %>% sum()
    
  }


###################################################################################################
### 2. EMPHT 
###################################################################################################

# type of distributions could be chosen in EMPHT program
dists <- data.frame(typeofdist = c("General phase-type",
  "Hyperexponential", "Sum of exponentials", "Coxian", "Coxian general"))


## EMpht fitting algorithm by Asumessen et al. 
#.
#. this function externally execute EMpht.exe compiled program of Asmussen's C script.
#.  (be aware of the file path)
#. for our purpose, log transformation and shifting are conducted in this function. 
#. arguements automatically used for EMpht.exe are...
#.  : sample as input, no censoring, unweigthed, random initiaion, 
#.    step length for Runge-Kutta procedure(default value)
#. except 'vec' and 'logPH', other inputs are going to be inserted via R console.  
#.
EMpht <- 
  
  function(vec, logPH = FALSE, type = NULL, phases = NULL, iter = NULL,
    show.output.on.console = FALSE, curve.fit = TRUE, ...) {
    
    if (sum(vec < 0) != 0) stop("vec should be postive numeric vector")
    
    ## Step 1. preprocessing target vector
    
    # log tranform
    if(logPH) y <- log(vec)
    else y <- vec
    # shifting constant
    if(min(y) < 0) sh <- abs(min(y)) 
    else sh <- 0 
    # target data to fit PH distribution
    v <- y + sh 
    
    # write a file to be used in EMpht.exe, put -1 to identify the end of the vector 
    write(c(v, -1), file = "unweighted")
    
    ## Step 2. get input from R console terminal
    
    if (is.null(type)) {
      
      print(dists)
      type <- readline("Select distiribution type : ")
      
    } else {
      
      type <- as.character(type)
      
    }
    
    if (is.null(phases))  phases <- readline("Number of phases : ")
    else phases <- as.character(phases)
    if (is.null(iter)) iter <- readline("Number of iteration : ")
    else iter <- as.character(iter)
    
    
    ## Step 3. run EMpht
    suppressWarnings({
      
      # path of EMpht.exe 
      empht <- "EMpht/EMpht.exe"
      
      system(empht, 
        input = c("1", "1", "3", phases, type,
          "4", iter, "1", "3"), show.output.on.console = show.output.on.console)
      # get estimated parameters 
      par_est <- as.matrix(read.table(file = "phases", header = F))
      pi <- par_est[,1]
      Q <- par_est[,2:ncol(par_est)]
      
    })
    
    # clean memory 
    unlink(c("phases", "inputdistr", "sample", "unweighted"))
    
    
    
    ## Step 4. report the result
    
    # return parameters 
    fit <- list(data = vec, dist.type = dists[type,], n.iter = as.numeric(iter),
      alpha = pi, T = Q, c = sh, log.trans = logPH) 
    
    if(curve.fit) ph.curvefit(fit, ...)
    
    return(fit)
    
  }


## find optimal number of phases via LRT 
#.
#. input x should be a numeric vector 
#.
phase.select <-
  
  function(x, start.phase = 2, max.phase = 7, logPH = FALSE, type = NULL, iter = NULL,
    show.output.on.console = FALSE) {
    
    ## Step 1. get input from R console terminal
    
    if (is.null(type)) {
      
      print(dists)
      type <- readline("Select distiribution type : ")
      
    }
    
    if (is.null(iter)) iter <- readline("Number of iteration : ")
    
    ## Step 2. initial fitting 
    p = start.phase
    model.0 <- EMpht(vec = x, phases = p, logPH = logPH, type = type, iter = iter,
      show.output.on.console = show.output.on.console)
    
    ## Step 3. loop
    p.value = 0
    while(p.value < 0.05 & p <= max.phase) {
      
      model.1 <- EMpht(vec = x, phases = p + 1, logPH = logPH, type = type, iter = iter,
        show.output.on.console = show.output.on.console)
      
      # log-liklihoods
      if(logPH) {
        
        l0 <- lscaled.logph(x = x, prob = model.0$alpha, rates = model.0$T, c = model.0$c)
        l1 <- lscaled.logph(x = x, prob = model.1$alpha, rates = model.1$T, c = model.1$c)
        
      } else {
        
        l0 <- lphtype(x = x, prob = model.0$alpha, rates = model.0$T)
        l1 <- lphtype(x = x, prob = model.1$alpha, rates = model.1$T)
        
      }
      
      p0 <- length(model.0$alpha) ; p1 <- length(model.1$alpha)
      # degrees of freedom
      if (type == 1 ) df <- p1^2 + p1 - p0^2 - p0
      else if (type == 2) df <- 2 * (p1 - p0)
      else if (type == 4) df <- 2 * (p1 - p0)
      else if (type == 5) df <- 3 * (p1 - p0)
      else stop("not available distribution type")
      
      # lambda statistic
      lambda <- 2 * (l1 - l0)
      
      # p.value
      p.value <- 1 - pchisq(q = lambda, df = df)
      
      # set up for next loop
      if(p.value <= 0.05) {
        p = p + 1
        model.0 = model.1 
      }
      
    }
    
    ph.curvefit(model.0)
    
    # PH parameter estimates of optimal phases 
    if(p.value > 0.05) 
      cat("   optimal number of phases : ", p, "(p.value =", p.value, ")\n") 
    else
      cat("number of phase have reached to maximum : ", max.phase, "\n")  
    
    return(model.0)    
    
  }

## drawing empirical histogram and fitted density curve of PH distribtion
ph.curvefit <- 
  
  function(fit, ...) {
    
    # sample histogram, estimated density, checking model fitting
    # assign legend on plot 
    leg.box <- c(as.character(fit$dist.type),
      paste("n.phase : ", length(fit$alpha), sep = ""),
      paste("iteration : ", fit$n.iter, sep =""))
    
    if (fit$log.trans) {
      
      par(mfrow = c(1,2))
      v <- log(fit$data) + fit$c
      v %>% hist(main = "Log-shifted data(PH)", probability = T, breaks = 50, xlab = "", ...) 
      x <- seq(from = 0, to = max(v), length.out = 1000)
      curve(dphtype(x, prob = fit$alpha, rates = fit$T), add = TRUE, lwd = 2)
      
      fit$data %>% hist(main = "Original data(Scaled LogPH)", probability = T, breaks = 50, ...)
      x <- seq(from = 0, to = max(v), length.out = 1000)
      curve(dscaled.logph(x, prob = fit$alpha, rates = fit$T, c = fit$c),  add = T, lwd = 2)
      legend('topright', legend = leg.box, box.col = "white")
      
      par(mfrow = c(1,1))
      
    } else {
      
      v <- fit$data
      v %>% hist(main = "", probability = T, breaks = 50, xlab = "", ...)
      x <- seq(from = 0, to = max(v), length.out = 1000)
      curve(dphtype(x, prob = fit$alpha, rates = fit$T), add = TRUE, lwd = 2)
      legend('topright', legend = leg.box, box.col = "white")
      
    }
    
    
  }


## drawing qqplot with fitted PH-family distribution 
#. insert fitted result of EMpht() or phase.select()
#. calculating parametric qunatiles consume significant time...
#. print AIC and BIC
ph.qqplot <-
  
  function(fit, ...) {
    
    y <- fit$data ;  n <- length(y) 
    type <- fit$dist.type ; p <- length(fit$alpha)
    x0 <- c(1:n)/(n+1)
    
    # save qqplot
    if (fit$log.trans) {
      x <- do.call("qscaled.logph", list(x0, prob = fit$alpha, rates = fit$T, c = fit$c))
      loglik <- sum(log(dscaled.logph(y, prob = fit$alpha, rates = fit$T, c = fit$c)))
    } else {
      x <- do.call("qphtype", list(x0, prob = fit$alpha, rates = fit$T))
      loglik <- sum(log(dphtype(y, prob = fit$alpha, rates = fit$T)))
    } 
    
    # number of parameters 
    if (type == dists[1,]) k <- p + p^2
    else if (type == dists[2,])  k <- p * 2
    else if (type == dists[4,]) k <- p * 2 - 1
    else if (type == dists[5,]) k <- p * 3 - 1
    
    # AIC & BIC
    AIC <- 2*k - 2*loglik
    BIC <- k * log(n) - 2*loglik
    print(c(AIC = AIC, BIC = BIC))
    
    qqplot(x, y, xlab = "Theoretical values", ylab = "Observed samples",
      main = "Model : Phase-type", xlim = c(0, max(y)), ...)
    abline(0, 1, col = 2, lty = 2)
    
  }


## import data of each station with its station id(GHCN data)
get.station.data <- 
  
  function(stationid) {
    
    # multiply 0.1 to trasform tenth of mm to mm
    fread(paste(FILE.DIR, stationid, ".dly", sep = "")) %>% unlist(., use.names = FALSE) * 0.1
    
  }


###################################################################################################
### 3. Functions for models in Papalexiou et al. [2012] and Li et al. [2013]
###################################################################################################


## general qqplot function with given data and model (quantile function and parameter estimation)
#.
#.  parameter estimation
#.  model quantile function
#.
#.

#. general-purpose estimation function for every considered model in analysis
generalEst <- 
  
  function(y, model, ...) {
  
    # sample
    y <- y[!is.na(y)]
  
    if(!(model %in% c("exp", "gamma", "kappa", "GG", "GB2", "GGP", "EGP", "PH")))
      # GGP := FK08, EGP := HEG in Li et al. 2012
      stop("invalid model name")
  
    if (model == "exp") {
      
      est = 1/mean(y) ; k <- 1
      
    } else if (model == "gamma") {
      
      fit <- MASS::fitdistr(y, 'gamma') ; k <- 2
      est <- fit$estimate
      
    } else if (model == "kappa") {
      
      fit <- fit.kappa(y, ...) ; k <- 3
      est <- fit
      
    } else if (model == "GG") {
      
      fit <- gamlssML(y, family = 'GG', ...) ; k <- 3
      est <- c(mu = fit$mu, sigma = fit$sigma, nu = fit$nu)
      
    } else if (model == "GB2") {
      
      fit <- gamlss::gamlssML(y, family = "GB2", ...) ; k <- 4
      est <- c(mu = fit$mu, sigma = fit$sigma, nu = fit$nu, tau = fit$tau)
      
    } else if (model == "GGP") {
      
      # dynamic input for theta? 
      fit <- GGP(y, ...) ; k <- 4
      est <- fit
      
    } else if (model == "EGP") {
      
      fit <- EGP(y, ...) ; k <- 3
      est <- fit
      
    } else if (model == "PH") {
      
      # setting default input
      fit <- EMpht(y, phases = 3, type = 4, iter = 5000, curve.fit = FALSE)
      k <- length(fit$alpha)*2 - 1
      est <- fit 
      
    }
  
    list(est = est, k = k)
  
  }


#. general-purpose infomation(AIC, BIC) generating function
generalInfo <-
  
  function(y, model = "exp gamma kappa GG GB2 GGP EGP PH") {
    
    # sample
    y <- y[!is.na(y)]
    n <- length(y)
    # theoretical quantile
    x0 <- c(1:n)/(n+1)
    
    models <- unlist(strsplit(model, split = " "))
    
    if(sum(!(models %in% c("exp", "gamma", "kappa", "GG", "GB2", "GGP", "EGP", "PH"))) != 0)
      # GGP := FK08, EGP := HEG in Li et al. 2012
      stop("invalid model name")
    
    sapply(models, function(x) {
      
      ests <- generalEst(y, model = x) ; k <- ests$k ; est <- ests$est
      if (x == 'gamma')
        v <- log(do.call(paste("d", x, sep = ""), list(y, shape = est[1], rate = est[2])))
      else if (x == 'GG') {
        v <- log(do.call(paste("d", x, sep = ""),
          list(y, mu = est[1], sigma = est[2], nu = est[3])))
      } else if (x == 'GB2') {
        v <- log(do.call(paste("d", x, sep = ""),
          list(y, mu = est[1], sigma = est[2], nu = est[3], tau = est[4])))
      } else {
        v <- log(do.call(paste("d", x, sep = ""), list(y, est)))
      }
      
      v[is.infinite(v)] <- NA
      loglik <- sum(v, na.rm = T)
      
      c(AIC = 2*k - 2*loglik, BIC = k * log(n) - 2*loglik)
      
    })
    
  }


#. general-purpose quantile function
generalQQplot <-
  
  function(y, model = "exp gamma kappa GG GB2 GGP EGP PH", as.quantile = FALSE) {
    
    # sample
    y <- y[!is.na(y)]
    n <- length(y)
    # theoretical quantile
    x0 <- c(1:n)/(n+1)
    
    models <- unlist(strsplit(model, split = " "))
    p <- length(models)
    
    if(sum(!(models %in% c("exp", "gamma", "kappa", "GG", "GB2", "GGP", "EGP", "PH"))) != 0)
      # GGP := FK08, EGP := HEG in Li et al. 2012
      stop("invalid model name")
    
    qmat <- sapply(models, function(x) {
      
      est <- generalEst(y, model = x)$est
      if (x == 'gamma')
        do.call(paste("q", x, sep = ""), list(x0, shape = est[1], rate = est[2]))
      else if (x == "GG"){
        do.call(paste("q", x, sep = ""), list(x0, mu = est[1], sigma = est[2], nu = est[3]))
      } else if (x == "GB2") {
        do.call(paste("q", x, sep = ""), 
          list(x0, mu = est[1], sigma = est[2], nu = est[3], tau = est[4]))
      } else {
        do.call(paste("q", x, sep = ""), list(x0, est)) 
      }
       
    })
    
    qmat %<>% data.table()
    
    if (as.quantile) return(qmat)
    # else plotting 
    # the layout should be assigned before
    else {
      
      
      invisible(sapply(models, function(m) {
        
        x <- select_(qmat, m) %>% unlist(., use.names = FALSE)
        qqplot(x, y, xlab = "Theoretical Quantiles", ylab = "Ordered Sample",
          xlim = c(0, max(y)), main = paste("Model : ", m, sep = ""))
        abline(0, 1,  lty = 2, col = 'red')
        
      }))
      
    }
    
  }

#. generally drawing qqplot via qmat made by generaQQplot(..., as.quantile=TRUE)
generalQQviaQmat <-
  
  function(qmat, layout.matrix) {
    
    layout(mat = layout.matrix)
    y <- qmat$y
    models <- colnames(qmat)
    models <- models[models != "y"]
    
    invisible(sapply(models, function(m) {
      
      x <- select_(qmat, m) %>% unlist(., use.names = FALSE)
      qqplot(x, y, xlab = "Theoretical Quantiles", ylab = "Ordered Sample",
        xlim = c(0, max(y)), main = paste("Model : ", m, sep = ""))
      abline(0, 1,  lty = 2, col = 'red')
      
    }))
    
  }


## GG and GB2 of Papalexiou et al. [2012]
#.
#.  using packages {gamlss}, {gamlss.dist}
#.  just embeded in general.qqplot()
#.
#.

## Kappa, Mielke, 1973


kappa_nLL <- function(y) {
  
  n <- length(y)
  
  function(alpha, beta, theta) {
    - ( n * log((alpha*theta)/beta) + (theta - 1) * sum(log(y/beta)) - ((alpha + 1) / alpha) * 
        sum(log(alpha + (y/beta)^(alpha * theta))) )
  }
  
}


#. fitting 3-parameter kappa distribution via ML method
fit.kappa <- 
  
  function(y) {
  
    start <- list(alpha = runif(1, 0.5, 2), beta =  runif(1, 0.5, 2), theta =  runif(1, 0.5, 2))
    
    trial <- function(start)
      stats4::mle(kappa_nLL(y), start, method = "L-BFGS-B", lower = c(0,0,0))
    
    out <- try(trial(start), TRUE)
    
    if (class(out) == "try-error") {
      
      while(class(out) == "try-error") {
        
        # reset starting point 
        start <- list(alpha = runif(1, 0.5, 2), beta =  runif(1, 0.5, 2), theta =  runif(1, 0.5, 2))
        
        out <- try(trial(start), TRUE)
        
      }
      
    }
    
    return(out)
    
  }

#. density, probability, quantile function of kappa
dkappa <- 
  
  function(x, fit) {
  
    alpha <- fit@coef[1] ; beta <- fit@coef[2] ; theta <- fit@coef[3]
    (alpha * theta / beta) * (x / beta)^(theta - 1) * 1 /
      (alpha + (x / beta)^(alpha * theta))^(1 + 1/alpha)
  
  }

pkappa <- 
  
  function(q, fit) {
  
    alpha <- fit@coef[1] ; beta <- fit@coef[2] ; theta <- fit@coef[3]
    ((q / beta)^(alpha * theta) / (alpha + (q / beta)^(alpha * theta)))^(1/alpha)

  }  

qkappa <- 
  
  function(p, fit) {
  
    alpha <- fit@coef[1] ; beta <- fit@coef[2] ; theta <- fit@coef[3]
    beta / ( (1/p^alpha - 1) / alpha)^(1/(alpha * theta))
  
  }


## GGP 
#.  
#.  1) estimating the gamma parameters with all the data
#.  2) determining a reasonable threshold and estimating GP parameters above a threshold
#.  3) adjusting GP scale parameter by eq. (6b) in Li et al. [2013]
#.  4) default theta is 55% quantile derived from the best result of ID24 of Li et al. 
#.  

# parameter estimation 
GGP <- 
  
  function(y, theta = quantile(y, 0.55), show.tail.fit = FALSE) {
    
    y <- y[!is.na(y)]
    if(sum(y < 0) != 0) stop("x should be positive")
    
    # step 1.
    gam_fit <- MASS::fitdistr(y, 'gamma')
    gam_par <- gam_fit$estimate
    
    # step 2.
    gp_fit <- gpd(y, threshold = theta, method = 'ml')
    gp_par <- gp_fit$par.ests
    
    if(show.tail.fit)  summary(fExtremes::gpdFit(y, u = theta))
    
    # step 3.
    gp_par[2] <- 1 / flexsurv::hgamma(theta, shape = gam_par[1], rate = gam_par[2])
    
    GGP_par <- c(gam_par, gp_par)
    
    list(data = y, par.ests = GGP_par, threshold = theta)
    
  }

# denstiy function

dGGP <-
  
  function(x, fit) {
    
    theta <- fit$threshold ; GGP_par <- fit$par.ests
    
    c(dgamma(x[x <= theta], shape = GGP_par[1], rate = GGP_par[2]), 
      (1 - pgamma(theta, shape = GGP_par[1], rate = GGP_par[2])) * 
        dgpd(x[x > theta], xi = GGP_par[3], beta = GGP_par[4], mu = theta))
    
  }


# distribution function

pGGP <-
  
  function(q, fit) {
    
    theta <- fit$threshold 
    
    under_theta <- pgamma(q[which(q <= theta)], shape = fit$par.ests[1], rate = fit$par.ests[2])
    over_theta <- pgamma(theta, shape = fit$par.ests[1], rate = fit$par.ests[2]) +
      pgpd(q[which(q > theta)], xi = fit$par.ests[3], beta = fit$par.ests[4], mu = theta) * 
      (1 - pgamma(theta, shape = fit$par.ests[1], rate = fit$par.ests[2]))
    
    c(under_theta, over_theta)
    
  }



# quantile function of model GGP
qGGP <- 
  
  function(p, fit) {
    
    y <- fit$data
    GGP_par <- fit$par.ests
    theta <- fit$threshold
    F_gamma_theta <- pgamma(q = theta, shape = GGP_par[1], rate = GGP_par[2])
    
    under_theta <- qgamma(p[which(p <= F_gamma_theta)], shape = GGP_par[1], rate = GGP_par[2])
    over_theta <- qgpd((p[which(p > F_gamma_theta)] - F_gamma_theta)/(1 - F_gamma_theta),
      xi = GGP_par[3], beta = GGP_par[4], mu = theta) 
    
    
    c(under_theta, over_theta)
  }



## EGP
#.
#.  main model of Li et al. [2013]
#.  obtain estimates via ML method, in R stats4::mle() with negative log-likelihood
#.  quantile function in closed form...
#.

#. negative log-likelihodd
EGP_nLL <- function(y) {
  
  n <- length(y)
  function(mu, xi, sigma) {
  
    theta <- -mu * log(mu/sigma)
    Z <- pexp(theta, rate = 1/mu) + 1 
    l <- - n * log(Z) + sum(- log(mu) - y[which(y <= theta)] / mu) + 
      sum(- log(sigma) - (1/xi + 1) * log(1 + xi * (y[which(y > theta)] - theta) / sigma))
    -l  
  
  }
  
}


# parameter estimation 
EGP <- 
  
  function(y, without.bound = FALSE) {
    
    setting <- function() {
      
      list(mu = median(y), 
        xi = 0.2,
        sigma = runif(1, 10, 20))
    
    }
    
    start <- setting()
    
    trial <- function(start) {
      
      if (without.bound)
        stats4::mle(EGP_nLL(y), start)
      
      else
        stats4::mle(EGP_nLL(y), start, method = "L-BFGS-B",
          lower = c(0,0,0), upper = c(Inf, 1, Inf))  
    }
      
    
    out <- try(trial(start), TRUE)
    
    if (class(out) == "try-error") {
      
      while(class(out) == "try-error") {
        
        # reset starting point 
        start <- setting()
        
        out <- try(trial(start), TRUE)
        
      }
      
    }
    
    
    return(out)
      
  }

# density function
dEGP <-
  
  function(x, fit) {
    
    mu = fit@coef[1] ; xi = fit@coef[2] ; sigma = fit@coef[3]
    theta = - mu * log(mu/sigma) ; Z = pexp(theta, rate = 1/mu) + 1
    
    under_theta <- dexp(x[which(x <= theta)], rate = 1/mu)
    over_theta <- dgpd(x[which(x > theta)], xi = xi, beta = sigma, mu = theta)
    c(under_theta, over_theta) / Z
    
  }
    

# distribution function
pEGP <-
  
  function(q, fit) {
    
    mu = fit@coef[1] ; xi = fit@coef[2] ; sigma = fit@coef[3]
    theta = - mu * log(mu/sigma) ; Z = pexp(theta, rate = 1/mu) + 1
    
    c(pexp(q[which(q <= theta)], rate = 1/mu),
      pexp(theta, 1/mu) + pgpd(q[which(q > theta)], xi = xi, beta = sigma, mu = theta)) / Z
    
  }


# quantile function in closed form
qEGP <-
  
  function(p, fit) {
    
    mu = fit@coef[1] ; xi = fit@coef[2] ; sigma = fit@coef[3]
    theta = - mu * log(mu/sigma) ; Z = pexp(theta, rate = 1/mu) + 1
    F_exp_theta <- pexp(theta, rate = 1/mu) 
    
    under_theta <- qexp(p[which(p <= F_exp_theta / Z)] * Z, rate = 1/mu)
    over_theta <- qgpd((p[which(p > F_exp_theta / Z)] * Z - F_exp_theta),
      xi = xi, beta = sigma, mu = theta) 
    
    c(under_theta, over_theta)
    
  }



###################################################################################################
### 4. Functions especially used in USCHN data
###################################################################################################

# flattening function
flat <- 
  
  # input vector y should be integer 
  
  function(y) {
    
    # flattening data
    z <- table(y)
    zz <- c()
    for ( i in 1:length(z) ) {
      
      interval <- c(as.numeric(names(z[i])) - 0.5, as.numeric(names(z[i])) + 0.5)
      v <- seq(from = interval[1], to = interval[2], length.out = z[i] + 2)
      zz <- c(zz, v[-c(1, length(v))])
      
    }
    
    zz
    
  }


# a function extracting daily data by each station without missing value and 0's 
by_station <-
  
  # insert a positive integer from 1 to 49 corresponding to station ID number needed 
  # flattening 
  
  function(ID, in.mm = TRUE, flat = FALSE) {
    
    id.value <- 
      tidy.ushcn.tx %>% select(COOP_ID) %>% unique() %>% slice(ID) %>% as.integer
    
    # extract a vector
    y <- 
      tidy.ushcn.tx %>% filter(as.integer(COOP_ID) == id.value) %>% 
      select(-c(COOP_ID, YEAR, MONTH)) %>% c() %>% unlist(., use.names = F) %>%
      subset(!(is.na(.)) & . != 0)
    
    if (flat) v <- flat(y) 
    else v <- y 
    
    # change unit(1/100 inch to 1/10 mm)
    if (in.mm) v * 0.254
    else v
    
  }

