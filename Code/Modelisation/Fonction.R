#####################
### Fconction utilisé
#####################
Spl.ln.pg <- function(data, u, parInit, cont = F, deriv = F){
  
  if(cont == F & deriv == F){
    
    f <- function(par){
      par[5] * dlnorm(data, par[1], par[2])/plnorm(u, par[1], par[2]) * I(data <= u) +
        (1-par[5]) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(5), c(rep(0, 4), -1))
    ci <- c(rep(0, 5), -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- mle$par
  }
  
  if(cont == T & deriv == F){
    
    w <- function(par) {
      (par[4]/plnorm(u, par[1], par[2]) * dlnorm(u, par[1], par[2]) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dlnorm(data, par[1], par[2])/plnorm(u, par[1], par[2]) * I(data <= u) +
        (1-w(par)) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
    mle$param <- c(mle$par, w(mle$par))
  }
  
  if(cont == T & deriv == T){
    
    mu <- function(par){
      log(par[1]) - par[2]^2 * par[1] * (1 + par[3])/par[4]
    }
    
    w <- function(par) {
      (par[4]/plnorm(par[1], mu(par), par[2]) * dlnorm(par[1], mu(par), par[2]) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dlnorm(data, mu(par), par[2])/plnorm(par[1], mu(par), par[2]) * I(data <= par[1]) +
        (1-w(par)) * dgpd(data, par[3], par[1], par[4]) * I(data > par[1])
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
    mle$lim <- mle$par[1]
    mle$param <- c(mu(mle$par), mle$par[2],  mle$par[3:4], w(mle$par))
  }
  mle
}
Spl.we.pg <- function(data, u, parInit, cont = F, deriv = F){
  if(cont == F & deriv == F){
    
    f <- function(par){
      par[5] * dweibull(data, par[1], par[2])/pweibull(u, par[1], par[2]) * I(data <= u) +
        (1-par[5]) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(5), c(-1, rep(0, 4)) , c(rep(0, 4), -1))
    ci <- c(rep(0, 5), -5, -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- mle$par
  }
  
  if(cont == T & deriv == F){
    
    w <- function(par) {
      (par[4]/pweibull(u, par[1], par[2]) * dweibull(u, par[1], par[2]) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dweibull(data, par[1], par[2])/pweibull(u, par[1], par[2]) * I(data <= u) +
        (1-w(par)) * dgpd(data, par[3], u, par[4]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    ui <- rbind(diag(4), c(-1, 0, 0, 0))
    ci <- c(rep(0, 4), -5)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par, w(mle$par))
  }
  
  if(cont == T & deriv == T){
    
    b <- function(par){
      (1/(par[2] * par[1]^par[2]) * (par[1] * (par[3] + 1)/par[4] - (par[2] - 1)))^(-1/par[2])
    }
    
    w <- function(par) {
      (par[4] * dweibull(par[1], par[2], b(par))/pweibull(par[1], par[2], b(par)) + 1)^-1
    }
    
    f <- function(par){
      w(par) * dweibull(data, par[2], b(par))/pweibull(par[1], par[2], b(par)) * I(data <= par[1]) +
        (1-w(par)) * dgpd(data, par[3], par[1], par[4]) * I(data > par[1])
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
    mle$lim <- mle$par[1]
    mle$param <- c(mle$par[2], b(mle$par),  mle$par[3:4], w(mle$par))
  }
  mle
}
spl.cox2.pg <- function(data, u, parInit, cont = F, deriv = F){
  
  if(cont==F & deriv==F){
    
    f <- function(par){
      par[6] * dcox2(data, par[1], par[2], par[3])/pcox2(u, par[1], par[2], par[3]) * I(data <= u) +
        (1-par[6]) * dgpd(data, par[4], u, par[5]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(6), c(-1, rep(0, 5)), c(rep(0, 5), -1))
    ci <- c(rep(0, 6), -1, -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- mle$par
    mle
  }
  
  if(cont == T & deriv == F){
    
    w <- function(par) {
      (par[5]/pcox2(u, par[1], par[2], par[3]) * dcox2(u, par[1], par[2], par[3]) + 1)^-1
    }
    
    f <- function(par){
      w(par) *  dcox2(data, par[1], par[2], par[3])/pcox2(u, par[1], par[2], par[3]) * I(data <= u) +
        (1-w(par)) * dgpd(data, par[4], u, par[5]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(5), c(-1, rep(0, 4)))
    ci <- c(rep(0, 5), -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par, w(mle$par))
  }
  mle
}
spl.gb2.pg <- function(data, u, parInit, cont = F){
  if(cont == F){
    
    f <- function(par) {
      par[7] * dgb2(data, par[1], par[2], par[3], par[4])/pgb2(u, par[1], par[2], par[3], par[4]) * I(data <= u) +
        (1 - par[7]) * dgpd(data, par[5], u, par[6]) * I(data > u)
    }
    
    logvrais <- function(par) {
      -sum(log(f(par)))
    }
    ui <- rbind(diag(7), c(rep(0, 6),-1))
    ci <- c(rep(0, 7),-1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui = ui, ci = ci)
    mle$param <- mle$par
  }
  
  if(cont == T){
    
    w <- function(par) {
      (par[6]/pgb2(u, par[1], par[2], par[3], par[4]) * dgb2(u, par[1], par[2], par[3], par[4]) + 1)^-1
    }
    
    f <- function(par) {
      w(par) * dgb2(data, par[1], par[2], par[3], par[4])/pgb2(u, par[1], par[2], par[3], par[4]) * I(data <= u) +
        (1 - w(par)) * dgpd(data, par[5], u, par[6]) * I(data > u)
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- diag(6)
    ci <- numeric(6)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par, w(mle$par))
  }
  mle
}

dln.pg <- function(x, u, par){ 
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * dlnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) +
    (1-p) * dgpd(x, xi, u, sigma) * I(x > u)
}
pln.pg <- function(x, u, par){ 
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; w <- par[5]
  
  w * plnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) + 
    (w + (1-w) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qln.pg <- function(k, u, par){
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; w <- par[5]
  ifelse(k <= w, qlnorm(k * plnorm(u, m, r)/w, m, r)*(k <= w), qgpd((k - w)/(1-w)*(k > w), xi, u, sigma))
}
tln.pg <- function(k, u, par){
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5] ; c <- qln.pg(k, u, par)
  
  kk <- function(i) (1 + xi/sigma * (i - u))
  
  (p/(plnorm(u, m, r) * (1-k)) * exp(m + r^2/2) * (pnorm((log(u) - m - r^2)/r) - pnorm((log(c) - m - r^2)/r)) + 
      (1-p)/(1-k) * (u + sigma/(1 - xi))) * I(c <= u) + 
    ((1-p)/(1 - k) * (c * kk(c)^(-1/xi) - sigma/(xi - 1) * kk(c)^(1 - 1/xi))) * I(c > u)
}
rln.pg <- function(n, u, par){
  U <- runif(n)
  qln.pg(U, u, par)
}
dwe.pg <- function(x, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * dweibull(x, t, b)/pweibull(u, t, b) * I(x <= u) +
    (1-p) * dgpd(x, xi, u, sigma) * I(x > u)
}
pwe.pg <- function(x, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * pweibull(x, t, b)/pweibull(u, t, b) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qwe.pg<- function(k, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  ifelse(k <= p, qweibull(pmin(k * pweibull(u, t, b)/p, 1), t, b), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
}
twe.pg <- function(k, u, par){
  t <- par[1] ; b <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5] ; c <- qwe.pg(k, u , par)
  
  kk <- function(i) (1 + xi/sigma * (i - u))
  
  (p/(pweibull(u, t, b) * (1-k)) * b * gamma(1 + 1/t) * (pgamma(u^t, 1+1/t, 1/b^t) - pgamma(c^t, 1+1/t, 1/b^t)) +
      (1-p)/(1-k) * (u + sigma/(1 - xi))) * I(c <= u) + 
    ((1-p)/(1 - k) * (c * kk(c)^(-1/xi) - sigma/(xi - 1) * kk(c)^(1 - 1/xi))) * I(c > u)
}
rwe.pg <- function(n, u, par){
  U <- runif(n)
  qwe.pg(U, u, par)
}

dcox2 <- function(x, p, b1, b2){
  p * dexp(x, b1) + (1-p) * (b2/(b2 - b1) * dexp(x, b1) + b1/(b1 - b2) * dexp(x, b2))
}
pcox2 <- function(x, p, b1, b2){
  p * pexp(x, b1) + (1-p) * (b2/(b2 - b1) * pexp(x, b1) + b1/(b1 - b2) * pexp(x, b2))
}
qcox2 <- function(k, p, b1, b2){
  
  f <- function(i) {optimize(function(x) abs(pcox2(x, p, b1, b2) - k[i]), c(0, 500000))$min}
  sapply(1:length(k), f)
}

pcox2.pg <- function(x, u, par){ 
  q <- par[1] ; b1 <- par[2] ; b2 <- par[3] ; 
  xi <- par[4] ; sigma <- par[5] ; p <- par[6]
  
  p * pcox2(x, q, b1, b2)/pcox2(u, q, b1, b2) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qcox2.pg <- function(k, u, par){
  q <- par[1] ; b1 <- par[2] ; b2 <- par[3]
  xi <- par[4] ; sigma <- par[5] ; p <- par[6]
  ifelse(k <= p, qcox2(pmin(k * pcox2(u, q, b1, b2)/p, 0.999), q, b1, b2), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
}
rcox2.pg <- function(n, u, par){
  U <- runif(n)
  qcox2.pg(U, u, par)
}

mgb2 <- function(k, a, b, p, q){
  
  b^k * beta(p + k/a, q - k/a)/beta(p, q)
  
}
pgb2.pg <- function(x, u, par){ 
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4]
  xi <- par[5] ; sigma <- par[6] ; w <- par[7]
  
  w * pGB2(x, b, a, p, q)/pGB2(u, b, a, p, q) * I(x <= u) + 
    (w + (1-w) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qgb2.pg <- function(k, u, par){ 
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4]
  xi <- par[5] ; sigma <- par[6] ; w <- par[7]
  
  qGB2(pmin(k * pGB2(u, b, a, p, q)/w, 0.9999), b, a, p, q) * I(k <= w) + 
    qgpd(pmax((k - w)/(1-w), 0), xi, u, sigma) * I(k > w)
}
rgb2.pg <- function(n, u, par){
  U <- runif(n)
  qgb2.pg(U, u, par)
}
tgb2.pg <- function(k, u, par){
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4] ;
  xi <- par[5] ; sigma <- par[6] ; w <- par[7] ; 
  c <- qgb2.pg(k, u, par) ; 
  
  kk <- function(i) (1 + xi/sigma * (i - u))
  h <- function(i) (i/b)^(a) / (1 + (i/b)^(a))
  
  (w/(pgb2(u, a, b, p, q) * (1-k)) * mgb2(1, a, b, p, q) * (pbeta(h(u), p + 1/a, q - 1/a) - pbeta(h(c), p + 1/a, q - 1/a)) + 
      (1-w)/(1-k) * (u + sigma/(1 - xi))) * I(c <= u) + 
    ((1-w)/(1 - k) * (c * kk(c)^(-1/xi) - sigma/(xi - 1) * kk(c)^(1 - 1/xi))) * I(c > u)
}

tr.ks <- function(data, u){
  
  param <- matrix(numeric(0), ncol = 2, nrow = length(u))
  stat <- numeric(length(u))
  p <- numeric(length(u))
  
  for(i in 1:length(u)){
    param[i, ] <- gpdFit(perte, u[i], method = "mle")$par.e
    test <- ks.test(data[data > u[i]], function(x) ReIns::pgpd(x, param[i, 2], u[i], param[i, 1]))
    stat[i] <- test[[1]]
    p[i] <- test[[2]]
  }
  list("lim - Stat min" = u[which.min(stat)])
}

VaR.PaG <- function(k, u){
  (s/xi * (((1 - k)/(1 - Fn(u)))^(-xi) - 1) + u) * I(k >= Fn(u))
} 
TVaR.PaG <- function(k, u){
  (s/xi * (((1 - k)/(1 - Fn(u)))^(-xi) * 1/(1-xi) - 1) + u) * I(k >= Fn(u))
} 

fond <- function(t, lam, prime, MI = 0, lim = Inf, info=F){
  freq <- rpois(t, lam)
  sev <- -sapply(freq, function(i) sum(pmin(rln.pg(i, u, param_ln.pg), lim)))
  balance <- c(MI, sev + prime)
  ruine <- ifelse(sum(cumsum(balance) < 0) == 0, 0, 1)
  if(info == T){
    list("Fréquence"=freq, "Sévérité"=sev, "Balance"=balance, "Nb.Perte"=ruine)
  }
  else {
    ruine
  }
}

## Ajout
Spl.ln.pg_an <- function(data, parInit){
  
    f <- function(par){
      par[5] * dlnorm(data, par[1], par[2])/plnorm(par[6], par[1], par[2]) * I(data <= par[6]) +
        (1-par[5]) * dgpd(data, par[3], par[6], par[4]) * I(data > par[6])
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    
    ui <- rbind(diag(6), c(rep(0, 4), -1, 0))
    ci <- c(rep(0, 6), -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par[1:2], 1/mle$par[3], mle$par[4]/mle$par[3], mle$par[5:6])
    mle
}
Spl.we.pg_an <- function(data, parInit){
    f <- function(par){
      par[5] * dweibull(data, par[1], par[2])/pweibull(par[6], par[1], par[2]) * I(data <= par[6]) +
        (1-par[5]) * dgpd(data, par[3], par[6], par[4]) * I(data > par[6])
    }
    
    logvrais <- function(par){
      -sum(log(f(par)))
    }
    ui <- rbind(diag(6), c(-1, numeric(5)), c(numeric(4), -1, 0))
    ci <- c(numeric(6), -5, -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=ui, ci=ci)
    mle$param <- c(mle$par[1:2], 1/mle$par[3], mle$par[4]/mle$par[3], mle$par[5:6])
    mle
}
spl.gb2.pg_an <- function(data, parInit){
    f <- function(par) {
      par[7] * dgb2(data, par[1], par[2], par[3], par[4])/pgb2(par[8], par[1], par[2], par[3], par[4]) * I(data <= par[8]) +
        (1 - par[7]) * dgpd(data, par[5], par[8], par[6]) * I(data > par[8])
    }
    
    logvrais <- function(par) {
      -sum(log(f(par)))
    }
    ui <- rbind(diag(8), c(rep(0, 6),-1, 0))
    ci <- c(rep(0, 8), -1)
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui = ui, ci = ci)
    mle$param <- c(mle$par[1:4], 1/mle$par[5], mle$par[6]/mle$par[5], mle$par[7:8])
    mle
}


pred_mod <- function(mod){
  th <- 6
  if(mod == 'gb2.pg'){
    th <- 8
  }
  p <- get(paste0("p", mod))
  par <- get(paste("par", mod, 'an', sep='_'))
  pred <- sapply(1:(length(perte_anne)-1), function(i) 1 - p(1e5, par[th, i], par[-th, i]))
  actual <- sapply(2:length(perte_anne), function(i) mean(perte_anne[[i]] > 1e5))

  
  unname(rbind(pred, actual))
}


modele <- c('ln.pg', 'gb2.pg')
#pred <- pred_mod('gb2.pg')
#round(pred, 2)


AIC_modele <- function(mod){
  th <- 6
  if(mod == 'gb2.pg'){
    th <- 8
  }
  2*get(paste('LV', mod, "an", sep='_')) + 2 * length(get(paste('par', mod, "an", sep='_'))[, 1])
}

mod <- 'ln.pg'
BIC_modele <- function(mod){
  th <- 6
  if(mod == 'gb2.pg'){
    th <- 8
  }
  2*get(paste('LV', mod, "an", sep='_')) + length(get(paste('par', mod, "an", sep='_'))[, 1]) * log(m)
}














