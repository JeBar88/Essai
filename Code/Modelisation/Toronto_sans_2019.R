## Package
library(ggplot2)
library(actuar)
library(formatR)
library(ggplot2)
library(gridExtra)
library(knitr)
library(evir)
library(qrmtools)
library(ReIns)
library(tea)
library(goftest)
library(DescTools)
library(GB2)
library(MASS)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd

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
  sigma <- par[4] ; p <- par[5]
  
  p * plnorm(x, m, r)/plnorm(u, m, r) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qln.pg <- function(k, u, par){
  m <- par[1] ; r <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  ifelse(k <= p, qlnorm(pmin(k * plnorm(u, m, r)/p, 0.999), m, r), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
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

pgb2.pg <- function(x, u, par){ 
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4]
  xi <- par[5] ; sigma <- par[6] ; w <- par[7]
  
  w * pgb2(x, a, b, p, q)/pgb2(u, a, b, p, q) * I(x <= u) + 
    (w + (1-w) * pgpd(x, xi, u, sigma)) * I(x > u)
}

qgb2.pg <- function(k, u, par){ 
  a <- par[1] ; b <- par[2] ; p <- par[3] ; q <- par[4]
  xi <- par[5] ; sigma <- par[6] ; w <- par[7]
  
  ifelse(k <= w, qgb2(pmin(k * pgb2(u, a, b, p, q)/w, 0.999), a, b, p, q), qgpd(pmax((k - w)/(1-w), 0), xi, u, sigma))
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

####################
### Analyse initiale -----------------------------------------------------------
####################

### Base de données ###
don <- read.csv("Fire Incidents Data.csv", stringsAsFactors = T)
patron <- "%Y-%m-%dT%H:%M:%S"
don$TFS_Alarm_Time <- as.Date(don$TFS_Alarm_Time, patron, tz = "UTC")
don$years <- as.factor(format(don$TFS_Alarm_Time, format = "%Y"))
perte <- don$Estimated_Dollar_Loss[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0 & don$years != "2019"]

### Statistique descriptive ###
summary(perte)
Esp <- mean(perte)
Var <- var(perte)
n <- length(perte)
cbind(Esp, Var, n)

########################
### PRocessus de Poisson--------------------------------------------------------
########################
annee <- don$years[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]
date <- don$TFS_Alarm_Time[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]

X <- sapply(levels(annee), function(i) sum(annee == i))
nb_acc <- unname(table(date))

total_jour <- as.numeric(as.Date("2019-06-30") - as.Date("2010-12-31"))

lam <- sum(nb_acc)/total_jour

nb_jour <- as.numeric(c(as.Date("2010-12-31"), sort(unique(date))) - as.Date("2010-12-31"))

matplot(c(as.Date("2010-12-31"), sort(unique(date)), as.Date("2019-06-30")), c(0, cumsum(nb_acc), sum(nb_acc)), type = "s", 
        ylab = "Nombre de sinistres cumulés", xlab = "Date", col = "red", lwd = 2)
title("Trajectoire du processus de Poisson")
lines(c(as.Date("2010-12-31"), sort(unique(date)), as.Date("2019-06-30")), c(nb_jour, total_jour) * lam, col = "blue", lwd = "2")



t <- seq(as.Date("2011-01-01"), as.Date("2019-06-30"), by="days")
simul <- replicate(5, rpois(total_jour, lam))
for (i in 1:5) {
  matplot(t, cumsum(simul[,i]), type = "s", add=T)
}

## Fond
# Montant total des sinsitre pour chaque jour
sini <- -c(0, (sapply(1:length(unique(date)),
                      function(i) sum(don$Estimated_Dollar_Loss[don$Estimated_Dollar_Loss > 0 & 
                                                                  !is.na(don$Estimated_Dollar_Loss) &
                                                                  don$TFS_Alarm_Time == sort(unique(date))[i]]))))

names(sini) <- c(min(date) -1, sort(unique(date)))


prime.p <- Esp * lam # Esp = 36394.08
sini <- c(0, sini[-1] + prime.p)

plot(c(min(date) - 1, sort(unique(date))), cumsum(sini), type = "s", xlab = "Date", ylab = "Montant dans le fond")






