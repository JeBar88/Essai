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
library(CASdatasets)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd


#####################
### Fconction utilisé
#####################
Spl.ln.pg <- function(data, u, parInit, cont = FALSE, deriv = FALSE){
  
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
Spl.we.pg <- function(data, u, parInit, cont = FALSE, deriv = FALSE){
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
    
    mle <- constrOptim(parInit, logvrais, grad = NULL, ui=diag(4), ci=rep(0, 4))
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
Spl.pa.pg <- function(data, u, parInit){
  
  f <- function(par){
    par[5] * actuar::dpareto(data, par[1], par[2])/
      (actuar::ppareto(u, par[1], par[2]) - actuar::ppareto(500, par[1], par[2])) * I(data <= u) +
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
  ifelse(k <= p, qlnorm(min(k * plnorm(u, m, r)/p, 0.999), m, r), qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
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

dpa.pg <- function(x, u, par){ 
  a <- par[1] ; lam <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * actuar::dpareto(x, a, lam)/actuar::ppareto(u, a, lam) * I(x <= u) +
    (1-p) * dgpd(x, xi, u, sigma) * I(x > u)
}
ppa.pg <- function(x, u, par){ 
  a <- par[1] ; lam <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  p * actuar::ppareto(x, a, lam)/actuar::ppareto(u, a, lam) * I(x <= u) + 
    (p + (1-p) * pgpd(x, xi, u, sigma)) * I(x > u)
}
qpa.pg<- function(k, u, par){
  a <- par[1] ; lam <- par[2] ; xi <- par[3]
  sigma <- par[4] ; p <- par[5]
  
  ifelse(k <= p, actuar::qpareto(pmin(k * actuar::ppareto(u, a, lam)/p, 1), a, lam), 
         qgpd(pmax((k - p)/(1-p), 0), xi, u, sigma))
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
### Analyse initiale
####################

install.packages("CASdatasets", repos = "http://cas.uqam.ca/pub/", type="source")
library(CASdatasets)

### Base de données ###
data("norfire")
perte <- norfire$Loss 

norfire

#don <- as.numeric(read.table("Data Fire Swedish 1982.txt")$V1)
#perte <- don[don > 0]

### Statistique descriptive ###
summary(perte)
Esp <- mean(perte)
Var <- var(perte)
n <- length(perte)
cbind(Esp, Var, n)

### Densité et log-densité ###
hist(perte, breaks = 1000, freq = F, main = "Histogramme des montants de sinsitres estimés",
     xlab = "Montants de sinsitres", ylab = "Densité")
hist(log(perte), breaks = 50, freq = F, main = "Histogramme du log des montants de sinsitres estimés",
     xlab = "Log des montants de sinsitres", ylab = "Densité")

### Fonction de réparition ###
Fn <- ecdf(perte)
plot(Fn, main = "Fonction de répartition empirique de X", lwd = 2)

### QQplot exponentielle, lognormale et Pareto ###
# Exponentielle QQ-Plot
ExpQQ(perte)

# QQ plot Lognormale
LognormalQQ(perte)

# QQ plot Pareto
ParetoQQ(perte)

### Fonction d'excès moyen ###
MeanExcess(perte, xlim=c(0, 50000), ylim=c(0, 20000),
           main = "Fonction d'excès moyen pour les montants de sinistres")

Hill(perte, plot = T, k=T)

###  Modèle Pareto sur toutes les données ###

# Méthode des moments
a1 <- -2*Var/(Esp^2 - Var)
lam1 <- Esp * (a1 - 1)

# Maximume de vraisemblance
logvrais.Pa <- function(par){
  -(sum(log(actuar::dpareto(perte, par[1], par[2])))  - 
      length(perte) * log(actuar::ppareto(500, par[1], par[2], low=F)))
}
mle.Pa <- constrOptim(c(a1, lam1), logvrais.Pa, grad = NULL, ui = diag(2), ci = c(0, 0))
a <- mle.Pa$par[1]
lam <- mle.Pa$par[2]
cbind(a, lam)

# Comparaison avec fonction de réparition empirique
plot(Fn, 
     main = "Fonction de répartition empirique\n comparé avec Pareto",
     cex.lab=1.3, cex.axis=1.3, cex.main=1.3)
curve(ppareto(x, a, lam), col = "red", lwd = 2, add = T)

###  Modèle lognormale sur toutes les données ###

# Méthode des moments
r1 <- sqrt(log(Var/Esp^2 + 1))
m1 <- log(Esp) - (r1^2)/2

# Maximume de vraisemblance
logvrais.LN <- function(par){
  -(sum(log(dlnorm(perte, par[1], par[2])) - length(perte) * log(plnorm(perte, par[1], par[2], low=F))))
}
mle.LN <- constrOptim(c(m1, r1), logvrais.LN, grad = NULL, ui = diag(2), ci = c(0, 0))
m <- mle.LN$par[1]
r <- mle.LN$par[2]
cbind(m, r)

# Comparaison avec fonction de réparition empirique
plot(Fn, main = "Fonction de répartition empirique comparé avec lognormale")
curve(plnorm(x, m, r), col = "red", lwd = 2, add = T)

### Méthode POT ###

# Minimiser la distance entre la statistique d'ordre la plus élevé et la pareto généralisé
mindist(sort(perte), method = "ks")$thres # 1E6

# Pas très constant avec b=10
danielsson(sort(perte), B=5)

# Constant mais trop bas
gomes(perte, B=2000)

## Minimiser la statistique KS
tr.ks(perte, seq(quantile(perte, 0.5), quantile(perte, 0.98), length.out=500))

# Minimiser la distance entre la fonction répartition empirique et celle du modèle
u <- 1672 # 1672
Fn(u)

# Trouver les paramètres de la portion Pareto généralisée
mle.PaG <- unname(gpdFit(perte, u, method = "mle")$par.e)
xi <- mle.PaG[2]
s <- mle.PaG[1]
cbind("alpha"=1/xi, "lambda"=s/xi, xi, s) # alpha = 1.37907

# Comparaison avec fonction de réparition empirique
Fx.PaG <- function(x, u) Fn(u) + (1-Fn(u)) * pgpd(x, xi, u, s)

plot(Fn, ylim=c(0.5, 1), xlim = c(0, 1e5),
     lwd = 2, main = "Fonction de répartition empirique comparé avec Pareto généralisée")
curve(Fx.PaG(x, u), col = "red", lwd = 2, add = T)
curve(plnorm(x, m, r), col = "blue", lwd = 2, add = T)
curve(ppareto(x, a, lam), col = "green", lwd = 2, add = T)

### VaR et TVaR de la portion Pareto généralisée ###
k <- c(0.95, 0.99, 0.995, 0.999)
VaR.PaG(k, u)
TVaR.PaG(k, u)


Info <- list("Stat_desc"=cbind(Esp, Var, n), "Param_Pa"=cbind(a, lam), "Param_LN"=cbind(m, r), "Param_PaG"=cbind(xi, s),
             "VaR.PaG_emp"=VaR.PaG(k, u), "TVaR.PaG_emp"=TVaR.PaG(k, u))


#############################
### Raccordement de deux lois
#############################

### Lognormale - Pareto généralisé ###
parInit <- c(m, r, xi, s, Fn(u))

mle.ln.pg <- Spl.ln.pg(perte, u, parInit, cont = F, deriv = F)
param_ln.pg <- mle.ln.pg$par
cbind("m"=param_ln.pg[1], "r"=param_ln.pg[2], "xi"=param_ln.pg[3], "sigma"=param_ln.pg[4], "p"=param_ln.pg[5])

plot(Fn, ylim=c(0, 1), 
     main = "Fonction de réparition du modèle comparée à celle empirique")
curve(pln.pg(x, u, param_ln.pg), add = T, col = "red", lwd = 2)

plot(Fn, ylim=c(0.99, 1), 
     main = "Fonction de réparition du modèle comparée à celle empirique")
curve(pln.pg(x, u, param_ln.pg), add = T, col = "red", lwd = 2)

plot(Fn, ylim=c(0.95, 1), 
     main = "Fonction de réparition du modèle comparée à celle empirique")
curve(pln.pg(x, u, param_ln.pg), add = T, col = "red", lwd = 2)

plot(Fn, ylim=c(0.90, 1), 
     main = "Fonction de réparition du modèle comparée à celle empirique")
curve(pln.pg(x, u, param_ln.pg), add = T, col = "red", lwd = 2)

### Pareto - Pareto généralisé ###

# Méthode des moments pour les paramètres de la Weibull

# Optimisation splicing
parInit <- c(a, lam, xi, s, Fn(u))
mle_pa.pg <- Spl.pa.pg(perte, u, parInit)
param_pa.pg <- mle_pa.pg
cbind("a"=param_pa.pg[1], "lam"=param_pa.pg[2], "xi"=param_pa.pg[3], "sigma"=param_pa.pg[4], "p"=param_pa.pg[5])

# Graphique
plot(Fn, ylim=c(0.0, 1))
curve(ppa.pg(x, u, param_pa.pg), add = T, col = "red", lwd = 2)

plot(Fn, ylim=c(0.99, 1))
curve(ppa.pg(x, u, param_pa.pg), add = T, col = "red", lwd = 2)


plot(Fn, ylim=c(0.95, 1))
curve(ppa.pg(x, u, param_pa.pg), add = T, col = "red", lwd = 2)

plot(Fn, ylim=c(0.90, 1))
curve(ppa.pg(x, u, param_pa.pg), add = T, col = "red", lwd = 2)


### Test quantitatif sur les modèles ###

ad.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
ad.test(perte, function(x) ppa.pg(x, u, param_pa.pg))$s

cvm.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
cvm.test(perte, function(x) ppa.pg(x, u, param_pa.pg))$s


#mean(replicate(20, wilcox.test(perte, rln.pg(1e6, u, param_ln.pg))$p.val))
#mean(replicate(20, wilcox.test(perte, rwe.pg(1e6, u, param_we.pg))$p.val))

### Espérance du modèle ###
m <- param_ln.pg[1]
r <- param_ln.pg[2]
xi <- param_ln.pg[3]
s <- param_ln.pg[4]
w <- param_ln.pg[5]
w * exp(m + r^2/2) * pnorm((log(u) - m - r^2)/r)/pnorm((log(u) - m)/r) + (1 - w) * (u + s/(1 - xi))

### Mesure de risque ###
k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qln.pg(k, u, param_ln.pg)
qwe.pg(k, u, param_we.pg)

tln.pg(k, u, param_ln.pg)
twe.pg(k, u, param_we.pg)

### Intervalle de confiance ###
Function_PaG <- function(input, index){
  u <- 75000
  Input <- input[index]
  Result <- gpdFit(Input, u, method = "mle")$par.e[2]
  return(Result)}
#Boot <- boot(perte, Function_PaG, R=5000)
#hist(Boot$t[,1])
#boot.ci(Boot, conf = 0.95, type = "bca") Long

### Simulation ###
X <- rln.pg(1e6, u, param_ln.pg)
mean(X)                         # Moyenne
mean(X[X > quantile(X, 0.95)])  # TVaR








