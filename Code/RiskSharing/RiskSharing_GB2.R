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
library(MASS)
library(gamlss.dist)
library(GB2)
library(QRM)
library(ismev)
library(xtable)
library(randomcoloR)
library(latex2exp)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd

####################
### Fonction utiles -----------------------------------------------------------
####################
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

####################
### Analyse initiale -----------------------------------------------------------
####################

##### Base de données #####
don <- read.csv("cafc_nfid_incident_new.csv", stringsAsFactors = T)
perte <- subset(don, DOLLOSSA > 0, select = DOLLOSSA)$DOLLOSSA
don_share <- subset(don,
                    DOLLOSSA > 0 & !is.na(GENCONST), 
                    select = c("DOLLOSSA", "PROPCLAS", "PROPGRP", "GENCONST"))

### Ajuster les données en format date
don$temps <- paste(don$YEAR, don$MONTH, don$DATE, sep = "-")
don$temps <- paste0(don$temps, 'T', don$TIME)
patron <- "%Y-%m-%dT%H:%M:%S"
don$temps <- as.POSIXct(don$temps, patron, tz = "UTC")
min <- 0
date <- sort(na.omit(don$temps[don$DOLLOSSA > min]))
don$YEAR <- as.factor(don$YEAR)

### Classe par type de construction ###
table(don_share$GENCONST)[2:6]

### Classification du batiment ###
table(don_share$PROPGRP)

### Liste avec les pertes pour les 5 catégorie
div <- 1e-5
perte_shar <- vector(mode = 'list', 5)
for (i in 1:5) {
  perte_shar[[i]] <- subset(don_share, GENCONST == i, select = DOLLOSSA)$DOLLOSSA
}

perte_shar <- lapply(perte_shar, '*', div)

# Nombre de jour
min(date)
max(date)
minDate <- as.Date("2004-12-31")
maxDate <- as.Date("2015-12-31")
total_jour <- as.numeric(maxDate - minDate)


lam <- sapply(perte_shar, function(i) length(i)/total_jour)
lam_an <- lam*365.25
lam <- readRDS('lam.RDS')
sum(sapply(perte_shar, sum))

meplot(perte_shar[[1]])


ParetoQQ(perte_shar[[1]], xlim=c(0, 60))
abline(h=log(3e5*div))
th <- 3e5*div

#### Répartition empirique ####
Fn1 <- ecdf(perte_shar[[1]])
Fn2 <- ecdf(perte_shar[[2]])
Fn3 <- ecdf(perte_shar[[3]])
Fn4 <- ecdf(perte_shar[[4]])
Fn5 <- ecdf(perte_shar[[5]])

##### Modélisation #####
### Portion Pareto généralisée
# mle_pg <- lapply(perte_shar, gpdFit, threshold = th)
# par_pg <- sapply(seq_along(mle_pg), function(i) unname(mle_pg[[i]]$par.est))
# xi <- par_pg[2,]
# b <- par_pg[1,]
# a <- 1/xi
# a

# ### Parametre de départ
# parInit <- lapply(seq_along(mle_pg),
#                   function(i) c(ml.gb2(perte_shar[[i]])$opt1$par, xi[i], b[i], get(paste0('Fn', i))(th)))
#
# ### Raccordement
# mle_gb2.pg <- lapply(seq_along(perte_shar),
#                      function(i) spl.gb2.pg(perte_shar[[i]], th, parInit[[i]]))
# par_gb2.pg <- sapply(seq_along(mle_gb2.pg), function(i) unname(mle_gb2.pg[[i]]$par))
par_gb2.pg <- readRDS('par_gb2.pg.RDS')

1/par_gb2.pg[7, ]

### Vérification graphique
plot(Fn5, xlim = c(0, 1e2), ylim = c(0.9, 1), lwd = 2)
curve(pgb2.pg(x, th, par_gb2.pg[, 5]), add = TRUE, col = 'darkblue', lwd=2)

##### Risk Sharing #####
### Avec discrétisation
nfft <- 2^25
h <- 0.01
m <- 2e5
k <- (0:(nfft - 1))*h
u <- 0.9
th <- 3
# # Discrétisé les v.a. B (sévérité)
# fb <- vector(mode = 'list', length(perte_shar))
# for (i in seq_along(perte_shar)) {
#   fb[[i]] <- discretize(pgb2.pg(x, th, par_gb2.pg[, i]), 0, m, h, method = 'lower')
#   fb[[i]] <- c(fb[[i]], numeric(nfft - length(fb[[i]])))
# }
fb <- readRDS('fb.RDS')

lapply(fb, sum)

## Calculer la FFT de B
# fbt <- lapply(fb, fft)

## Calculer la FFT de X
# fgp_po <- function(t, lam) exp(lam*(t - 1))
# fxt <- sapply(seq_along(fbt), function(i) fgp_po(fbt[[i]], lam[i]))

## Calculer fs et Fs
# fst <- apply(fxt, 1, prod)
fst <- readRDS('fst.RDS')
# fs <- Re(fft(fst, TRUE))/nfft
fs <- readRDS('fs.RDS')
Fs <- cumsum(fs)
VaRS <- k[min(which(Fs >= u))]
TVaRS <- VaRS + sum(pmax(k - VaRS, 0)*fs)/(1 - u)
cbind(VaRS, TVaRS)

sum(fs*k)

## Calculer la FFT de phi et e1
# phit <- sapply(seq_along(fb), function(i) fft((k + 1)*c(fb[[i]][-1], 0)))
phit <- readRDS('phit.RDS')
e1 <- exp((-2i*k*pi)/nfft)

## Calculer la conditional expectation
mut <- lapply(seq_along(fb), function(i) lam[i] * e1 * phit[, i] * fst)
mu <- lapply(seq_along(fb), function(i) Re(fft(mut[[i]], TRUE))/nfft)
cm <- lapply(mu, '*', 1/fs)



CVaRX <- sapply(seq_along(fb), function(i) cm[[i]][VaRS + 1])
CTVaRX <- sapply(seq_along(fb),
                 function(i) (sum(mu[[i]][k > VaRS]) + CVaRX[i] * (Fs[max(which(k <= VaRS))] - u))/(1-u))
sum(CTVaRX)

### Avec simulation
ms <- 5e6
# X <- sapply(seq_along(perte_shar),
#             function(i) rcomppois(ms, lam[i], rgb2.pg(th, par_gb2.pg[, i])))
X <- readRDS('X.RDS')

mean(X[, 3])

S <- rowSums(X)
Ss <- sort(S)

VaR <- function(u) Ss[u*ms]
TVaR <- function(u) mean(S[S > VaR(u)])

VaR(u)
TVaR(u)

CTVaR <- function(u, port) mean(port * I(S > VaR(u)))/(1-u)
CTVaR(u, X[, 1])

CTVaRXi <- sapply(1:5, function(i) CTVaR(u, X[, i]))
names(CTVaRXi) <- paste0('X', 1:5)
CTVaRXi

sum(CTVaRXi)

##### Sharing risque non identique #####

EX <- apply(X, 2, mean)
ES <- mean(S)

w <- EX/ES

C <- w*S











