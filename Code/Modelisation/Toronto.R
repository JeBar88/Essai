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
library(gamlss.dist)
library(GB2)
library(MASS)
library(xtable)
library(randomcoloR)
library(latex2exp)
library(Distributacalcul)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd


####################
### Analyse initiale -----------------------------------------------------------
####################

### Base de données ###
don <- read.csv("Fire Incidents Data.csv", stringsAsFactors = T)
perte <- don$Estimated_Dollar_Loss[!is.na(don$Estimated_Dollar_Loss) & don$Estimated_Dollar_Loss > 0]

## Ajuster les données en format date
patron <- "%Y-%m-%dT%H:%M:%S"
don$TFS_Alarm_Time <- as.POSIXct(don$TFS_Alarm_Time, patron, tz = "UTC")
don$Date <- as.Date(don$TFS_Alarm_Time, format = "%Y-%m-%dT", tz = "UTC")
don$years <- as.factor(format(don$TFS_Alarm_Time, format = "%Y"))

## Perte en fonction de l'annne
perte_don <- subset(don, years != 2019 & Estimated_Dollar_Loss > 0, select = c(Estimated_Dollar_Loss, years))
perte_anne <- unstack(perte_don[, c("Estimated_Dollar_Loss", "years")])[-9]

### Statistique descriptive ###
summary(perte)
Esp <- mean(perte)
Var <- var(perte)
n <- length(perte)
m <- sapply(perte_anne, length)
cbind(Esp, Var, n)

### Densité et log-densité ###
hist(perte, breaks = 10000, freq = T, main = "Histogramme des montants de sinistres estimés",
     xlab = "Montants de sinistres", ylab = "Nombre",
     xlim = c(0, 3e5))
hist(log(perte), breaks = 50, freq = T, main = "",
     xlab = "Log des montants de sinistres", ylab = "Nombre")

### Fonction de réparition ###
Fn <- ecdf(perte)
plot(Fn, main = "Fonction de répartition empirique de X", xlim = c(0, 6e7))

Fm <- ecdf(perte[perte > 1e6])


### QQplot exponentielle, lognormale et Pareto ###
# Exponentielle QQ-Plot
ExpQQ(perte)

# QQ plot Lognormale
LognormalQQ(perte)

# QQ plot Pareto
ParetoQQ(perte, xlim = c(0, 20))
abline(h=log(75000), col='red')


### Fonction d'excès moyen ###
MeanExcess(perte, main = "Fonction d'excès moyen pour les montants de sinistres")

Hill(perte, plot = T, k=T)

###  Modèle lognormale sur toutes les données ###

# Méthode des moments
r1 <- sqrt(log(Var/Esp^2 + 1))
m1 <- log(Esp) - (r1^2)/2

# Maximume de vraisemblance
logvrais.LN <- function(par){
  -sum(log(dlnorm(perte, par[1], par[2])))
}

# Paramètres estimés
mle.LN <- constrOptim(c(m1, r1), logvrais.LN, grad = NULL, ui = diag(2), ci = c(0, 0))
m <- mle.LN$par[1]
r <- mle.LN$par[2]
cbind(m, r)

# Comparaison avec fonction de réparition empirique
plot(Fn, xlim = c(0, 1e7), ylim=c(0.9, 1), main = "")
curve(plnorm(x, m, r), col = "red", lwd = 2, add = T)

k <- c(0.9, 0.95, 0.99)
varLnorm(m, r)
VaR <- VatRLnorm(k, m, r)
TVaR <- TVatRLnorm(k, m, r)

var(perte)



xtable(rbind(sapply(k, function(u) mean(perte[perte > quantile(perte, u)])), TVaR), digits = 0)


### Méthode POT ### 

# Minimiser la distance entre la statistique d'ordre la plus élevé et la pareto généralisé
#u <- mindist(sort(perte), method = "ks")$thres # 1E6

# Pas très constant avec b=10
#danielsson(sort(perte), B=1000)

# Constant mais trop bas
#gomes(perte, B=1000) # 150000

# Test la normalité pas constant
#u <- TH(sort(perte), seq(1E5, 1E7,,150))

# Test de seuil
u <- 75000
l <- 150000
# Trouver les paramètres de la portion Pareto généralisée
mle.PaG <- unname(gpdFit(perte, u, method = "mle")$par.e)
xi <- mle.PaG[2]
s <- mle.PaG[1]
cbind("alpha"=1/xi, "lambda"=s/xi, xi, s)

mle.PaG2 <- unname(gpdFit(perte, l, method = "mle")$par.e)
xi2 <- mle.PaG2[2]
s2 <- mle.PaG2[1]
cbind("alpha2"=1/xi2, "lambda2"=s2/xi2, xi2, s2)

# Comparaison avec fonction de réparition empirique
Fx.PaG <- function(x, u) ifelse(x >= u, (Fn(u) + (1-Fn(u)) * pgpd(x, xi, u, s)), NA)
Fx.PaG2 <- function(x, u) ifelse(x >= u, (Fn(u) + (1-Fn(u)) * pgpd(x, xi2, u, s2)), NA)

Fx.PaG(10e6, u)

plot(Fn, xlim=c(1e6, 5e7), ylim=c(0, 1),
     lwd = 2, main = "", ylab=TeX('$F_n(x)$'))
curve(Fx.PaG(x, u), col = "red", lwd = 2, add = T)
curve(Fx.PaG2(x, l), col = "green", lwd = 2, add = T)
curve(plnorm(x, m, r), col = "blue", lwd = 2, add = T)
legend(3e6, 0.99, legend = c("LN", "PaG(u = 75000)", "PaG(u=150000)"), col = c("blue", "red", "green"), lwd = 2)
legend(3e7, 0.95, legend = "PaG(u = 75000)", col = "red", lwd = 2)



### VaR et TVaR de la portion Pareto généralisée ###
k <- c(0.95, 0.99, 0.995, 0.999)
VaR.PaG(k, u)
TVaR.PaG(k, u)

## PHT
parR <- c(0.75, 0.9, 1)
u + s/(parR - xi)

## GS
k <- c(0.9, 0.95, 0.99)
del <- 0.25
s/xi*((1/(1 - xi) * (1 - k)^(-xi) - 1)) + 1/(1 - k)^xi * 2*del*s/((1 - xi)*(2 - xi)) + u


Info <- list("Stat_desc"=cbind(Esp, Var, n), "Param_LN"=cbind(m, r), "Param_PaG"=cbind(xi, s),
             "VaR.PaG_emp"=VaR.PaG(k, u), "TVaR.PaG_emp"=TVaR.PaG(k, u))

#############################
### Raccordement de deux lois --------------------------------------------------
#############################
########## Lognormale - Pareto généralisée ###########
  
parInit <- c(m, r, xi, s, Fn(u))

mle_ln.pg <- Spl.ln.pg(perte, u, parInit, cont = F, deriv = F)
param_ln.pg <- mle_ln.pg$par
cbind("m"=param_ln.pg[1], "r"=param_ln.pg[2], "a"= 1/param_ln.pg[3], "l"=param_ln.pg[4]/param_ln.pg[3], "p"=param_ln.pg[5])

# Graphique
plot(Fn, ylim=c(0, 1), xlim=c(0, 1e6), 
     main = "Comparaison de la fonction de répartition de la loi LN-PaG \navec la fonction de répartition empirique")
curve(pln.pg(x, u, param_ln.pg), add = T, col = "red", lwd = 2)

# Graphique
#plot(ecdf(perte_anne[[3]]), ylim=c(0, 1), xlim=c(0, 1e6))
#curve(pln.pg(x, param_ln_annee[6], param_ln_annee[-6]), add = T, col = "red", lwd = 2)

########## weibull - Pareto généralisée ########### 
  
# Méthode des moments pour les paramètres de la Weibull
f <- function(par) Esp^2 * (gamma(1 + 2/par)/gamma(1 + 1/par)^2 - 1)
t <- optimize(function(x) abs(f(x) - Var), c(0, 50))$min
b <- Esp/gamma(1 + 1/t)

# Optimisation splicing
parInit <- c(t, b, xi, s, Fn(u))
mle_we.pg <- Spl.we.pg(perte, u, parInit, cont = F)
mle_we.pg.d <- Spl.we.pg(perte, u, parInit[-5], cont = T)
param_we.pg <- mle_we.pg$param
param_we.pg.d <- mle_we.pg.d$param
cbind("t"=param_we.pg[1], "b"=param_we.pg[2], "a"=1/param_we.pg[3], "l"=param_we.pg[4]/param_we.pg[3], "p"=param_we.pg[5])
cbind("t"=param_we.pg.d[1], "b"=param_we.pg.d[2], "a"=1/param_we.pg.d[3], "l"=param_we.pg.d[4]/param_we.pg.d[3], "p"=param_we.pg.d[5])

# Graphique
plot(Fn, ylim=c(0, 1), xlim = c(0, 1e6), 
     main = "Comparaison de la fonction de répartition de la loi We-PaG \navec la fonction de répartition empirique")
curve(pwe.pg(x, u, param_we.pg), add = T, col = "red", lwd = 2)
curve(pwe.pg(x, u, param_we.pg.d), add = T, col = "red", lwd = 2)

########### coxienne-2 - Pareto généralisée ###########
# Optimisation splicing
MOM <- function(par){
  p <- 1 - 1/par[2] * (Esp - par[1])
  2*(p*par[1]^2 + (1 - p) * (par[1]^2 + par[2]^2 + (par[1] * par[2]))) - Esp^2
}
mle.MOM <- constrOptim(c(1e4, 1e6), function(x) abs(MOM(x) - Var), grad = NULL, ui=diag(2), ci=c(0, 0))

b1 <- 1/mle.MOM$par[1]
b2 <- 1/mle.MOM$par[2]
p <- 1 - b2 * (Esp - 1/b1)

lg <- function(par){
  -sum(log(dcox2(perte, par[1], par[2], par[3])))
}
ui <- rbind(diag(3), c(-1, rep(0, 2)))
ci <- c(rep(0, 3), -1)
mle.cox2 <- constrOptim(c(p, b1, b2), lg, grad = NULL, ui=ui, ci=ci)

parInit <- c(mle.cox2$par[1], mle.cox2$par[2], mle.cox2$par[3], xi, s, Fn(u))
mle_cox2.pg <- spl.cox2.pg(perte, u, parInit, cont = F)
mle_cox2.pg.d <- spl.cox2.pg(perte, u, parInit[-6], cont = T)
param_cox2.pg <- mle_cox2.pg$param
param_cox2.pg.d <- mle_cox2.pg.d$param
cbind("q"=param_cox2.pg[1], "b1"=param_cox2.pg[2], "b2"=param_cox2.pg[3],
      "a"=1/param_cox2.pg[4], "l"=param_cox2.pg[5]/param_cox2.pg[4], "p"=param_cox2.pg[6])
cbind("q"=param_cox2.pg.d[1], "b1"=param_cox2.pg.d[2], "b2"=param_cox2.pg.d[3],
      "a"=1/param_cox2.pg.d[4], "l"=param_cox2.pg.d[5]/param_cox2.pg.d[4], "p"=param_cox2.pg.d[6])

# Graphique
plot(Fn, ylim=c(0.0, 1), xlim=c(0, 1e6),
     main = "Comparaison de la fonction de répartition de la loi Cox2-PaG \navec la fonction de répartition empirique")
curve(pcox2.pg(x, u, param_cox2.pg.d), add = T, col = "red", lwd = 2)


########## bêta de type 2 - Pareto généralisée ##########

parInit <- c(ml.gb2(perte)$opt1$par, xi, s, Fn(u))

mle_gb2.pg <- spl.gb2.pg(perte, u, parInit, cont = F)
mle_gb2.pg.d <- spl.gb2.pg(perte, u, parInit[-7], cont = T)

param_gb2.pg <- mle_gb2.pg$param
param_gb2.pg.d <- mle_gb2.pg.d$param
rbind("a"=param_gb2.pg[1], "b"=param_gb2.pg[2], "p"=param_gb2.pg[3], "q"=param_gb2.pg[4],
      "alpha"=1/param_gb2.pg[5], "l"=param_gb2.pg[6]/param_gb2.pg[5], "w"=param_gb2.pg[7])
rbind("a"=param_gb2.pg.d[1], "b"=param_gb2.pg.d[2], "p"=param_gb2.pg.d[3], "q"=param_gb2.pg.d[4],
      "alpha"=1/param_gb2.pg.d[5], "l"=param_gb2.pg.d[6]/param_gb2.pg.d[5], "w"=param_gb2.pg.d[7])

# Graphique
plot(Fn, ylim=c(0, 1), xlim=c(0, 1e6), 
     main = "")
curve(pgb2.pg(x, u, param_gb2.pg.d), add = T, col = "red", lwd = 2)

# ## QQplots
# QQplot_an <- function(data, an){
#   uk <- ppoints(length(data))
#   x <- qln.pg(uk, param_ln_annee[6, an], param_ln_annee[-6, an])
#   y <- qwe.pg(uk, param_we_annee[6, an], param_we_annee[-6, an])
#   z <- qgb2.pg(uk, param_gb2_annee[8, an], param_gb2_annee[-8, an])
#   plot(log(x), log(sort(data)), xlab="", ylab="", main = "LN-PaG", axes=F, frame.plot=T)
#   matplot(log(sort(data)), log(sort(data)), col = "red", type = "l", add=T)
#   plot(log(y), log(sort(data)), xlab="", ylab="", main = "We-PaG", axes=FALSE, frame.plot=TRUE)
#   matplot(log(sort(data)), log(sort(data)), col = "red", type = "l", add=T)
#   plot(log(z), log(sort(data)), xlab="", ylab="", main = "GB2-PaG", axes=FALSE, frame.plot=TRUE)
#   matplot(log(sort(data)), log(sort(data)), col = "red", type = "l", add=T)
# }
# 
# par(mfrow=c(1, 3))
# QQplot_an(perte_anne[[1]])
# 
# par(mfrow=c(4, 3))
# par(mar = c(2, 2, 2, 2))
# sapply(seq_along(perte_anne)[1:4], function(i) QQplot_an(perte_anne[[i]], i))
# sapply(seq_along(perte_anne[5:8]), function(i) QQplot_an(perte_anne[[i]], i))
# 
# 
# length(names(perte_anne))


########## Test quantitatif sur les modèles ##########
## Tous les donnee
ad.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
ad.test(perte, function(x) pwe.pg(x, u, param_we.pg))$s
ad.test(perte, function(x) pwe.pg(x, u, param_we.pg.d))$s
ad.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg))$s
ad.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg.d))$s
ad.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg.d))$s

cvm.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
cvm.test(perte, function(x) pwe.pg(x, u, param_we.pg))$s
cvm.test(perte, function(x) pwe.pg(x, u, param_we.pg.d))$s
cvm.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg))$s
cvm.test(perte, function(x) pcox2.pg(x, u, param_cox2.pg.d))$s
cvm.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg.d))$s

### AIC ###
2*mle_ln.pg$value - 2 * length(mle_ln.pg$par) 
2*mle_we.pg$value - 2 * length(mle_we.pg$par) 
2*mle_we.pg.d$value - 2 * length(mle_we.pg.d$par)  
2*mle_cox2.pg$value - 2 * length(mle_cox2.pg$par)  
2*mle_cox2.pg.d$value - 2 * length(mle_cox2.pg.d$par)
2*mle_gb2.pg.d$value - 2 * length(mle_gb2.pg.d$par)

### BIC ###
2*mle_ln.pg$value - length(mle_ln.pg$par)  * log(n)
2*mle_we.pg$value - length(mle_we.pg$par) * log(n)
2*mle_we.pg.d$value - length(mle_we.pg.d$par)  * log(n)
2*mle_cox2.pg$value - length(mle_cox2.pg$par) * log(n)
2*mle_cox2.pg.d$value - length(mle_cox2.pg.d$par) * log(n)
2*mle_gb2.pg.d$value - length(mle_gb2.pg.d$par) * log(n)

####################################
### Information sur le modèle LN-PaG -------------------------------------------
####################################
### Espérance du modèle ###
m <- param_ln.pg[1]
r <- param_ln.pg[2]
xi <- param_ln.pg[3]
s <- param_ln.pg[4]
w <- param_ln.pg[5]
E <- w * exp(m + r^2/2) * pnorm((log(u) - m - r^2)/r)/pnorm((log(u) - m)/r) + (1 - w) * (u + s/(1 - xi))

### Mesure de risque ###
k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qln.pg(k, u, param_ln.pg)
tln.pg(k, u, param_ln.pg)

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
mean(X[X > quantile(X, 0.995)])  # TVaR

qemp <- sapply(k, function(k) quantile(perte, k))
temp <- sapply(k, function(k) mean(perte[perte > quantile(perte, k)]))

####################################
### Information sur le modèle GB2-Pg -------------------------------------------
####################################
### Espérance du modèle ###
a <- param_gb2.pg.d[1]
b <- param_gb2.pg.d[2]
p <- param_gb2.pg.d[3]
q <- param_gb2.pg.d[4]
xi <- param_gb2.pg.d[5]
s <- param_gb2.pg.d[6]
w <- param_gb2.pg.d[7]

E.gb2.pg <- function(par){
  a <- par[1]
  b <- par[2]
  p <- par[3]
  q <- par[4]
  xi <- par[5]
  s <- par[6]
  w <- par[7]
  h <- function(x){
    (x/b)^(a)/(1 + (x/b)^(a))
  }
  w/pbeta(h(u), p, q) * mgb2(1, a, b, p, q) * pbeta(h(u), p + 1/a, q - 1/a) + (1 - w) * (u + s/(1 - xi))
}

E.gb2.pg(param_gb2.pg.d)

E.gb2.pg_tr <- function(d, par){
  a <- par[1]
  b <- par[2]
  p <- par[3]
  q <- par[4]
  xi <- par[5]
  s <- par[6]
  w <- par[7]
  H <- 1 + xi/s*(d - u)
  w/pbeta(h(u), p, q) * mgb2(1, a, b, p, q) * pbeta(h(u), p + 1/a, q - 1/a) +
    (1 - w) * (u - d*H^(-1/xi) - s/(1 - xi)*(H^(1 - 1/xi) - 1)) * I(d > u)
}

alpha <- c(seq(1.2, 2, 0.2), 2.5)
#alpha <- 1.2
xi <- 1/alpha

resu1 <- numeric(6)
for (i in seq_along(xi)) {
  resu1[i] <- (E.gb2.pg_tr(8e6, c(a, b, p, q, xi[i], s, w)) - E.gb2.pg_tr(4e6, c(a, b, p, q, xi[i], s, w)))
}
resu1

resu2 <- numeric(6)
for (i in seq_along(xi)) {
  resu2[i] <- (E.gb2.pg_tr(12e6, c(a, b, p, q, xi[i], s, w)) - E.gb2.pg_tr(8e6, c(a, b, p, q, xi[i], s, w)))
}
resu2

resu3 <- numeric(6)
for (i in seq_along(xi)) {
  resu3[i] <- (E.gb2.pg_tr(16e6, c(a, b, p, q, xi[i], s, w)) - E.gb2.pg_tr(12e6, c(a, b, p, q, xi[i], s, w)))
}
resu3

resu4 <- numeric(6)
for (i in seq_along(xi)) {
  resu4[i] <- (E.gb2.pg_tr(20e6, c(a, b, p, q, xi[i], s, w)) - E.gb2.pg_tr(16e6, c(a, b, p, q, xi[i], s, w)))
}
resu4

res <- rbind(resu1, resu2, resu3, resu4)



xtable(res, digits = 0)

mean(perte * I(perte >= 4e6 & perte < 8e6))/(Fn(8e6) - Fn(4e6))



mean(perte* I(perte > u))/mean(perte)

### Mesure de risque ###
k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qgb2.pg(k, u, param_gb2.pg.d)
tgb2.pg(k, u, param_gb2.pg.d)

pgb2(30, 1, 2, 3, 4)
pGB2(30, 2, 1, 3, 4)

### Simulation ###
X <- rgb2.pg(1e6, u, param_gb2.pg.d)
mean(X)                         # Moyenne
mean(X[X > quantile(X, 0.95)])  # TVaR


####################################
### Analyse annuel                   -------------------------------------------
####################################

### Trouver les paramètres ###
parInit <- c(m, r, xi, s, Fn(u), u)
par_ln.pg_an <- sapply(perte_anne, function(i) Spl.ln.pg_an(i, parInit)$par)

parInit <- c(t, b, xi, s, Fn(u), u)
par_we.pg_an <- sapply(perte_anne, function(i) Spl.we.pg_an(i, parInit)$par)

parInit <- sapply(perte_anne, function(i) c(ml.gb2(i)$opt1$par, xi, s, Fn(u), u))
par_gb2.pg_an <- sapply(seq_along(perte_anne), function(i) spl.gb2.pg_an(perte_anne[[i]], parInit[, i])$par)


xtable(rbind(par_ln.pg_an, par_we.pg_an, par_gb2.pg_an), digits = 2)

### Nombre d'observations par années
n_an <- sapply(perte_anne, length)

### Mesure de risque ###
VaR_ln.pg_an <- sapply(seq_along(perte_anne), function(i) qln.pg(0.95, par_ln.pg_an[6, i], par_ln.pg_an[-6, i]))
VaR_we.pg_an <- sapply(seq_along(perte_anne), function(i) qwe.pg(0.95, par_we.pg_an[6, i], par_we.pg_an[-6, i]))
VaR_gb2.pg_an <- sapply(seq_along(perte_anne), function(i) qgb2.pg(0.95, par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]))

TVaR_ln.pg_an <- sapply(seq_along(perte_anne), function(i) tln.pg(0.95, par_ln.pg_an[6, i], par_ln.pg_an[-6, i]))
TVaR_we.pg_an <- sapply(seq_along(perte_anne), function(i) twe.pg(0.95, par_we.pg_an[6, i], par_we.pg_an[-6, i]))
TVaR_gb2.pg_an <- sapply(seq_along(perte_anne), function(i) tgb2.pg(0.95, par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]))

VaR_TVaR_modele <- rbind(VaR_ln.pg_an, TVaR_ln.pg_an,
                         VaR_we.pg_an, TVaR_we.pg_an,
                         VaR_gb2.pg_an, TVaR_gb2.pg_an)
colnames(VaR_TVaR_modele) <- names(perte_anne)

print(xtable(Tableau), format.args = list(big.mark = " ", decimal.mark = ","))


### Prédiction ###
modele <- c('ln.pg', "we.pg", 'gb2.pg')
pred <- pred_mod('ln.pg')
for (i in modele[-1]) {
  pred <- rbind(pred, pred_mod(i))
}
xtable(pred)

### Critère d'information ###
### logvraisemblance
parInit <- c(m, r, xi, s, Fn(u), u)
LV_ln.pg_an <- sapply(perte_anne, function(i) Spl.ln.pg_an(i, parInit)$val)

parInit <- c(t, b, xi, s, Fn(u), u)
LV_we.pg_an <- sapply(perte_anne, function(i) Spl.we.pg_an(i, parInit)$val)

parInit <- sapply(perte_anne, function(i) c(ml.gb2(i)$opt1$par, xi, s, Fn(u), u))
LV_gb2.pg_an <- sapply(seq_along(perte_anne), function(i) spl.gb2.pg_an(perte_anne[[i]], parInit[, i])$val)

### AIC
AIC_mod <- AIC_modele("ln.pg")
for (i in modele[-1]) {
  AIC_mod <- rbind(AIC_mod, AIC_modele(i))
}
AIC_mod

### BIC
BIC_mod <- BIC_modele("ln.pg")
for (i in modele[-1]) {
  BIC_mod <- rbind(BIC_mod, BIC_modele(i))
}
BIC_mod

CI <- rbind(LV_ln.pg_an, LV_we.pg_an, LV_gb2.pg_an, AIC_mod, BIC_mod)
xtable(CI, digits = 0)

### Test statistique ###
ad_ln.pg <- sapply(seq_along(perte_anne), function(i)
  ad.test(perte_anne[[i]], function(x) pln.pg(x, par_ln.pg_an[6, i], par_ln.pg_an[-6, i]))[1])

ad_we.pg <- sapply(seq_along(perte_anne), function(i)
  ad.test(perte_anne[[i]], function(x) pwe.pg(x, par_we.pg_an[6, i], par_we.pg_an[-6, i]))[1])

ad_gb2.pg <- sapply(seq_along(perte_anne), function(i)
  ad.test(perte_anne[[i]], function(x) pgb2.pg(x, par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]))[1])


ks_ln.pg <- sapply(seq_along(perte_anne), function(i)
  ks.test(perte_anne[[i]] + rnorm(length(perte_anne[[i]]), 0, 0.0001),
          function(x) pln.pg(x, par_ln.pg_an[6, i], par_ln.pg_an[-6, i]))[1])

ks_we.pg <- sapply(seq_along(perte_anne), function(i)
  ks.test(perte_anne[[i]] + rnorm(length(perte_anne[[i]]), 0, 0.0001),
          function(x) pwe.pg(x, par_we.pg_an[6, i], par_we.pg_an[-6, i]))[1])

ks_gb2.pg <- sapply(seq_along(perte_anne), function(i)
  ks.test(perte_anne[[i]] + rnorm(length(perte_anne[[i]]), 0, 0.0001),
          function(x) pgb2.pg(x, par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]))[1])


xtable(rbind(ks_ln.pg, ks_we.pg, ks_gb2.pg, ad_ln.pg, ad_we.pg, ad_gb2.pg), digits = 2)

### Test d'adéquation graphique ###
### QQplots

qq_ln.pg_an <- function(data){
  titre <- names(data)
  for (i in seq_along(data)) {
    qu <- function(u) 1:n_an[i]/(n_an[i] + 1)
    plot(qln.pg(qu(u), par_ln.pg_an[6, i], par_ln.pg_an[-6, i]), sort(data[[i]]), ylab="Quantiles observes",
         xlab="Quantiles theoriques", xlim=c(0, u), ylim=c(0, u),
         axes=F, frame=T,
         main=paste0('LN-PaG', ' ', '(', titre[i], ")"))
    abline(0, 1)
  }
  for (i in seq_along(data)) {
    qu <- function(u) 1:n_an[i]/(n_an[i] + 1)
    plot(qln.pg(qu(u), par_ln.pg_an[6, i], par_ln.pg_an[-6, i]), sort(data[[i]]), ylab="Quantiles observes",
         xlab="Quantiles theoriques",
         xlim=c(u, 1e6), ylim=c(u, 1e6),
         axes=F, frame=T)
    abline(0, 1)
  }
}
qq_we.pg_an <- function(data){
  titre <- names(data)
  for (i in seq_along(data)) {
    qu <- function(u) 1:n_an[i]/(n_an[i] + 1)
    plot(qwe.pg(qu(u), par_we.pg_an[6, i], par_we.pg_an[-6, i]), sort(data[[i]]), ylab="Quantiles observes",
         xlab="Quantiles theoriques", xlim=c(0, u), ylim=c(0, u),
         axes=F, frame=T,
         main=paste0('we-PaG', ' ', '(', titre[i], ")"))
    abline(0, 1)
  }
  for (i in seq_along(data)) {
    qu <- function(u) 1:n_an[i]/(n_an[i] + 1)
    plot(qwe.pg(qu(u), par_we.pg_an[6, i], par_we.pg_an[-6, i]), sort(data[[i]]), ylab="Quantiles observes",
         xlab="Quantiles theoriques",
         xlim=c(u, 1e6), ylim=c(u, 1e6),
         axes=F, frame=T)
    abline(0, 1)
  }
}
qq_gb2.pg_an <- function(data){
  titre <- names(data)
  for (i in seq_along(data)) {
    qu <- function(u) 1:n_an[i]/(n_an[i] + 1)
    plot(qgb2.pg(qu(u), par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]), sort(data[[i]]), ylab="Quantiles observes",
         xlab="Quantiles theoriques", xlim=c(0, u), ylim=c(0, u),
         axes=F, frame=T,
         main=paste0('GB2-PaG', ' ', '(', titre[i], ")"))
    abline(0, 1)
  }
  for (i in seq_along(data)) {
    qu <- function(u) 1:n_an[i]/(n_an[i] + 1)
    plot(qgb2.pg(qu(u), par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]), sort(data[[i]]), ylab="Quantiles observes",
         xlab="Quantiles theoriques",
         xlim=c(u, 1e6), ylim=c(u, 1e6),
         axes=F, frame=T)
    abline(0, 1)
  }
}

par(mfrow=c(2, 8))
par(mar = rep(1, 4))
qq_ln.pg_an(perte_anne)
qq_we.pg_an(perte_anne)
qq_gb2.pg_an(perte_anne)


########################
### Processus de Poisson--------------------------------------------------------
########################
## Sélectionne le minimum de montant de sinistre
min <- 0
annee <- subset(don, Estimated_Dollar_Loss > min, select = years)$years
date <- sort(subset(don, Estimated_Dollar_Loss > min, select = Date)$Date)

## Maximum de vraisemblance pour l'ensemble des données
X <- sapply(levels(annee), function(i) sum(annee == i))

### Construction du nombre d'accidents
date_comp <- seq(min(date), max(date), "days")
nb_acc <- data.frame(date_comp, 'nb'=0)
nb_acc[nb_acc$date_comp %in% unique(date), 'nb'] <- unname(table(date))
nb_acc$nb

total_jour <- as.numeric(as.Date("2019-06-30") - as.Date("2010-12-31"))
nb_jour <- as.numeric(c(as.Date("2010-12-31"), unique(date)) - as.Date("2010-12-31"))
lam <- sum(nb_acc$nb)/total_jour

## Graphique pour nombre de sinistres cumulé 
matplot(c(as.Date("2010-12-30"), date_comp), cumsum(c(0, nb_acc$nb)), type = "s", 
        ylab = "Nombre de sinistres cumulés", xlab = "Date", col = 2, lwd = 2)
title("Trajectoire du processus de Poisson pour \nles sinistres supérieurs à 1 million")
lines(c(as.Date("2010-12-31"), sort(unique(date))), nb_jour * lam, col = 1, lwd = "2")

# Simulation de parcours du processus de Poisson homogène
set.seed(2024)
couleur <- distinctColorPalette(7)
t <- seq(as.Date("2011-01-01"), as.Date("2019-06-30"), by="days")
simul <- replicate(7, rpois(total_jour, lam))
for (i in 1:7) {
  matplot(t, cumsum(simul[,i]), type = "s", col = couleur[i], add=T)
}

## Temps entre les sinistres
dif <- diff(sort(don$TFS_Alarm_Time))
tis <- as.numeric(dif)
Fc <- ecdf(tis)

## Maximum de vraisemblance pour une exponentielle
param_tis <- length(tis)/sum(tis)

## graphique comparaison
plot(Fc, ylim=c(0, 1), lwd = 2, xlab= "Heures", main="")
#title(main = "Comparaison de la fonction de répartition de la loi exponentielle \navec la fonction de répartition empirique pour le\n temps entre les sinsitres supérieur à 1 million")
curve(pexp(x, param_tis), xlim = c(0, 2e5), col = 4, lwd = 2, add = T)

ks.test(tis, function(x) pexp(x, param_tis))$p.val
ad.test(tis, function(x) pexp(x, param_tis))$p.val
cvm.test(tis, function(x) pexp(x, param_tis))$p.val

########################
### Fonds mutuelle fictif ------------------------------------------------------
########################

### Montant total des sinsitre pour chaque jour ###
don_sini <- subset(don, Estimated_Dollar_Loss > 0 &
         !is.na(don$Estimated_Dollar_Loss) &
         don$Date %in% unique(date),
         select = c('Estimated_Dollar_Loss', 'Date'))
sini <- -c(0, unname(tapply(don_sini$Estimated_Dollar_Loss, don_sini$Date, FUN = sum)))

### Par jour ###
E_gb2 <- E.gb2.pg(param_gb2.pg)

prime <- E_gb2 * lam
sini <- c(0, sini[-1] + prime)

plot(c(min(date) - 1, unique(date)), cumsum(sini)/1e6, type = "s", 
     xlab = "Date", ylab = "Montant cumulé dans le fonds (millions)",
     main = "Évolution d'un fonds commun avec les sinistres \nde la base de données de la ville de Toronto")
abline(h=0, col = 2, lty = 2, lwd = 1.5)

1.12*prime

### Simulation autres parcourts ###
fond <- function(t, lam, prime, MI = 0, lim = Inf, info=F){
  freq <- rpois(t, lam)
  sev <- -sapply(freq, function(i) sum(pmin(rgb2.pg(i, u, param_gb2.pg), lim)))
  balance <- cumsum(c(MI, sev + prime))
  ruine <- any(balance < 0)
  if(info == T){
    list("Fréquence"=freq, "Sévérité"=sev, "Balance"=balance, "ruine"=ruine)
  }
  else {
    ruine
  }
}

test <- fond(365*5, lam, 551570, info = T, MI = 0)
plot(0:(365*5), test$Balance, type = 's')
abline(h=0)
test$ruine

test <- -fond(1e5, lam, prime, info = TRUE)$Sév
quantile(test, 0.95)

pp <- mean(replicate(1e4, fond(365, lam, 551570 , MI=2e6, lim=Inf)))
pp1.1 <- mean(replicate(1e4, fond(365.25, lam, prime*1.1, 2e6)))
cbind(pp, pp1.1) # prime pure 0.81968 / prime * 1.1 + 1e6 0.5027

1 - param_gb2.pg[7]

1 - pgb2.pg(u, u, param_gb2.pg)


set.seed(2022)
m <- 10 # Nombre de simulation voulue
X <- replicate(m, fond(365, lam, 551570, 2e6, 5e5, info = T)$Balance)

plot(0:(365.25), X[, 1]/1e6, type = "s", xlab = "Date", ylab = "Montant cumulé dans le fond (millions)",
     main = "Évolution d'un fonds commun avec sinistres \nincendies simulés avec la loi GB2-PaG")

couleur <- distinctColorPalette(m-1)
for (i in 2:m) {
  matplot(0:(365), X[, i]/1e6, type = "s", col = couleur[i], add = TRUE)
}
abline(h=0, col = 1, lty=2, lwd=1.5)


########### Analyse pour l'article Forme ###########

don_F <- subset(don, Estimated_Dollar_Loss >= 4e6, select = c('Estimated_Dollar_Loss',
                                                           "Property_Use",
                                                           "Sprinkler_System_Presence",
                                                           "Building_Status"))
xtable(don_F, digits = 0)

levels(don$Sprinkler_System_Presence)[2:3] <- 'sprinkler'

perte_F <- subset(don,
                  Estimated_Dollar_Loss > 1e6 & (Sprinkler_System_Presence == 'sprinkler' | Sprinkler_System_Presence == '3 - No sprinkler system'),
                  select = c('Estimated_Dollar_Loss',"Sprinkler_System_Presence"))


prop.table(table(perte_F$Sprinkler_System_Presence))

aggregate(Estimated_Dollar_Loss ~ Sprinkler_System_Presence, data = perte_F, mean)


mean(perte_F$Estimated_Dollar_Loss)

table(don$Sprinkler_System_Presence)







