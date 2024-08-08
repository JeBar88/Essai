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
### Analyse initiale -----------------------------------------------------------
####################

### Base de données ###
don <- read.csv("cafc_nfid_incident_new.csv", stringsAsFactors = T)
perte <- subset(don, don$DOLLOSSA > 0, select = DOLLOSSA)$DOLLOSSA

## Ajuster les données en format date
don$temps <- paste(don$YEAR, don$MONTH, don$DATE, sep = "-")
don$temps <- paste0(don$temps, 'T', don$TIME)
patron <- "%Y-%m-%dT%H:%M:%S"
don$temps <- as.POSIXct(don$temps, patron, tz = "UTC")
min <- 0
date <- sort(na.omit(don$temps[don$DOLLOSSA > min]))
don$YEAR <- as.factor(don$YEAR)

## Perte en fonction de l'annne
perte_don <- subset(don, DOLLOSSA > 0, select = c(DOLLOSSA, YEAR))
perte_anne <- unstack(perte_don[, c("DOLLOSSA", "YEAR")])[5:11]
perte_anne <- sapply(perte_anne, sort)


# Fréquence des sinistres
summary_annee <- function(x){
  n_t <- length(x)
  n_1 <- sum(I(x >= 1 & x < 5e4))
  n_2 <- sum(I(x >= 5e4 & x < 5e5))
  n_3 <- sum(I(x >= 5e5 & x < 1e6))
  n_4 <- sum(I(x >= 1e6 & x < 1e7))
  n_5 <- sum(I(x >= 1e7 & x < 5e7))
  n_6 <- sum(I(x >= 5e7))
  a <- rbind(n_1, n_2, n_3, n_4, n_5, n_6)
  a
}


xtable(sapply(perte_anne, summary_annee), digits = 1)

# Graphique des pertes en fonction du temps
plot(don$temps[don$DOLLOSSA > 5e6], don$DOLLOSSA[don$DOLLOSSA > 5e6]/1e6,
     type='h', ylab='Perte (millions)', xlab='Temps')

length(don$DOLLOSSA[don$DOLLOSSA > 5e6 & !is.na(don$DOLLOSSA)])



### Statistique descriptive ###
summary(perte)
Esp <- mean(perte)
Var <- var(perte)
n <- length(perte)
cbind(Esp, Var, n)

### Densité et log-densité ###
hist(perte, breaks = 10000, freq = T, main = "Histogramme des montants de sinistres estimés",
     xlab = "Montants de sinistres", ylab = "Nombre",
     xlim = c(0, 3e5))
hist(log(perte), breaks = 50, freq = T, main = "",
     xlab = "Log des montants de sinistres", ylab = "Nombre")

### Fonction de réparition ###
Fn <- ecdf(perte)
plot(Fn, main = "Fonction de répartition empirique de X", xlim = c(0, 5e6))

### QQplot exponentielle, lognormale et Pareto ###
# Exponentielle QQ-Plot
ExpQQ(perte)

# QQ plot Lognormale
LognormalQQ(perte)

# QQ plot Pareto
ParetoQQ(perte[ perte >= 0], xlim = c(0, 30))
abline(h=log(3e5), col="red")

# MEplot(perte, omit=25)

### Fonction d'excès moyen ###
# MeanExcess(perte, main = "Fonction d'excès moyen pour les montants de sinistres")

# Hill(perte, plot = T, k=T)

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

### Méthode POT ### 

# Minimiser la distance entre la statistique d'ordre la plus élevé et la pareto généralisé
#mindist(sort(perte), method = "ks")$thres # 1E6

# Pas très constant avec b=10
#danielsson(sort(perte), B=1000)

# Constant mais trop bas
#gomes(perte, B=500) # 150000

# Test la normalité pas constant
#u <- TH(sort(perte), seq(1E5, 1E7,,150))

# Test de seuil
u <- 3.5e5
l <- 500000
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
Fx.PaG2 <- function(x, u) ifelse(x >= u, (Fn(u) + (1-Fn(u)) * pgpd(x, xi2, l, s2)), NA)

Fx.PaG(5e6, u)

plot(Fn, xlim=c(0, 5e6), ylim=c(0.9, 1),
     lwd = 2, main = "", ylab=TeX('$F_n(x)$'))
curve(Fx.PaG(x, u), col = "red", lwd = 2, add = T)
curve(Fx.PaG2(x, l), col = "green", lwd = 2, add = T)
curve(plnorm(x, m, r), col = "blue", lwd = 2, add = T)
legend(3e6, 0.99, legend = c("LN", "PaG(u = 75000)", "PaG(u=150000)"), col = c("blue", "red", "green"), lwd = 2)

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

(u + param_ln.pg[4]/(0.9 - param_ln.pg[3])) * (1 - param_ln.pg[5]) +
param_ln.pg[5]*45811.418402/plnorm(u, param_ln.pg[1], param_ln.pg[2])

# Graphique
plot(Fn, ylim=c(0.9, 1), xlim=c(0, 5e6), 
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
param_we.pg <- mle_we.pg$param
cbind("t"=param_we.pg[1], "b"=param_we.pg[2], "a"=1/param_we.pg[3], "l"=param_we.pg[4]/param_we.pg[3], "p"=param_we.pg[5])

# Graphique
plot(Fn, ylim=c(0, 1), xlim = c(0, 1e6), 
     main = "Comparaison de la fonction de répartition de la loi We-PaG \navec la fonction de répartition empirique")
curve(pwe.pg(x, u, param_we.pg), add = T, col = "red", lwd = 2)

########## bêta de type 2 - Pareto généralisée ##########
parInit <- c(ml.gb2(perte)$opt1$par, xi, s, Fn(u))

mle_gb2.pg <- spl.gb2.pg(perte, u, parInit, cont = F)
# mle_gb2.pg.d <- spl.gb2.pg(perte, u, parInit[-7], cont = T)

param_gb2.pg <- mle_gb2.pg$param
# param_gb2.pg.d <- mle_gb2.pg.d$param
round(cbind("a"=param_gb2.pg[1], "b"=param_gb2.pg[2], "p"=param_gb2.pg[3], "q"=param_gb2.pg[4],
      "alpha"=1/param_gb2.pg[5], "l"=param_gb2.pg[6]/param_gb2.pg[5], "w"=param_gb2.pg[7]), 2)
# rbind("a"=param_gb2.pg.d[1], "b"=param_gb2.pg.d[2], "p"=param_gb2.pg.d[3], "q"=param_gb2.pg.d[4],
#      "alpha"=1/param_gb2.pg.d[5], "l"=param_gb2.pg.d[6]/param_gb2.pg.d[5], "w"=param_gb2.pg.d[7])

# Graphique
plot(Fn, ylim=c(0.9, 1), xlim=c(0, 5e6), main = "")
curve(pgb2.pg(x, u, param_gb2.pg), add = T, col = "red", lwd = 2)

1/param_gb2.pg[7]
param_ln.pg[4]/param_ln.pg[3]

########## Test quantitatif sur les modèles ##########
## Tous les donnee
ad_ln.pg <- ad.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
ad_we.pg <- ad.test(perte, function(x) pwe.pg(x, u, param_we.pg))$s
ad_gb2.pg <- ad.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg))$s
ad <- rbind(ad_ln.pg, ad_we.pg, ad_gb2.pg)

cvm_ln.pg <- cvm.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
cvm_we.pg <- cvm.test(perte, function(x) pwe.pg(x, u, param_we.pg))$s
cvm_gb2.pg <- cvm.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg))$s
cvm <- rbind(cvm_ln.pg, cvm_we.pg, cvm_gb2.pg)

### AIC ###
aic_ln.pg <- 2*mle_ln.pg$value - 2 * length(mle_ln.pg$par) 
aic_we.pg <- 2*mle_we.pg$value - 2 * length(mle_we.pg$par) 
aic_gb2.pg <- 2*mle_gb2.pg$value - 2 * length(mle_gb2.pg$par)
aic <- rbind(aic_ln.pg, aic_we.pg, aic_gb2.pg)

### BIC ###
bic_ln.pg <- 2*mle_ln.pg$value - length(mle_ln.pg$par)  * log(n)
bic_we.pg <- 2*mle_we.pg$value - length(mle_we.pg$par) * log(n)
bic_gb2.pg <- 2*mle_gb2.pg$value - length(mle_gb2.pg$par) * log(n)
bic <- rbind(bic_ln.pg, bic_we.pg, bic_gb2.pg)

res_test <- cbind(cvm, ad, aic, bic)

print(xtable(res_test, digits = 0), format.args = list(big.mark = " ", decimal.mark = "."))


####################################
### Information sur le modèle GB2-Pg -------------------------------------------
####################################
### Espérance du modèle ###
a <- param_gb2.pg[1]
b <- param_gb2.pg[2]
p <- param_gb2.pg[3]
q <- param_gb2.pg[4]
xi <- param_gb2.pg[5]
s <- param_gb2.pg[6]
w <- param_gb2.pg[7]

h <- function(x){
  (x/b)^(a)/(1 + (x/b)^(a))
}

mean(perte)

E.gb2.pg <- w/pbeta(h(u), p, q) * mgb2(1, a, b, p, q) * pbeta(h(u), p + 1/a, q - 1/a) + (1 - w) * (u + s/(1 - xi))
E.gb2.pg_tr <- function(d){
  w/pbeta(h(u), p, q) * mgb2(1, a, b, p, q) * pbeta(h(u), p + 1/a, q - 1/a)
}

ratio <- (E.gb2.pg - E.gb2.pg_tr(u))/E.gb2.pg
ratio

mean(perte* I(perte > u))/mean(perte)

k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qgb2.pg(k, u, param_gb2.pg)
tgb2.pg(k, u, param_gb2.pg)

xtable(rbind(quantile(perte, k), qgb2.pg(k, u, param_gb2.pg)), digits = 0)

quantile(perte, 0.995)
qln.pg(0.995, u, param_ln.pg)
qgb2.pg(0.995, u, param_gb2.pg)



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

param_ln.pg

### Mesure de risque ###
k <- c(0.90, 0.95, 0.99, 0.995, 0.999)
qln.pg(k, u, param_ln.pg)
tln.pg(k, u, param_ln.pg)

### Simulation ###
X <- rln.pg(1e7, u, param_ln.pg)
mean(X)                         # Moyenne
mean(X[X > quantile(X, 0.99)])  # TVaR

qemp <- sapply(k, function(k) quantile(perte, k))
temp <- sapply(k, function(k) mean(perte[perte > quantile(perte, k)]))

tln.pg(0.99, u, param_ln.pg)

tab <- rbind(tln.pg(k, u, param_ln.pg), tgb2.pg(k, u, param_gb2.pg))

####################################
### Analyse annuel                   -------------------------------------------
####################################

### Trouver les paramètres ###
parInit <- c(m, r, xi, s, Fn(u), u)
par_ln.pg_an <- sapply(perte_anne, function(i) Spl.ln.pg_an(i, parInit)$par)

parInit <- sapply(perte_anne, function(i) c(ml.gb2(i)$opt1$par, xi, s, Fn(u), u))
par_gb2.pg_an <- sapply(seq_along(perte_anne), function(i) spl.gb2.pg_an(perte_anne[[i]], parInit[, i])$par)


### Nombre d'observation par années
n_an <- sapply(perte_anne, length)

xtable(rbind(par_ln.pg_an, par_gb2.pg_an), digits = 2)

### Mesure de risque ###
VaR_ln.pg_an <- sapply(seq_along(perte_anne), function(i) qln.pg(c(0.9, 0.95), par_ln.pg_an[6, i], par_ln.pg_an[-6, i]))
VaR_gb2.pg_an <- sapply(seq_along(perte_anne), function(i) qgb2.pg(c(0.9, 0.95), par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]))

VaR_modele <- rbind(VaR_ln.pg_an, VaR_gb2.pg_an )
colnames(VaR_modele) <- names(perte_anne)

xtable(VaR_modele, digits = 0)

### Prédiction ###
modele <- c('ln.pg', 'gb2.pg')
pred <- pred_mod('gb2.pg')
round(pred, 2)

pred <- rbind(pred_mod('ln.pg'), pred_mod('gb2.pg'))
xtable(pred)

### Statistique ###
### Kolmogorov_smirnov
KS_modele <- function(mod){
  n <- sapply(perte_anne, length)
  p <- get(paste0('p', mod))
  par <- get(paste0('par', "_", mod, '_an'))
  nbar <- dim(par)[1]
  sapply(seq_along(perte_anne),
         function(i) max(abs(p(perte_anne[[i]], par[nbar, i], par[-nbar, i]) - (1:n[i] - 1)/n[i]),
                         abs(p(perte_anne[[i]], par[nbar, i], par[-nbar, i]) - 1:n[i]/n[i])))
}

KS_modele('ln.pg')

KS <- rbind(KS_modele('ln.pg'), KS_modele('gb2.pg'))

AD_ln.pg <- sapply(seq_along(perte_anne), function(i) 
  ad.test(perte_anne[[i]], function(x) pln.pg(x, par_ln.pg_an[6, i], par_ln.pg_an[-6, i]))$stat)
AD_gb2.pg <- sapply(seq_along(perte_anne), function(i) 
  ad.test(perte_anne[[i]], function(x) pgb2.pg(x, par_gb2.pg_an[8, i], par_gb2.pg_an[-8, i]))$stat)

xtable(rbind(KS, AD_ln.pg, AD_gb2.pg))



### Critère d'information ###
### logvraisemblance
parInit <- c(m, r, xi, s, Fn(u), u)
LV_ln.pg_an <- sapply(perte_anne, function(i) Spl.ln.pg_an(i, parInit)$val)

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

CI <- rbind(LV_ln.pg_an, LV_gb2.pg_an, AIC_mod, BIC_mod)
xtable(CI, digits = 0)

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
         xlim=c(u, 1e7), ylim=c(u, 1e7),
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
         xlim=c(u, 1e7), ylim=c(u, 1e7),
         axes=F, frame=T)
    abline(0, 1)
  }
}

par(mfrow=c(2, 7))
par(mar = rep(1, 4))
qq_ln.pg_an(perte_anne)
qq_gb2.pg_an(perte_anne)

(1+0.92)/2


k <- 0.95
b <- 1 - k
mu <- 2
s <- 0.6
1/(1-k)*exp(2 + s^2/2)*(1 - pnorm(qnorm(k) - s))


1/(1-k)*exp(2 + s^2/2)*(pnorm(s - qnorm(k)))

qnorm(0.4)




##################################
### Analyse des plus gros sinistre ---------------------------------------------
##################################
### Perte supérieur à 5 millions
perte_elever <- subset(don, DOLLOSSA > 5e6, 
                       select = c(YEAR, MONTH, DATE,DOLLOSSA, MAJOCC, MAJOCGRP, PROPCLAS, GENCONST))
perte_elever <- perte_elever[order(perte_elever$DOLLOSSA, decreasing=T),]
table(perte_elever$MAJOCC)
table(perte_elever$GENCONST)


perte_don <- subset(don, DOLLOSSA > 0, 
                    select = c(DOLLOSSA, YEAR, MONTH, DATE, TIME, MAJOCC, MAJOCGRP, PROPCLAS, GENCONST))
perte_xx <- perte_don[grepl("62", substr(perte_don$PROPCLAS, 1,2) ,fixed=TRUE), ]
perte_xx <- perte_xx[order(perte_xx$DOLLOSSA, decreasing=T),]

perte_bois <- perte_don[perte_don$GENCONST %in% 3, ]
perte_bois <- perte_bois[order(perte_bois$DOLLOSSA, decreasing=T),]
perte_bois$temps <- paste(perte_bois$YEAR, perte_bois$MONTH, perte_bois$DATE, sep = "-")
perte_bois$temps <- paste0(perte_bois$temps, 'T', perte_bois$TIME)
patron <- "%Y-%m-%dT%H:%M:%S"
perte_bois$temps <- as.POSIXct(perte_bois$temps, patron, tz = "UTC")


plot(perte_bois$temps, perte_bois$DOLLOSSA, type = 'h', 
     xlab = "Temps", ylab = 'Perte estimée')


hist(perte_bois$DOLLOSSA, breaks = 10000, xlim = c(0, 5e5))

var(perte_bois$DOLLOSSA)

table(perte_bois$PROPCLAS)
table(perte_bois$PROPCLAS)[table(perte_bois$PROPCLAS) > 0]
summary(perte_bois$DOLLOSSA)

sapply(perte_anne, mean)

########################
### Processus de Poisson -------------------------------------------------------
########################

table_NA <- sapply(levels(don$YEAR), function(i) sum(is.na(don$temps[don$YEAR == i])))

## Maximum de vraisemblance pour l'ensemble des données
X <- sapply(levels(don$YEAR), function(i) sum(don$YEAR == i))
nb_acc <- unname(table(date))



table(date)[which.max(table(date))]


# Nombre de jour
min(date)
max(date)
minDate <- as.Date("2004-12-31")
maxDate <- as.Date("2015-12-31")
total_jour <- as.numeric(maxDate - minDate)

# Paramètre lambda journalier
lam <- sum(nb_acc)/total_jour
nb_jour <- as.numeric(c(minDate, sort(unique(date))) - minDate)

## Graphique pour nombre de sinistres cumulé 
matplot(c(minDate, sort(unique(date)), maxDate), c(0, cumsum(nb_acc), sum(nb_acc)), type = "s", 
        ylab = "Nombre de sinistres cumulés", xlab = "Date", col = "red", lwd = 2)
#title("Trajectoire du processus de Poisson pour \nles sinistres supérieurs à 1 million")
lines(c(minDate, sort(unique(date)), maxDate), c(nb_jour, total_jour) * lam, col = "blue", lwd = "2")



# Simulation de parcours du processus de Poisson homogène
set.seed(2023)
couleur <- distinctColorPalette(7)
t <- seq(as.Date("2005-01-01"), maxDate, by="days")
simul <- replicate(7, rpois(total_jour, lam))
for (i in 1:7) {
  matplot(t, cumsum(simul[,i]), type = "s", col = couleur[i], add=T)
}


## Date ou 690 sisnistre son arrivé

table(date)[which.max(table(date))]

test <- subset(don, DOLLOSSA > 0 & temps >= '2011-05-15' & temps < '2011-05-16', 
       select = c(DOLLOSSA, MAJOCC, MAJOCGRP, PROPCLAS, GENCONST, JURIS, temps))

dim(test)
table(test$PROPCLAS)[table(test$PROPCLAS) > 0]
table(test$JURIS)[table(test$JURIS) > 0]

## Temps entre les sinistres
dif <- diff(date)
tis <- as.numeric(dif)
Fc <- ecdf(tis)

## Maximum de vraisemblance pour une exponentielle
param_tis <- length(tis)/sum(tis)

## graphique comparaison
plot(Fc, ylim=c(0, 1), xlim=c(0, 4e6), lwd = 2, xlab= "Heures", main="")
#title(main = "Comparaison de la fonction de répartition de la loi exponentielle \navec la fonction de répartition empirique pour le\n temps entre les sinsitres supérieur à 1 million")
curve(pexp(x, param_tis), col = "blue", lwd = 2, add = T)

ExpQQ(tis)

tis <- tis + rnorm(length(tis), 0, 0.00001)
ks.test(tis, pexp, rate=param_tis)



########################
### Analyse style Swiss Re -----------------------------------------------------
########################

don_type_facility <- subset(don, DOLLOSSA > 0, select = c("PROPCLAS", "DOLLOSSA"))
levels(don_type_facility$PROPCLAS) <- substr(levels(don_type_facility$PROPCLAS), 1, 2)
prop <- sort(prop.table(table(don_type_facility$PROPCLAS)))
barplot(tail(prop, 10))

sum(tail(prop))


don_type_facility <- subset(don, DOLLOSSA > 0, select = c("PROPCLAS", "DOLLOSSA"))
levels(don_type_facility$PROPCLAS) <- substr(levels(don_type_facility$PROPCLAS), 1, 1)
prop <- tail(prop.table(table(don_type_facility$PROPCLAS)), 7)
names(prop) <- c("Résidentiel",
                 "Services",
                 "Commerciale",
                 "Industrielle",
                 "Stockage",
                 "test",
                 "Divers")

par(mar=c(4, 6, 4, 4))
barplot(sort(tail(prop, 7)),
        horiz = T,
        las=1,
        xlim = c(0, 0.45),
        xaxt="n",
        col="#F8766D",
        border="#F8766D",
        space = 0.15)
axis(1, at = seq(0, 0.45, 0.05))
grid(nx=9, NA, lty=1) # grid only in y-direction








montant <- aggregate(DOLLOSSA ~ PROPCLAS, data = don_type_facility, sum)[-c(1:3),]
montant <- droplevels(montant)
montant$DOLLOSSA <- sort(montant$DOLLOSSA)
montant$DOLLOSSA <- montant$DOLLOSSA/sum(montant$DOLLOSSA)


barplot(DOLLOSSA ~ PROPCLAS, data=montant,
        horiz = T,
        las=1,
        names.arg = c("Résidentiel",
                                "Services",
                                "Commerciale",
                                "Industrielle",
                                "Stockage",
                                "test",
                                "Divers"),
        ylab = "",
        xlab = "",
        xlim=c(0, 0.6),
        xaxt="n",
        col="#F8766D",
        border="#F8766D",
        space = 0.15)
axis(1, at = seq(0, 0.6, 0.05))
grid(nx=12, NA, lty=1) # grid only in y-direction


cat('\\item', letters,'\n')



for (i in names(don)) {
  cat('\\item', i,'\n')
}














