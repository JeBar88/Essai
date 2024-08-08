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
library(randomcoloR)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd

####################
### Analyse initiale -----------------------------------------------------------
####################

### Base de données ###
don <- read.csv("Fire_Incidents_SF.csv", stringsAsFactors = T)
perte <- subset(don, Estimated.Property.Loss > 0, select = Estimated.Property.Loss)$Estimated.Property.Loss

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
ParetoQQ(perte[perte >= 75000], xlim = c(0, 20))
abline(h=log(150000), col="red")


### Fonction d'excès moyen ###
#MeanExcess(perte, main = "Fonction d'excès moyen pour les montants de sinistres")

#Hill(perte, plot = T, k=T)

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
#danielsson(sort(perte), B=50)

# Constant mais trop bas
#gomes(perte, B=2000)

# Test la normalité pas constant

# Test de seuil
u <- 150000
l <- 1e5
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

Info <- list("Stat_desc"=cbind(Esp, Var, n), "Param_LN"=cbind(m, r), "Param_PaG"=cbind(xi, s),
             "VaR.PaG_emp"=VaR.PaG(k, u), "TVaR.PaG_emp"=TVaR.PaG(k, u))


sum(perte > u)


sum(sort(perte)[sort(perte) > u])/sum(perte)

########## bêta de type 2 - Pareto généralisée ##########
parInit <- c(ml.gb2(perte)$opt1$par, xi, s, Fn(u))
mle_gb2.pg <- spl.gb2.pg(perte, u, parInit, cont = F)
# mle_gb2.pg.d <- spl.gb2.pg(perte, u, parInit[-7], cont = T)

param_gb2.pg <- mle_gb2.pg$param
# param_gb2.pg.d <- mle_gb2.pg.d$param
cbind("a"=param_gb2.pg[1], "b"=param_gb2.pg[2], "p"=param_gb2.pg[3], "q"=param_gb2.pg[4],
      "alpha"=1/param_gb2.pg[5], "l"=param_gb2.pg[6]/param_gb2.pg[5], "w"=param_gb2.pg[7])
# rbind("a"=param_gb2.pg.d[1], "b"=param_gb2.pg.d[2], "p"=param_gb2.pg.d[3], "q"=param_gb2.pg.d[4],
#      "alpha"=1/param_gb2.pg.d[5], "l"=param_gb2.pg.d[6]/param_gb2.pg.d[5], "w"=param_gb2.pg.d[7])

# Graphique
plot(Fn, ylim=c(0, 1), xlim=c(0, 1e6), main = "")
curve(pgb2.pg(x, u, param_gb2.pg), add = T, col = "red", lwd = 2)

ad.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg))$s
cvm.test(perte, function(x) pgb2.pg(x, u, param_gb2.pg))$s
2*mle_gb2.pg$value - length(mle_gb2.pg$par) * log(n)

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

E.gb2.pg <- w/pbeta(h(u), p, q) * mgb2(1, a, b, p, q) * pbeta(h(u), p + 1/a, q - 1/a) + (1 - w) * (u + s/(1 - xi))

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

########################
### Processus de Poisson--------------------------------------------------------
########################
## Ajuster les données en format date
patron <- "%Y-%m-%dT%H:%M:%S"
don$Arrival.DtTm<- as.POSIXct(don$Arrival.DtTm, patron, tz = "UTC")
don$years <- as.factor(format(don$Arrival.DtTm, format = "%Y"))

## Sélectionne le minimum de montant de sinistre
min <- 0
annee <- subset(don, Estimated.Property.Loss > min, select = years)$years
date <- sort(subset(don, Estimated.Property.Loss > min, select = Arrival.DtTm)$Arrival.DtTm)

## Maximum de vraisemblance pour l'ensemble des données
X <- sapply(levels(annee), function(i) sum(annee == i))
nb_acc <- unname(table(date))

# Paramètre lambda journalier
lam <- sum(nb_acc)/total_jour

min(date)
max(date)
minDate <- as.Date("2002-12-31")
maxDate <- as.Date("2023-05-15")
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
t <- seq(as.Date("2003-01-01"), maxDate, by="days")
simul <- replicate(7, rpois(total_jour, lam))
for (i in 1:7) {
  matplot(t, cumsum(simul[,i]), type = "s", col = couleur[i], add=T)
}

## Temps entre les sinistres
dif <- diff(date)
tis <- as.numeric(dif)
Fc <- ecdf(tis)

## Maximum de vraisemblance pour une exponentielle
param_tis <- length(tis)/sum(tis)
unname(dif)

## graphique comparaison
plot(Fc, ylim=c(0, 1), xlim=c(0, 5e5), lwd = 2, xlab= "Heures", main="")
#title(main = "Comparaison de la fonction de répartition de la loi exponentielle \navec la fonction de répartition empirique pour le\n temps entre les sinsitres supérieur à 1 million")
curve(pexp(x, param_tis), col = "blue", lwd = 2, add = T)

ExpQQ(tis)

tis <- tis + rnorm(length(tis), 0, 0.00001)
ks.test(tis, function(x) pexp(x, param_tis))



