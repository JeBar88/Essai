## Package
library(ggplot2)
library(ggpubr)
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
library(forcats)
dgpd <- ReIns::dgpd
pgpd <- ReIns::pgpd
qgpd <- ReIns::qgpd
rgpd <- ReIns::rgpd

####################
### Fonction utiles -----------------------------------------------------------
####################
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
  ifelse(k <= w, qlnorm(pmin(k * plnorm(u, m, r)/w, 0.999), m, r), qgpd(pmax((k - w)/(1-w), 0), xi, u, sigma))
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

####################
### Analyse initiale -----------------------------------------------------------
####################
##### Base de données #####
don <- read.csv("cafc_nfid_incident_new.csv", stringsAsFactors = T)
perte <- subset(don, DOLLOSSA > 0, select = DOLLOSSA)$DOLLOSSA
don_share <- subset(don,
                    DOLLOSSA > 0 & !is.na(GENCONST), 
                    select = c("DOLLOSSA", "PROPCLAS", "PROPGRP", "GENCONST", 'RISKVALA'))


summary(don_share$DOLLOSSA)

which(don_share$DOLLOSSA == 2e8)

summary(don_share[, 5])[3]

sum(don$DOLLOSSA > 5e6, na.rm = TRUE)

### Ajuster les données en format date
don$temps <- paste(don$YEAR, don$MONTH, don$DATE, sep = "-")
don$temps <- paste0(don$temps, 'T', don$TIME)
patron <- "%Y-%m-%dT%H:%M:%S"
don$temps <- as.POSIXct(don$temps, patron, tz = "UTC")
min <- 0
date <- sort(na.omit(don$temps[don$DOLLOSSA > min]))
don$YEAR <- as.factor(don$YEAR)

don_graph <- na.omit(subset(don,
                    DOLLOSSA > 0 & !is.na(GENCONST), 
                    select = c("DOLLOSSA", "PROPCLAS", "PROPGRP", "GENCONST", 'temps')))
don_graph <- subset(don_graph, DOLLOSSA > 5e6)
don_graph$DOLLOSSA <- don_graph$DOLLOSSA/1e5



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

length(date)



lam <- sapply(perte_shar, function(i) length(i)/total_jour)
lam_an <- lam*365.25

sum(sapply(perte_shar, sum))

#### Répartition empirique ####
Fn1 <- ecdf(perte_shar[[1]])
Fn2 <- ecdf(perte_shar[[2]])
Fn3 <- ecdf(perte_shar[[3]])
Fn4 <- ecdf(perte_shar[[4]])
Fn5 <- ecdf(perte_shar[[5]])

##### Graphiques sinsitres #####
### Perte estimées selon la date
matplot(don_graph[don_graph$GENCONST == 1, 5],
        don_graph[don_graph$GENCONST == 1, 1], type = 'h', ylim = c(0, 1000),
        ylab = "Pertes estimées (100 milles $)", xlab = 'Date',
        main = 'Pertes estimées supérieures à 5 millions selon la date du sinistre', lwd = 2)
matplot(don_graph[don_graph$GENCONST == 2, 5],
        don_graph[don_graph$GENCONST == 2, 1], type = 'h', col = 'red', add = TRUE, lwd = 2)
matplot(don_graph[don_graph$GENCONST == 3, 5],
        don_graph[don_graph$GENCONST == 3, 1], type = 'h', col = 'blue', add = TRUE, lwd=2)
matplot(don_graph[don_graph$GENCONST == 4, 5],
        don_graph[don_graph$GENCONST == 4, 1], type = 'h', col = 'orange', add = TRUE, lwd = 2)
matplot(don_graph[don_graph$GENCONST == 5, 5],
        don_graph[don_graph$GENCONST == 5, 1], type = 'h', col = 'darkgreen', add = TRUE, lwd=2)
legend('topright', legend = c('Bois non-protéger', 'Bois protéger', 'Bois massif', 'Acier non-protéger', 'Acier protéger'),
       col = c('black', 'red', 'blue', 'orange', 'darkgreen'), lwd = 2)

##### Graphique contribution
theme_set(theme_pubr())
vk <- c(0.01, 0.025, 0.05, 0.075, 0.1)
g <- sapply(vk, function(i) sum(sort(perte, dec = TRUE)[1:(floor(i*length(perte)))])/sum(perte))
contrib <- data.frame(vk, g)
ggplot2::ggplot(contrib, aes(x=factor(vk), y=g)) + 
  geom_col(fill = 'royalblue') + 
  ylim(0, 1) +
  geom_abline(intercept = 0.5, slope = 0, linetype = 'dashed') +
  labs(y="Proportion of the total losses", x="Percentage of observations",
       title="Proportion of the largest losses to the total losses") +
  theme_pubclean()

### Pertes estimées totals selon le portefeuille
perte_tot_port <- aggregate(DOLLOSSA~GENCONST, FUN = sum, data = don_share)[2:6, ]

ggplot2::ggplot(perte_tot_port, aes(x=factor(GENCONST), y=DOLLOSSA/1e5)) +
  geom_bar(stat="identity", fill = 'royalblue') +
  geom_text(aes(label= round(DOLLOSSA/(1e5 * sum(DOLLOSSA/1e5)), 2)), vjust=-0.5) +
  labs(title = 'Pertes estimées totales selon le portefeuille', x = 'Portefeuille', y = 'Pertes estimées total (100 milles $)') + 
  theme_pubclean()

### Grosses pertes
don_perte_shar <- subset(don_share, GENCONST %in% 1:5, select = c('DOLLOSSA', 'GENCONST'))
don_perte_shar <- don_perte_shar[order(don_perte_shar$DOLLOSSA, decreasing = TRUE), ]

perte_gros <- don_perte_shar[1:100, ]

xtable(t(table(don_perte_shar$GENCONST)))

# ggplot2::ggplot(perte_gros, aes(GENCONST)) + 
#   geom_bar(fill = 'royalblue') + 
#   labs(title = 'Nombre de gros sinistres dans chaque portefeuille',
#        x='Portefeuille', y='Nombre') +
#   theme_pubclean()

perte_tot_gros <- aggregate(DOLLOSSA~GENCONST, FUN = sum, data=perte_gros)

ggplot2::ggplot(perte_tot_gros, aes(x=factor(GENCONST), y=DOLLOSSA/1e5)) +
  geom_bar(stat="identity", fill = 'royalblue') +
  geom_text(aes(label= round(DOLLOSSA/sum(DOLLOSSA), 2)), vjust=-0.5) +
  labs(title = 'Pertes estimées totales selon le portefeuille', subtitle = 'Les 100 plus grosses pertes',
       x = 'Portefeuille', y = 'Pertes estimées total (100 milles $)') + 
  theme_pubclean()

meplot(perte_shar[[1]])
MEplot(perte_shar[[1]], xlim=c(0, 100), ylim=c(1, 100))
abline(v=3)

round(sapply(perte_shar, length)/length(unlist(perte_shar)), 3)

aggregate(DOLLOSSA~GENCONST, FUN = length, data=perte_gros)$DOLLO/100

### Effet detecteur de fumé et système d'arroseur
don_DET_SPRIN <- subset(don, DOLLOSSA > 1e6 & RISKVALA > 0 & SPRINPRO != 0 & FIREDET != 0 & SPRINPRO != 8 & SPRINPRO != 9,
                        select = c("DOLLOSSA", 'SPRINPRO', 'FIREDET', 'RISKVALA'))
dim(don_DET_SPRIN)

don_DET_SPRIN$SPRINPRO <- as.factor(don_DET_SPRIN$SPRINPRO)
don_DET_SPRIN$FIREDET <- as.factor(don_DET_SPRIN$FIREDET)

levels(don_DET_SPRIN$FIREDET)[1:9] <- 1
levels(don_DET_SPRIN$FIREDET)[2] <- 0

levels(don_DET_SPRIN$SPRINPRO)[1:6] <- 1
levels(don_DET_SPRIN$SPRINPRO)[2] <- 0

sprin_mean <- aggregate(DOLLOSSA ~ SPRINPRO, FUN = mean, data=don_DET_SPRIN)
sprin_len <- table(don_DET_SPRIN$SPRINPRO)

round(sprin_len/sum(sprin_len), 3)

xtable(cbind(sprin_mean[, 2], sprin_len), digit = 0)

test <- aggregate(RISKVALA ~ SPRINPRO, FUN = mean, data=don_DET_SPRIN)

sprin_mean[ , 2]/test[, 2]

### Analyse des type des batiment
don_type <- subset(don, DOLLOSSA > 0 & RISKVALA > 0 & GENCONST %in% 1:6 & PROPGRP != 0,
                   select = c('DOLLOSSA', 'RISKVALA', 'GENCONST', 'PROPGRP', 'PROPCLAS'))

dim(don_type)

PROPGRP_pm <- aggregate(DOLLOSSA ~ PROPGRP, FUN = sum, data=don_type)
PROPGRP_vm <- aggregate(RISKVALA ~ PROPGRP, FUN = mean, data=don_type)
GENCONST_pm <- aggregate(DOLLOSSA ~ GENCONST, FUN = sum, data=don_type)
GENCONST_vm <- aggregate(RISKVALA ~ GENCONST, FUN = mean, data=don_type)

## Selon le secteur d'utilité
ggplot2::ggplot(PROPGRP_pm , aes(x=factor(PROPGRP), y=DOLLOSSA/1e9)) +
  geom_bar(stat="identity", fill = 'royalblue') +
  geom_text(aes(label= round(DOLLOSSA/sum(DOLLOSSA), 2)), vjust=-0.5) +
  labs(title = "Pertes totales selon le secteur d'utilité du bâtiment",
       x = "Secteur d'utilité du bâtiment", y = 'Pertes totales (milliards $)') + 
  theme_pubclean()

ggplot2::ggplot(don_type, aes(x=factor(PROPGRP))) +
  geom_bar(fill = 'royalblue') +
  # geom_text(aes(label= ..count..), vjust=-0.5) +
  labs(title = "Nombre d'incendie selon le secteur d'utilité du bâtiment",
       x = "Secteur d'utilité du bâtiment", y = "Nombre d'incendie") + 
  theme_pubclean()

ggplot2::ggplot(GENCONST_pm, aes(x=factor(GENCONST), y=DOLLOSSA/1e9)) +
  geom_col(fill = 'royalblue') +
  geom_text(aes(label= round(DOLLOSSA/sum(DOLLOSSA), 2)), vjust=-0.5) +
  labs(title = "Pertes totale selon le type de construction",
       x = "Type de construction", y = 'Pertes moyennes (milliards $)') + 
  theme_pubclean()

ggplot2::ggplot(GENCONST_vm, aes(x=factor(GENCONST), y=RISKVALA/1e6)) +
  geom_col(fill = 'royalblue') +
  geom_text(aes(label= round(RISKVALA/sum(RISKVALA), 2)), vjust=-0.5) +
  labs(title = "Nombre d'incendie selon le type de construction",
       x = "Type de construction", y = "Nombre d'incendie") + 
  theme_pubclean()

## Analyse gros sinistres
don_gros <- subset(don, DOLLOSSA > 0 & GENCONST %in% 5 & PROPGRP != 0,
                   select = c('DOLLOSSA', 'RISKVALA', 'GENCONST', 'PROPGRP', 'PROPCLAS'))

don_gros <- don_gros[order(don_gros$DOLLOSSA, decreasing = TRUE), ]
don_gros[, 1:2] <- don_gros[, 1:2]/1e6
print(xtable(don_gros[1:10, -3], digits = 1), include.row = FALSE)



### Analyse construction
don_const <- subset(don, DOLLOSSA > 0 & GENCONST %in% 1:6 & PROPCLAS == 8310,
                   select = c('DOLLOSSA', 'RISKVALA', 'GENCONST', 'PROPCLAS'))

const_tot <- aggregate(RISKVALA~GENCONST, FUN = mean, data = don_const)

ggplot2::ggplot(const_tot, aes(x=factor(GENCONST), y=RISKVALA/1e6))+
  geom_col(fill = 'royalblue') +
  # geom_text(aes(label= round(DOLLOSSA/sum(DOLLOSSA), 3)), vjust=-0.5) + 
  labs(title = "Valeur moyenne du bâtiments selon le type de construction",
       subtitle = 'Bâtiments en construction',
       x = "Type de construction", y = "Valeurs moyennes (millions $)") + 
  theme_pubclean()

### Analyse répartition des sinistres selon PROPGRP
don_type_facility <- subset(don, DOLLOSSA > 0 & PROPGRP %in% 1000:8000 & GENCONST %in% 1:5,
                            select = c("PROPGRP", "DOLLOSSA", 'GENCONST'))
don_type_facility$PROPGRP <- as.factor(don_type_facility$PROPGRP)
don_type_facility$GENCONST <- as.factor(don_type_facility$GENCONST)

ggplot2::ggplot(don_type_facility, aes(fct_infreq(PROPGRP)), fill = GENCONST) +
  geom_bar(fill = 'royalblue') +
  geom_text(stat='count', aes(label= round(..count../sum(..count..), 2)), vjust=-0.5) +
  labs(title = "Nombre d'incendies selon le secteur d'utilité du bâtiment",
       x="Secteur d'utilité du bâtiment",  y = "Nombre d'incendies") +
  scale_x_discrete(labels=c('Résidentiel', 'Divers', 'Stockage',
                            'Assamblage', 'Commercial', 'Industrielle',
                            'Service', 'Institutionelle')) +
  theme_pubclean()


#### Modélisation ######### MGENCONSTodélisation #####
### POT
ParetoQQ(perte_shar[[1]], xlim=c(0, 60))
abline(h=log(5e5*div))
th <- 3e5*div

## Mean excess plot
par(mfrow=c(2, 3), mar=c(3, 2, 3, 2))
for (i in seq_along(perte_shar)) {
  MEplot(perte_shar[[i]], xlim=c(0, 30), ylim=c(1, 30))
  abline(v=3)
}

## Portion Pareto généralisée
mle_pg <- lapply(perte_shar, gpdFit, threshold = th)
par_pg <- sapply(seq_along(mle_pg), function(i) unname(mle_pg[[i]]$par.est))
xi <- par_pg[2,]
b <- par_pg[1,]
a <- 1/xi
a

sapply(perte_shar, function(i) sum(i > th))

### Parametre de départ
# Paramètres estimés
parInit <- c(2, 8, xi[1], b[1], 0.9)

### Raccordement
mle_ln.pg <- lapply(seq_along(perte_shar),
                     function(i) Spl.ln.pg(perte_shar[[i]], th, parInit))

par_ln.pg <- sapply(seq_along(mle_ln.pg), function(i) unname(mle_ln.pg[[i]]$par))

1/par_ln.pg[3, ]

### Vérification graphique
par(mfrow=c(2, 3), mar=c(3, 2, 3, 2))
for (i in seq_along(perte_shar)) {
  plot(get(paste0('Fn', i)), xlim = c(0, 1e2), ylim = c(0.90, 1), lwd = 2, 
       main = paste('Portefeuille', i), xlab = '', ylab = '')
  curve(pln.pg(x, th, par_ln.pg[, i]), add = TRUE, col = 'royalblue', lwd=2)
}


#################
### Risk Sharing ---------------------------------------------------------------
#################
##### Coditional mean #####
### Avec discrétisation
nfft <- 2^24
h <- 0.01
m <- 1e5
k <- (0:(nfft - 1))
u <- c(0.9, 0.95, 0.99, 0.995, 0.999)
th <- 3
vk <- seq_along(lam)

# Discrétisé les v.a. B (sévérité)
fb <- vector(mode = 'list', length(perte_shar))
for (i in seq_along(perte_shar)) {
  fb[[i]] <- discretize(pln.pg(x, th, par_ln.pg[, i]), 0, m, h, method = 'lower')
  fb[[i]] <- c(fb[[i]], numeric(nfft - length(fb[[i]])))
}
# fb <- readRDS('fb.RDS')

EX <- sapply(vk, function(i) sum(k*fb[[i]]))*h*lam

Fb <- lapply(fb, cumsum)

## Calculer la FFT de B
fbt <- lapply(fb, fft)

## Calculer la FFT de X
fgp_po <- function(t, lam) exp(lam*(t - 1))
fxt <- lapply(vk, function(i) fgp_po(fbt[[i]], lam[i]))

## Calculer fs et Fs
fst <- Reduce('*', fxt)
# fst <- readRDS('fst.RDS')
fs <- Re(fft(fst, TRUE))/nfft
#fs <- readRDS('fs.RDS')
Fs <- cumsum(fs)
VaRS <- sapply(u, function(i) k[min(which(Fs >= i))])
TVaRS <- sapply(seq_along(u), function(i) VaRS[i] + sum(pmax(k - VaRS[i], 0)*fs)/(1 - u[i]))
cbind(VaRS, TVaRS)*h

ES <- sum(fs*k)*h

## Calculer la FFT de phi et e1
phit <- lapply(vk, function(i) fft((k + 1)*c(fb[[i]][-1], 0)))
# phit <- readRDS('phit.RDS')
e1 <- exp((-2i*k*pi)/nfft)


## Calculer la conditional expectation
mut <- lapply(vk, function(i) lam[i] * e1 * phit[[i]] * fst)
mu <- lapply(vk, function(i) Re(fft(mut[[i]], TRUE))/nfft)
cm <- lapply(mu, '*', h/fs)
cm_tout <- Reduce('+', cm)

##### risque non identique #####
s <- 5e5*div/h
1 - Fs[s + 1]

min(which(Fs >= 0.995))

k[11756]*h/div

EX <- sapply(fb, function(f) sum(f * k)*h)*lam
cm_res <- sapply(vk, function(i) cm[[i]][s + 1])
pm_res <- sapply(vk, function(i) EX[i]/ES * s*h)

lol <- seq(5e5, 1e9, 5e5)
test1 <- sapply(lol*div/h, function(s) cm[[1]][s + 1])/(lol*div)
test2 <- sapply(lol*div/h, function(s) cm[[2]][s + 1])/(lol*div)
test3 <- sapply(lol*div/h, function(s) cm[[3]][s + 1])/(lol*div)
test4 <- sapply(lol*div/h, function(s) cm[[4]][s + 1])/(lol*div)
test5 <- sapply(lol*div/h, function(s) cm[[5]][s + 1])/(lol*div)

par(mar = c(5.1, 5.1, 1.1, 2.1))
plot(lol/1e6, test1*100, type='l', lwd=2, ylim=c(0, 60), 
     xlab='s (millions)', ylab=TeX('$C^{CM}_j$ (%)'))
for (i in 2:5) {
  matplot(lol/1e6, get(paste0('test', i))*100, type='l', lwd=2, col=i, add = TRUE)
}
legend('topright', legend = c('CC', 'PCC', 'HTC', 'NCC', 'PNCC'), col = 1:5, lwd = 2)

xtable(rbind(cm_res, pm_res), digits = 2)

summary(don_share$DOLLOSSA)

median(na.omit(don_share$RISKVALA))

aggregate(RISKVALA~GENCONST, FUN=median, don_share)

### Grapgique des contributions
vk <- as.integer(seq(1e6, 1.5e9, 1e5)*div/h)
cm_val <- t(sapply(vk, function(s) sapply(1:5, function(i) cm[[i]][s + 1])))

matplot(vk*h, cm_val[, 1], type = 'l', lwd = 2, ylim = c(0, 1300), xlim = c(0, 3000),
        xlab = "Perte totale (100 milles $)", ylab = "Contribution (100 milles $)")
for (i in 2:5) {
  matplot(vk*h, cm_val[, i], type = 'l', lwd = 2, col = i, add = TRUE)
}
legend('topleft', legend = paste('Port.', 1:5),
       col = 1:5, lwd = 2)
# abline(v=ES, lty = 2, lwd = 2)

### Calcule de la pente (linéarité ?)
bu <- round(k[min(which(Fs >= 0.9999))]*h, -2)
bl <- 300
delta <- (cm_val[which(vk == bu/h), ] - cm_val[which(vk == bl/h), ])/(bu - bl)

cm_delta_val <- t(sapply(vk, function(s) sapply(1:5, function(i) (delta[i]*s*h))))
pm_val <- t(sapply(vk, function(s) sapply(1:5, function(i) (EX[i]/ES*s*h))))

j <- 1
matplot(vk*h, cm_val[, j], type = 'l', lwd = 2, xlim = c(0, 1500), ylim = c(0, 700),
        xlab = "Perte totale (100 milles $)", ylab = "Contribution (100 milles $)")
matplot(vk*h, cm_delta_val[, j], add = TRUE, type = 'l', col = 2, lwd = 2)
matplot(vk*h, pm_val[, j], add = TRUE, type = 'l', col = 3, lwd = 2)
legend('topleft', legend = c('cm', 'cm_delta', 'pm'),
       col = 1:3, lwd = 2)

matplot(vk*h, cm_val[, 5], type = 'l', lwd = 2, ylim = c(0, 300), xlim = c(0, 1500),
        xlab = "Perte totale (100 milles $)", ylab = "Contribution (100 milles $)")
matplot(vk*h, cm_delta_val[, 5], add = TRUE, type = 'l', col = 2, lwd = 2)
matplot(vk*h, pm_val[, 5], add = TRUE, type = 'l', col = 3, lwd = 2)

838.36/div

### Graphique répartition
vk <- seq(1, 1e5, 1)
matplot(k[vk]*h, Fs[vk], xlim = c(0, 100), type = 'l', lwd = 2, ylim = c(0, 1))
for (i in seq_along(fb)) {
  matplot(k[vk]*h, Fb[[i]][vk], type = 'l', lwd = 2, add = TRUE, col = i + 1)
}
abline(v=ES, lwd = 2, lty = 2)


# lim <- 1:1e6
# plot(k[lim], cm_tout[lim], type= 'l', xlim = c(0, 1e6), ylim=c(-1000, 2e5))
# 

### Contribution à la TVaR
CVaRX <- sapply(seq_along(fb), function(i) cm[[i]][VaRS + 1])

CTVaRX <- sapply(seq_along(u), 
                 function(j) sapply(seq_along(fb),
                 function(i) (sum(mu[[i]][k > VaRS[j]]) + CVaRX[j, i] * (Fs[max(which(k <= VaRS[j]))] - u[j]))/(1 - u[j])))*h


xtable(round(t(CTVaRX)/(TVaRS*h), 2))



# matplot(u, CTVaRX[1, ], type = 'l', lwd=2, ylim = c(5, 150), xlim = c(0.95, 1),
#         ylab = TeX('$CTVaR(X_i,S)$'), xlab = TeX('$\\kappa$'))
# matplot(u, CTVaRX[2, ], type = 'l', lwd=2, col = 'red', add = TRUE)
# matplot(u, CTVaRX[3, ], type = 'l', lwd=2, col = 'blue', add = TRUE)
# matplot(u, CTVaRX[4, ], type = 'l', lwd=2, col = 'orange', add = TRUE)
# matplot(u, CTVaRX[5, ], type = 'l', lwd=2, col = 'darkgreen', add = TRUE)
# legend('topleft', legend = c('Bois non-proteger', 'Bois proteger', 'Bois massif', 'Acier non-proteger', 'Acier proteger'), 
#        col = c('black', 'red', 'blue', 'orange', 'darkgreen'), lwd = 2)


E_ln.pg <- function(par){
  mu <- par[1] ; r <- par[2] ; xi <- par[3]
  s <- par[4] ; w <- par[5] 
  w*mlnorm(1, mu, r)*pnorm((log(th) - mu - r^2)/r)/pnorm((log(th) - mu)/r) + (1-w)*(s/(1 - xi) + th)
}

sum(lam*sapply(1:5, function(i) E_ln.pg(par_ln.pg[, i])))


### Avec simulation
set.seed(2024)
m <- 1e6
X <- sapply(seq_along(perte_shar),
             function(i) rcomppois(m, lam[i], rln.pg(th, par_ln.pg[, i])))

S <- rowSums(X)
Ss <- sort(S)

rea <- cbind(round(X, 2), rowSums(round(X, 2)))
head(rea)

cm_simul <- t(sapply(rea[, 6], function(s) sapply(vk, function(i) cm[[i]][s/h + 1])))
cm_simul <- round(cm_simul, 2)

head(rea[, -6])
head(cm_simul)

xtable(cbind(head(rea[, 5], 5), head(cm_simul[, 5], 5)))

colMeans(head(rea[, -6], 5))

colMeans(rea[, -6])
colMeans(cm_simul)

xtable(rbind(colMeans(rea[, -6]), colMeans(cm_simul)), digit = 2)


matplot(1:m, cumsum(rea[, 1]), type = 'l', lwd = 2, pch = 1)
matplot(1:m, cumsum(cm_simul[, 1]), type='l', add = TRUE, col = 2, lwd = 2, pch = 1)


### Contribution par simulation
VaR <- function(u) Ss[u*ms]
TVaR <- function(u) mean(S[S > VaR(u)])

VaR(u)
sapply(u, TVaR)

CTVaR <- function(u, port) mean(port * I(S > VaR(u)))/(1-u)
sapply(u, function(i) CTVaR(i, X[, 1]))

CTVaRXi <- sapply(u, function(j) sapply(1:5, function(i) CTVaR(j, X[, i])))
colSums(CTVaRXi)

CTVaRXi[, 3]

##### Caractéristique Esp, median, TVaR #####
u <- c(0.95, 0.99, 0.999)
VaR_shar <- function(u, don) quantile(don, u)
TVaR_shar <- function(u, don){
  mean(don[don > VaR_shar(u, don)])
}

### Summary
Me_Med <- sapply(perte_shar, summary)[3:4, ]
max_shar <- sapply(perte_shar, function(i) tail(sort(i), 3))

### VaR
VaR <- sapply(perte_shar, VaR_shar, u=u)

info <- rbind(Me_Med, VaR, max_shar)

xtable(info)

##### Analyse gros sinistres #####

# Sample dataframe
df <- data.frame(
  Name = c("John", "Alice", "Bob", "Emily"),
  Age = c(25, 30, 22, 35),
  Score = c(80, 90, 75, 85)
)

# Sort dataframe by the 'Age' column
sorted_df <- df[order(df$Age), ]

# Print sorted dataframe
print(sorted_df)



