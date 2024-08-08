

test <- MEfit(perte)

mle.me <- SpliceFitGPD(perte, tsplice = u)

curve(pSplice(x, mle.me), xlim = c(0, 1e6))


plot(Fn, main = "Fonction de répartition empirique comparé avec lognormale",
     xlim = c(0, 5e6), ylim=c(0, 1))
curve(pSplice(x, mle.me), col = "red", lwd = 2, add = T)


ad.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
ad.test(perte, function(x) pSplice(x, mle.me))$s

cvm.test(perte, function(x) pln.pg(x, u, param_ln.pg))$s
cvm.test(perte, function(x) pSplice(x, mle.me))$s
