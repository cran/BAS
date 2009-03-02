data(Hald)
hald.gprior =  bas.lm(Y~ ., data=Hald, n.models=2^4, alpha=13,
                      prior="g-prior", initprobs="eplogp")

hald.gprior
summary(hald.gprior)
image(hald.gprior, subset=-1, vlas=0)

hald.coef = coefficients(hald.gprior)
hald.coef
plot(hald.coef)
predict(hald.gprior, hald.gprior$X, top=5)
fitted(hald.gprior, type="HPM")
hald.EB = update(hald.gprior, newprior="EB-global")
hald.bic = update(hald.gprior,newprior="BIC")
hald.zs = update(hald.bic, newprior="ZS-null")
