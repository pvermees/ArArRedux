ns <- 10
relerr <- 0.02
Ri <- c(0.9847, 1.0033, 0.9996, 1.0335, 0.9953,
        1.0072, 0.9625, 1.0046, 1.0006, 1.0314)
#Ri <- rnorm(ns,mean=1,sd=err)
sRi <- Ri*relerr
J <- 1
sJ <- J*relerr
lambda <- 0.00055305

ti <- log(1+J*Ri)/lambda
sti <- sqrt((Ri*sJ)^2 + (J*sRi)^2)/(lambda*(1+J*Ri))

ERJ <- diag(c(sRi,sJ)^2)
dtdRi <- (1/lambda)*(J/(1+J*Ri))
dtdJ <- (1/lambda)*(Ri/(1+J*Ri))
Jac <- matrix(0,ns,ns+1)
Jac[cbind(1:ns,1:ns)] <- dtdRi
Jac[,ns+1] <- dtdJ
Et <- Jac %*% ERJ %*% t(Jac)
sti2 <- sqrt(diag(Et))

tbarwrong <- sum(ti/sti^2)/sum(1/sti^2)
stbarwrong <- sqrt(1/sum(1/sti^2))
MSWDwrong <- sum(((ti-tbarwrong)/sti)^2)/(ns-1)

Rbar <- sum(Ri/sRi^2)/sum(1/sRi^2)
sRbar <- 1/sum(1/sRi^2)
tbar <- log(1+J*Rbar)/lambda
stbar <- sqrt((Rbar*sJ)^2 + (J*sRbar)^2)/(lambda*(1+J*Rbar))

vartbar <- solve(matrix(1,1,ns) %*% solve(Et) %*% matrix(1,ns,1))
stbar2 <- sqrt(vartbar)
tbar2 <- vartbar*( matrix(1,1,ns) %*% solve(Et) %*% matrix(ti,ns,1) )

MSWD <- sum(((Ri-Rbar)/sRi)^2)/(ns-1)
dt <- matrix(ti,1,ns) - matrix(tbar2,1,ns)
MSWD2 <- (dt %*% solve(Et) %*% t(dt))/(ns-1)

print(MSWD)
