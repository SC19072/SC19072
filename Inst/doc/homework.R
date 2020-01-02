## -----------------------------------------------------------------------------
set.seed(1025)
par(mar=c(1,1,1,1))
par(mfrow=c(2,2))
n<-1.025e4;
sigma<-c(1,3,4,6);
for (i in 1:4)
{F<-runif(n)
x<-sqrt(-2*(sigma[i])^2*log(F))     ## F=exp(-x^2/(2*sigma^2))
hist(x,prob=TRUE,breaks=20,col="red",main=expression(f(x)==x/sigma^2*exp(-x^2/(2*sigma^2))))
y<-seq(0,20,0.001)
## Fit the curve of the theoretical model.
lines(y,y/(sigma[i])^2*exp(-y^2/(2*(sigma[i])^2)),col="blue") }


## -----------------------------------------------------------------------------
set.seed(1225)
par(mar=c(1,1,1,1))
par(mfrow=c(3,3))
n<-1000
p1<-c(0.75,0.7,0.6,0.5,0.45,0.4,0.35,0.2,0.1)
p2<-1-p1
for (i in 1:9)
{
n1<-rbinom(1,n,p1[i]) 
##Generate n1 samples based on binomial distribution.
n2<-n-n1
x1<-rnorm(n1,0,1)
x2<-rnorm(n2,3,1)
x<-c(x1,x2)                ## combine x1 with x2
hist(x,prob=TRUE,breaks=15,col="brown",main=p1[i])
y<-seq(-4,6,0.01)
lines(y,p1[i]*(sqrt(2*pi))^(-1)*exp(-y^2/2)+p2[i]*(sqrt(2*pi)^(-1))*exp(-(y-3)^2/2),col="blue")
}

## -----------------------------------------------------------------------------
set.seed(1234)
miracle<-function(n,sigma)
{
d<-dim(sigma)
##Generate  a zero matrix T based on the dimesion of matrix sigma.
T<-matrix(0,nrow=d[1],ncol=d[2]) 
## Cholesky decomposition of matirx sigma produces an upper triangular matrix, so we need to transpose it to get the lower triangular matrix L. 
L<-t(chol(sigma))
## Construct a matrix T.
 for (i in 1:d[1] )
   {for (j in 1:d[2]) 
         {if (i>j)
            {T[i,j]<-rnorm(1,0,1)
             }
             else if (i<j)
               {
                   T[i,j]<-0
                 }
           else
                {
                    T[i,j]<-sqrt(rchisq(1,n-i+1))
                  }
           
   }

 }
A<-T%*%t(T)         
S<-L%*%A%*%t(L)
S
}
n<-10
sigma<-matrix(c(1,0.1,0.3,0.1,1,0.2,0.3,0.2,1),nrow=3)
x<-miracle(n,sigma)
x

## -----------------------------------------------------------------------------
set.seed(12345678)
n<-3e6;
x<-runif(n,min=0,max=pi/3);
gamma.hat<-mean(sin(x))*pi/3            ## the value of the estimation
print(c(gamma.hat,-cos(pi/3)+cos(0)))   ## the value of the integration

## -----------------------------------------------------------------------------
set.seed(123456)
theta <- function(n =3e4, antithetic = TRUE) {
  a<- runif(n/2,0,1)
  if (!antithetic) b<-runif(n/2) else b<-1-a
  u <- c(a, b)
  theta<-numeric(1)
  g <-exp(-u)/(1+u^2)
  theta<- mean(g)
  theta
}

m <- 2e3
theta1 <- theta2 <- numeric(m)
for (i in 1:m) 
{
  theta1[i] <-theta(n =3e3, anti=FALSE)
  theta2[i] <-theta(n =3e3, anti=TRUE)
}
theta1.hat<-mean(theta1);        ## the estimation of integration without antithetic variables
theta2.hat<-mean(theta2);        ## the estimation of integration with antithetic variables
var1<-var(theta1);
var2<-var(theta2);
var12<-var2-var1;                 ## the difference of two estimations' variance
data.frame(theta1.hat,theta2.hat,var1,var2,var12)

## -----------------------------------------------------------------------------
set.seed(1025)
N<-1e4;
K<-5;                ## the number of stratum
re<-N/K;
vector<-rep(0,K)
mir<- matrix(0, 50, 2)
g<-function(x)
  {exp(-x-log(1+x^2))}   ## g(x)
fg<-function(x)  
{
 g(x)/((exp(-x)/(1-exp(-1))))*(x<1)*(x>0)  ## fg(x)=g(x)/f(x)
}

for (i in 1:50) {
  for(j in 1:K)
    {
    Y<-runif(re,j-1,j)        ## sample from five equal interval between 0 and 5
    x<--log(1-Y/5*(1-exp(-1)))        ## Inverse transformation for X
    vector[j]<-mean(fg(x))          ## the estimation of per staturm
    }
    mir[i, 1]<-mean(vector)
    mir[i, 2]<-sd(vector)
}
apply(mir,2,mean)

## -----------------------------------------------------------------------------
set.seed(1025)
n<-20
k<-3.65e3
N<-numeric(1)
XMLL<-XMUL<-numeric(k)
alpha<-0.05
for (i in 1:k)
{
x<-rchisq(n,2);
XMLL[i]<-mean(x)-qt(1-alpha/2,n-1)*sd(x)/sqrt(n)  
## the lower confidence limit of estimated parameter
XMUL[i]<-mean(x)+qt(1-alpha/2,n-1)*sd(x)/sqrt(n)
## the upper coonfidence limit of estimated parameter
if (XMLL[i]<=2 & XMUL[i]>=2)
 { N<-N+1}              ## the number of times that A falls into the confidence interval 
else
 { N<-N}
}
pvalue<-N/k             ## coverage probability 
pvalue

## -----------------------------------------------------------------------------
set.seed(1025)
n<-1.3e3;
m<-1.3e3;
c<-c(0.025,0.05,0.95,0.975);
a<-a1<-a2<-a3<-sk<-numeric(n)
rq<-simq<-b<-sdsimq<-f<-numeric(4)
for (j in 1:n)
{
  x<-rnorm(m,0,sqrt(6/n))
  a1[j]<-mean(x)
  a2[j]<-mean((x-a1[j])^3)
  a3[j]<-(mean((x-a1[j])^2))^1.5
  sk[j]<-a2[j]/a3[j]
}
a<-sk[order(sk)]      
## order the estimated skewnesses we get from the lowest to the highest
for (i in 1:4)
{
  rq[i]<-qnorm(c[i],0,sqrt(6/n))
  b[i]<-c[i]*n                                    
  simq[i]<-a[b[i]]                               ## the quantiles of given normal distribution
  f[i]<-1/sqrt(2*pi*6/n)*exp(-(rq[i])^2/(2*6/n)) ## estimated quantiles
  sdsimq[i]<-sqrt((c[i]*(1-c[i]))/(n*(f[i])^2))  ## the standard error of estimated quantiles
}
data.frame(rq,simq,sdsimq)

## -----------------------------------------------------------------------------
set.seed(1101)
n<-350
m<-1200
alpha<-c(seq(0.1,2,0.1),seq(1,73),2)
num1<-betaskew<-matrix(0,length(alpha),m)
p.reject<-cv<-numeric(length(alpha))
for (i in 1:length(alpha))
{
  for (j in 1:m)
{
  x<-rbeta(n,alpha[i],alpha[i])
  betaskew[i,j]<-mean((x-mean(x))^3)/((mean((x-mean(x))^2))^1.5)    
  ## the skewness of random statistic from beta distribution with alpha and beta equal 0.05
  cv[i]<-qnorm(0.975, 0,sqrt(6*(n-2)/((n+1)*(n+3))))                     
  ## critical value
  num1[i,j]<-as.integer(abs(betaskew[i,j])>=cv[i])       #test decision is 1 (reject) or 0)
 
}
p.reject[i]<-mean(num1[i,])
## the probability that the skewness falls in the reject domain for one loop
}
plot(alpha,p.reject,type="b",lty=5,pch=6,col="blue")
abline(h=.05,lty=4,col="red")

## -----------------------------------------------------------------------------
set.seed(1101)
n<-450
m<-1200
df<-seq(2,50,1)
num1<-betaskew<-matrix(0,length(df),m)
p.reject<-cv<-numeric(length(df))
for (i in 1:length(df))
{
  for (j in 1:m)
{
  x<-rt(n,df[i])
  betaskew[i,j]<-mean((x-mean(x))^3)/((mean((x-mean(x))^2))^1.5)    
  ## the skewness of random statistic from beta distribution with alpha and beta equal 0.05
  cv[i]<-qnorm(0.975, 0,sqrt(6*(n-2)/((n+1)*(n+3))))                     
  ## critical value
  num1[i,j]<-as.integer(abs(betaskew[i,j])>=cv[i])       #test decision is 1 (reject) or 0)
 
}
p.reject[i]<-mean(num1[i,])
## the probability that the skewness falls in the reject domain for one loop
}
plot(df,p.reject,type="b",lty=6,pch=8,col="red")
abline(h=.05,lty=3)

## -----------------------------------------------------------------------------
set.seed(1025)
n<-800
m<-1200
epsilon<-c(seq(0,0.1,0.01),seq(0.1,1,0.05))
num1<-betaskew<-matrix(0,length(epsilon),m)
pw<-cv<-numeric(length(epsilon))
for (i in 1:length(epsilon))
{
  for (j in 1:m)
{
  alpha<-sample(c(2,25),replace = TRUE,size=n,prob=c(1-epsilon[i],epsilon[i]))
  x<-rbeta(n,alpha,alpha)
  betaskew[i,j]<-mean((x-mean(x))^3)/((mean((x-mean(x))^2))^1.5)    
  ## the skewness of random statistic from beta distribution with alpha and beta equal 0.05
  cv[i]<-qnorm(0.975, 0,sqrt(6/n))                    
  ## critical value
  num1[i,j]<-as.integer(abs(betaskew[i,j])>=cv[i])       #test decision is 1 (reject) or 0)
 
}
pw[i]<-mean(num1[i,])
## the probability that the skewness falls in the reject domain for one loop
}
plot(epsilon,pw,type="b",lty=12,pch=    )
se <- sqrt(pw*(1-pw)/m)                #add standard errors
lines(epsilon, pw+se, lty = 6,col="blue")
lines(epsilon, pw-se, lty = 6,col="blue")
abline(h=.05,lty=2,col="brown")

## -----------------------------------------------------------------------------
set.seed(1101)
n<-250
m<-1200
epsilon<-c(seq(0,0.1,0.01),seq(0.1,1,0.05))
num1<-betaskew<-matrix(0,length(epsilon),m)
pw<-cv<-numeric(length(epsilon))
for (i in 1:length(epsilon))
{
  for (j in 1:m)
{
  df<-sample(c(2,50),replace = TRUE,size=n,prob=c(1-epsilon[i],epsilon[i]))
  x<-rt(n,df)
  betaskew[i,j]<-mean((x-mean(x))^3)/((mean((x-mean(x))^2))^1.5)    
  ## the skewness of random statistic from beta distribution with alpha and beta equal 0.05
  cv[i]<-qnorm(0.975, 0,sqrt(6*(n-2)/((n+1)*(n+3))))                    
  ## critical value
  num1[i,j]<-as.integer(abs(betaskew[i,j])>=cv[i])       #test decision is 1 (reject) or 0)
 
}
pw[i]<-mean(num1[i,])
## the probability that the skewness falls in the reject domain for one loop
}
plot(epsilon,pw,type="b",lty=4,pch=8,col="red",ylim =c(0,1))
se <- sqrt(pw*(1-pw)/m) #add standard errors
lines(epsilon, pw+se, lty = 6)
lines(epsilon, pw-se, lty = 6)
abline(h=.05,lty=2,col="blue")

## -----------------------------------------------------------------------------
set.seed(1025)
alpha<-c(0.01,0.025,0.05,0.1)
n<-1200                         ## the number of sample for one replicate
m<-1.35e3;                       ## the number of repilcates  
p1<-p2<-p3<-numeric(m)
p.hat1<-p.hat2<-p.hat3<-0
for (i in 1:m)
{
 x1<-rchisq(n,1) 
 x2<-runif(n,0,2)
 x3<-rexp(n,1)
 ## the means of above three distributions all are one.
 ttest1<-t.test(x1,alternative="two.side",mu=1)          
 ## the t-test of chi-square distribution
 ttest2<-t.test(x2,alternative="two.side",mu=1)  
 ## the t-test of uniform distribution
 ttest3<-t.test(x3,alternative="two.side",mu=1)
 ## the t-test of exponetion distribution
 p1[i]<-ttest1$p.value
 ## the simulation result of chi-square distribution for one alpha
 p2[i]<-ttest2$p.value
 ## the simulation result of uniform distributon for one alpha
 p3[i]<-ttest3$p.value
 # the simulation result of exponetion distribution for one alpha
 
}
for (j in 1:4)
{
p.hat1[j]<-round(mean(p1<=alpha[j]),3)
## the simulation result of chi-square distribution for four different alphas
p.hat2[j]<-round(mean(p2<=alpha[j]),3)
## the simulation result of uniform distributon for four different alphas
p.hat3[j]<-round(mean(p3<=alpha[j]),3)
# the simulation result of exponetion distribution for four different alphas
}
data.frame(alpha,chi.phat=p.hat1,unif.phat=p.hat2,exp.phat=p.hat3)

## -----------------------------------------------------------------------------
set.seed(1025)
alpha<-c(0.01,0.025,0.05,0.1)
n<-80                         ## the number of sample for one replicate
m<-1.35e3;                       ## the number of repilcates  
p1<-p2<-p3<-numeric(m)
p.hat1<-p.hat2<-p.hat3<-0
for (i in 1:m)
{
 x1<-rchisq(n,1) 
 x2<-runif(n,0,2)
 x3<-rexp(n,1)
 ## the means of above three distributions all are one.
 ttest1<-t.test(x1,alternative="two.side",mu=1)          
 ## the t-test of chi-square distribution
 ttest2<-t.test(x2,alternative="two.side",mu=1)  
 ## the t-test of uniform distribution
 ttest3<-t.test(x3,alternative="two.side",mu=1)
 ## the t-test of exponetion distribution
 p1[i]<-ttest1$p.value
 ## the simulation result of chi-square distribution for one alpha
 p2[i]<-ttest2$p.value
 ## the simulation result of uniform distributon for one alpha
 p3[i]<-ttest3$p.value
 # the simulation result of exponetion distribution for one alpha
 
}
for (j in 1:4)
{
p.hat1[j]<-round(mean(p1<=alpha[j]),3)
## the simulation result of chi-square distribution for four different alphas
p.hat2[j]<-round(mean(p2<=alpha[j]),3)
## the simulation result of uniform distributon for four different alphas
p.hat3[j]<-round(mean(p3<=alpha[j]),3)
# the simulation result of exponetion distribution for four different alphas
}
data.frame(alpha,chi.phat=p.hat1,unif.phat=p.hat2,exp.phat=p.hat3)

## -----------------------------------------------------------------------------
library(bootstrap);
library(boot);
library(MASS);
set.seed(1025)
cor<-matrix(0,5,5)
opar<-par(no.readonly=TRUE)
par(mar=c(1,1,1,1))
par(mfrow=c(3,4))
for (i in 1:4)
{
  for (j in (i+1):5)
  {
    plot(t(scor[i]),t(scor[j]),pch=8,col="red",xlab=colnames(scor)[i],ylab=colnames(scor)[j])
  }
}
par(opar)

cor<-cor(scor)     ## coefficient mareix 
cor

set.seed(1025)
k<-45                 ## the size of one sample
a<-c(1,3,3,4)
b<-c(2,4,5,5)
original<-bias<-se<-numeric(4)
b.cor <-function(x1,i)
  { cor(x1[i,1],x1[i,2])}
for (m in 1:4)
{
  {
  x1<-mvrnorm(k,rep(0,2),matrix(c(cor[a[m],a[m]],cor[a[m],b[m]],cor[a[m],b[m]],cor[b[m],b[m]]),2))
  ## Generate random numbers from a binary normal distribution with a mean vector of (0,0)'and       a covariance matrix required.
  randnums <- boot(data=x1,statistic=b.cor,R=2500)
  original[m]<-round(randnums$t0,4)
  bias[m]<-round(mean(randnums$t)-randnums$t0,3)
  se[m]<-round(sd(randnums$t),3)
  }
}
estimated.parameters<-c("rho12","rho34","rho35","rho45")
data.frame(estimated.parameters,original,bias,se)

## -----------------------------------------------------------------------------
set.seed(1025)
library(boot)
mu<-0
n<-20
interval.norm<-interval.basic<-interval.perc<-matrix(0,1000,2)
skew.f<- function(x,i) {
  #computes the sample skewness coeff.
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3)
  m2 <- mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}

for (i in 1:1000) {
  x<-rnorm(n,0,2)
  obj<-boot(x,statistic=skew.f,R=1000)
  inter<- boot.ci(obj,type=c("norm","basic","perc"))
  interval.norm[i,]<-inter$norm[2:3];
  interval.basic[i,]<-inter$basic[4:5];
  interval.perc[i,]<-inter$percent[4:5];
}
# the coverage probability of three types of interval
norm<-mean(interval.norm[,1]<=mu & interval.norm[,2]>=mu)
basic<-mean(interval.basic[,1]<=mu & interval.basic[,2]>=mu)
perc<-mean(interval.perc[,1]<=mu & interval.perc[,2]>=mu)
# the left side probability of three types of interval
norm.left<-mean(interval.norm[,1]>=mu)
basic.left<-mean(interval.basic[,1]>=mu)
perc.left<-mean(interval.perc[,1]>=mu)
# the right side probability of three types of interval
norm.right<-mean(interval.norm[,2]<=mu)
basic.right<-mean(interval.basic[,2]<=mu)
perc.right<-mean(interval.perc[,2]<=mu)
coverage.p<-c(norm,basic,perc)
left.p<-c(norm.left,basic.left,perc.right)
right.p<-c(norm.right,basic.right,perc.right)
distribution<-c("N(0,4)")
data.frame(distribution,coverage.p,left.p,right.p)

## -----------------------------------------------------------------------------
set.seed(1025)
library(boot)
mu<-sqrt(8/5)
n<-20
interval.norm<-interval.basic<-interval.perc<-matrix(0,1000,2)
skew.f<- function(x,i) {
  #computes the sample skewness coeff.
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3)
  m2 <- mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}

for (i in 1:1000) {
  x<-rchisq(n,5)
  obj<-boot(x,statistic=skew.f,R=1200)
  inter<- boot.ci(obj,type=c("norm","basic","perc"))
  interval.norm[i,]<-inter$norm[2:3];
  interval.basic[i,]<-inter$basic[4:5];
  interval.perc[i,]<-inter$percent[4:5];
}
# the coverage probability of three types of interval
norm<-mean(interval.norm[,1]<=mu & interval.norm[,2]>=mu)
basic<-mean(interval.basic[,1]<=mu & interval.basic[,2]>=mu)
perc<-mean(interval.perc[,1]<=mu & interval.perc[,2]>=mu)
# the left side probability of three types of interval
norm.left<-mean(interval.norm[,1]>=mu)
basic.left<-mean(interval.basic[,1]>=mu)
perc.left<-mean(interval.perc[,1]>=mu)
# the right side probability of three types of interval
norm.right<-mean(interval.norm[,2]<=mu)
basic.right<-mean(interval.basic[,2]<=mu)
perc.right<-mean(interval.perc[,2]<=mu)
coverage.p<-c(norm,basic,perc)
left.p<-c(norm.left,basic.left,perc.right)
right.p<-c(norm.right,basic.right,perc.right)
distribution<-c("chi-square(5)")
data.frame(distribution,coverage.p,left.p,right.p)

## -----------------------------------------------------------------------------
library(bootstrap)
n <- 88
theta.hat <- numeric(n)
cov <- cov(scor)
 ## covariance matrix
eig<- (eigen(cov))$values
 ## eigenvalues computed from the covariance matrix
eigen<- eig[order(-(eig))]
 ## order the eigenvalues from the largest to the lowest
theta <- eigen[1]/sum(eigen)
for (i in 1:n) 
{
  scor.jack <- scor[-i, ]
  cov.hat <- cov(scor.jack)
  ## the jackknife estimates of covariance matrix 
  eig.jack<- (eigen(cov.hat))$values
  ## the jackknife estimates of eigenvalues computed from the covariance matrix
  eigen.jack<- eig.jack[order(-(eig.jack))]
   ## order the eigenvalues from the largest to the lowest
  theta.hat[i] <- eigen.jack[1]/sum(eigen.jack)
}
bias <- (n-1)*(mean(theta.hat)-theta)
se <- sqrt((n-1)*mean((theta.hat-theta)^2))
round(c(original = theta, bias = bias, se = se), 5)


## -----------------------------------------------------------------------------
library(DAAG,quietly=TRUE);
attach(ironslag)
n <- length(magnetic)   #in DAAG ironslag
error1 <- error2 <- error3 <- error4 <- numeric(n)
y.hat1 <- y.hat2 <- y.hat3 <- y.hat4 <- numeric(n)
R1 <- R2 <- R3 <- R4 <- numeric(n)
T1 <- T2 <- T3 <- T4 <- numeric(n)
ymean <-  mean(magnetic)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  ## linear model
  model1 <- lm(y ~ x)
  y.hat1[k] <- model1$coef[1] + model1$coef[2] * chemical[k]
  error1[k] <- magnetic[k] - y.hat1[k]
  R1[k] <- (y.hat1[k]-ymean)^2
  T1[k] <- (magnetic[k]-ymean)^2
  
  ## quadratic model
  model2 <- lm(y ~ x + I(x^2))
  y.hat2[k] <- model2$coef[1] + model2$coef[2] * chemical[k] +model2$coef[3] *    chemical[k]^2
  error2[k] <- magnetic[k] - y.hat2[k]
  R2[k] <- (y.hat2[k]-ymean)^2
  T2[k] <- (magnetic[k]-ymean)^2
  
  ## log-log model
  model3 <- lm(log(y) ~ x)
  y.hat3[k] <- exp(model3$coef[1] + model3$coef[2] * chemical[k])
  error3[k] <- magnetic[k] - y.hat3[k]
  R3[k] <- (y.hat3[k]-ymean)^2
  T3[k] <- (magnetic[k]-ymean)^2
  
  ## cubic polynomial model
  model4 <- lm(y ~ x + I(x^2) + I(x^3))
  y.hat4[k] <- model4$coef[1] + model4$coef[2] * chemical[k] +model4$coef[3] * chemical[k]^2 + model4$coef[4] * chemical[k]^3
  error4[k] <- magnetic[k] - y.hat4[k]
  R4[k] <- (y.hat4[k]-ymean)^2
  T4[k] <- (magnetic[k]-ymean)^2
}
error.square<-c(mean(error1^2), mean(error2^2), mean(error3^2), mean(error4^2))
adjusted.Rsquare<-c(sum(R1)/sum(T1)*51/50, sum(R2)/sum(T2)*51/49, sum(R3)/sum(T3)*51/50, sum(R4)/sum(T4)*51/48)
## adjusted R.square=(SSE/n-p-1)/(SST/n-1) ,n=52
models<-c("linear","quadratic","log-log","cubic")
data.frame(models,error.square,adjusted.Rsquare)


## -----------------------------------------------------------------------------
set.seed(1205)
n<-85
m<-90
N<-1000
extnum<-seq(5,20,1)
p0<-numeric(length(extnum))
p<-numeric(length(extnum))
extp<-numeric(N)
x<-rnorm(n,0,3)
y<-rnorm(m,5,3)
z <- c(x, y);
K <- 1:(n+m);

count5test <- function(x, y,extnum)
  {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) >extnum ))
}


 for (i in 1:length(extnum)) 
 {
 p0[i]<-count5test(x, y,extnum[i])
 for (j in 1:N)
 {
 k <- sample(K, size = n, replace = FALSE)
 x1 <- z[k]
 y1 <- z[-k] #complement of x1
 x1<-x1-mean(x1)
 y1<-y1-mean(y1)
 extp[j]<-count5test(x1, y1,extnum[i])
 }
 p[i]<-mean(abs(c(p0[i],extp)))
 }
data.frame(extnum,p)

plot(extnum,p,type="b",pch=8,col="red",xlab = "the number of the extreme points",ylab="p-value")


## -----------------------------------------------------------------------------
library(boot)
library(MASS)
library(Ball)
alpha<-0.05
dCov <- function(x, y) { 
  x <- as.matrix(x);  y <- as.matrix(y)
  n <- nrow(x); m <- nrow(y)
  if (n != m || n < 2) stop("Sample sizes must agree") 
  if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
  Akl <- function(x) {
    d <- as.matrix(dist(x))
    m <- rowMeans(d); M <- mean(d)
    a <- sweep(d, 1, m); b <- sweep(a, 2, m) 
    b + M
  }
  A <- Akl(x);  B <- Akl(y)
  sqrt(mean(A * B)) 
}

## -----------------------------------------------------------------------------
ndCov2 <- function(z, ix, dims) {
  #dims contains dimensions of x and y 
  p <- dims[1]
  q <- dims[2]
  d <- p + q
  x <- z[ , 1:p] #leave x as is
  y <- z[ix, -(1:p)] #permute rows of y 
  return(nrow(z) * dCov(x, y)^2)
}

## -----------------------------------------------------------------------------
num<-seq(20,200,20)
p.cor1<-p.cor2<-p.ball1<-p.ball2<-numeric(100)
pc.pw1<-pc.pw2<-pb.pw1<-pb.pw2<-numeric(length(num))
for (i in 1:length(num))
{
 for (j in 1:50)
 {
  x0<-mvrnorm(num[i],c(0,0),matrix(c(1,0,0,1),ncol=2))
  e<-mvrnorm(num[i],c(0,0),matrix(c(1,0,0,1),ncol=2))
  y1<-x0/4+e
  y2<-x0/4*e
  z1<- cbind(x0,y1)
  z2<- cbind(x0,y2)
  boot.obj1 <- boot(data = z1, statistic = ndCov2, R = 99,
                 sim = "permutation", dims = c(2, 2))
  boot.obj2 <- boot(data = z2, statistic = ndCov2, R = 99,
                 sim = "permutation", dims = c(2, 2))
  tb1 <- c(boot.obj1$t0, boot.obj1$t)
  tb2 <- c(boot.obj2$t0, boot.obj2$t)
  p.cor1[j] <- mean(tb1>=tb1[1])
  p.cor2[j] <- mean(tb2>=tb2[1])
  p.ball1[j] <- bcov.test(z1[,1:2],z1[,3:4],R=99,seed=i*j*2)$p.value
  p.ball2[j] <- bcov.test(z2[,1:2],z2[,3:4],R=99,seed=i*j*4)$p.value
 }
  pc.pw1[i]<-mean(p.cor1<=alpha)
  ## model1:the power of distance correlation test
  pb.pw1[i]<-mean(p.ball1<=alpha)
  ## model2:the power of ball covariance test
  pc.pw2[i]<-mean(p.cor2<=alpha)
   ## model2:the power of distance correlation test
  pb.pw2[i]<-mean(p.ball2<=alpha)
  ## model2:the power of ball covariance test
}
plot(num,pc.pw1,type="b",pch=8,col="red",xlab="sample size",ylab="power",main="model1:distance versus ball")
lines(num,pb.pw1,type="b",pch=6,col="blue")
legend("bottomright",c("distance","ball"),lty=c(1,2),pch=c(8,6),col=c("red","blue"))
plot(num,pc.pw2,type="b",pch=5,col="brown",xlab="sample size",ylab="power",main="model2:distance versus ball")
lines(num,pb.pw2,type="b",pch=4,col="green")
legend("bottomright",c("distance","ball"),lty=c(1,2),pch=c(5,4),col=c("brown","green"))

## -----------------------------------------------------------------------------
set.seed(1025)
    rw.Metropolis <- function(sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (
                    u[i] <= ((1/2*exp(-abs(y)))/(1/2*exp(-abs(x[i-1])))))
        # the density of Laplace function f(x)=1/2*exp(-abs(x)
                    x[i] <- y  
                else {
                    x[i] <- x[i-1]
                    k <- k + 1
                }
            }
        return(list(x=x, k=k))
    }

    n <- 4  #degrees of freedom for target Student t dist.
    N <- 4000
    sigma <- c(.05, .5, 6,  20)

    x0 <- 10
    rw1 <- rw.Metropolis(sigma[1], x0, N)
    rw2 <- rw.Metropolis(sigma[2], x0, N)
    rw3 <- rw.Metropolis(sigma[3], x0, N)
    rw4 <- rw.Metropolis(sigma[4], x0, N)

    #number of candidate points rejected
    reject=c(rw1$k, rw2$k, rw3$k, rw4$k)
    reject.rate<-reject/N
    accept<-N-reject
    accept.rate<-accept/N
    ABC <- data.frame(sigma=sigma,accept,accept.rate,reject,reject.rate)
    knitr::kable(ABC)

## -----------------------------------------------------------------------------

    par(mar=c(1,1,1,1))
    par(mfrow=c(2,2))  #display 4 graphs together
    refline <- c(log(0.05),-log(0.05))
    #  Use inverse function to find quantiles
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
    par(mar=c(1,1,1,1))
    par(mfrow=c(1,1)) #reset to default
    

## -----------------------------------------------------------------------------
a <- c(.05, seq(.1, .9, .1), .95)
    Q1<-log(2*a[1:6])      
    #  Use inverse function to find quantiles
    Q2<--log(2-2*a[7:11])
    Q <- c(Q1,Q2)
    rw <- cbind(rw1$x, rw2$x, rw3$x, rw4$x)
    mc <- rw[501:N, ]
    Qrw <- apply(mc, 2, function(x) quantile(x, a))
    ABC <- data.frame(round(cbind(Q, Qrw), 3))
    names(ABC) <- c('True','sigma=0.05','sigma=0.5','sigma=2','sigma=16')
    knitr::kable(ABC) 

## -----------------------------------------------------------------------------
x<-seq(2,16,2)
a<-log(exp(x))
b<-exp(log(x))
a==b
a==x
b==x
isTRUE(all.equal(a,b))

## -----------------------------------------------------------------------------
k<-c(5:25,100) ##The root of equation hasn't been solved when k equals 4,so the loop starts when k equals 5
root<-numeric(length(k))
for (i in 1:length(k))
{
  coef<-function(k)
  {
   return(gamma((k+1)/2)/(gamma(k/2)*sqrt(k*pi)))
    ## the coefficient of integral
  }
 f<-function(a)
{
(coef(k[i])*(integrate(function(x)(1+(x^2)/k[i])^(-(1+k[i])/2),lower=sqrt((a^2)*k[i]/(k[i]+1-a^2)),upper=Inf)$value)- coef(k[i]-1)*(integrate(function(x) (1+(x^2)/(k[i]-1))^(-k[i]/2),lower=sqrt((a^2)*(k[i]-1)/(k[i]-a^2)),upper=Inf)$value))
}
root[i]<-uniroot(f,lower =0.01,upper = 1+sqrt(k[i])/2)$root
}
data.frame(k,root)

## -----------------------------------------------------------------------------
k<-c(4:25,100)
## ##The root of equation hasn't been solved when k equals 500 and 1000,so the loop starts from k equals 4 to k equals 100.
root<-numeric(length(k))
for (i in 1:length(k))
{
  coef<-function(k)
  {
   return(gamma((k+1)/2)/(gamma(k/2)*sqrt(k*pi)))
    ## the coefficient of integral
  }
 f<-function(a)
{
(coef(k[i])*(integrate(function(x)(1+(x^2)/k[i])^(-(1+k[i])/2),lower=0,upper=sqrt((a^2)*k[i]/(k[i]+1-a^2)))$value)- coef(k[i]-1)*(integrate(function(x) (1+(x^2)/(k[i]-1))^(-k[i]/2),lower=0,upper=sqrt((a^2)*(k[i]-1)/(k[i]-a^2)))$value))
}
root[i]<-uniroot(f,lower =0.01,upper = 1+sqrt(k[i])/2)$root
}
data.frame(k,root)

## -----------------------------------------------------------------------------
p <- q <- r <- numeric(1000)
nA. <- 28; nB. <- 24; nOO <- 41; nAB <- 70
p[1] <- 0.09; q[1] <- 0.23;r[1] <- (1- p[1]- q[1])
threshold <- 2e-6


f <- function(a,b) {
  return((nB.*b/(2-b-2*a)+nB.+nAB)/(nA.*a/(2-a-2*b)+nA.+nAB))
}
g <- function(a,b) {
 return(((1-a/(2-a-2*b))*nA.+(1-b/(2-b-2*a))*nB.+2*nOO)/((nB.*b/(2-b-2*a)+
                          nB.+nAB)))
}

for (i in 2:1000)
   {
   p[i] <- 1/(1+f(p[i-1],q[i-1])*(1+g(p[i-1],q[i-1])))
   ## the maximum likelihood estimation of p
   q[i] <- f(p[i-1],q[i-1])/(1+f(p[i-1],q[i-1])*(1+g(p[i-1],q[i-1])))
   ## the maximum likelihood estimation of q
   r[i] <- 1- p[i] - q[i]
  
   if((p[i]-p[i-1] <= threshold) & (q[i]-q[i-1] <= threshold) &
      (r[i]-r[i-1] <= threshold))
   ##  the loop will stop if the p's and q's differece of values between two   iterations are both samller the given threshold.
       {print(c(i, p[i], q[i],r[i]))
       break
    }
}
x <-c(1:i)
plot(x, p[1:i], type="b",pch=8, col = "purple",ylim=c(0,0.7), main = "The log-maximum likelihood values in M-steps" ,lty = 4, xlab = "The number of iterations", ylab = "The parameters' value of iterations")
lines(x, q[1:i], type= "b", pch=11,col = "red",lty = 5)
lines(x, r[1:i], type= "b", pch=6, col = "blue",lty = 6)
legend("bottomright", legend = c("p", "q", "r"),pch=c(8,11,6),lty =c(4,5,6), col = c("purple", "red", "blue"))

## -----------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
a<-numeric(4)
for (i in seq_along(formulas))
{
a[i]<-lapply(i,function(i) {lm(formulas[[i]],data=mtcars)})
}
a


## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
lapply(bootstraps, lm, formula = mpg ~ disp)


## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
b<-numeric(10)
for (i in 1:10)
{
  b[i]<-lapply(bootstraps[i], lm, formula = mpg ~ disp)

}
b

## -----------------------------------------------------------------------------
formulas <- list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
a<-numeric(4)
rsq <- function(mod) summary(mod)$r.squared
for (i in 1:4){
  a[i]<-lapply(i,function(i){rsq(lm(formulas[[i]],data=mtcars))})
}
a

## -----------------------------------------------------------------------------

bootstraps <- lapply(1:10, function(i) {
rows <- sample(1:nrow(mtcars), rep = TRUE)
mtcars[rows, ]
})
rsq <- function(mod) summary(mod)$r.squared

coe<-numeric(10)
for (i in 1:10)
{
 coe[i]<-lapply(i,function(i){rsq(lm(mpg~disp,data=bootstraps[[i]]))})
}
coe


## -----------------------------------------------------------------------------
trials <- replicate(100,
           t.test(rpois(10, 10), rpois(7, 10)),
           simplify = FALSE
           )
sapply(1:100, function(i) {trials[[i]]$p.value})   

## -----------------------------------------------------------------------------
sapply(trials, "[[",3)

## -----------------------------------------------------------------------------
library(parallel)
# mcsapply()
mcsapply<-function(k,f){
cl <- makeCluster(getOption("cl.cores", 4))
result<-parLapply(cl,k,f) 
stopCluster(cl) 
return(unlist(result))
} 
trials <- replicate(
         2000,
         t.test(rnorm(11,1,3), rnorm(12,2,4)),
         simplify = FALSE
       )
system.time(mcsapply(trials,function(x) unlist(x)[3]))

## -----------------------------------------------------------------------------
system.time(sapply(trials,function(x) unlist(x)[3]))

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

set.seed(1025)
rw_MetropolisR <- function(sigma, x0, N) 
{
  #Metropolis Randomwalk using R
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-(abs(y) - abs(x[i-1]))))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
library(Rcpp)
#// This is the rw_MetropolisC.cpp
#include <Rcpp.h>
#using namespace Rcpp;
#// [[Rcpp::export]]
cppFunction('NumericVector rw_MetropolisC(double sigma, double x0, int N) 
{
  //Metropolis Randomwalk using C
  NumericVector x(N);
  x[0] = x0;
  double u, y;
  int k = 0;
  for (int i = 1; i < N; i++) 
  {
    y = rnorm(1, x[i-1], sigma)[0];
    u = runif(1)[0];
    if (u <= exp(-(abs(y) - abs(x[i-1])))) 
    {
      x[i] = y; 
    }
    else 
    {
      x[i] = x[i-1];
      k++;
    }
  }
  return x;
}')

## -----------------------------------------------------------------------------
N = 4000
sigma <- c(0.5,12,60)
x0 = 12
for (i in 1:length(sigma)){
ts = microbenchmark(rwR = rw_MetropolisR(sigma[i], x0, N)$x, 
                    rwC = rw_MetropolisC(sigma[i], x0, N))
print(summary(ts)[, c(1,3,5,6)])

rwR = rw_MetropolisR(sigma[i], x0, N)$x
rwC = rw_MetropolisC(sigma[i], x0, N)
par(mar=c(1,1,1,1))
par(mfrow = c(1, 3))
b <- 1001 #discard the burnin sample
y <- (rwR)[b:N]
a <- ppoints(100)
QR <- ifelse(a <= 1/2, log(2*a), -log(2-2*a)) #quantiles of Laplace
QQ1 <- quantile(rwR, a)
qqplot(QR, QQ1, main=paste("R sigma=",sigma[i],seq=NULL), xlab="Laplace Quantiles", ylab="Sample Quantiles",col="blue")
abline(a=0, b=1)

y <- (rwC)[b:N]
a <- ppoints(100)
QR <- ifelse(a <= 1/2, log(2*a), -log(2-2*a)) #quantiles of Laplace
QQ2 <- quantile(rwC, a)
qqplot(QR, QQ2, main=paste("C sigma=",sigma[i],seq=NULL), xlab="Laplace Quantiles", ylab="Sample Quantiles",col="purple")
abline(a=0, b=1)

qqplot(QQ1, QQ2, main=paste("C-R sigma=",sigma[i],seq=NULL), xlab="R Samplee Quantiles", ylab="C Sample Quantiles",col="red")
abline(a=0, b=1)
}

