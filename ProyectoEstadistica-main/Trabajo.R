library(mixtools)
library(mclust)
library(ks)
library(plot3D)
library(tidyverse)
library(MASS)
library(plotly)
library(mvnTest)
library(nortest)
library(readxl)
library(data.table)


#######################
testAD<-function(data,mus,sigmapob,lambdapob){
  if (!is.data.frame(data) && !is.matrix(data)) 
    stop('data supplied must be either of class \"data frame\" or \"matrix\"')
  if (dim(data)[2] < 2 || is.null(dim(data))) 
  {stop("data dimesion has to be more than 1")}
  if (dim(data)[1] < 3) {stop("not enough data for assessing mvn")}
  data.name <- deparse(substitute(data))
  xp <- as.matrix(data)
  p <- dim(xp)[2]
  n <- dim(xp)[1]
  ## getting MLEs...
  s.mean <- colMeans(xp)
  s.cov <- (n-1)/n*cov(xp)
  s.cov.inv <- solve(s.cov) # inverse matrix of S (matrix of sample covariances)
  D <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
  for (j in 1:n)
    D[j] <- t(xp[j,]-s.mean)%*%(s.cov.inv%*%(xp[j,]-s.mean))
  D.or <- sort(D) ## get ordered statistics
  Gp <- pchisq(D.or,df=p)
  ## getting the value of A-D test...
  ind <- c(1:n)
  an <- (2*ind-1)*(log(Gp[ind])+log(1 - Gp[n+1-ind]))
  AD <- -n - sum(an) / n
  ## getting the p-value...
  N <- 1e4
  U <- rep(0,N) ## initializing values of the AD test
  for (i in 1:N) { ## loop through N reps
    dat<-rmvnorm.mixt(1000, mus=mus, Sigmas=sigmapob, props=lambdapob)
    mean1 <- colMeans(dat)
    cov1 <- (n-1)/n*cov(dat)
    cov.inv <- solve(cov1) # inverse matrix of S (matrix of sample covariances)
    D <- rep(NA,n) # vector of (Xi-mu)'S^-1(Xi-mu)...
    for (j in 1:n)
      D[j] <- t(data[j,]-mean1)%*%(cov.inv%*%(data[j,]-mean1))
    Gp <- pchisq(sort(D),df=p)
    ## getting the value of A-D test...
    an <- (2*ind-1)*(log(Gp[ind])+log(1 - Gp[n+1-ind]))
    U[i] <- -n - sum(an) / n
  }
  p.value <- (sum(U >= AD)+1)/(N+1)
  result<-new('ad',AD=AD,p.value=p.value,data.name=data.name)
  result
}


#############################

###DATA BMI
data<-read.csv("bmi.csv")
altura<-data$Height
peso<-data$Weight
ind<-peso/altura^2
plot(density(ind))
data1<-data %>% filter(BmiClass=='Normal Weight')
altura<-data1$Height
peso<-data1$Weight
ind1<-peso/altura^2

plot(density(peso))
plot(density(altura))
plot(density(peso/altura^2))


###MULTIVARIANTE
den3d <- kde2d(peso, altura)
persp(den3d, box=FALSE)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()

matriz<-matrix(c(peso,altura),nrow=342,ncol=2)
colnames(matriz)<-c('Peso','Altura')
em<-mvnormalmixEM(matriz,k=3)

mus<-rbind(em$mu[[1]],em$mu[[2]],em$mu[[3]])
sigmas<-rbind(em$sigma[[1]],em$sigma[[2]],em$sigma[[3]])
lambdas<-as.vector(em$lambda)
dat<-rmvnorm.mixt(5000, mus=mus, Sigmas=sigmas, props=as.vector(lambdas))
indicador_sim<-(dat[,1]/dat[,2]^2)
plot(density(indicador_sim))
den3d <- kde2d(dat[,1],dat[,2])
persp(den3d, box=FALSE)
plot_ly(x=den3d$x, y=den3d$y, z=den3d$z) %>% add_surface()


####TEST
test<-testAD(matriz,mus,sigmas,lambdas)
test
#Si sigue una distribución de mixturas multinormales

#MONTECARLO
nsim<-5000
nx<-100
media<-numeric(nsim)
varianza<-numeric(nsim)
desv<-numeric(nsim)
for(i in 1:5000){
  datos_sim<-rmvnorm.mixt(nx, mus=mus, Sigmas=sigmas, props=as.vector(lambdas))
  indicador_sim<-datos_sim[,1]/datos_sim[,2]^2
  media[i]<-mean(indicador_sim)
  varianza[i]<-var(indicador_sim)
  desv[i]<-sd(indicador_sim)
}
#Estimacion media, varianza, desvest(precision)
estimediaMont<-mean(media)
estivarMont<-mean(varianza)
estidesMont<-mean(desv)

#Sesgo
sesgoMontmedia<-estimediaMont-mean(ind1)
sesgoMont
sesgoMontvar<-estivarMont-var(ind1)
sesgoMontvar
sesgoMontdes<-estidesMont-sd(ind1)
sesgoMontdes
#Intervalos confianza 
#media
alfa <- 0.05
ic<-quantile(media,c(alfa/2, 1 - alfa/2))
ic

#CONSTRASTE DE HIPOTESIS
#Ho: u=22.64
#Ha: u!=22.64
#Constatar normalidad de la media del indicador
nsim<-5000
nx<-1000
pvalor<-numeric(nsim)
#Constrastes
for (i in 1:nsim) {
  datos_sim<-rmvnorm.mixt(nx, mus=mus, Sigmas=sigmas, props=as.vector(lambdas))
  indicador_sim<-datos_sim[,1]/datos_sim[,2]^2
  t<-(mean(indicador_sim)-22.64)/(sd(indicador_sim)/sqrt(nx))
  p.value<-1-pt(abs(t),nx-1)+pt(-abs(t),nx-1)
  pvalor[i]<-p.value
}
{
  cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
  cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
  cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
}



###BOOTSTRAP
x <- ind1
n<-length(x)
h <- bw.SJ(x)
npden <- density(x, bw = h)
plot(npden)
range_x<-range(npden$x)

#MEDIA, VARIANZA, PRECISION

B <- 1000
n<-length(ind1)
estadistico_boot <- numeric(B)
var_boot <- numeric(B)
desv_boot <- numeric(B)
for (k in 1:B) {
  remuestra <- sample(ind1, n, replace = TRUE)
  estadistico_boot[k] <- mean(remuestra)
  var_boot[k]<-var(remuestra)
  desv_boot[k]<-sd(remuestra)
  
}
estimedia_boot <- mean(estadistico_boot)  
estimedia_boot
estivar_boot <- mean(var_boot)  
estivar_boot
estisd_boot <- mean(desv_boot)  
estisd_boot 
#SESGO
sesgobootmedia<-estimedia_boot-mean(ind1)
sesgobootmedia
sesgobootvar<-estivar_boot-var(ind1)
sesgobootvar
sesgobootdes<-estisd_boot-sd(ind1)
sesgobootdes




#INTERVALOS - METODO PERCENTIL T
alfa<-0.05
B <- 1000
n<-length(ind1)
remuestra<-numeric(n)
estadistico_boot <- numeric(B)
for (k in 1:B) {
  remuestra <- sample(ind1, n, replace = TRUE)
  x_barra_boot<-mean(remuestra)
  cuasi_dt_boot <- sd(remuestra)
  estadistico_boot[k] <- sqrt(n) * (x_barra_boot - mean(ind1))/cuasi_dt_boot  
}
pto_crit <- quantile(estadistico_boot, c(alfa/2, 1 - alfa/2))
# Construcción del IC
ic_inf_boot <- mean(ind1) - pto_crit[2] * sd(ind1)/sqrt(n)
ic_sup_boot <- mean(ind1) - pto_crit[1] * sd(ind1)/sqrt(n)
IC_boot <- c(ic_inf_boot, ic_sup_boot) 
names(IC_boot) <- paste0(100*c(alfa/2, 1-alfa/2), "%")
IC_boot


#CONTRASTE DE HIPOTESIS 

#CONSTRASTE DE HIPOTESIS
#Ho: u=22.64
#Ha: u!=22.64
n<-length(ind1)
nsim<-5000
pvalor<-numeric(nsim)
#Constrastes
for (i in 1:nsim) {
  remuestra <- sample(ind1, n, replace = TRUE)
  t<-(mean(remuestra)-22.64)/(sd(remuestra)/sqrt(n))
  p.value<-1-pt(abs(t),n-1)+pt(-abs(t),n-1)
  pvalor[i]<-p.value
}
{
  cat("\nProporción de rechazos al 1% =", mean(pvalor < 0.01), "\n")
  cat("Proporción de rechazos al 5% =", mean(pvalor < 0.05), "\n")
  cat("Proporción de rechazos al 10% =", mean(pvalor < 0.1), "\n")
}
