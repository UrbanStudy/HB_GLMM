# HB_GLMM
Hierarchical Byes for Generalized Linear Mixed Effect Model

## The First Example

### NHTS data

[The data source](https://nhts.ornl.gov/)


```{r, eval=F}
NHTS2017 <- (read.csv("~/trippub.csv"))[,c(1,30,62,64,69,72,85)]
# NHTS2017 <- NHTS2017[complete.cases(NHTS2017),]
NHTS2017 <- NHTS2017[NHTS2017$VMT_MILE!=-1&NHTS2017$HHFAMINC>=0&NHTS2017$HH_CBSA!="XXXXX", ]
nhts2017 <- NHTS2017[sample(nrow(NHTS2017), 10000,replace =F), ]
save(nhts2017, file="nhts2017.RData")
```

Sampling 10000 observations from the original data excluded the zero-miles VMT, negative household income, and unknown CBSA id (XXXXX)

```{r, eval=T}
load("nhts2017.RData")
str(nhts2017)
summary(nhts2017)
table(nhts2017$HH_CBSA)
```

```{r}
ids<-sort(unique(nhts2017$HH_CBSA)) 
m<-length(ids)
Y<-list() ; X<-list() ; N<-NULL
for(j in 1:m) 
{
  Y[[j]]<-nhts2017[nhts2017$HH_CBSA==ids[j],2] 
  N[j]<- sum(nhts2017$HH_CBSA==ids[j])
  xj<-nhts2017[nhts2017$HH_CBSA==ids[j], 4] 
  xj<-(xj-mean(xj))
  X[[j]]<-cbind( rep(1,N[j]), xj  )
}
```

```{r}
#### OLS fits
S2.LS<-BETA.LS<-NULL
for(j in 1:m) {
  fit<-lm(Y[[j]]~-1+X[[j]] )
  BETA.LS<-rbind(BETA.LS,c(fit$coef)) 
  S2.LS<-c(S2.LS, summary(fit)$sigma^2) 
}
```

```{r}
plot( range(nhts2017[,4]),c(0,40),type="n",xlab="HHFAMINC", ylab="VMT")  # range(NHTS2017[,2])
for(j in 1:m) {    abline(BETA.LS[j,1],BETA.LS[j,2],col="gray")  }

BETA.MLS<-apply(BETA.LS,2,mean)
abline(BETA.MLS[1],BETA.MLS[2],lwd=2)

plot(N,BETA.LS[,1],xlab="sample size",ylab="intercept")
abline(h= BETA.MLS[1],col="black",lwd=2)
plot(N,BETA.LS[,2],xlab="sample size",ylab="slope")
abline(h= BETA.MLS[2],col="black",lwd=2)
```
