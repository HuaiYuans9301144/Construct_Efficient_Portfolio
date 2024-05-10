
if (!require(Rsolnp)) { #載入 Rsolnp 函數庫
  install.packages("Rsolnp") #安裝 Rsolnp (會自動下載，故需連接網路)
  require(Rsolnp)# 載入 Rsolnp 函數庫
}

#拉氏乘數法
# W: 投資組合權重
# MU.i: 各資產期望報酬率
# SIGMA.i: 資產報酬率共變數矩陣
#The return of portfolio
mu.p <- function(W, MU.i, SIGMA.i){
  return (sum(W*MU.i)) #use the sum function to return real value
}

#The variance of portfolio
sigma.2.p <- function(W, MU.i, SIGMA.i){
  return (sum((W%*%SIGMA.i)%*%W)) #use the sum function to return real value
}

#The sum of weights
eqf1 <- function(W, MU.i, SIGMA.i){
  return (sum(W))
}
#The sum of weights and the expected return
eqf2 <- function(W, MU.i, SIGMA.i){
  SumWeight=sum(W)
  ExpectedReturn=mu.p(W, MU.i, SIGMA.i)
  return (c(SumWeight, ExpectedReturn))
}
#######################################################
#在要求的期望報酬之下，極小化投資組合變異數
#######################################################

#找三檔個股 分別是群創 長榮 欣興
#rf使用央行定存利率
setwd("C:\\Users\\user\\Desktop")
data<-read.csv("stockdata1.csv")

RF<-data[,4]
STOCK<-data[,1:3]
MKT<-data[,5]

W0 <-c(1/3, 1/3, 1/3)
MU.i0<-colMeans((STOCK-RF/12))*12
SIGMA.i0 <- cov((STOCK-RF/12))*12

mu.p(W0, MU.i0, SIGMA.i0)
sigma.2.p(W0, MU.i0, SIGMA.i0)

lb<-rep(0,length(MU.i0))
ub<-NULL

if(1){
  #MKT<- rowMeans(STOCK)
  PR<- cbind(STOCK,MKT)-(RF/12)
  Ys <- colnames(STOCK)
  Rs<-list()
  for(i in (1:length(Ys))){
    Rs[[i]] <- lm(paste(Ys[i],"~MKT"),PR)
    if(i>1){
      Yhat<- cbind(Yhat,Rs[[i]]$fitted.values)
    }else{
      Yhat<- Rs[[i]]$fitted.values
    }
    
  }
  MU.i0 <- colMeans(Yhat)*12
  SIGMA.i0 <- cov(Yhat)*12
  for(i in(1:length(Ys))){
    SIGMA.i0[i,i]<-SIGMA.i0[i,i]+summary(Rs[[i]])$sigma^2*12
  }
}

rf<-mean(RF)
#The required expected return
TargetReturn = 0.1
minimum.variance.portfolio<-solnp(pars=W0,
                                  fun=sigma.2.p,
                                  MU.i = MU.i0, SIGMA.i = SIGMA.i0, 
                                  eqfun=eqf2, eqB=c(1, TargetReturn),
                                  control=list(trace=FALSE)
)


minimum.variance.portfolio$pars
mu.p(minimum.variance.portfolio$pars, MU.i0, SIGMA.i0)
sigma.2.p(minimum.variance.portfolio$pars, MU.i0, SIGMA.i0)
sum(minimum.variance.portfolio$pars)
#######################################################
#The efficient frontier.
#######################################################
TargetReturns <- seq(0, 1, 0.01)
Min.STD <- rep(0, length(TargetReturns))
i <- 0 ;
for (TargetReturn in TargetReturns){
  i <- i+1;
  W <- solnp(pars=W0,
             fun=sigma.2.p,
             MU.i = MU.i0, SIGMA.i = SIGMA.i0,
             eqfun=eqf2, eqB=c(1, TargetReturn),
             control=list(trace=FALSE)
  )
  Min.STD[i] <- sigma.2.p(W$pars, MU.i0, SIGMA.i0)^0.5
}
plot( Min.STD, TargetReturns, type="l", lwd=2) #Plot the efficient frontier.
#######################################################
#Maximize the sharpe ratio.
#The Sharpe Ratio for maximization
sharpe.ratio <- function(W, MU.i, SIGMA.i, RF){
  sp <- (mu.p(W, MU.i, SIGMA.i)-RF)/sigma.2.p(W, MU.i, SIGMA.i)^.5
  return(sp)
}
#The negative Sharpe Ratio for minimization
neg.sharpe.ratio <- function(W, MU.i, SIGMA.i, RF){
  nsp <- -sharpe.ratio(W, MU.i, SIGMA.i, RF)
  return(nsp)
}

#The sum of weights for the maximixatio of sharpe ratio
sharpe.ratio.eqf1 <- function(W, MU.i, SIGMA.i, RF){
  return (eqf1(W, MU.i, SIGMA.i))
}

neg.sharpe.ratio(W0, MU.i = MU.i0, SIGMA.i = SIGMA.i0, RF=rf)
sharpe.ratio(W0, MU.i = MU.i0, SIGMA.i = SIGMA.i0, RF=rf)
sp.W <- solnp(pars=W0,
              fun=neg.sharpe.ratio,
              MU.i = MU.i0, SIGMA.i = SIGMA.i0, RF=rf,
              eqfun=sharpe.ratio.eqf1, eqB=1,
              control=list(trace=FALSE)
)
sp.W
sp<-sharpe.ratio(sp.W$pars, MU.i0, SIGMA.i0, rf) #計算最大sharpe ratio
sp.W.M <- mu.p(sp.W$pars, MU.i0, SIGMA.i0) #最大sharpe ratio的expected return 
sp.W.D <- sigma.2.p(sp.W$pars, MU.i0, SIGMA.i0)^.5 #最大sharpe ratio的 STD
sp.W.M #最大sharpe ratio的expected return 
sp.W.D #最大sharpe ratio的 STD

ER = rf + sp*Min.STD; #計算CAL線函數值
CAL <- lm(ER~Min.STD) #建構CAL線函數模型
points(sp.W.D,sp.W.M,pch=16,col="red") #標出最大sharpe ratio的位置(紅)
abline(CAL,lwd=1, col="blue") #劃出CAL線(藍)

