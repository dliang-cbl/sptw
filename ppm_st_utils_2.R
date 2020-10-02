## loglik function
## x: the quasi-likelihood fit
## data: the observed and quadrature points (default the input data)
## wcoord: name for weight variable (default wt)
## value: log-likelihood of the fitted poisson process
## adjust the total area so it match
ppm.st.lik <- function(x,data=NULL,wcoord=NULL){
  
  if(is.null(data)){
    data.nona <- subset(x$model,!is.na(z))
  }else{
    data.nona <- subset(data,!is.na(z))
  }
  
  if(!is.null(wcoord)){
    data.nona$wt <- data.nona[,wcoord]
  }
  
  if(is.null(data)){
    lambda <- fitted(x)
  }else{
    lambda <- predict(x,newdata=data.nona,type="response")
  }
  
  l.lik <- with(data.nona,{
    wt * (z*log(lambda)-lambda)
  })
  sum(l.lik) - sum(log(1:sum(data.nona$z>0)))
}

## value: residual measure of each quadrature point
ppm.st.residuals <- function(x,data){
  l.data <- subset(data,!is.na(z))
  lambda <- predict(x,newdata=l.data,type="response")
  l.res <- with(l.data,{
    wt*(z-lambda)*1/sqrt(lambda)
  })
  l.res
}
test <- function()
{
  ppm.st.lik(gam0,glm0)
  ## -200.7945
}

qaic <- function(gamobj,wcoord="(weights)"){
  D <- -2*ppm.st.lik(gamobj,wcoord = wcoord)
  llik <- logLik(gamobj)
  df <- attr(llik,"df")
  data.frame(D=D,df=df,AIC=D+2*df)
}

## x: predicted data with residuals
## formula: define residual on left and area on right and time period on right
## value:
## a residual plot
## plus the aggregated time
plot.ppm.st.residuals <- function(formula,x,...)
{
  tmp1 <- aggregate(formula,data=x,FUN=sum)
  vars <- all.vars(formula)
  wt <- tmp1[,vars[1]]
  upr <- 2*sqrt(wt)
  lwr <- -2*sqrt(wt)
  e <- tmp1[,vars[2]]
  t_ <- tmp1[,vars[3]]
  plot(t_,e,ylim=range(e,upr,lwr),
       xlab="Time period",ylab="Pearson Residual",...)
  polygon(x=c(t_,rev(t_)),y=c(lwr,rev(upr)),col=grey(0.9))
  points(t_,e,...)
  lines(lowess(t_,e),col="red")
  tmp1
}

test <- function()
{
  
}