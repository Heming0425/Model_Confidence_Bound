warnings('off')

#Transfrom a list of variable name denoted selection results to 0-1 matrix result
f01<-function(object){

  matrix.01<-data.frame(matrix(0,length(object),p));colnames(matrix.01)<-full.var
  for(i in 1:length(object)){matrix.01[i,(full.var %in% object[[i]])]<-1}
  return(matrix.01)
}

#Transform 0-1 matrix variable selection result ot a list of variable name denoted selection resuls
f02<-function(object)
{

  var.matrix<-vector(mode="list",length=dim(object)[1])
  for(i in 1:dim(object)[1])
  {
    var.matrix[[i]]<-full.var[object[i,]>0]
  }
  return(var.matrix)
}

#big_matrix_illustration function:numerical vector to name list transformation
f03 <- function(object){

  name_vector <- full.var[object > 0]
}

#function3
inter.union<-function(object){
  j=3;
  inter.sect<-intersect(object[[1]],object[[2]])
  union.sect<-union(object[[1]],object[[2]])
  while(j<=k){
    inter.sect<-intersect(inter.sect,object[[j]])
    union.sect<-union(union.sect,object[[j]])
    j<-j+1
  }
  return(list(inter<-inter.sect,union<-union.sect))
}

#online adaptive lasso code:http://www4.stat.ncsu.edu/~boos/var.select/lasso.adaptive.html
lasso.adapt.bic2<-function(x,y){

  # adaptive lasso from lars with BIC stopping rule
  # this one uses the "known variance" version of BIC with RSS/(full model mse)
  # must use a recent version of R so that normalize=FALSE can be used in lars

  require(lars)
  ok<-complete.cases(x,y)
  x<-x[ok,]                            # get rid of na's
  y<-y[ok]                             # since regsubsets can't handle na's
  m<-ncol(x)
  n<-nrow(x)
  x<-as.matrix(x)                      # in case x is not a matrix

  #  standardize variables like lars does
  one <- rep(1, n)
  meanx <- drop(one %*% x)/n
  xc <- scale(x, meanx, FALSE)         # first subtracts mean
  normx <- sqrt(drop(one %*% (xc^2)))
  names(normx) <- NULL
  xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)

  out.ls=lm(y~xs)                      # ols fit on standardized
  beta.ols=out.ls$coeff[2:(m+1)]       # ols except for intercept
  w=abs(beta.ols)                      # weights for adaptive lasso
  xs=scale(xs,center=FALSE,scale=1/w)  # xs times the weights
  object=lars(xs,y,type="lasso",normalize=FALSE)

  # get min BIC
  # bic=log(n)*object$df+n*log(as.vector(object$RSS)/n)   # rss/n version
  sig2f=summary(out.ls)$sigma^2        # full model mse
  bic2=log(n)*object$df+as.vector(object$RSS)/sig2f       # Cp version
  step.bic2=which.min(bic2)            # step with min BIC

  fit=predict.lars(object,xs,s=step.bic2,type="fit",mode="step")$fit
  coeff=predict.lars(object,xs,s=step.bic2,type="coef",mode="step")$coefficients
  coeff=coeff*w/normx                  # get back in right scale
  st=sum(coeff !=0)                    # number nonzero
  mse=sum((y-fit)^2)/(n-st-1)          # 1 for the intercept

  # this next line just finds the variable id of coeff. not equal 0
  if(st>0) x.ind<-as.vector(which(coeff !=0)) else x.ind<-0
  intercept=as.numeric(mean(y)-meanx%*%coeff)
  return(list(fit=fit,st=st,mse=mse,x.ind=x.ind,coeff=coeff,intercept=intercept,object=object,
              bic2=bic2,step.bic2=step.bic2))
}

#Get the bootsrap variable selection result
# Note you can use this function for different variable selection method, but the default one is adaptive lasso
BOOT.CI<-function(x, dep.var.index, r){
  var.instances<-vector(mode="list",length=r)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    boot.data<-x[ind,]
    # fit<-regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
    # var.instances[[j]]<-names(coef(fit,which.min(summary(fit)$bic)))[-1]
    fit<-lasso.adapt.bic2(x=boot.data[,1:p],y=boot.data$y)
    var.instances[[j]]<-names(fit$coeff)[(fit$coeff)!=0]
    # adalasso<-adalasso(X=as.matrix(boot.data[,1:p]),y=boot.data$y,k=10)
    # var.instances[[j]]<-full.var[adalasso$coefficients.adalasso!=0]
    # opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    # lasso.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    # beta4<-coef(lasso.fit,s=(opt.lambda+0.2))[,1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta4)[beta4!=0]]
    # fit <- cv.lqa(y.train = boot.data[,dep.var.index],x.train = as.matrix(boot.data[,-dep.var.index]),intercept = TRUE,lambda.candidates = list(c(1)),family = gaussian(),penalty.family = lasso,n.fold = 10,loss.func = "aic.loss")
    # beta<- coef(fit$best.obj)[-1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var.instances)
}

# GET the residual bootsrap variable selection results by using adaptive lasso and
# you can give a specified lambda
RES.BOOT.CI2<-function(x, dep.var.index, r,lmbd,seed){
  var.instances<-vector(mode="list",length=r)
  ## = adaLASSO = ##
  tau=1
  lasso_init=glmnet(as.matrix(x[,1:p]),x$y) #nlambda = 100
  first.step.coef=lasso_init$beta[,which.min(abs(lasso_init$lambda-lmbd))]
  penalty.factor=abs(first.step.coef+1/sqrt(nrow(x)))^(-tau)
  adalasso=glmnet(as.matrix(x[,1:p]),x$y,penalty.factor=penalty.factor)
  beta_est=adalasso$beta[,which.min(abs(adalasso$lambda-lmbd))] #lmbd找最近的
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta_est
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta_est
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    lasso_init_boot=glmnet(as.matrix(boot.data[,1:p]),boot.data$y)
    first.step.coef_boot=lasso_init_boot$beta[,which.min(abs(lasso_init_boot$lambda-lmbd))]
    penalty.factor_boot=abs(first.step.coef_boot+1/sqrt(nrow(x)))^(-tau)
    adalasso_boot=glmnet(as.matrix(boot.data[,1:p]),boot.data$y,penalty.factor=penalty.factor_boot)
    beta_est_boot=adalasso_boot$beta[,which.min(abs(adalasso_boot$lambda-lmbd))]
    var.instances[[j]]<-full.var[beta_est_boot!=0]
  }
  return(var.instances)
}

#Get the residual bootsrap variable selection result by using adaptive lasso
# Note you can use this function for different variable selection method, but the default one is adaptive lasso
RES.BOOT.CI<-function(x, dep.var.index, r,seed=seed){
  var.instances<-vector(mode="list",length=r)
  adalasso<-adalasso(X=as.matrix(x[,1:p]),y=x$y,k=10) #k=10
  res_original <- x$y - as.matrix(x[,1:p]) %*% adalasso$coefficients.adalasso
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% adalasso$coefficients.adalasso
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    # fit<-regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
    # var.instances[[j]]<-names(coef(fit,which.min(summary(fit)$bic)))[-1]
    # fit<-lasso.adapt.bic2(x=boot.data[,1:p],y=boot.data$y)
    # var.instances[[j]]<-names(fit$coeff)[(fit$coeff)!=0]
    adalasso_tem<-adalasso(X=as.matrix(boot.data[,1:p]),y=boot.data$y,k=10)
    var.instances[[j]]<-full.var[adalasso_tem$coefficients.adalasso!=0]
    # opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    # lasso.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    # beta4<-coef(lasso.fit,s=(opt.lambda+0.2))[,1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta4)[beta4!=0]]
    # fit <- cv.lqa(y.train = boot.data[,dep.var.index],x.train = as.matrix(boot.data[,-dep.var.index]),intercept = TRUE,lambda.candidates = list(c(1)),family = gaussian(),penalty.family = lasso,n.fold = 10,loss.func = "aic.loss")
    # beta<- coef(fit$best.obj)[-1]
    # var.instances[[j]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var.instances)
}


# get residual bootstrap variable selection models from SCAD and MCP
RES.BOOT.CI3 <- function(x, dep.var.index, r, pnlt, lmbd,seed){
  var<-vector(mode="list",length=r)
  fit<-cv.ncvreg(X=x[,1:p],y=x$y,penalty=pnlt)
  fit2<-fit$fit
  #若无lmbd,选最小,若有,选择最接近的
  if(lmbd == ''){
    beta<-fit2$beta[,fit$min]
  }else{
    beta<-fit2$beta[,which.min(abs(fit$lambda-lmbd))]
  }
  beta<-beta[2:(p+1)]
  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- cv.ncvreg(X=as.matrix(boot.data[,1:p]),y=boot.data$y,penalty=pnlt)
    fit2_tem<-tem$fit
    if(lmbd==''){
      beta_tem<-fit2_tem$beta[,tem$min]
    }else{
      beta_tem<-fit2_tem$beta[,which.min(abs(tem$lambda-lmbd))]
    }
    beta_tem <- beta_tem[2:(p+1)]
    var[[j]]<-names(beta_tem)[which(beta_tem!=0)]
  }
  return(var)
}

# Use LAD and Sqrt method to get residual bootstrap variable models
RES.BOOT.CI4 <- function(x, dep.var.index, r, q, lmbd,seed){
  var<-vector(mode="list",length=r)
  fit<-slim(X=as.matrix(x[,1:p]),Y=x$y,method="lq",q=q)
  if(lmbd==''){
    beta <- fit$beta[,which.min(abs(fit$lambda-1))]
  }else{
    beta <- fit$beta[,which.min(abs(fit$lambda-lmbd))]
  }

  res_original <- x$y - as.matrix(x[,1:p]) %*% beta
  res_after_center <- res_original - mean(res_original)
  constant <- as.matrix(x[,1:p]) %*% beta
  set.seed(seed=seed)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x[,1:p], new_response)
    colnames(boot.data)[p+1] <- "y"
    tem <- slim(X=as.matrix(boot.data[,1:p]),Y=boot.data$y,method = 'lq', q = q,verbose = FALSE)
    if(lmbd==''){
      var[[j]]<-full.var[tem$beta[,which.min(abs(fit$lambda-1))]!=0]
    }else{
      var[[j]]<-full.var[tem$beta[,which.min(abs(fit$lambda-lmbd))]!=0]
    }

  }
  return(var)
}

# GET residual boostrap variable selection models with stepwise BIC
RES.BOOT.CI5 <- function(x, p, r,lmbd,seed){
  if(lmbd==''){
    var<-vector(mode="list",length=r)
    fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
    if (names(coef(fit,which.min(summary(fit)$bic)))[1] == '(Intercept)'){
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))[-1]
    } else {
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))
    }
    res_original <- x$y - as.matrix(x[,1:p]) %*% beta
    res_after_center <- res_original - mean(res_original)
    constant <- as.matrix(x[,1:p]) %*% beta
    set.seed(seed=seed)
    for(j in 1:r) {
      ind=sample(1:nrow(x),nrow(x),replace=T)
      new_response <- constant + res_after_center[ind]
      boot.data <- cbind(x[,1:p], new_response)
      colnames(boot.data)[p+1] <- "y"
      tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
      if (names(coef(tem,which.min(summary(tem)$bic)))[1] == '(Intercept)'){
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))[-1]
      } else {
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))
      }
    }
    return(var)
  }else{
    var<-vector(mode="list",length=r)
    fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
    if (names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))[1] == '(Intercept)'){
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))] <- coef(fit,which.min(summary(fit)$bic))[-1]
    } else {
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))] <- coef(fit,which.min(summary(fit)$bic))
    }
    res_original <- x$y - as.matrix(x[,1:p]) %*% beta
    res_after_center <- res_original - mean(res_original)
    constant <- as.matrix(x[,1:p]) %*% beta
    for(j in 1:r) {
      ind=sample(1:nrow(x),nrow(x),replace=T)
      new_response <- constant + res_after_center[ind]
      boot.data <- cbind(x[,1:p], new_response)
      colnames(boot.data)[p+1] <- "y"
      tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
      if (names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))[1] == '(Intercept)'){
        var[[j]]<-names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))[-1]
      } else {
        var[[j]]<-names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))
      }
    }
    return(var)
  }
}

# Get boostrap models from the modified lasso
BOOT.MODI.LASSO <- function(x, dep.var.index, r, lmbd,seed){
  var_lasso<-vector(mode="list",length=r)
  set.seed(seed=seed)
  for (i in 1:r){
    ind=sample(1:nrow(x),nrow(x),replace=T)
    boot.data<-x[ind,]
    opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    lasso.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=1)
    beta<-coef(lasso.fit,s=opt.lambda)[,1][which(abs(coef(lasso.fit,s=opt.lambda)[,1]) > lmbd)]
    var_lasso[[i]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var_lasso)
}

BOOT.MODI.ELASTIC <- function(x, dep.var.index, r, lmbd,seed){
  var_lasso<-vector(mode="list",length=r)
  set.seed(seed=seed)
  for (i in 1:r){
    ind=sample(1:nrow(x),nrow(x),replace=T)
    boot.data<-x[ind,]
    opt.lambda<-cv.glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],alpha=1)$lambda.min
    elastic.fit<-glmnet(x=as.matrix(boot.data[,-dep.var.index]),y=boot.data[,dep.var.index],family='gaussian',alpha=0.5)
    beta<-coef(elastic.fit,s=opt.lambda)[,1][which(abs(coef(elastic.fit,s=opt.lambda)[,1]) > lmbd)]
    var_lasso[[i]]<-full.var[full.var%in%names(beta)[beta!=0]]
  }
  return(var_lasso)
}

#Calculate the Model confidences bounds and their corresponding freq rate
CI<-function(var.list,var.matrix,p)
{
  full.var<-colnames(var.matrix)
  colsum<-apply(var.matrix,2,sum)
  order<-order(colsum,decreasing = T)
  freq<-vector(length=p+1);freq[1]<-0
  lower<-vector(mode="list",length=p+1)
  upper<-vector(mode="list",length=p+1)
  for(i in 0:p)
  {
    cap<-vector(length=p+1)
    cap[1]<-0
    for(j in 0:(p-i))
    {
      if (j==0 & i!=0)
      {
        uppertest<-full.var[order[1:i]]
        for(r in 1:length(var.list))
        {
          if(all(var.list[[r]]%in%uppertest)) cap[j+1]<-cap[j+1]+1
        }
      }else{
        if(j!=0){
          lowtest<-full.var[order[1:j]]
          uppertest<-full.var[order[1:(j+i)]]
          for(r in 1:length(var.list))
          {
            if(all(all(lowtest%in%var.list[[r]]),all(var.list[[r]]%in%uppertest))) cap[j+1]<-cap[j+1]+1
          }
        }
        if (j == 0){
          for (r in 1:length(var.list)){
            if (identical(var.list[[r]],character(0))) cap[j+1] <- cap[j+1] + 1
          }
        }
      }
    }
    freq[i+1]<-max(cap)/r
    maxlocation<-which.max(cap)
    if(maxlocation==1)
    {
      if (i != 0){
        lower[[i+1]]<-''
        upper[[i+1]]<-full.var[order[1:i]]
      } else if (i == 0){
        lower[[1]]<-''
        upper[[1]]<-''
      }
    }else{
      lower[[i+1]]<-full.var[order[1:(maxlocation-1)]]
      upper[[i+1]]<-full.var[order[1:(maxlocation-1+i)]]
    }
  }
  result<-list(freq=freq,lower=lower,upper=upper)
  return(result)
}

# FULL MODEL - CI
CI.Ufull <- function(var.list,var.matrix,p)
{
  full.var<-colnames(var.matrix)
  colsum<-apply(var.matrix,2,sum)
  order<-order(colsum,decreasing = T)
  freq<-vector(length=p+1);freq[1]<-0
  lower<-vector(mode="list",length=p+1)
  upper<-vector(mode="list",length=p+1)
  for(i in 0:p)
  {
    cap<-vector(length=p+1)
    cap[1]<-0
    j <- p - i
    if (j==0 & i!=0){
      uppertest<-full.var[order[1:p]]
      for(r in 1:length(var.list)){
        if(all(var.list[[r]]%in%uppertest)) cap[j+1]<-cap[j+1]+1
      }
    }else{
      if(j!=0){
        lowtest<-full.var[order[1:j]]
        uppertest<-full.var[order[1:p]]
        for(r in 1:length(var.list))
        {
          if(all(all(lowtest%in%var.list[[r]]),all(var.list[[r]]%in%uppertest))) cap[j+1]<-cap[j+1]+1
        }
      }
      if (j == 0){
        for (r in 1:length(var.list)){
          if (identical(var.list[[r]],character(0))) cap[j+1] <- cap[j+1] + 1
        }
      }
    }
    freq[i+1]<-max(cap)/r

    if (i == 0){
      lower[[i+1]]<-full.var[order[1:p]]
      upper[[i+1]]<-full.var[order[1:p]]
    }else if(i == p){
      lower[[i+1]]<-""
      upper[[i+1]]<-full.var[order[1:p]]
    }else{
      lower[[i+1]]<-full.var[order[1:(p-i)]]
      upper[[i+1]]<-full.var[order[1:p]]
    }
  }
  result<-list(freq=freq,lower=lower,upper=upper)
  return(result)
}

# A function used in calculating optimal MCB
var.f<-function(k){
  full.combine<-full.var
  var.include<-vector("list",length=choose(p,k))
  for(i in 1:length(var.include)){
    var.include[[i]]<-full.combine[t(combn(p,k))[i,]]
  }
  return(var.include)
}

# A function used in calculating optimal MCB
f1<-function(varmatrix,varinstance)
{
  freq<-rep(0,length(varmatrix))
  for(i in 1:length(varmatrix))
  {
    for(j in 1:length(varinstance))
    {
      freq[i]<-freq[i]+all(varmatrix[[i]] %in% varinstance[[j]])
    }
    freq[i]<-freq[i]/length(varinstance)
  }
  return(max(freq))
}

# A function used in calculating optimal MCB
f2<-function(varmatrix,varinstance)
{
  freq<-rep(0,length(varmatrix))
  for(i in 1:length(varmatrix))
  {
    for(j in 1:length(varinstance))
    {
      freq[i]<-freq[i]+all(varinstance[[j]] %in% varmatrix[[i]])
    }
    freq[i]<-freq[i]/length(varinstance)
  }
  return(max(freq))
}

# This function is to all possiable lower bound with length x, and there is no repeat with choose((p-low.len),width)
low.var.norepeat<-function(p,low.len){
  var.low<-vector("list")
  var<-as.matrix(t(combn(full.var,low.len)))
  var<-apply(var,1,paste,collapse="#")
  return(sapply(var,strsplit,split="#",USE.NAMES = F))
}

# A function used in calculating optimal MCB
low_up.var<-function(p,low.len,width){
  # This function is to all possiable lower bound with length x, and there is repeat with choose((p-low.len),width)
  low.var.repeat<-function(p,low.len,width){
    var.low<-vector("list")
    var<-as.matrix(t(combn(full.var,low.len)))
    var<-apply(var,1,paste,collapse="#")
    object<-as.matrix(rep(var,each=choose((p-low.len),width)))
    return(sapply(object,strsplit,split="#",USE.NAMES = F))
  }
  low.var.object<-low.var.repeat(p,low.len,width)

  add.var<-function(p,low.len,width){
    var<-as.matrix(t(combn(full.var,low.len)))
    ind<-t(apply(var,1,function(var){!full.var%in%var}))
    object<-list()
    for(i in 1:nrow(ind)){
      object[[i]]<-full.var[ind[i,]]
    }
    add<-lapply(object,combn,m=width)
    add1<-lapply(add,t);len<-nrow(add1[[1]])
    temp<-matrix(,nr=length(add1)*len,nc=ncol(add1[[1]]))
    for(i in 1:length(add1)){
      temp[((i-1)*len+1):(i*len),]<-add1[[i]]
    }
    return(temp)
  }
  add.var.object<-add.var(p,low.len,width)

  up.var<-list()
  for(i in 1:length(low.var.object)){
    up.var[[i]]<-c(low.var.object[[i]],add.var.object[i,])
  }
  return(list(low=low.var.object,up=up.var))
}

# The function to calculate optimal MCB
##i is the  length of lower bound
##w is the width of the confidence bound
cioptimal<-function(var.instances, p,r){
  res<-data.frame(matrix(,nr=p,nc=p))
  colnames(res)<-paste("width",0:(p-1),sep="")
  for(w in 0:(p-1)){
    for(i in 1:(p-w)){
      var.set<-low_up.var(p,i,w)
      count<-rep(0,length(var.set$low))
      cap<-vector(length=length(count))
      for(j in 1:length(var.set$low)){

        for(k in 1:r){
          count[j]<-count[j]+all(all(var.set$low[[j]] %in% var.instances[[k]]),all(var.instances[[k]] %in% var.set$up[[j]]))
        }
        cap[j]<-count[j]/r
      }
      res[i,(w+1)]<-max(cap)
    }
  }
  cap_res<-apply(res,2,max,na.rm=T)
  return(cap_res)
}

# Given the data set, to see whether MCBs could cover the true model
captrue_model<-function(x, dep.var.index, r, lambda_list)
{
  freq_matrix <- matrix(0, nrow = length(lambda_list), ncol = dep.var.index)
  truecapture_matrix <- matrix(0, nrow = length(lambda_list), ncol = dep.var.index)
  for (i in 1:length(lambda_list)){
    var.instances <- RES.BOOT.CI2(x, dep.var.index, r, lambda_list[i])
    var.01 <- f01(var.instances)
    result<-CI(var.instances,var.01,p)
    lower<-result$lower
    upper<-result$upper
    freq_matrix[i,] <- result$freq
    truecapture_matrix[i,] <- cap_true(lower,upper,p)
  }
  results<-list(freq=freq_matrix,truecapture=truecapture_matrix)
  return(results)
}

# A functional function for the function "captrue_model"
cap_true<-function(low_list,upper_list,p){
  truecapture<-vector()
  for(i in 1:(p+1))
  {
    if (low_list[[i]][1] == ""){
      if (all(all(true.model %in% upper_list[[i]]))) truecapture[i]<-1 else truecapture[i]<-0
    } else{
      if(all(all(low_list[[i]]%in%true.model),all(true.model %in% upper_list[[i]]))) truecapture[i]<-1 else truecapture[i]<-0
    }
  }
  return(truecapture)
}

# A function a analysis_of_bootstrap_validity: used to quantify the selection method consistency rate
judge_equal <- function(data, compared_model){
  adalasso_tem <- adalasso(X=as.matrix(data[,1:p]),y=data$y,k=10)
  var.instances <- full.var[adalasso_tem$coefficients.adalasso!=0]
  if (length(var.instances) == length(compared_model) && sum(var.instances%in%compared_model) == length(compared_model)) equal <-1 else equal <- 0
  return(equal)
}



#MODEL RANKING
model_ranking <- function(input_matrix, input_vector){
  row_dim <- dim(input_matrix)[1]
  if (row_dim == 0){
    input_matrix = rbind(input_matrix, input_vector)
  } else{
    insert_index = 1
    for (i in 1:row_dim){
      k = 0
      if (sum(input_vector) > sum(input_matrix[i,])){
        insert_index <- insert_index + 1
        k = 1
      } else if (sum(input_vector) == sum(input_matrix[i,])){
        model_dif <- input_vector - input_matrix[i,]
        for (j in 1:p){
          if (model_dif[j] == -1){
            insert_index <- insert_index + 1
            k = 1
            break
          } else if (model_dif[j] == 1) {
            k = 0
            break
          }
        }
      }
      if (k == 0){
        break
      }
    }
    if (insert_index == (row_dim + 1)){
      input_matrix <- rbind(input_matrix, input_vector)
    } else{
      input_matrix <- insertRow(input_matrix, input_vector, insert_index)
    }
  }
  return(input_matrix)
}

# A function in big_matrix_illustration.R
insertRow <- function(existingDF, newrow, r){
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}



# A revised version of MCB_compare_other, so you can set the tuning parameter of adaptive lasso by yourself
MCB_compare_other2 <- function(x, dep.var.index, r, candidate_lambda){
  if (calculate_MCB == 1){
    MCB_result <- MCB_acc_cardi_calculate(x, r, candidate_lambda)
    cardinality_MCB <- tryCatch(MCB_result$cardinality_MCB_sub, error = function(e){print('error3'); problem = MCB_result})
    capture_MCB <- tryCatch(MCB_result$capture_MCB_sub, error = function(e){print('error3'); problem = MCB_result})
    score_MCB <- tryCatch(MCB_result$score_MCB_sub,error = function(e){print('error3'); problem = MCB_result})
  }
  if (calculate_VSCS == 1){
    VSCS_result <- VSCS_acc_cardi_calculate(x)
    cardinality_VSCS <- tryCatch(VSCS_result$cardinality_VSCS_sub,error = function(e){print('error4'); problem = VSCS_result})
    capture_VSCS <- tryCatch(VSCS_result$capture_VSCS_sub,error = function(e){print('error4'); problem = VSCS_result})
  }
  if (calculate_MCB == 1 && calculate_VSCS == 1) return_result<-list(cardinality_MCB=cardinality_MCB,cardinality_VSCS=cardinality_VSCS,capture_MCB=capture_MCB,capture_VSCS=capture_VSCS,score_MCB= score_MCB)
  if (calculate_MCB == 1 && calculate_VSCS == 0) return_result<-list(cardinality_MCB=cardinality_MCB,capture_MCB=capture_MCB,score_MCB= score_MCB)
  if (calculate_MCB == 0 && calculate_VSCS == 1) return_result<-list(cardinality_VSCS=cardinality_VSCS,capture_VSCS=capture_VSCS)
  return(return_result)
}

# Calcuate the true model capture and cardinality of VSCS
VSCS_acc_cardi_calculate <- function(dataset){
  cardinality_VSCS_sub <- vector('numeric', length = length(alpha))
  capture_VSCS_sub <- vector('numeric', length = length(alpha))
  x_transform <- cbind(dataset[,p+1], dataset[,1:p])
  colnames(x_transform)[1] <- 'y'
  result_VSCS <- ecs(full.data=x_transform, full.family='gaussian')
  for (zz in 1:length(alpha)){
    capture_VSCS_sub[zz] <- VSCS_capture_true_model(result_VSCS[[zz+1]], true.model_vector)
    if (is.null(result_VSCS[[zz+1]])) cardinality_VSCS_sub[zz] <- 0 else if (is.vector(result_VSCS[[(zz+1)]])) {
      cardinality_VSCS_sub[zz] <- 1
    } else {cardinality_VSCS_sub[zz] <- dim(result_VSCS[[zz+1]])[1]}
  }
  return_value <- list(cardinality_VSCS_sub=cardinality_VSCS_sub, capture_VSCS_sub= capture_VSCS_sub)
  return(return_value)
}

# Calculate the true model capture and cardinality and score of MCB
MCB_acc_cardi_calculate <- function(dataset, r, lambda_list){
  cardinality_MCB_sub <- matrix(0, nrow = length(lambda_list), ncol = length(alpha))
  capture_MCB_sub <- matrix(0, nrow = length(lambda_list), ncol = length(alpha))
  score_MCB_sub <- matrix(0, nrow = length(lambda_list), ncol = length(alpha))
  for (i in 1:length(lambda_list)){
    lambda <- lambda_list[i]
    var.instances<-RES.BOOT.CI2(dataset,p + 1,r,lambda)
    var.01<-f01(var.instances)
    result_MCB<-CI(var.instances,var.01,p)
    ci_lower<-result_MCB$lower
    ci_upper<-result_MCB$upper
    # print("print out MCB")
    # print(lambda_list[i])
    # print("print lower bound")
    # print(ci_lower)
    # print("upper bound")
    # print(ci_upper)
    # print(result_MCB$freq)
    for (j in 1:length(alpha)){
      talpha <- alpha[j]
      lower_alpha<-ci_lower[which(result_MCB$freq > (1 - talpha))[1]]
      upper_alpha<-ci_upper[which(result_MCB$freq > (1 - talpha))[1]]
      if (lower_alpha[[1]][1] == ""){
        if (all(all(true.model %in% upper_alpha[[1]]))) {
          capture_MCB_sub[i, j] <- 1
          score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha)
        }else {
          capture_MCB_sub[i, j] <- 0
          score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha) + 2/talpha * (length(true.model) - length(upper_alpha))
        }
      } else{
        if(all(all(lower_alpha[[1]]%in%true.model),all(true.model %in% upper_alpha[[1]]))) {
          capture_MCB_sub[i, j] <- 1
          score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha)
        }else {
          capture_MCB_sub[i, j] <- 0
          score_MCB_sub[i,j] <- length(upper_alpha) - length(lower_alpha) + 2/talpha * max((length(true.model) - length(upper_alpha)),0) + 2/talpha * max((length(lower_alpha) - length(true.model)),0)
        }
      }
      lower_alpha<-f01(lower_alpha)
      upper_alpha<-f01(upper_alpha)
      cardinality_MCB_sub[i,j] <- 2^(sum(upper_alpha) - sum(lower_alpha))
    }
  }
  return_value = list(cardinality_MCB_sub=cardinality_MCB_sub, capture_MCB_sub= capture_MCB_sub, score_MCB_sub=score_MCB_sub)
  return(return_value)
}


# A function used in MCB and VSCS comparison
VSCS_capture_true_model <- function(model_matrix, real_model_vector){
  capture <- 0
  if (is.null(dim(model_matrix))){return(capture)} else{
    if (is.vector(model_matrix)) {
      if(sum(abs(model_matrix - real_model_vector)) == 0){
        capture <- 1
        return(capture)
      } else {
        return(capture)
      }
    } else if (dim(model_matrix)[1] == 0) return(capture) else{
      for (i in 1:dim(model_matrix)[1]){
        if (sum(abs(model_matrix[i,] - real_model_vector)) == 0){
          capture <- 1
          break
        }
      }
      return(capture)
    }
  }
}


getmcb <- function(result,c){
  n = length(result$freq)
  fit_freq = c(result$freq)[result$freq - c >= 0]
  best_fit = which.min(fit_freq - c) + n - sum(result$freq - c >= 0)
  mc = list()
  mc$lbm <- result$lower[[best_fit]]
  mc$ubm <- result$upper[[best_fit]]
  mc$bcr <- result$freq[best_fit]
  return(mc)
}

getmcbt <- function(result){
  mcbf <- data.frame(lbm = matrix(result$lower), bcr = result$freq, ubm = matrix(result$upper))
  n = length(mcbf$lbm)
  mcbf$width <- c(0:(n-1))
  mcbf <- mcbf[,c('width','lbm','bcr','ubm')]
  return(mcbf)
}

# 主函数

mcb.compare <- function(x, y, B=200, lambdas=NULL, methods = NULL, level=0.95,seed=122, ufull=FALSE){
  #var.x 传入可供选择的变量,var.x = c("var1","var2".....)
  #pred.y 被预测变量y
  #data 数据集
  #r bootstrap次数
  #methods 变量选择方法,methods = c('adalasso', 'lasso', 'SCAD', 'MCP','stepwise','LAD','SQRT'),默认为所有
  #lmbd 只寻找最近的lmbd
  warnings('off')

  x = as.data.frame(x)
  n = length(x)
  for (i in 1:n) {
    name = paste('x',i,sep='')
    colnames(x)[i] = name
  }
  y = as.data.frame(y)

  #数据变换
  data = cbind(x,y)
  full.var <<- paste(colnames(x),sep="")
  colnames(data)[dim(data)[2]] = "y"
  p <<- dim(data)[2]-1
  r = B

  final_result = ''

  if(is.null(lambdas)){
    lambdas = ''
  }

  if(is.null(methods)){
    methods <- c('aLasso', 'Lasso','Elastic', 'SCAD', 'MCP','stepwise','LAD','SQRT')
  }

  mcbfit = list()
  mcbframe = list()

  #adaptivelasso
  if('aLasso' %in% methods){
    if(lambdas == ''){
      var_adaptivelasso <- RES.BOOT.CI2(data, p + 1, r,lmbd=1,seed=seed)
    }else{
      var_adaptivelasso <- RES.BOOT.CI2(data, p + 1, r,lmbd=lambdas[which(methods=='adaptiveLasso')],seed=seed)
    }
    var_01_ada_lasso<-f01(var_adaptivelasso)
    if(ufull){
      result_aLasso <- CI.Ufull(var_adaptivelasso, var_01_ada_lasso, p)
    }else{
      result_aLasso <- CI(var_adaptivelasso, var_01_ada_lasso, p)
    }
    mcbfit$aLasso <- getmcb(result_aLasso,c=level)
    mcbframe$aLasso <- getmcbt(result_aLasso)
  }

  #lasso
  if('Lasso' %in% methods){
    if(lambdas == ''){
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=0.05,seed=seed)
    }else{
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=lambdas[which(methods=='Lasso')],seed=seed)
    }
    var_01_lasso<-f01(var_lasso)
    if(ufull){
      result_Lasso <- CI.Ufull(var_lasso, var_01_lasso, p)
    }else{
      result_Lasso <- CI(var_lasso, var_01_lasso, p)
    }
    mcbfit$Lasso<-getmcb(result_Lasso,c=level)
    mcbframe$Lasso<-getmcbt(result_Lasso)
  }

  #Elastic
  if('Elastic' %in% methods){
    if(lambdas == ''){
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=0.05,seed=seed)
    }else{
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=lambdas[which(methods=='Lasso')],seed=seed)
    }
    var_01_elastic<-f01(var_elastic)
    if(ufull){
      result_Elastic <- CI.Ufull(var_elastic, var_01_elastic, p)
    }else{
      result_Elastic <- CI(var_elastic, var_01_elastic, p)
    }
    mcbfit$Elastic<-getmcb(result_Elastic,c=level)
    mcbframe$Elastic<-getmcbt(result_Elastic)
  }

  #SCAD与MCP的lmbd默认都是找最小
  #SCAD
  if('SCAD' %in% methods){
    if(lambdas==''){
      var_SCAD <- RES.BOOT.CI3(data, p + 1, r, pnlt='SCAD',lmbd=lambdas,seed=seed)
    }else{
      var_SCAD <- RES.BOOT.CI3(data, p + 1, r, pnlt='SCAD',lmbd=lambdas[which(methods=='SCAD')],seed=seed)
    }
    var_01_SCAD<-f01(var_SCAD)
    if(ufull){
      result_SCAD <- CI.Ufull(var_SCAD, var_01_SCAD, p)
    }else{
      result_SCAD <- CI(var_SCAD, var_01_SCAD, p)
    }
    mcbfit$SCAD <- getmcb(result_SCAD,c=level)
    mcbframe$SCAD <- getmcbt(result_SCAD)
  }

  #MCP
  if('MCP' %in% methods){
    if (lambdas==''){
      var_MCP <- RES.BOOT.CI3(data, p + 1, r, pnlt='MCP',lmbd=lambdas,seed=seed)
    }else{
      var_MCP <- RES.BOOT.CI3(data, p + 1, r, pnlt='MCP',lmbd=lambdas[which(methods=='MCP')],seed=seed)
    }
    var_01_MCP<-f01(var_MCP)
    if(ufull){
      result_MCP <- CI.Ufull(var_MCP, var_01_MCP, p)
    }else{
      result_MCP <- CI(var_MCP, var_01_MCP, p)
    }
    mcbfit$MCP <- getmcb(result_MCP,c=level)
    mcbframe$MCP <- getmcbt(result_MCP)
  }

  #lmbd默认最小
  #stepwise
  if('stepwise' %in% methods){
    if(lambdas==''){
      var_stepwise <- RES.BOOT.CI5(data, p, r,lmbd=lambdas,seed=seed)
    }else{
      var_stepwise <- RES.BOOT.CI5(data, p, r,lmbd=lambdas[which(methods=='stepwise')],seed=seed)
    }
    var_01_stepwise <- f01(var_stepwise)
    if(ufull){
      result_stepwise <- CI.Ufull(var_stepwise, var_01_stepwise, p)
    }else{
      result_stepwise <- CI(var_stepwise, var_01_stepwise, p)
    }
    mcbfit$stepwise <- getmcb(result_stepwise,c=level)
    mcbframe$stepwise <- getmcbt(result_stepwise)
  }

  #LAD与SQRT的lmbd默认为最接近1
  #LAD
  if('LAD' %in% methods){
    if(lambdas==''){
      var_LAD <- RES.BOOT.CI4(data, p+1, r, q = 1, lmbd = lambdas,seed=seed)
    }else{
      var_LAD <- RES.BOOT.CI4(data, p+1, r, q = 1, lmbd = lambdas[which(methods=='LAD')],seed=seed)
    }
    var_01_LAD <- f01(var_LAD)
    if(ufull){
      result_LAD <- CI.Ufull(var_LAD, var_01_LAD, p)
    }else{
      result_LAD <- CI(var_LAD, var_01_LAD, p)
    }
    mcbfit$LAD <- getmcb(result_LAD,c=level)
    mcbframe$LAD <- getmcbt(result_LAD)
  }


  #SQRT
  if('SQRT' %in% methods){
    if(lambdas==''){
      var_SQRT <- RES.BOOT.CI4(data, p+1, r, q = 2, lmbd = lambdas,seed=seed)
    }else{
      var_SQRT <- RES.BOOT.CI4(data, p+1, r, q = 2, lmbd = lambdas[which(methods=='SQRT')],seed=seed)
    }
    var_01_SQRT <- f01(var_SQRT)
    if(ufull){
      result_SQRT <- CI.Ufull(var_SQRT, var_01_SQRT, p)
    }else{
      result_SQRT <- CI(var_SQRT, var_01_SQRT, p)
    }
    mcbfit$SQRT <- getmcb(result_SQRT,c=level)
    mcbframe$SQRT <- getmcbt(result_SQRT)
  }

  final_result = ''
  select_name = paste('result_',methods,sep = '')
  for(i in select_name){
    final_result <- cbind(final_result,get(i)$freq)
  }

  final_result <- final_result[,c(2:(length(methods)+1))]
  colnames(final_result) <- methods

  all_result = list()

  all_result$mcb <- mcbfit
  all_result$mcbframe <- mcbframe

  return(all_result)
}


####修改完成
mcb <- function(x, y, B=200, lambda=NULL, method = 'Lasso', level=0.95, seed = 122,ufull=FALSE){
  #var.x 传入可供选择的变量,var.x = c("var1","var2".....)
  #pred.y 被预测变量y
  #data 数据集
  #r bootstrap次数
  #methods 变量选择方法,methods = c('adalasso', 'lasso', 'SCAD', 'MCP','stepwise','LAD','SQRT'),默认为所有
  #lmbd 只寻找最近的lmbd
  #数据变换
  warnings('off')

  x = as.data.frame(x)
  n = length(x)
  for (i in 1:n) {
    name = paste('x',i,sep='')
    colnames(x)[i] = name
  }
  y = as.data.frame(y)


  data = cbind(x,y)
  full.var <<- paste(colnames(x),sep="")
  colnames(data)[dim(data)[2]] = "y"
  p <<- dim(data)[2]-1
  r = B

  if(is.null(lambda)){
    lambda = ''
  }

  #adaptivelasso
  if(method == 'aLasso'){
    if(lambda == ''){
      var_adaptivelasso <- RES.BOOT.CI2(data, p + 1, r,lmbd=1,seed=seed)
    }else{
      var_adaptivelasso <- RES.BOOT.CI2(data, p + 1, r,lmbd=lambda,seed=seed)
    }
    var_01_ada_lasso<-f01(var_adaptivelasso)
    if(ufull){
      result <- CI.Ufull(var_adaptivelasso, var_01_ada_lasso, p)
    }else{
      result <- CI(var_adaptivelasso, var_01_ada_lasso, p)
    }
  }

  #lasso
  if(method == 'Lasso'){
    if(lambda == ''){
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=0.05,seed=seed)
    }else{
      var_lasso <- BOOT.MODI.LASSO(data, p + 1, r, lmbd=lambda,seed=seed)
    }
    var_01_lasso<-f01(var_lasso)
    if(ufull){
      result <- CI.Ufull(var_lasso, var_01_lasso, p)
    }else{
      result <- CI(var_lasso, var_01_lasso, p)
    }
  }

  #Elastic
  if(method == 'Elastic'){
    if(lambda == ''){
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=0.05,seed=seed)
    }else{
      var_elastic <- BOOT.MODI.ELASTIC(data, p + 1, r, lmbd=lambda,seed=seed)
    }
    var_01_elastic<-f01(var_elastic)
    if(ufull){
      result <- CI.Ufull(var_elastic, var_01_elastic, p)
    }else{
      result <- CI(var_elastic, var_01_elastic, p)
    }
  }

  #SCAD与MCP的lmbd默认都是找最小
  #SCAD
  if(method == 'SCAD'){
    var_SCAD <- RES.BOOT.CI3(data, p + 1, r, pnlt='SCAD',lmbd=lambda,seed=seed)
    var_01_SCAD<-f01(var_SCAD)
    if(ufull){
      result <- CI.Ufull(var_SCAD, var_01_SCAD, p)
    }else{
      result <- CI(var_SCAD, var_01_SCAD, p)
    }
  }

  #MCP
  if(method=='MCP'){
    var_MCP <- RES.BOOT.CI3(data, p + 1, r, pnlt='MCP',lmbd=lambda,seed=seed)
    var_01_MCP<-f01(var_MCP)
    if(ufull){
      result <- CI.Ufull(var_MCP, var_01_MCP, p)
    }
    else{
      result <- CI(var_MCP, var_01_MCP, p)
    }
  }

  #lmbd默认最小
  #stepwise
  if(method=='stepwise'){
    var_stepwise <- RES.BOOT.CI5(data, p, r,lmbd=lambda,seed=seed)
    var_01_stepwise <- f01(var_stepwise)
    if(ufull){
      result <- CI.Ufull(var_stepwise, var_01_stepwise, p)
    }else{
      result <- CI(var_stepwise, var_01_stepwise, p)
    }
  }

  #LAD与SQRT的lmbd默认为最接近1
  #LAD
  if(method=='LAD'){
    var_LAD <- RES.BOOT.CI4(data, p+1, r, q = 1, lmbd = lambda,seed=seed)
    var_01_LAD <- f01(var_LAD)
    if(ufull){
      result <- CI.Ufull(var_LAD, var_01_LAD, p)
    }else{
      result <- CI(var_LAD, var_01_LAD, p)
    }
  }


  #SQRT
  if(method=='SQRT'){
    var_SQRT <- RES.BOOT.CI4(data, p+1, r, q = 2, lmbd = lambda,seed=seed)
    var_01_SQRT <- f01(var_SQRT)
    if(ufull){
      result <- CI.Ufull(var_SQRT, var_01_SQRT, p)
    }else{
      result <- CI(var_SQRT, var_01_SQRT, p)
    }
  }

  all_result <- list()

  fit_freq = c(result$freq)[result$freq - level >= 0]
  best_fit = which.min(fit_freq - level) + n+1 - sum(result$freq - level >= 0)
  mcb = list()
  mcb$lbm <- result$lower[[best_fit]]
  mcb$ubm <- result$upper[[best_fit]]
  mcb$bcr <- result$freq[best_fit]
  all_result$mcb <- mcb

  mcbframe <- data.frame(lbm = matrix(result$lower), bcr = result$freq, ubm = matrix(result$upper))
  mcbframe$width <- c(0:n)
  mcbframe <- mcbframe[,c('width','lbm','bcr','ubm')]
  all_result$mcbframe <- mcbframe

  return(all_result)
}

mucplot <- function(obj){

  obj = obj$mcbframe

  method = c()
  width = c()
  bcr = c()

  if(is.null(obj$width)){
    tag = 1
  }else{
    width = c(width, obj$width)
    bcr = c(bcr, obj$bcr)
    tag = 0
  }

  if(tag==1){
    method = c()
    width = c()
    bcr = c()
    for(i in 1:length(names(obj))){
      method = c(method, rep(c(names(obj)[i]),length(obj[[1]]$width)))
      width = c(width, obj[[i]]$width)
      bcr = c(bcr, obj[[i]]$bcr)
    }
  }

  if(tag==1){
    data <- data.frame(method=method, width=width, bcr=bcr)
    data$width <- data$width / max(data$width)
  }else{
    data <- data.frame(width=width, bcr=bcr)
    data$width <- data$width / max(data$width)
  }

  if(tag==1){
    gp <- ggplot(data, aes(width, bcr, colour = as.factor(method))) +
      geom_line(stat = "identity") +
      labs(x = "w/p", y = "bcr", title = "MUC") +
      scale_colour_discrete(name = "Method") +
      ylim(limits=c(0,1)) +
      xlim(limits=c(0,1))
    return(gp)
  }else{
    gp <- ggplot(data, aes(width, bcr)) +
      geom_line(stat = "identity",colour="red") +
      labs(x = "w/p", y = "bcr", title = "MUC") +
      ylim(limits=c(0,1)) +
      xlim(limits=c(0,1))
    return(gp)
  }

}
