########################################################################
#####Cheng's estimators for estimating causal effects ##################
####################################  ###############################

est_reg_con<-function(adjset,obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {form <- formula("Y~Tr")} else {
      form <- formula(paste("Y~Tr", paste(adjcov, collapse="+"), sep="+"))}
    fit <- lm(form, data=obs)
    est_ATE <- summary(fit)$coefficients["Tr", "Estimate"]
    est_sd <- summary(fit)$coefficients["Tr", "Std. Error"]
  }
  return(c(est_ATE, est_sd))
}
est_reg_bin <- function(adjset, obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    return(NA)
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {form <- formula("Y~Tr")} else {
      form <- formula(paste("Y~Tr", paste(adjcov, collapse="+"), sep="+"))}
    fit <- glm(form, data=obs, family="binomial")
    stdreg <- tryCatch(stdGlm(fit, obs, X="Tr"), error=function(x){NA})
    if (anyNA(stdreg)) {return(NA)} else {
      p_0 <- stdreg$est[1]
      p_1 <- stdreg$est[2]
      est_logMCOR <- log(p_1/(1-p_1) / (p_0/(1-p_0)))
    }
    return(est_logMCOR)
  }
}


est_matchit_con<-function(adjset,obs, method = "subclasss", subclass=6){
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else if (sum(adjset)==0) { 
    est_ATE <- mean(obs$Y[obs$Tr==1]) - mean(obs$Y[obs$Tr==0])
    n_1 <- sum(obs$Tr==1)
    n_0 <- sum(obs$Tr==0)
    s2_pooled <- (var(obs$Y[obs$Tr==1])*(n_1-1) +
                    var(obs$Y[obs$Tr==0])*(n_1-0)) / (n_1+n_0-2)
    est_sd <- sqrt(( 1/n_1 + 1/n_0 ) * s2_pooled)
  } else {
    require(Zelig)
    if(method == "subclass"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "subclass", subclass = subclass)
      # learest square regression for con outcome
      zmodel.out<-zelig(fmla2, data=match.data(m.out), model = "ls", cite = FALSE, by="subclass")
      control.out<-setx(zmodel.out,Tr=0) #control
      treat.out<-setx(zmodel.out,Tr=1) #Tread
      s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)          
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
    if(method == "nearest"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "nearest")
      z.out <- zelig(fmla2, data =match.data(m.out), model = "ls",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd)) 
    }
    if(method == "cem"){
      require(cem)
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset], collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      glm1<-glm(fmla1,family = binomial,  data = obs)
      m.out<-matchit(glm1, method = "cem",data = obs)
      z.out <- zelig(fmla2, data =match.data(m.out), model = "ls",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
  }
  
  return(c(est_ATE, est_sd))
}
est_matchit_bin<-function(adjset,obs, method = "subclasss", subclass=6){
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else if (sum(adjset)==0) { 
    est_ATE <- mean(obs$Y[obs$Tr==1]) - mean(obs$Y[obs$Tr==0])
    n_1 <- sum(obs$Tr==1)
    n_0 <- sum(obs$Tr==0)
    s2_pooled <- (var(obs$Y[obs$Tr==1])*(n_1-1) +
                    var(obs$Y[obs$Tr==0])*(n_1-0)) / (n_1+n_0-2)
    est_sd <- sqrt(( 1/n_1 + 1/n_0 ) * s2_pooled)
  } else {
    require(Zelig)
    if(method == "subclass"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "subclass", subclass = subclass)
      # learest square regression for con outcome
      zmodel.out<-zelig(fmla2, data=match.data(m.out), model = "logit", cite = FALSE, by="subclass")
      control.out<-setx(zmodel.out,Tr=0) #control
      treat.out<-setx(zmodel.out,Tr=1) #Tread
      s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)          
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
    if(method == "nearest"){
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset],collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      m.out <- matchit(fmla1, data = obs, method = "nearest")
      z.out <- zelig(fmla2, data =match.data(m.out), model = "logit",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd)) 
    }
    if(method == "cem"){
      require(cem)
      fmla1<-as.formula(paste("Tr~", paste(covs[adjset], collapse = "+")))
      txnam<-c("Tr", covs[adjset])
      fmla2<-as.formula(paste("Y~", paste(txnam, collapse = "+")))
      glm1<-glm(fmla1,family = binomial,  data = obs)
      m.out<-matchit(glm1, method = "cem",data = obs)
      z.out <- zelig(fmla2, data =match.data(m.out), model = "logit",cite = FALSE)
      x.out <- setx(z.out, Tr = 0)
      x1.out <- setx(z.out, Tr = 1)
      s.out <- sim(z.out, x = x.out, x1 = x1.out)
      
      est_ATE <- mean(unlist(s.out$sim.out$x1$fd))
      est_sd <- sd(unlist(s.out$sim.out$x1$fd))  
    }
  }
  
  return(c(est_ATE, est_sd))
}



est_match_con <- function(adjset, obs, type){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  # type      "iv" (for inverse variance) or "PS" (for propensity score)
  
  # if(!(exists("run"))){run <<- 0}
  # run <<- run+1
  # cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else if (sum(adjset)==0) { 
    est_ATE <- mean(obs$Y[obs$Tr==1]) - mean(obs$Y[obs$Tr==0])
    n_1 <- sum(obs$Tr==1)
    n_0 <- sum(obs$Tr==0)
    s2_pooled <- (var(obs$Y[obs$Tr==1])*(n_1-1) +
                    var(obs$Y[obs$Tr==0])*(n_1-0)) / (n_1+n_0-2)
    est_sd <- sqrt(( 1/n_1 + 1/n_0 ) * s2_pooled)
  } else {
    if (type=="PS") {
      if (sum(adjset)==0) {PSform <- formula("Tr~1")} else {
        PSform <- formula(paste("Tr~", paste(covs[adjset], collapse="+"),
                                sep="")) }
      PSmodel <- glm(PSform, family="binomial", data=obs)
      converged_PS <- PSmodel$converged
      PS <- predict(PSmodel, type="response")  
      covadj <- PS
    } else {
      all <- obs[ ,1:(ncol(obs)-2)]
      covadj <- all[ ,adjset]
    }
    Mres <- Match(Tr=obs$Tr, Y=obs$Y, X=covadj, estimand="ATE", M=1,replace=FALSE)
    est_ATE <- Mres$est
    est_sd <- Mres$se
  }
  return(c(est_ATE, est_sd))
}


est_match_bin <- function(adjset, obs, type){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  # type      "iv" (for inverse variance) or "PS" (for propensity score)
  
  # if(!(exists("run"))){run <<- 0}
  # run <<- run+1
  # cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_lOR <- NA
  } else if (sum(adjset)==0) {
    mu_0 <- mean(obs$Y[obs$Tr==0])
    mu_1 <- mean(obs$Y[obs$Tr==1])
    est_lOR <- log( mu_1/(1-mu_1) / (mu_0/(1-mu_0)) )
  } else {
    if (type=="PS") {
      PSform <- formula(paste("Tr~", paste(covs[adjset], collapse="+"),
                              sep=""))
      PSmodel <- glm(PSform, family="binomial", data=obs)
      converged_PS <- PSmodel$converged
      PS <- predict(PSmodel, type="response")
      names(PS) <- NULL
      covadj <- as.matrix(PS)
    } else {
      all <- obs[ ,1:(ncol(obs)-2)]
      covadj <- all[ ,adjset]
    }
    
    Mres <- tryCatch({Match(Tr=obs$Tr, Y=obs$Y, X=covadj, estimand="ATE",
                            Weight=1)},
                     error=function(e){return(NA)}
    )
    if (anyNA(Mres)) {est_lOR <- NA} else {
      wei <- Mres$weights
      m_Y <- Mres$mdata$Y*wei
      m_Tr <- Mres$mdata$Tr
      m_n <- Mres$orig.nobs
      mu_0 <- sum(m_Y[m_Tr==0])/m_n
      mu_1 <- sum(m_Y[m_Tr==1])/m_n
      est_lOR <- log( mu_1/(1-mu_1) / (mu_0/(1-mu_0)) )
    }
  }
  return(est_lOR)
}

est_dr_con <- function(adjset, obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  # 
  #   if(!(exists("run"))){run <<- 0}
  #   run <<- run+1
  #   cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_ATE <- NA
    est_sd <- NA
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {
      form_Tr <- formula("Tr~1")
      form_Y <- formula("Y~1")
    } else {
      form_Tr <- formula(paste("Tr~", paste(adjcov, collapse="+"), sep="+"))
      form_Y <- formula(paste("Y~", paste(adjcov, collapse="+"), sep="+"))
    }
    ppi.glm <- glm(form_Tr, data=obs, family=binomial)
    X <- model.matrix(ppi.glm)
    ppi.hat <- ppi.glm$fitted
    
    eta1.glm <- glm(form_Y, subset=Tr==1, data=obs)
    eta1.hat <- predict.glm(eta1.glm, type="response", newdata=obs)
    eta0.glm <- glm(form_Y, subset=Tr==0, data=obs)
    eta0.hat <- predict.glm(eta0.glm, type="response", newdata=obs)
    
    #ppi.hat treated as known
    out.lik <- tryCatch({ate.clik(obs$Y, obs$Tr, ppi.hat, g0=cbind(1,eta0.hat),
                                  g1=cbind(1,eta1.hat))}, error=function(e){return(NA)}
    )
    if (anyNA(out.lik)) {est_ATE <- est_sd <- NA} else {
      est_ATE <- out.lik$diff
      est_sd <- sqrt(out.lik$v.diff)
    }
  }
  return(c(est_ATE, est_sd))
}


est_dr_bin <- function(adjset, obs){
  # adjset    selected adjustment set (logical vector)
  # obs       data frame with observation data
  
  # if(!(exists("run"))){run <<- 0}
  # run <<- run+1
  # cat(run, "\n")
  
  covs <- names(obs)[1:(ncol(obs)-2)]
  if (anyNA(adjset)) {
    est_logMCOR <- NA
  } else { 
    adjcov <- covs[adjset]
    if (length(adjcov)==0) {
      form_Tr <- formula("Tr~1")
      form_Y <- formula("Y~1")
    } else {
      form_Tr <- formula(paste("Tr~", paste(adjcov, collapse="+"), sep="+"))
      form_Y <- formula(paste("Y~", paste(adjcov, collapse="+"), sep="+"))
    }
    ppi.glm <- glm(form_Tr, data=obs, family=binomial)
    X <- model.matrix(ppi.glm)
    ppi.hat <- ppi.glm$fitted
    
    eta1.glm <- glm(form_Y, subset=Tr==1, data=obs, family=binomial)
    eta1.hat <- predict.glm(eta1.glm, type="response", newdata=obs)
    eta0.glm <- glm(form_Y, subset=Tr==0, data=obs, family=binomial)
    eta0.hat <- predict.glm(eta0.glm, type="response", newdata=obs)
    
    #ppi.hat treated as known
    out.lik <- tryCatch({ate.clik(obs$Y, obs$Tr, ppi.hat, g0=cbind(1,eta0.hat),
                                  g1=cbind(1,eta1.hat))}, error=function(e){return(NA)}
    )
    if (anyNA(out.lik)) {est_logMCOR <- NA} else {
      p_0 <- out.lik$mu[2]
      p_1 <- out.lik$mu[1]
      est_logMCOR <- log(p_1/(1-p_1) / (p_0/(1-p_0)))
    }
  }
}