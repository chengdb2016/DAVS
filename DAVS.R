############################################################################################
# DAVS.R
#
# The main functions for DAVS
#
# Code by
#
#  - Debo cheng (chedy055@mymail.unisa.edu.au)
#
# Current version: 18/10/2019
# First version: 12/11/2021

#' @title Towards precise causal effect estimation from data with hidden variables
#'
#' @description
#'
#'
#' @param dat       dataset
#' @param Q         the coso variable 
#' @param pag       the learned PAG 
#' @param alpha     significance level
#' @param models    match or logitreg or linearreg or matchit
#' @param method    subclass, cem, nearest
#' @param subclass  for the method of subclass.
#'
#'
#' @return A list containing several outcomes of interest:
#'   \item{\code{DAVS_ACE}}{a suggested average causal effect}
#'   \item{\code{DAVS_SD}}{standard deviation}
#'   \item{\code{tetrad_scores}}{a list of corresponding scores for each candidate causal effect}
#'   \item{\code{valid Z}}{the valid adjustment sets}
#'   \item{\code{runtime]}{runing time}
#'
#' @details
#' This provides a toy example on how to run a simulation study using the main method
#' developed in this package. See the description of \code{Run.toy.experiment} for more details.
#'
#' @export


Davs.con.causaleffect<-function(dat,Q,pag,alpha=0.05,models="match",
                          method = "subclass", subclass=6, type="ps"){
  ###  dat:observational data; PCv: PCwy and candidate s
  ###  type      "iv" (for inverse variance) or "PS" (for propensity score)
  # require(pcalg)
  starttime<-proc.time()
  pag <- matrix(pag,nrow(pag),ncol(pag),dimnames = NULL)
  w.pos <- ncol(expdat) - 1 
  y.pos <- ncol(expdat)
  pdes <- possibleDe(pag,w.pos) ## trying this out 11/27/2016
  if (!(y.pos %in% pdes)) return(0)
  # possible d-sep
  if(possdsep=="small"){
    pdsep <- union(pdsepset.reach(w.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
  }
  if(possdsep=="big"){
    pdsep <- union(reach(w.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
  }
  pdsep <- sort(setdiff(pdsep,c(w.pos,y.pos)))
  
  pdsepset <- as.vector(setdiff(pdsep,c(Q,pdes)))
  
  valid_Z<-list()
  DAVS.ACE<-list()
  DAVS.sd<-list()
  
  tempdata<-expdat[pdsepset]
  VarQ <- expdat[Q]
  Tr<-expdat$Tr;Y<-expdat$Y
  tempdat<-cbind(VarQ,tempdata, Tr,Y)
  suffStat<-list(C=cor(tempdat),n=nrow(tempdat))
  ordnum <- 1; Q <- 1;nvar<-2:(length(pdsepset)+1)
  
  ###    Rule 2  #### 
  Pvaluewy<-gaussCItest(ncol(tempdat)-1,ncol(tempdat),c(),suffStat)
  Pvaluesw<-gaussCItest(Q,ncol(tempdat)-1,c(),suffStat)
  Pvaluesy<-gaussCItest(Q,ncol(tempdat),c(),suffStat)
  if(Pvaluewy > alpha || (Pvaluesw < alpha && Pvaluesy > alpha)){
    DAVS.ACE[[ordnum]]<-NULL;DAVS.sd[[ordnum]]<-NULL
    ordnum <- ordnum+1
    # cat("The causal effect of W on Y does not be estimated from this dataset.","\n")
  }else{
    ##Rule 1  ####
    # create the size is equal to 1 Called C1(candidate set)  --- L1
    L<-combn(nvar,1)
    L.tmp <- list() 
    for (ii in 1:ncol(L)) {
      Z<-L[,ii]
      Pvalue01<-gaussCItest(Q,ncol(tempdat),Z,suffStat)
      # cat("P-VALUE01:",Pvalue01,"\n")
      if(Pvalue01<alpha){
        M<-c(ncol(tempdat)-1,Z) # W \cup Z
        Pvalue02<-gaussCItest(Q,ncol(tempdat),M,suffStat)
        if(Pvalue02>alpha){
          if(!(list(pdsepset[Z-1]) %fin% valid_Z)){
            cat("Q:",Q,"\n");cat("Z:",pdsepset[Z-1],"\n")
            valid_Z[[ordnum]]<-pdsepset[Z-1]
            L.tmp[[ordnum]]<-L[,ii]  
            if(models == "match"){
              rest<-est_match_con(Z,tempdat,type) 
            } 
            if(models == "logitreg"){
              rest<-est_reg_bin(Z,tempdat) 
            }
            if(models == "linearreg"){
              rest<-est_reg_con(Z,tempdat) 
            }
            if(models == "matchit" && method == "subclass"){
              rest<-est_matchit_con(Z, tempdat, method = "subclass", subclass=subclass)
            }
            if(models == "matchit" && method == "nearest"){
              rest<-est_matchit_con(Z, tempdat, method = "nearest")
            }
            if(models == "matchit" && method == "cem"){
              rest<-est_matchit_con(Z, tempdat, method = "cem", subclass=subclass)
            }
            DAVS.ACE[[ordnum]]<-rest[1]
            DAVS.sd[[ordnum]]<-rest[2]
            ordnum<-ordnum+1 
          }
        }
      }
    }
    Fk<-setdiff(L, L.tmp)
    if(length(Fk)==2){
      Pvalue01<-gaussCItest(Q,ncol(tempdat),Fk,suffStat)
      if(Pvalue01<alpha){
        M<-c(ncol(tempdat)-1,Fk) # W \cup Z
        Pvalue02<-gaussCItest(Q,ncol(tempdat),Fk,suffStat) 
        if(Pvalue02>alpha){
          if(!(list(pdsepset[Z]) %fin% valid_Z)){
            valid_Z[[ordnum]]<-pdsepset[Fk] 
            if(models == "match"){
              rest<-est_match_con(Z,tempdat,type) 
            } 
            if(models == "logitreg"){
              rest<-est_reg_bin(Z,tempdat) 
            }
            if(models == "linearreg"){
              rest<-est_reg_con(Z,tempdat) 
            }
            if(models == "matchit" && method == "subclass"){
              rest<-est_matchit_con(Z,tempdat, method="subclass",subclass=subclass)
            }
            if(models == "matchit" && method == "nearest"){
              rest<-est_matchit_con(Z, tempdat, method = "nearest")
            }
            if(models == "matchit" && method == "cem"){
              rest<-est_matchit_con(Z,tempdat, method="cem",subclass=subclass)
            }
            DAVS.ACE[[ordnum]]<-rest[1]
            DAVS.sd[[ordnum]]<-rest[2]
            ordnum<-ordnum+1 
          }
        }
      }
    }
    if(length(Fk)>2){
      for (k in 2:(length(Fk))) {
        Ck<-create_Ck(Fk,k-1)   # Candidate generation functiOn
        if(length(Ck)==0){
          break
        }
        Ltmp <- list()
        for (jj in 1:length(Ck)) {
          Z<-Ck[[jj]]
          Pvalue01<-gaussCItest(Q,ncol(tempdat),Z,suffStat)
          if(Pvalue01<alpha){
            M<-c(ncol(tempdat)-1,Z) # W \cup Z
            Pvalue02<-gaussCItest(Q,ncol(tempdat),M,suffStat) 
            if(Pvalue02>alpha){
              if(!(list(pdsepset[Z]) %fin% valid_Z)){
                valid_Z[[ordnum]]<-pdsepset[Z] 
                if(models == "match"){
                  rest<-est_match_con(Z,tempdat,type) 
                } 
                if(models == "logitreg"){
                  rest<-est_reg_bin(Z,tempdat) 
                }
                if(models == "linearreg"){
                  rest<-est_reg_con(Z,tempdat) 
                }
                if(models == "matchit" && method == "subclass"){
                  rest<-est_matchit_con(Z,tempdat,method="subclass",subclass=subclass)
                }
                if(models == "matchit" && method == "nearest"){
                  rest<-est_matchit_con(Z,tempdat,method="nearest")
                }
                if(models == "matchit" && method == "cem"){
                  rest<-est_matchit_con(Z,tempdat,method="cem",subclass=subclass)
                }
                DAVS.ACE[[ordnum]]<-rest[1]
                DAVS.sd[[ordnum]]<-rest[2]
                ordnum<-ordnum+1 
              }
            }else{
              Ltmp[[jj]]<-Ck[[jj]]
            }
          }
        }
        if(length(Ltmp)==0){
          break
        }
        Ltmp1 = Ltmp[-which(sapply(Ltmp, is.null))]
        if(length(Ltmp1)==0){
          Fk<-Ck 
        }else{
          Fk<-Ltmp1  
        }
        if(length(Fk)==1){
          break
        }
      }
    }
  }
  
  if(length(DAVS.ACE)==0 || length(DAVS.sd)==0){
    DAVS_ACE<-list(); ACE_DCE<-vector();SD_DCE<-vector() 
  }else{
    DAVS_ACE<-unlist(DAVS.ACE)
    ACE_DCE<-mean(DAVS_ACE)
    DAVS_sd<-unlist(DAVS.sd)
    SD_DCE<-mean(DAVS_sd)
  }
  runtime <- proc.time()-starttime
  Runtime<-runtime[1]
  
  retu<-list(DAVS_ACE,ACE_DCE,SD_DCE,valid_Z,Runtime)
  names(retu) <- c("DAVS_ACE","ACE_DCE","SD_DCE","valid_Z","Runtime")
  return(retu)
}

Davs.bin.causaleffect<-function(expdat,Q,pag,alpha=0.05,possdsep="small",models="match", type="ps"){
  ###  dat:  observational binary data; adjv: adjwy and candidate s
  ### models: match, regression, 
  ###  type      "iv" (for inverse variance) or "PS" (for propensity score)
  # require(pcalg)
  starttime<-proc.time()
  pag <- matrix(pag,nrow(pag),ncol(pag),dimnames = NULL)
  w.pos <- ncol(expdat) - 1 
  y.pos <- ncol(expdat)
  pdes <- possibleDe(pag,w.pos) ## trying this out 11/27/2016
  if (!(y.pos %in% pdes)) return(0)
  # possible d-sep
  if(possdsep=="small"){
    pdsep <- union(pdsepset.reach(w.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
  }
  if(possdsep=="big"){
    pdsep <- union(reach(w.pos,y.pos,-1,pag),y.pos) ##added y.pos to this set 8/29
  }
  pdsep <- sort(setdiff(pdsep,c(w.pos,y.pos)))
  
  pdsepset <- as.vector(setdiff(pdsep,c(Q,pdes)))
  
  valid_Z<-list()
  DAVS.ACE<-list()
  DAVS.sd<-list()
  tempdata<-expdat[pdsepset]
  VarQ <- expdat[Q]
  Tr<-expdat$Tr;Y<-expdat$Y
  tempdat<-cbind(VarQ,tempdata, Tr,Y)
  suffStat<-list(dm=tempdat,adaptDF=FALSE)
  ordnum <- 1; Q <- 1;nvar<-2:(length(pdsepset)+1)
  
  ###    Rule 2  #### 
  Pvaluewy<-binCItest(ncol(tempdat)-1,ncol(tempdat),c(),suffStat)
  Pvaluesw<-binCItest(Q,ncol(tempdat)-1,c(),suffStat)
  Pvaluesy<-binCItest(Q,ncol(tempdat),c(),suffStat)
  if(Pvaluewy > alpha || (Pvaluesw < alpha && Pvaluesy > alpha)){
    DAVS.ACE<-NULL;DAVS.sd<-NULL 
    # cat("The causal effect of W on Y does not be estimated from this dataset.","\n")
  }else{
    ##Rule 1  ####
    # create the size is equal to 1 Called C1(candidate set)  --- L1
    if(length(nvar)==1){
      Z<-nvar
      L.tmp <- list() 
      Pvalue01<-binCItest(Q,ncol(tempdat),Z,suffStat)
      # cat("P-VALUE01:",Pvalue01,"\n")
      if(Pvalue01<alpha){
        M<-c(ncol(tempdat)-1,Z) # W \cup Z
        Pvalue02<-binCItest(Q,ncol(tempdat),M,suffStat) 
        if(Pvalue02>alpha){
          valid_Z[[ordnum]]<-pdsepset[Z-1]
          L.tmp[[ordnum]]<-Z
          rest<-est_reg_bin(Z,tempdat) 
          DAVS.ACE[[ordnum]]<-rest[1]
          DAVS.sd[[ordnum]]<-rest[2] 
          ordnum<-ordnum+1 
        }
      }
    }else{
      L<-combn(nvar,1) 
      L.tmp <- list() 
      for (ii in 1:ncol(L)) {
        Z<-L[,ii]
        Pvalue01<-binCItest(Q,ncol(tempdat),Z,suffStat)
        # cat("P-VALUE01:",Pvalue01,"\n")
        if(Pvalue01<alpha){
          M<-c(ncol(tempdat)-1,Z) # W \cup Z
          Pvalue02<-binCItest(Q,ncol(tempdat),M,suffStat) 
          if(Pvalue02>alpha){
            if(!(list(pdsepset[Z-1]) %fin% valid_Z)){
              valid_Z[[ordnum]]<-pdsepset[Z-1]
              L.tmp[[ordnum]]<-L[,ii] 
              if(models == "match"){
                rest<-est_match_bin(Z,tempdat,type) 
                DAVS.ACE[[ordnum]]<-rest[1]
                DAVS.sd[[ordnum]]<-rest[2]
              } 
              if(models == "regression"){
                rest<-est_reg_bin(Z,tempdat) 
                DAVS.ACE[[ordnum]]<-rest[1]
                DAVS.sd[[ordnum]]<-rest[2] 
              }
              if(models == "matchit"){
                fmla1<-as.formula(paste("Tr~",paste(pdsepset[Z-1],collapse = "+")))
                txnam<-c("Tr",pdsepset[Z-1])
                fmla2<-as.formula(paste("Y~",paste(txnam,collapse = "+")))
                m.out <- matchit(fmla1,data = tempdat, method = "subclass")
                zmodel.out<-zelig(fmla2,data=match.data(m.out),model = "logit", cite = FALSE)
                control.out<-setx(zmodel.out,Tr=0) #control
                treat.out<-setx(zmodel.out,Tr=1) #Ted
                s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)
                DAVS.ACE[[ordnum]]<-mean(sim(zmodel.out,x=treat.out)$get_qi())-mean(sim(zmodel.out,x=control.out)$get_qi())
                DAVS.sd[[ordnum]]<-abs(sd(sim(zmodel.out,x=treat.out)$get_qi()-sim(zmodel.out,x=control.out)$get_qi()))
              }
              ordnum<-ordnum+1 
            }
          }
        }
      } 
      Fk<-setdiff(L, L.tmp)
      if(length(Fk)==2){
        Pvalue01<-binCItest(Q,ncol(tempdat),Fk,suffStat)
        if(Pvalue01<alpha){
          M<-c(ncol(tempdat)-1,Fk) # W \cup Z
          Pvalue02<-binCItest(Q,ncol(tempdat),Fk,suffStat) 
          if(Pvalue02>alpha){
            if(!(list(pdsepset[Z-1]) %fin% valid_Z)){
              valid_Z[[ordnum]]<-pdsepset[Fk] 
              if(models == "match"){
                rest<-est_match_bin(Z,tempdat,type) 
                DAVS.ACE[[ordnum]]<-rest[1]
                DAVS.sd[[ordnum]]<-rest[2]
              } 
              if(models == "regression"){
                rest<-est_reg_bin(Z,tempdat) 
                DAVS.ACE[[ordnum]]<-rest[1]
                DAVS.sd[[ordnum]]<-rest[2] 
              }
              if(models == "matchit"){
                fmla1<-as.formula(paste("Tr~",paste(AdjWY[Z],collapse = "+")))
                txnam<-c("Tr",pdsepset[Z])
                fmla2<-as.formula(paste("Y~",paste(txnam,collapse = "+")))
                m.out <- matchit(fmla1,data = tempdat, method = "subclass")
                zmodel.out<-zelig(fmla2,data=match.data(m.out),model = "logit", cite = FALSE)
                control.out<-setx(zmodel.out,Tr=0) #control
                treat.out<-setx(zmodel.out,Tr=1) #Ted
                s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)
                DAVS.ACE[[ordnum]]<-mean(sim(zmodel.out,x=treat.out)$get_qi())-mean(sim(zmodel.out,x=control.out)$get_qi())
                DAVS.sd[[ordnum]]<-abs(sd(sim(zmodel.out,x=treat.out)$get_qi()-sim(zmodel.out,x=control.out)$get_qi()))
              }
              ordnum<-ordnum+1 
            }
          }
        }
      }
      if(length(Fk)>2){
        for (k in 2:(length(Fk))) {
          Ck<-create_Ck(Fk,k-1)   # Candidate generation functiOn
          if(length(Ck)==0){
            break
          }
          Ltmp <- list()
          for (jj in 1:length(Ck)) {
            Z<-Ck[[jj]]
            Pvalue01<-binCItest(Q,ncol(tempdat),Z,suffStat)
            if(Pvalue01<alpha){
              M<-c(ncol(tempdat)-1,Z) # W \cup Z
              Pvalue02<-binCItest(Q,ncol(tempdat),M,suffStat) 
              if(Pvalue02>alpha){
                if(!(list(pdsepset[Z]) %fin% valid_Z)){
                  valid_Z[[ordnum]]<-pdsepset[Z] 
                  if(models == "match"){
                    rest<-est_match_bin(Z,tempdat,type) 
                    DAVS.ACE[[ordnum]]<-rest[1]
                    DAVS.sd[[ordnum]]<-rest[2]
                  } 
                  if(models == "regression"){
                    rest<-est_reg_bin(Z,tempdat) 
                    DAVS.ACE[[ordnum]]<-rest[1]
                    DAVS.sd[[ordnum]]<-rest[2] 
                  }
                  if(models == "matchit"){
                    fmla1<-as.formula(paste("Tr~",paste(AdjWY[Z],collapse = "+")))
                    txnam<-c("Tr",pdsepset[Z])
                    fmla2<-as.formula(paste("Y~",paste(txnam,collapse = "+")))
                    m.out <- matchit(fmla1,data = tempdat, method = "subclass")
                    zmodel.out<-zelig(fmla2,data=match.data(m.out),model = "logit", cite = FALSE)
                    control.out<-setx(zmodel.out,Tr=0) #control
                    treat.out<-setx(zmodel.out,Tr=1) #Ted
                    s.out<-sim(zmodel.out, x = control.out, x1 = treat.out)
                    DAVS.ACE[[ordnum]]<-mean(sim(zmodel.out,x=treat.out)$get_qi())-mean(sim(zmodel.out,x=control.out)$get_qi())
                    DAVS.sd[[ordnum]]<-abs(sd(sim(zmodel.out,x=treat.out)$get_qi()-sim(zmodel.out,x=control.out)$get_qi()))
                  }
                  ordnum<-ordnum+1 
                }
              }else{
                Ltmp[[jj]]<-Ck[[jj]]
              }
            }
          }
          if(length(Ltmp)==0){
            break
          }
          Ltmp1 = Ltmp[-which(sapply(Ltmp, is.null))]
          if(length(Ltmp1)==0){
            Fk<-Ck 
          }else{
            Fk<-Ltmp1  
          }
          if(length(Fk)==1){
            break
          }
        }
      }
    }
  }
  if(length(DAVS.ACE)==0 || length(DAVS.sd)==0){
    DAVS_ACE<-list();ACE_DCE<-vector();SD_DCE<-vector() 
  }else{
    DAVS_ACE<-unlist(DAVS.ACE)
    ACE_DCE<-mean(DAVS_ACE)
    DAVS_sd<-unlist(DAVS.sd)
    SD_DCE<-mean(DAVS_sd)
  }
  runtime<-proc.time() -starttime
  Runtime<-runtime[1]  
  
  retu<-list(DAVS_ACE,ACE_DCE,SD_DCE,valid_Z,Runtime)
  names(retu) <- c("DAVS_ACE","ACE_DCE","SD_DCE","valid_Z","Runtime")
  
  return(retu)
}
