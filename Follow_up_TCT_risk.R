library(phyloseq)
library(microbiome)
library(metafor)
library(metap)
library(tidyverse)
library(foreach)
library(crayon)
library(vegan)
library(MASS)
library(lme4)
library(ggsci)
library(ggpubr)
library(ggVennDiagram)
library(car)
library(pheatmap)
library(ggrepel)
library(picante)

### ===============================================================================================
### set cox regression functions
Assess_cox<-function(dat,time,outcome,exp,cov,strat){
  tmp<-dat[,c(time,outcome,exp,cov,strat)]
  tmp<-na.omit(tmp)
  covid<-paste0("Cov",c(1:length(cov)))
  
  if(is.null(strat)){
    colnames(tmp)<-c("time","out","exp",covid)
    formula<-as.formula(paste0("Surv(time,out==1)~exp+",paste(covid,collapse = "+")))
  }else{
    colnames(tmp)<-c("time","out","exp",covid,"strat")
    formula<-as.formula(paste0("Surv(time,out==1)~exp+",paste(covid,collapse = "+"),"+strata(strat)"))
  }
  
  cox.hpv<-coxph(formula = formula,data = tmp)
  coef<-summary(cox.hpv)
  if(is.null(strat)){
    res<-data.frame(exp=exp,out=outcome,strat=NA,
                    HR=coef$conf.int[1,1],
                    CI_low=coef$conf.int[1,3],
                    CI_high=coef$conf.int[1,4],
                    pval=coef$coefficients[1,5],
                    N=dim(tmp)[1],Ncase.out=length(which(tmp$out==1)))
  }else{
    res<-data.frame(exp=exp,out=outcome,strat=strat,
                    HR=coef$conf.int[1,1],
                    CI_low=coef$conf.int[1,3],
                    CI_high=coef$conf.int[1,4],
                    pval=coef$coefficients[1,5],
                    N=dim(tmp)[1],Ncase.out=length(which(tmp$out==1)))
  }
  
  return(res)
  
}

### ==================================================
### (1) baseline MH -> HPV progress (00-01,10-11)
### ==================================================
tmp0<-dat.full[which(dat.full$HPV==0),]
res.mh.hpv0<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_HPV",exp = "MH",
                        cov = c("AGE","BMI","center"),strat = NULL)
res.mh.hpv0.sc<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_HPV",exp = "MH_SC",
                        cov = c("AGE","BMI","center"),strat = NULL)

tmp1<-dat.full[which(dat.full$HPV==1),]
res.mh.hpv1<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_HPV",exp = "MH",
                        cov = c("AGE","BMI","center"),strat = NULL)

res.mh.hpv<-rbind(res.mh.hpv0,res.mh.hpv1)
res.mh.hpv$out<-c("HPV(00-01)","HPV(10-11)")
res.mh.hpv$group<-"all samples"

### ==================================================
### (2) baseline shannon -> MH/HPV (00-01,10-11)
### ==================================================

### a. HPV
tmp0<-dat.full[which(dat.full$HPV==0),]
res.alpha.hpv0<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_HPV",exp = "Shannon",
                        cov = c("AGE","BMI","center"),strat = NULL)

tmp1<-dat.full[which(dat.full$HPV==1),]
res.alpha.hpv1<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_HPV",exp = "Shannon",
                           cov = c("AGE","BMI","center"),strat = NULL)

### b. MH
tmp0<-dat.full[which(dat.full$MH==0),]
res.alpha.mh0<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_MH",exp = "Shannon",
                           cov = c("AGE","BMI","center"),strat = NULL)

tmp1<-dat.full[which(dat.full$MH==1),]
res.alpha.mh1<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_MH",exp = "Shannon",
                           cov = c("AGE","BMI","center"),strat = NULL)

res.alpha.sti<-rbind(res.alpha.hpv0,res.alpha.hpv1,res.alpha.mh0,res.alpha.mh1)
res.alpha.sti$out<-c("HPV(00-01)","HPV(10-11)","MH(00-01)","MH(10-11)")
res.alpha.sti$group<-"all samples"

### =========================================================================
### (3) baseline MH/HPV/Shannon -> TCT progress (in HPV+/HPV- samples)
### =========================================================================

### a. all samples TCT00-01,10-11
tmp0<-dat.full[which(dat.full$TCT_CAT==0),]
res.hpv.tct0<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "HPV",
                           cov = c("AGE","BMI","center"),strat = NULL)
res.mh.tct0<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "MH",
                         cov = c("AGE","BMI","center"),strat = NULL)
res.alpha.tct0<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "Shannon",
                         cov = c("AGE","BMI","center"),strat = NULL)

tmp1<-dat.full[which(dat.full$TCT_CAT==1),]
res.hpv.tct1<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "HPV",
                         cov = c("AGE","BMI","center"),strat = NULL)
res.mh.tct1<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "MH",
                        cov = c("AGE","BMI","center"),strat = NULL)
res.alpha.tct1<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "Shannon",
                           cov = c("AGE","BMI","center"),strat = NULL)

res.tct.all<-rbind(res.hpv.tct0,res.mh.tct0,res.alpha.tct0,
                   res.hpv.tct1,res.mh.tct1,res.alpha.tct1)
res.tct.all$out<-c("TCT(00-01)","TCT(00-01)","TCT(00-01)","TCT(10-11)","TCT(10-11)","TCT(10-11)")
res.tct.all$group<-"all samples"

### b. HPV- samples TCT00-01,10-11
tmp0<-dat.full[which(dat.full$TCT_CAT==0 & dat.full$HPV==0),]
res.mh.tct0.hn<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "MH",
                        cov = c("AGE","BMI","center"),strat = NULL)
res.alpha.tct0.hn<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "Shannon",
                           cov = c("AGE","BMI","center"),strat = NULL)

tmp1<-dat.full[which(dat.full$TCT_CAT==1 & dat.full$HPV==0),]
res.mh.tct1.hn<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "MH",
                        cov = c("AGE","BMI","center"),strat = NULL)
res.alpha.tct1.hn<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "Shannon",
                           cov = c("AGE","BMI","center"),strat = NULL)

res.tct.hn<-rbind(res.mh.tct0.hn,res.alpha.tct0.hn,res.mh.tct1.hn,res.alpha.tct1.hn)
res.tct.hn$out<-c("TCT(00-01)","TCT(00-01)","TCT(10-11)","TCT(10-11)")
res.tct.hn$group<-"HPV- samples"

### c. HPV+ samples TCT00-01,10-11
tmp0<-dat.full[which(dat.full$TCT_CAT==0 & dat.full$HPV==1),]
res.mh.tct0.hp<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "MH",
                           cov = c("AGE","BMI","center"),strat = NULL)
res.alpha.tct0.hp<-Assess_cox(dat = tmp0,time = "TIME",outcome = "F_TCT_CAT",exp = "Shannon",
                              cov = c("AGE","BMI","center"),strat = NULL)

tmp1<-dat.full[which(dat.full$TCT_CAT==1 & dat.full$HPV==1),]
res.mh.tct1.hp<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "MH",
                           cov = c("AGE","BMI","center"),strat = NULL)
res.alpha.tct1.hp<-Assess_cox(dat = tmp1,time = "TIME",outcome = "F_TCT_CAT",exp = "Shannon",
                              cov = c("AGE","BMI","center"),strat = NULL)

res.tct.hp<-rbind(res.mh.tct0.hp,res.alpha.tct0.hp,res.mh.tct1.hp,res.alpha.tct1.hp)
res.tct.hp$out<-c("TCT(00-01)","TCT(00-01)","TCT(10-11)","TCT(10-11)")
res.tct.hp$group<-"HPV+ samples"

res.tct<-rbind(res.tct.all,res.tct.hn,res.tct.hp)

