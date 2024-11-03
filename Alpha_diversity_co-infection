library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggrepel)
library(matrixStats)
library(ggsci)
library(ggpubr)
library(vegan)
library(ade4)
library(foreach)
library(crayon)
library(ggnewscale)
library(picante)

## ========================================================================================================================================
## association Infection-free vs sti, and Infection-free vs sti-mono

alpha_stat<-function(sti,dat,col.alpha,col.STIn,col.HPVn,cov,mono){
  
  tmp<-dat[,c(col.alpha,sti,col.STIn,col.HPVn,cov)]
  colid<-LETTERS[seq(from=1,to=length(cov))]
  colnames(tmp)<-c("alpha","test","N1","N2",colid)
  
  if(mono==T){
    if(sti %in% hpv.subtype){
      tmp$pheno<-NA; tmp$pheno[which(tmp$test==1 & tmp$N1==1 & tmp$N2==1)]<-1; tmp$pheno[which(tmp$N1==0)]<-0
    }else{
      tmp$pheno<-NA; tmp$pheno[which(tmp$test==1 & tmp$N1==1)]<-1; tmp$pheno[which(tmp$N1==0)]<-0
    }
  }else{
    tmp$pheno<-NA; tmp$pheno[which(tmp$test==1)]<-1; tmp$pheno[which(tmp$N1==0)]<-0
  }

  ## mean of alpha
  a1<-mean(tmp$alpha[which(tmp$pheno==1)])
  ## associate alpha to infection phenotype
  tmp.data<-na.omit(tmp[,c("alpha","pheno",colid)])
  formula<-as.formula(paste("pheno~", "alpha+",paste(colid, collapse="+")))
  tmp.mod1=glm(formula = formula,data = tmp.data,family = "binomial")
  coef<-summary(tmp.mod1)
  
  tmp.result=data.frame(STI1=sti,mean.alpha=a1,
                        estimate=coef$coefficients[2,1],se=coef$coefficients[2,2],
                        z=coef$coefficients[2,3],Pval=coef$coefficients[2,4],N=length(which(tmp$pheno==1)))
  
  return(tmp.result)
}


## ========================================================================================================================================
## box plot of alpha diversity

## general infection data
tmp1<-df.pheno[which(df.pheno$N_infect==0),c("sampleid","diversity_shannon","pd")]
tmp1$group<-"Infection-free"

df.plot1<-tmp1
for(i in col.sti){
  tmp<-df.pheno[which(df.pheno[,i]==1),c("sampleid","diversity_shannon","pd")]
  tmp$group<-i
  df.plot1<-rbind(df.plot1,tmp)
}

df.plot1$group<-factor(df.plot1$group,levels = unique(df.plot1$group))

## mono infection data
df.plot2<-tmp1
for(i in col.sti){
  if(i %in% c("HPV16","HPV18","HPV52","HPV58")){
    tmp<-df.pheno[which(df.pheno[,i]==1 & df.pheno$Mono_infection=="Mono-infection" & df.pheno$N_HPV==1),c("sampleid","diversity_shannon","pd")]
  }else{
    tmp<-df.pheno[which(df.pheno[,i]==1 & df.pheno$Mono_infection=="Mono-infection"),c("sampleid","diversity_shannon","pd")]
  }
  tmp$group<-i
  df.plot2<-rbind(df.plot2,tmp)
}

df.plot2$group<-factor(df.plot2$group,levels = unique(df.plot2$group))

### (1) generate new data frame for plot
## shannon
df.plot10<-df.plot1[,c("sampleid","diversity_shannon","group")]
df.plot10$sig1<-result.shannon$Pval[match(df.plot10$group,result.shannon$STI1)]
df.plot10$type="Infection"

df.plot11<-df.plot2[,c("sampleid","diversity_shannon","group")]
df.plot11$sig1<-result.shannon_momo$Pval[match(df.plot11$group,result.shannon_momo$STI1)]
df.plot11$type="Mono-Infection"

df.plot.a<-rbind(df.plot10,df.plot11)
df.plot.a$sig1[which(df.plot.a$sig1>=0.05)]<-1
df.plot.a$sig1[which(df.plot.a$sig1<0.05)]<-2
df.plot.a$sig1[which(df.plot.a$group=="Infection-free")]<-0
df.plot.a$sig1<-as.factor(df.plot.a$sig1)
levels(df.plot.a$sig1)<-c("Reference","NS","P<0.05")
df.plot.a$diversity<-"Shannon diversity"

## plot Phylogenetic diversity
## pd
df.plot20<-df.plot1[,c("sampleid","pd","group")]
df.plot20$sig1<-result.pd$Pval[match(df.plot20$group,result.pd$STI1)]
df.plot20$type="Infection"

df.plot21<-df.plot2[,c("sampleid","pd","group")]
df.plot21$sig1<-result.pd_momo$Pval[match(df.plot21$group,result.pd_momo$STI1)]
df.plot21$type="Mono-Infection"

df.plot.b<-rbind(df.plot20,df.plot21)
df.plot.b$sig1[which(df.plot.b$sig1>=0.05)]<-1
df.plot.b$sig1[which(df.plot.b$sig1<0.05)]<-2
df.plot.b$sig1[which(df.plot.b$group=="Infection-free")]<-0
df.plot.b$sig1<-as.factor(df.plot.b$sig1)
levels(df.plot.b$sig1)<-c("Reference","NS","P<0.05")
df.plot.b$diversity<-"Phylogenetic diversity"

colnames(df.plot.a)[2]<-"alpha"
colnames(df.plot.b)[2]<-"alpha"

give.n <- function(x){
  return(c(y = median(x)*1.1, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

p1<-ggplot(df.plot.a,aes(x=group,y=alpha,fill=sig1))+
  geom_boxplot(alpha=0.8)+
  #geom_jitter(alpha=0.4)+
  geom_hline(yintercept = median(tmp1$diversity_shannon), linetype="dashed", color = "gray3")+
  theme_bw()+xlab("")+ylab("Alpha diversity")+ggtitle("Shannon diversity")+
  facet_grid(type~.,scales = "free",space = "free")+
  scale_fill_manual(values = c("Reference"="#6F99ADFF","NS"="#FFDC91FF","P<0.05"="#BC3C29FF"))+
  theme(legend.position = "bottom")+labs(fill="")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75))

p2<-ggplot(df.plot.b,aes(x=group,y=alpha,fill=sig1))+
  geom_boxplot(alpha=0.8)+
  #geom_jitter(alpha=0.4)+
  geom_hline(yintercept = median(tmp1$pd), linetype="dashed", color = "gray3")+
  theme_bw()+xlab("")+ylab("Alpha diversity")+ggtitle("Phylogenetic diversity")+
  facet_grid(type~.,scales = "free",space = "free")+
  scale_fill_manual(values = c("Reference"="#6F99ADFF","NS"="#FFDC91FF","P<0.05"="#BC3C29FF"))+
  theme(legend.position = "bottom")+labs(fill="")+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median, 
               position = position_dodge(width = 0.75))




