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


### function1 general assocaition between taxa and infection
Assess_assoc_taxa=function(mtx.taxa,dat.pheno,tmp.taxa,infection,cov.fix){
  ## merge data
  tmp.data=merge(mtx.taxa[,tmp.taxa,drop=F],dat.pheno[,c(infection,cov.fix),drop=F],by="row.names",all=F)
  tmp.data=na.omit(tmp.data)
  ## generate formular
  if(is.null(cov.fix)){
    colnames(tmp.data)<-c("ID","Bacteria","Infection")
    formula<-as.formula(paste("Infection~", "Bacteria"))
  }else{
    colid<-LETTERS[seq(from=1,to=length(cov.fix))]
    colnames(tmp.data)<-c("ID","Bacteria","Infection",colid)
    formula<-as.formula(paste("Infection~", "Bacteria+",paste(colid, collapse="+")))
  }
  ## run model
  tmp.mod1=glm(formula = formula,data = tmp.data,family = "binomial")
  coef<-summary(tmp.mod1)
  ## summary
  tb<-table(tmp.data$Infection)
  tmp.result=data.frame(Bacteria=tmp.taxa,Infection=infection,
                        estimate=coef$coefficients["Bacteria",1],se=coef$coefficients["Bacteria",2],
                        z=coef$coefficients["Bacteria",3],Pval=coef$coefficients["Bacteria",4],
                        case=tb[[2]],N=tb[[1]]+tb[[2]])
  return(tmp.result)
}



### ----------------------------------------------------------------
### run model1 Infection ~ Taxa (yse/no) + AGE + BMI + center

result.model1 = foreach(i=1:length(testlist1),.combine = rbind) %do%  {
  infection=testlist1[i]
  cov.fix=c("AGE","BMI","center")
  mtx.taxa<-t(taxa)
  
  nn=foreach(n=1:ncol(mtx.taxa),.combine = rbind) %do%  {
    tmp.taxa=colnames(mtx.taxa)[n]
    result=Assess_assoc_taxa(mtx.taxa = mtx.taxa,dat.pheno = df.pheno,tmp.taxa = tmp.taxa,
                             infection = infection,cov.fix = cov.fix)
    cat(yellow(tmp.taxa,"+++++",infection,"\n"))
    
    return.string=result
  }
  cat(black(infection,"\n"))
  return.string=nn
}

result.model1$model<-"model1"
result.model1$FDR<-p.adjust(result.model1$Pval,method = "BH")


### -------------------------------------------------------------------
### run model2 Infection ~ Taxa (Infection-free vs yes) + AGE + BMI + center 

result.model2 = foreach(i=1:length(testlist2),.combine = rbind) %do%  {
  infection=testlist2[i]
  cov.fix=c("AGE","BMI","center")
  mtx.taxa<-t(taxa)
  
  nn=foreach(n=1:ncol(mtx.taxa),.combine = rbind) %do%  {
    tmp.taxa=colnames(mtx.taxa)[n]
    result=Assess_assoc_taxa(mtx.taxa = mtx.taxa,dat.pheno = df.pheno,tmp.taxa = tmp.taxa,
                             infection = infection,cov.fix = cov.fix)
    cat(yellow(tmp.taxa,"+++++",infection,"\n"))
    
    return.string=result
  }
  cat(black(infection,"\n"))
  return.string=nn
}

result.model2$model<-"model2"
result.model2$FDR<-p.adjust(result.model2$Pval,method = "BH")


### -------------------------------------------------------------------
### run model3 Infection ~ Taxa (Infection-free vs mono-infection) + AGE + BMI + center 

result.model3 = foreach(i=1:length(testlist3),.combine = rbind) %do%  {
  infection=testlist3[i]
  cov.fix=c("AGE","BMI","center")
  mtx.taxa<-t(taxa)
  
  nn=foreach(n=1:ncol(mtx.taxa),.combine = rbind) %do%  {
    tmp.taxa=colnames(mtx.taxa)[n]
    result=Assess_assoc_taxa(mtx.taxa = mtx.taxa,dat.pheno = df.pheno,tmp.taxa = tmp.taxa,
                             infection = infection,cov.fix = cov.fix)
    cat(yellow(tmp.taxa,"+++++",infection,"\n"))
    
    return.string=result
  }
  cat(black(infection,"\n"))
  return.string=nn
}

result.model3$model<-"model3"
result.model3$FDR<-p.adjust(result.model3$Pval,method = "BH")

### -------------------------------------------------------------------
### run model4 Infection ~ Taxa (yes vs no) + AGE + BMI + center + MH


result.model4 = foreach(i=1:length(testlist4),.combine = rbind) %do%  {
  infection=testlist4[i]
  cov.fix=c("AGE","BMI","center","MH")
  mtx.taxa<-t(taxa)
  
  nn=foreach(n=1:ncol(mtx.taxa),.combine = rbind) %do%  {
    tmp.taxa=colnames(mtx.taxa)[n]
    result=Assess_assoc_taxa(mtx.taxa = mtx.taxa,dat.pheno = df.pheno,tmp.taxa = tmp.taxa,
                             infection = infection,cov.fix = cov.fix)
    cat(yellow(tmp.taxa,"+++++",infection,"\n"))
    
    return.string=result
  }
  cat(black(infection,"\n"))
  return.string=nn
}

result.model4$model<-"model4"
result.model4$FDR<-p.adjust(result.model4$Pval,method = "BH")


### -------------------------------------------------------------------
### run model5 co-Infection ~ Taxa (infection-free vs co-infection) + AGE + BMI + center
testlist5<-colnames(df.pheno)[105:113]

result.model5 = foreach(i=1:length(testlist5),.combine = rbind) %do%  {
  infection=testlist5[i]
  cov.fix=c("AGE","BMI","center")
  mtx.taxa<-t(taxa)
  
  nn=foreach(n=1:ncol(mtx.taxa),.combine = rbind) %do%  {
    tmp.taxa=colnames(mtx.taxa)[n]
    result=Assess_assoc_taxa(mtx.taxa = mtx.taxa,dat.pheno = df.pheno,tmp.taxa = tmp.taxa,
                             infection = infection,cov.fix = cov.fix)
    cat(yellow(tmp.taxa,"+++++",infection,"\n"))
    
    return.string=result
  }
  cat(black(infection,"\n"))
  return.string=nn
}

result.model5$model<-"model5"
result.model5$model<-"model5"
result.model5$FDR<-p.adjust(result.model5$Pval,method = "BH")



### ===========================================
### plot resutls 

p1<-ggplot(result.sig,aes(x=STI,y=Bacteria,fill=z))+
  geom_tile(colour="gray")+
  scale_fill_gradient2(low="#0052A5", high = "#d11141", mid = "white", midpoint = 0)+
  #scale_fill_manual(values=c("#000099","#0033CC","#0000FF","#6666FF","#9999FF","white","#FFF3B2","#FEB24C","#FC4E2A","#E31A1C","#990000"))+
  geom_text(aes(label=label), size=4,color="black") +
  theme_classic()+
  scale_y_discrete(limits = taxa_order1)+
  labs(fill="Zscore")+
  facet_grid(~model,scales = "free_x",space = "free")+
  #coord_fixed(ratio = 1)+
  theme(#axis.text.y = element_text(hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
    axis.text.y = element_text(face = "italic",size = 12),
    panel.border = element_rect(colour = "black", fill=NA, size=0.8),
    legend.position = "right")+ylab("")+xlab("")
