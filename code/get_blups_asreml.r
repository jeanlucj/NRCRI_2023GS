fitASfunc<-function(fixedFormula,randFormula, MultiTrialTraitData,...){
  # test arguments for function
  # ----------------------
  # MultiTrialTraitData<-dbdata$MultiTrialTraitData[[7]]
  # #Trait<-dbdata$Trait[[7]]
  # fixedFormula<-dbdata$fixedFormula[[7]]
  # randFormula<-dbdata$randFormula[[7]]
  #test<-fitASfunc(fixedFormula,randFormula,MultiTrialTraitData)
  # ----------------------
  MultiTrialTraitData %<>%
    mutate(across(c(GID,yearInLoc,
                    CompleteBlocks,
                    IncompleteBlocks,
                    trialInLocYr,
                    repInTrial,
                    blockInRep),as.factor)) %>% 
    droplevels
  
  require(asreml); 
  fixedFormula<-as.formula(fixedFormula)
  randFormula<-as.formula(randFormula)
  # fit asreml 
  out<-asreml(fixed = fixedFormula,
              random = randFormula,
              data = MultiTrialTraitData, 
              maxiter = 40, workspace=1000e6, 
              na.action = na.method(x="omit"))
  #### extract residuals - Round 1
  
  outliers1<-which(abs(scale(out$residuals))>3.3)
  
  if(length(outliers1)>0){
    
    x<-MultiTrialTraitData[-outliers1,]
    # re-fit
    out<-asreml(fixed = fixedFormula,
                random = randFormula,
                data = x, 
                maxiter = 40, workspace=1000e6, 
                na.action = na.method(x="omit"))
    #### extract residuals - Round 2
    outliers2<-which(abs(scale(out$residuals))>3.3)
    if(length(outliers2)>0){
      #### remove outliers
      x<-x[-outliers2,]
      # final re-fit
      out<-asreml(fixed = fixedFormula,
                  random = randFormula,
                  data = x, maxiter = 40,workspace=1000e6, 
                  na.action = na.method(x="omit"))
    } else outliers2<-NULL
  } else outliers1<-outliers2<-NULL

  ll<-summary(out,all=T)$loglik
  varcomp<-summary(out,all=T)$varcomp
  Vg<-varcomp["GID!GID","component"]
  Ve<-varcomp["units!R","component"]
  H2=Vg/(Vg+Ve)
  blups<-summary(out,coef=T)$coef.random %>%
    as.data.frame %>%
    rownames_to_column(var = "GID") %>%
    dplyr::select(GID,solution,std.error) %>%
    filter(grepl("GID",GID)) %>%
    rename(BLUP=solution) %>%
    mutate(GID=gsub("GID_","",GID),
           PEV=std.error^2, # asreml specific
           REL=1-(PEV/Vg), # Reliability
           drgBLUP=BLUP/REL, # deregressed BLUP
           WT=(1-H2)/((0.1 + (1-REL)/REL)*H2)) # weight for use in Stage 2
  out<-tibble(loglik=ll,Vg,Ve,H2,
              blups=list(blups),
              varcomp=list(varcomp),
              outliers1=list(outliers1),
              outliers2=list(outliers2))
  gc()
  return(out) }
