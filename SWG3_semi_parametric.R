SWG4=function(chemin.in,yearclimcond,seas,nyear,andebut,ht_seuil){
  
  #Read station
  kandi_val<- read_csv(chemin.in, col_types = cols(year = col_character()))
  Clim_real=kandi_val
  
  
##########################Model VAR for others variables#####################
  
  #Zscore and temporal serie
  
  no_rain=Clim_real[,-c(4,5)]
  
  #log(max(rad)-rad) et sqrt(vent)
  maxi_insol=max(no_rain$RAD,na.rm = T)
  
  rain_rain=Clim_real%>%transmute(year=year,month=month,day=day,PRCP=PRCP,TMAX=TMAX, TMIN=TMIN,UMAX=UMAX,UMIN=UMIN, RAD=log(maxi_insol-RAD),WIND=sqrt(WIND))
  
  
  no_rain=no_rain%>%transmute(year=year,month=month,day=day,TMAX=TMAX, TMIN=TMIN,UMAX=UMAX,UMIN=UMIN, RAD=log(maxi_insol-RAD),WIND=sqrt(WIND))
  
  no_rain_wide=reshape2::melt(data=no_rain, measure.vars=names(no_rain)[-c(1,2,3)], id.vars=c("year","month","day"))
  
  no_rain_wide <- dcast(no_rain_wide, month + day+variable ~ year)
  t=no_rain_wide[,-c(1,2,3)]
  
  nb_cores <- detectCores() - 2
  cl <- makeCluster(nb_cores)
  registerDoParallel(cl)
  zscore = foreach(i =1:dim(t)[1], .combine=rbind)%dopar%t(scale(as.numeric(t[i,])))
  stopCluster(cl)
  
  colnames(zscore)=names(t)
  zscore=cbind(no_rain_wide[,1:3],zscore)
  zscore=reshape2::melt(data=zscore, measure.vars=names(zscore)[-c(1,2,3)], id.vars=c("month","day","variable"))
  colnames(zscore)=c("month","day","elmnt","year","value")
  
  ###Tri
  zscore=zscore[order( zscore[,3] ),]
  zscore<- dcast(zscore, year + month + day~ elmnt)
  zscore=zscore[complete.cases(zscore),]
  zscore=xts(zscore[,-c(1,2,3)],as.Date(paste(zscore$year,zscore$month,zscore$day, sep = "-")))
  
  
  
  # Set length of Month
  if(seas%in%c(3,5,8)){lmois=c(31,30,31)}
  if(seas%in%c(4,9)){lmois=c(30,31,30)}
  if(seas==6){lmois=c(30,31,31)}
  if(seas==7){lmois=c(31,31,30)}  
  
  
  #Compute SWGs parameters
  jjseas=0
  cum_moy_seas_real=0
  trans_mat=list(0)
  parametr1=list(0)
  parametr2=list(0)
  quantmois=list(0)
  coefA=list(0)
  coef0=list(0)
  coef1=list(0)
  coefS=list(0)
  u_sigma0=list(0)
  u_sigma1=list(0)
  moy_nonprcp=list(0)
  
  repeat{
    #modèle var par mois pour les autres variables autres que PRCP
    #Model VAR for month for others varaiables (no rain)
    model.var=VAR(zscore[month(zscore)==as.numeric(seas+jjseas),],1, type = "none")
    coeftampon=Acoef(model.var)[[1]]
    #matrix A
    coefA[[jjseas+1]]=coeftampon
    
    acftampon=acf(zscore[month(zscore)==as.numeric(seas+jjseas),],type = "correlation", plot = F,lag.max = 1)
    
    #Matrix Mo,M1 and residus (écart type)
    coef0[[jjseas+1]]=acftampon$acf[1,,]
    coef1[[jjseas+1]]=acftampon$acf[2,,]
    coefS[[jjseas+1]]=coef0[[jjseas+1]]-coefA[[jjseas+1]]%*%t(coef1[[jjseas+1]])
    
    #Mean and SD for rain days/no rain days
    
    u_tampon1=colMeans(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas)&rain_rain$PRCP>=ht_seuil,-c(1,2,3,4)], na.rm = T)
    
    sigma_tampon1=apply(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas)&rain_rain$PRCP>=ht_seuil,-c(1,2,3,4)],2,sd, na.rm = T)
    
    tampon_u_sig1=as.data.frame(cbind(u_tampon1,sigma_tampon1))
    colnames(tampon_u_sig1)=c("moy","sd")
    
    u_tampon0=colMeans(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas)&rain_rain$PRCP<ht_seuil,-c(1,2,3,4)], na.rm = T)
    
    sigma_tampon0=apply(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas)&rain_rain$PRCP<ht_seuil,-c(1,2,3,4)],2,sd, na.rm = T)
    
    tampon_u_sig0=as.data.frame(cbind(u_tampon0,sigma_tampon0))
    colnames(tampon_u_sig0)=c("moy","sd")
    
    u_sigma0[[jjseas+1]]= tampon_u_sig0
    
    u_sigma1[[jjseas+1]]= tampon_u_sig1
    
    #Mean for factor
    moy_nonprcp[[jjseas+1]]=colMeans(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas),-c(1,2,3,4)], na.rm = T)
    
    ### Initialization of month and year consistent with climatology
    A=matrix(nrow = 3, ncol = 3, 0)
    cum_real=0
    
    for(i in yearclimcond){
      #matrix of transition
      month1=Clim_real[Clim_real$year==i&Clim_real$mois==(seas+jjseas),]$PRCP
      cum_real=cum_real+sum(month1, na.rm = T)
      
      ####states
      q80=quantile(month1[which(month1>=ht_seuil)], probs = 0.8)
      isec=which(month1<ht_seuil)
      iH1=which(month1>=ht_seuil&month1<=q80)
      #iH2=which(month1>5&month1<=10)
      iH2=which(month1>q80)
      
      month1[isec]<-"sec"
      month1[iH1]<-"H1"
      month1[iH2]<-"H2"
      #month1[iH3]<-"H3"
      
      
      
      month1.fit<- markovchainFit(data = month1,method = "mle")
      
      A=A+month1.fit$estimate@transitionMatrix ##matrices de transition
    }
    A=A/length(yearclimcond)
    trans_mat[[jjseas+1]]=A
    cum_real=cum_real/length(yearclimcond)
    cum_moy_seas_real=cum_moy_seas_real+cum_real
    
    
    #rainfall amount
    month1_amount=NULL
    for(i in yearclimcond){
      monthh=filter(Clim_real,year==i,mois==(seas+jjseas))
      month1_amount=c(month1_amount,monthh$PRCP)
    }
    
    
    #   month1_amount= foreach(i =yearclimcond, .combine=c) %dopar%
    #     (subset(Clim_real,year=i,mois==(seas+jjseas),select =c(PRCP))$PRCP)
    
    q80=quantile(month1_amount[which(month1_amount>=ht_seuil)], probs = 0.8)
    
    month1_amount1=month1_amount[which(month1_amount>=ht_seuil&month1_amount<=q80)] 
    
    month1_amount3=month1_amount[which(month1_amount>q80)]
    
    
    
    #Ajustement distribution des pluies journalières
    #Fitting daily rain
    amount_pargm1<-kde(month1_amount1)
    #param1=amount_pargm1$estimate
    
    amount_pargm2<-kde(unique(month1_amount3))
    #param2=amount_pargm2$estimate
    
    parametr1[[jjseas+1]]=amount_pargm1
    parametr2[[jjseas+1]]=amount_pargm2
    quantmois[[jjseas+1]]=q80
    
    jjseas=jjseas+1
    if (jjseas==3){
      break()
    }
  }
  
  para_gene_pluie=list(matrices=trans_mat,para_qant=list(parametr1,parametr2),cum_real.clim=cum_moy_seas_real,quantmois)

  #Simulation   
  #initialisation toutes les années
  nbyear=0
  startyear=andebut
  prcpSimSeasonyear=NULL
  nonprcpSimSeasonyear=NULL
  
  repeat{ 
    #initialisation annuel
    prcpSimSeason=NULL
    nonprcpSimseason=NULL
    jseas=0
    repeat{
      #initialisation autoregression
      Z0=mvrnorm(n = 1, rep(0,6), coef0[[jseas+1]])
      
      pi1=(trans_mat[[jseas+1]][3,1]*(1+trans_mat[[jseas+1]][3,2]-trans_mat[[jseas+1]][2,2])-trans_mat[[jseas+1]][3,2]*(trans_mat[[jseas+1]][3,1]-trans_mat[[jseas+1]][2,1]))/((1+trans_mat[[jseas+1]][3,1]-trans_mat[[jseas+1]][1,1])*(1+trans_mat[[jseas+1]][3,2]-trans_mat[[jseas+1]][2,2])-(trans_mat[[jseas+1]][3,1]-trans_mat[[jseas+1]][2,1])*(trans_mat[[jseas+1]][3,2]-trans_mat[[jseas+1]][1,2]))
      
      pi2=(trans_mat[[jseas+1]][3,2]-pi1*(trans_mat[[jseas+1]][3,2]-trans_mat[[jseas+1]][1,2]))/(1+trans_mat[[jseas+1]][3,2]-trans_mat[[jseas+1]][2,2])
      
      pi0=1-pi1-pi2
      
      #calcul day 1
      u=runif(1,0,1)
      if(u<pi0){
        month11=0
        nonprcpSim=(Z0*u_sigma0[[jseas+1]][,2]+u_sigma0[[jseas+1]][,1])
      }else{
        
        nonprcpSim=(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1])
        
        #index1=which(trans_mat[[jseas+1]]==max(trans_mat[[jseas+1]][,-dim(trans_mat[[jseas+1]])[2]]),arr.ind = T)
        
        if(u>=pi0&u<(pi0+pi1)){
          month11=qkde(runif(1,0,1),parametr1[[jseas+1]])
        }else{
          month11=qkde(runif(1,0,1),parametr2[[jseas+1]])
          
        }
        
      }
      
      
      ms=(lmois[jseas+1]-1) #mois pour la sim
      
      for(i in 1:ms){
        Z0=coefA[[jseas+1]]%*%Z0+mvrnorm(n = 1, rep(0,6), coefS[[jseas+1]])
        
        if(month11[i]<=ht_seuil){tampon="sec"}
        if((month11[i]>ht_seuil&month11[i]<=quantmois[[jseas+1]])){tampon="H1"}
        if(month11[i]>quantmois[[jseas+1]]){tampon="H2"}
        
        u=runif(1,0,1) #generation de nombres aléatoires
        
        if(u<trans_mat[[jseas+1]][tampon,3]){
          
          month11=c(month11,0)
          #non precipitation variables 
          nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma0[[jseas+1]][,2]+u_sigma0[[jseas+1]][,1]))
          
        }else{###############
          #non precipitation var
          nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1]))
          
          # index1=which(u<=(trans_mat[[jseas+1]][tampon,-dim(trans_mat[[jseas+1]])[2]]),arr.ind = T)
          # min.index=min(trans_mat[[jseas+1]][tampon,index1])
          # index.bon=which(trans_mat[[jseas+1]][tampon,]==min.index,arr.ind = T)
          
          if(u>=trans_mat[[jseas+1]][tampon,3]&u<(trans_mat[[jseas+1]][tampon,3]+trans_mat[[jseas+1]][tampon,1])){
            month11=c(month11,month11=qkde(runif(1,0,1),parametr1[[jseas+1]]))
          }else{
            month11=c(month11,qkde(runif(1,0,1),parametr2[[jseas+1]]))
            
          }
          
        }
        #####################
      }
      
      nonprcpSim=as.data.frame(t(nonprcpSim))
      moy_sim_nonprcp=colMeans(nonprcpSim,na.rm = T)
      fact_nonprcp=moy_nonprcp[[jseas+1]]/moy_sim_nonprcp
      nonprcpSim=nonprcpSim%>%transmute(TMAX=TMAX*fact_nonprcp[1],TMIN=TMIN*fact_nonprcp[2],UMAX=ifelse(UMAX*fact_nonprcp[3]>100,100,UMAX*fact_nonprcp[3]),UMIN=abs(UMIN*fact_nonprcp[4]),RAD=ifelse((maxi_insol-exp(RAD*fact_nonprcp[5]))<0,0.1,maxi_insol-exp(RAD*fact_nonprcp[5])), WIND=(WIND*fact_nonprcp[6])**2)
      
      
      
      nonprcpSimseason=rbind(nonprcpSimseason,nonprcpSim)
      prcpSimSeason= c(prcpSimSeason,month11)
      jseas=jseas+1
      if (jseas==3){
        break()
      }
    }
    
    mm1=rep(seas,lmois[1])
    j1=1:lmois[1]
    mm2=rep(seas+1,lmois[2])
    j2=1:lmois[2]
    mm3=rep(seas+2,lmois[3])
    j3=1:lmois[3]
    sim_mois_jour=cbind(rep(startyear,length(prcpSimSeason)),c(mm1,mm2,mm3),c(j1,j2,j3),abs(round(prcpSimSeason,1)))
    remove(mm1,mm2,mm3,j1,j2,j3)
    
    nonprcpSimSeasonyear=rbind(nonprcpSimSeasonyear,nonprcpSimseason)
    prcpSimSeasonyear=rbind(prcpSimSeasonyear,sim_mois_jour)
    nbyear=nbyear+1
    startyear=startyear+1
    if (nbyear==nyear){
      break()
    }
    
  }
  
  colnames(prcpSimSeasonyear)=c("year","month","day","rain")
  cum_ss=as.data.frame(prcpSimSeasonyear)%>%mutate( season=rep(rep(1,dim(prcpSimSeasonyear)[1])))
  cum_sim_year=as.data.frame(cum_ss)%>%group_by(year, season)%>%summarise(rain=sum(rain))
  
 # Boughton(1999) Correction of Variance
  Clim_cond1=NULL
  sdtampon=NULL
  for(i in yearclimcond){
    Clim_cond11= filter(Clim_real,year%in%i,mois%in%c(seas,seas+1,seas+2))
    ggg=as.data.frame(Clim_cond11)%>%group_by(year)%>%summarise(PRCP=sum(PRCP))
    Clim_cond1=rbind(Clim_cond1,Clim_cond11)
    sdtampon=rbind(sdtampon,ggg)
  }
  
  sd0_adj=sd(sdtampon$PRCP,na.rm = T)
  
  
  sd1_adj=sd(cum_sim_year$rain,na.rm = T)
  f_boughton=sd0_adj/sd1_adj
  prcpSimSeasonyearCor=NULL
  anfin=andebut+nyear
  
  
  for(i in (andebut:anfin)){
    
    ratio.tampon=(cum_moy_seas_real+(cum_sim_year[cum_sim_year$year==i,]$rain-cum_moy_seas_real)*f_boughton)/cum_sim_year[cum_sim_year$year==i,]$rain
    transfoorm=as.data.frame(prcpSimSeasonyear)%>%filter(year==i)%>%mutate(rain_adj=round(rain*ratio.tampon,1))
    prcpSimSeasonyearCor=rbind(prcpSimSeasonyearCor,transfoorm)
    
  }

  cum_sim_year_norm_adj=as.data.frame(prcpSimSeasonyearCor)%>%group_by(year)%>%summarise(rain=sum(rain),rain_adj=sum(rain_adj))
  nonprcpSimSeasonyear=round(nonprcpSimSeasonyear,1)
  data_sim=as.data.frame(cbind(prcpSimSeasonyearCor,nonprcpSimSeasonyear))
  return(list(data_sim=data_sim,cm_seas_sim=cum_sim_year_norm_adj,parametre=para_gene_pluie))

}

