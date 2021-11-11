SWG2=function(chemin.in,yearclimcond,seas,nyear,andebut,ht_seuil){
  
  #read station files
  kandi_val<- read_csv(chemin.in, col_types = cols(year = col_character()))
  Clim_real=kandi_val
  
  #####VAR Model for other variables (TMAX,TMIN...)
  
  #Standardization 
  no_rain=Clim_real[,-c(4,5)]
  
  #Transform radiation and wind
  #log(max(rad)-rad) et sqrt(wind)
  maxi_insol=max(no_rain$RAD,na.rm = T)
  
  rain_rain=Clim_real%>%transmute(year=year,month=month,day=day,PRCP=PRCP,TMAX=TMAX, TMIN=TMIN,UMAX=UMAX,UMIN=UMIN, RAD=log(maxi_insol-RAD),WIND=sqrt(WIND))
  
  no_rain=no_rain%>%transmute(year=year,month=month,day=day,TMAX=TMAX, TMIN=TMIN,UMAX=UMAX,UMIN=UMIN, RAD=log(maxi_insol-RAD),WIND=sqrt(WIND))
  
  no_rain_wide=reshape2::melt(data=no_rain, measure.vars=names(no_rain)[-c(1,2,3)], id.vars=c("year","month","day"))
  
  no_rain_wide <- dcast(no_rain_wide, month + day+variable ~ year)
  
  t=no_rain_wide[,-c(1,2,3)]
  
  zscore = foreach(i =1:dim(t)[1], .combine=rbind)%do%t(scale(as.numeric(t[i,])))
  
  colnames(zscore)=names(t)
  
  zscore=cbind(no_rain_wide[,1:3],zscore)
  
  zscore=reshape2::melt(data=zscore, measure.vars=names(zscore)[-c(1,2,3)], id.vars=c("month","day","variable"))
  colnames(zscore)=c("month","day","elmnt","year","value")
  
  ###order and put in temporal series
  zscore=zscore[order( zscore[,3] ),]
  zscore<- dcast(zscore, year + month + day~ elmnt)
  zscore=zscore[complete.cases(zscore),]
  zscore=xts(zscore[,-c(1,2,3)],as.Date(paste(zscore$year,zscore$month,zscore$day, sep = "-")))
  
  
  # Set month length
  if(seas%in%c(3,5,8)){lmois=c(31,30,31)}
  if(seas%in%c(4,9)){lmois=c(30,31,30)}
  if(seas==6){lmois=c(30,31,31)}
  if(seas==7){lmois=c(31,31,30)}  
  
  
  #Computing SWGs Parameters for rain and other variables (TMAX,TMIN...)
  jjseas=0
  cum_moy_seas_real=0
  trans_mat=list(0)
  parametr=list(0)
  coefA=list(0)
  coef0=list(0)
  coef1=list(0)
  coefS=list(0)
  u_sigma0=list(0)
  u_sigma1=list(0)
  moy_nonprcp=list(0)
  repeat{
    #VAR model by month for non_rain variable 
    model.var=VAR(zscore[month(zscore)==as.numeric(seas+jjseas),],1, type = "none")
    coeftampon=Acoef(model.var)[[1]]
    #matrix A
    coefA[[jjseas+1]]=coeftampon
    
    acftampon=acf(zscore[month(zscore)==as.numeric(seas+jjseas),],type = "correlation", plot = F,lag.max = 1)
    
    #Matrix Mo,M1 et des résidus
    coef0[[jjseas+1]]=acftampon$acf[1,,]
    coef1[[jjseas+1]]=acftampon$acf[2,,]
    coefS[[jjseas+1]]=coef0[[jjseas+1]]-coefA[[jjseas+1]]%*%t(coef1[[jjseas+1]])
    
    #mean and sd for rain/no rain days
    
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
    
    #mean for Bougthon factor possibility
    moy_nonprcp[[jjseas+1]]=colMeans(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas),-c(1,2,3,4)], na.rm = T)
    
    ### Initialization month and year consistent with seasonal forecast prob mois et année climato cohérente
    ### transition matrix for rain occurence  
    A=matrix(nrow = 2, ncol = 2, 0)
    cum_real=0
    
    for(i in yearclimcond){
      #transition matrix
      month1=Clim_real[Clim_real$year==i&Clim_real$mois==(seas+jjseas),]$PRCP
      cum_real=cum_real+sum(month1, na.rm = T)
      
      isec=which(month1<ht_seuil)
      ipluie=which(month1>=ht_seuil)
      
      month1[isec]<-"sec"
      month1[ipluie]<-"pluie"
      
      month1.fit<- markovchainFit(data = month1,method = "mle")
      
      A=A+month1.fit$estimate@transitionMatrix ##matrices de transition
    }
    A=A/length(yearclimcond)
    trans_mat[[jjseas+1]]=A
    cum_real=cum_real/length(yearclimcond)
    cum_moy_seas_real=cum_moy_seas_real+cum_real
    
    
    #rain amount
    month1_amount=NULL
    for(i in yearclimcond){
      monthh=filter(Clim_real,year==i,mois==(seas+jjseas))
      month1_amount=c(month1_amount,monthh$PRCP)
    }
    
    
    
    month1_amount=month1_amount[which(month1_amount>=ht_seuil)]
    
    #Fitting daily rain amount 

    parametr[[jjseas+1]]=VGAM::Coef(VGAM::vglm(month1_amount~1,mix2exp(nsimEIM = 200),trace=F))

    jjseas=jjseas+1
    if (jjseas==3){
      break()
    }
  }
  
  para_gene_pluie=list(matrices=trans_mat,para_qant=parametr,cum_real.clim=cum_moy_seas_real)
  
  ################Simulation###################################  
  #initialization for all years
  nbyear=0
  startyear=andebut
  prcpSimSeasonyear=NULL
  nonprcpSimSeasonyear=NULL
  
  repeat{ 
    #initialization 
    prcpSimSeason=NULL
    nonprcpSimseason=NULL
    jseas=0
    repeat{
      #initialization autoregression
      Z0=mvrnorm(n = 1, rep(0,6), coef0[[jseas+1]])
      
      #Inconditional probability of rainy days
      pi=trans_mat[[jseas+1]][2,1]/(1+trans_mat[[jseas+1]][2,1]-trans_mat[[jseas+1]][1,1])
      
      month11=NULL
      u=runif(1,0,1)
      
      if(u<=pi){
        nonprcpSim=(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1])
      }else{
        nonprcpSim=(Z0*u_sigma0[[jseas+1]][,2]+u_sigma0[[jseas+1]][,1])
      }
      
      month11=ifelse(u<=pi,Renext::qmixexp2(p = runif(1,0,1), prob1 = parametr[[jseas+1]][1],
               rate1 = parametr[[jseas+1]][2], rate2 = parametr[[jseas+1]][3]),0)

      ms=(lmois[jseas+1]-1) #mois pour la sim
      
      for(i in 1:ms){
        Z0=coefA[[jseas+1]]%*%Z0+mvrnorm(n = 1, rep(0,6), coefS[[jseas+1]])
        
        u=runif(1,0,1)
        if(month11[i]<=ht_seuil&u<=trans_mat[[jseas+1]][2,1]){# voir hauteur limites ici 
          
          #non précipitation variables
          nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1]))
          
          #precipitation variables
          
          month11=c(month11,Renext::qmixexp2(p = runif(1,0,1), prob1 = parametr[[jseas+1]][1],
                                   rate1 = parametr[[jseas+1]][2], rate2 = parametr[[jseas+1]][3]))
  
        }else{
          if(month11[i]>ht_seuil&u<=trans_mat[[jseas+1]][1,1]){
            #non precip variables
            nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1]))
            
            #precipitation variables
            month11=c(month11,Renext::qmixexp2(p = runif(1,0,1), prob1 = parametr[[jseas+1]][1],
                                               rate1 = parametr[[jseas+1]][2], rate2 = parametr[[jseas+1]][3]))
          }else{
            month11=c(month11,0)
            nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma0[[jseas+1]][,2]+u_sigma0[[jseas+1]][,1]))
          }
          
        }
      }
      
      nonprcpSim=as.data.frame(t(nonprcpSim))
      moy_sim_nonprcp=colMeans(nonprcpSim,na.rm = T)
      fact_nonprcp=moy_nonprcp[[jseas+1]]/moy_sim_nonprcp
      nonprcpSim=nonprcpSim%>%transmute(TMAX=TMAX*fact_nonprcp[1],TMIN=TMIN*fact_nonprcp[2],UMAX=ifelse(UMAX*fact_nonprcp[3]>100,100,UMAX*fact_nonprcp[3]),UMIN=abs(UMIN*fact_nonprcp[4]),RAD=ifelse((maxi_insol-exp(RAD*fact_nonprcp[5]))<0,0.1,maxi_insol-exp(RAD*fact_nonprcp[5])), WIND=(WIND*fact_nonprcp[6])**2)
      
      
      
      nonprcpSimseason=rbind.data.frame(nonprcpSimseason,nonprcpSim)
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
    sim_mois_jour=cbind(rep(startyear,length(prcpSimSeason)),c(mm1,mm2,mm3),c(j1,j2,j3),round(prcpSimSeason,1))
    remove(mm1,mm2,mm3,j1,j2,j3)
    
    nonprcpSimSeasonyear=rbind.data.frame(nonprcpSimSeasonyear,nonprcpSimseason)
    prcpSimSeasonyear=rbind.data.frame(prcpSimSeasonyear,sim_mois_jour)
    nbyear=nbyear+1
    startyear=startyear+1
    if (nbyear==nyear){
      break()
    }
    
  }
  
  colnames(prcpSimSeasonyear)=c("year","month","day","rain")
  
  cum_ss=as.data.frame(prcpSimSeasonyear)%>%mutate( season=rep(rep(1,dim(prcpSimSeasonyear)[1])))
  
  cum_sim_year=as.data.frame(cum_ss)%>%group_by(year, season)%>%summarise(rain=sum(rain))
  cum_moy_seas_sim=mean(cum_sim_year$rain,na.rm = T)
  
  f_prep=cum_moy_seas_real/cum_moy_seas_sim
  prcpSimSeasonyear=as.data.frame(prcpSimSeasonyear)
  rainadj=round(prcpSimSeasonyear$rain*f_prep,1)
  
  prcpSimSeasonyear=as.data.frame(prcpSimSeasonyear) %>%
    transmute(year=year,month=month, day=day, rain=rain,rain_adj=round(rainadj,1))
  
  cum_sim_year_normal_adj=as.data.frame(prcpSimSeasonyear)%>%group_by(year)%>%summarise(rain=sum(rain),rain_adj=sum(rain_adj))
  
  nonprcpSimSeasonyear=round(nonprcpSimSeasonyear,1)
  data_sim=as.data.frame(cbind.data.frame(prcpSimSeasonyear,nonprcpSimSeasonyear))
  return(list(data_sim=data_sim,cm_seas_sim=cum_sim_year_normal_adj,parametre=para_gene_pluie))
}
