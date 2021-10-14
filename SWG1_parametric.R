two_parts_model=function(chemin.in,seas,prevu,n,distrubtion,nyear,andebut,ht_seuil){
  #library
  library(dplyr)
  library(foreach)
  library(readr)
  library(doParallel)
  library(parallel)
  library(markovchain)
  library(Rsolnp)
  library(fitdistrplus)
  library(vars)
  library(lubridate)
  library(matlib)
  library(reshape2)
  library(xts)
  
  
  nb_cores <- detectCores() - 1
  cl <- makeCluster(nb_cores)
  registerDoParallel(cl)  
  
  #read station files
  kandi_val<- read_csv(chemin.in, col_types = cols(year = col_character()))
  Clim_real=kandi_val

####################Briggs and Wilks Methodology################## 
  #Select season
  saison_tampon=subset(Clim_real,mois%in%seas:(seas+2))
  
  #Compute seasonal cumulative rainfall
  saison_tampon=saison_tampon %>%
    transmute(year=year,PRCP=PRCP,season=rep(1,dim(saison_tampon)[1]))
  saison_tampon=saison_tampon%>%group_by(year,season)%>%summarise(PRCP=sum(PRCP))
  #cleaning df
season=saison_tampon%>%dplyr::select(year=year,rain=PRCP)
moy=mean(season$rain,na.rm = T)
remove(saison_tampon)

#compute terciles
  tercile=quantile(season$rain,probs = seq(0,1,1/3))

# Find years in each categorial

  belowyear=as.vector(subset(season, 
                             rain<tercile[2],select = year)$year)
  aboveyear=as.vector(subset(season, 
                             rain>tercile[3],select = year)$year)  
  normalyear=as.vector(subset(season, 
                      rain>tercile[2]&rain<tercile[3],select = year)$year) 

  
  
  
#Bootsrapping
  sample_below=sample(belowyear,n*prevu[1],replace=T)
  sample_normal=sample(normalyear,n*prevu[2],replace=T)
  sample_above=sample(aboveyear,n*prevu[3],replace = T)
  
#Consitent climatology
  yearclimcond=c(sample_below,sample_normal,sample_above)
 

   ##############Modele VAR pour les autres paramètres#######
   
  #Mise en forme zscore et série temporelle
  no_rain=Clim_real[,-c(4,5)]
  
  #log(max(rad)-rad) et sqrt(vent)
  maxi_insol=max(no_rain$RAD)
   
  rain_rain=Clim_real%>%transmute(year=year,month=month,day=day,PRCP=PRCP,TMAX=TMAX, TMIN=TMIN,UMAX=UMAX,UMIN=UMIN, RAD=log(maxi_insol-RAD),WIND=sqrt(WIND))
  

  no_rain=no_rain%>%transmute(year=year,month=month,day=day,TMAX=TMAX, TMIN=TMIN,UMAX=UMAX,UMIN=UMIN, RAD=log(maxi_insol-RAD),WIND=sqrt(WIND))
  
  no_rain_wide=reshape2::melt(data=no_rain, measure.vars=names(no_rain)[-c(1,2,3)], id.vars=c("year","month","day"))
  
  no_rain_wide <- dcast(no_rain_wide, month + day+variable ~ year)
  t=no_rain_wide[,-c(1,2,3)]
  zscore = foreach(i =1:dim(t)[1], .combine=rbind)%dopar%t(scale(as.numeric(t[i,])))
  colnames(zscore)=names(t)
  zscore=cbind(no_rain_wide[,1:3],zscore)
  zscore=reshape2::melt(data=zscore, measure.vars=names(zscore)[-c(1,2,3)], id.vars=c("month","day","variable"))
  colnames(zscore)=c("month","day","elmnt","year","value")
  
  ###tri
  zscore=zscore[order( zscore[,3] ),]
  zscore<- dcast(zscore, year + month + day~ elmnt)
  zscore=zscore[complete.cases(zscore),]
  zscore=xts(zscore[,-c(1,2,3)],as.Date(paste(zscore$year,zscore$month,zscore$day, sep = "-")))
  

  

  
  
 # fixer les longueurs des mois
 if(seas%in%c(3,5,8)){lmois=c(31,30,31)}
 if(seas%in%c(4,9)){lmois=c(30,31,30)}
 if(seas==6){lmois=c(30,31,31)}
 if(seas==7){lmois=c(31,31,30)}  
  
  
  #calcul des Paramètres du générateurs
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
    #modèle var par mois pour les autres variables autres que PRCP
model.var=VAR(zscore[month(zscore)==as.numeric(seas+jjseas),],1, type = "none")
coeftampon=Acoef(model.var)[[1]]
#matrice A
coefA[[jjseas+1]]=coeftampon

acftampon=acf(zscore[month(zscore)==as.numeric(seas+jjseas),],type = "correlation", plot = F,lag.max = 1)

#Matrices Mo,M1 et des résidus(écart type)
coef0[[jjseas+1]]=acftampon$acf[1,,]
coef1[[jjseas+1]]=acftampon$acf[2,,]
coefS[[jjseas+1]]=coef0[[jjseas+1]]-coefA[[jjseas+1]]%*%t(coef1[[jjseas+1]])

#moyennes et écart types des jours pluvieux/non pluvieux

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

 #moyenne pour facteur
moy_nonprcp[[jjseas+1]]=colMeans(rain_rain[rain_rain$year%in%yearclimcond&as.numeric(rain_rain$month)==(jjseas+seas),-c(1,2,3,4)], na.rm = T)

 ### initialisation mois et année climato cohérente
  A=matrix(nrow = 2, ncol = 2, 0)
  cum_real=0
  
  for(i in yearclimcond){
  #matrices de transition
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
  
  
  #quantité de pluie
  month1_amount=NULL
  for(i in yearclimcond){
    monthh=filter(Clim_real,year==i,mois==(seas+jjseas))
    month1_amount=c(month1_amount,monthh$PRCP)
  }
  
  
  #month1_amount= foreach(i =yearclimcond, .combine=c) %dopar%
    #(subset(Clim_real,year=i,mois==(seas+jjseas),select =c(PRCP))$PRCP)
  
  month1_amount=month1_amount[which(month1_amount>=ht_seuil)]
  
  #Ajustement distribution des pluies journalières
  if(distrubtion==1){
    amount_pargm <- fitdist(month1_amount, "exp")
    param=amount_pargm$estimate #paramètres du modèle
  }else{ 
    amount_pargm <- fitdist(month1_amount, "gamma")
    param=amount_pargm$estimate #paramètres du modèle
  }
  parametr[[jjseas+1]]=param
  
  jjseas=jjseas+1
  if (jjseas==3){
    break()
  }
  }
  
para_gene_pluie=list(matrices=trans_mat,para_qant=parametr,cum_real.clim=cum_moy_seas_real)
  
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
#initialisation autoragression
Z0=mvrnorm(n = 1, rep(0,6), coef0[[jseas+1]])

#proba inconditionelle des jours pluvieux
pi=trans_mat[[jseas+1]][2,1]/(1+trans_mat[[jseas+1]][2,1]-trans_mat[[jseas+1]][1,1])

 month11=NULL
 u=runif(1,0,1)

 if(u<=pi){
nonprcpSim=(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1])
  }else{
nonprcpSim=(Z0*u_sigma0[[jseas+1]][,2]+u_sigma0[[jseas+1]][,1])
  }
 
   
 if(distrubtion==1){
 month11=ifelse(u<=pi,qexp(runif(1,0,1),parametr[[jseas+1]]),0)
  }else{
 month11=ifelse(u<=pi,qgamma(runif(1,0,1),parametr[[jseas+1]][1],parametr[[jseas+1]][2]),0)
  }
 
 ms=(lmois[jseas+1]-1) #mois pour la sim
 
for(i in 1:ms){
  Z0=coefA[[jseas+1]]%*%Z0+mvrnorm(n = 1, rep(0,6), coefS[[jseas+1]])
  
  u=runif(1,0,1)
  if(month11[i]<=ht_seuil&u<=trans_mat[[jseas+1]][2,1]){# voir hauteur limites ici 
    
  #non précipitation variables
 nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1]))
    
  #precipitation variables
 
    #choix de la loi
    if(distrubtion==1){
month11=c(month11,qexp(runif(1,0,1),parametr[[jseas+1]]))
    }else{
month11=c(month11,qgamma(runif(1,0,1),parametr[[jseas+1]][1],parametr[[jseas+1]][2]))
    }
    #fin choix
    
  }else{
    if(month11[i]>ht_seuil&u<=trans_mat[[jseas+1]][1,1]){
      #non precip variables
nonprcpSim=cbind(nonprcpSim,(Z0*u_sigma1[[jseas+1]][,2]+u_sigma1[[jseas+1]][,1]))

     #precipitation variables
      #choix de la loi

      if(distrubtion==1){
        month11=c(month11,qexp(runif(1,0,1),parametr[[jseas+1]]))
      }else{
        month11=c(month11,qgamma(runif(1,0,1),parametr[[jseas+1]][1],parametr[[jseas+1]][2]))
      }    
      #fin choix
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
sim_mois_jour=cbind(rep(startyear,length(prcpSimSeason)),c(mm1,mm2,mm3),c(j1,j2,j3),round(prcpSimSeason,1))
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
cum_moy_seas_sim=mean(cum_sim_year$rain,na.rm = T)

f_prep=cum_moy_seas_real/cum_moy_seas_sim
prcpSimSeasonyear=as.data.frame(prcpSimSeasonyear)
rainadj=round(prcpSimSeasonyear$rain*f_prep,1)

prcpSimSeasonyear=as.data.frame(prcpSimSeasonyear) %>%
  transmute(year=year,month=month, day=day, rain=rain,rain_adj=round(rainadj,1))

cum_sim_year_normal_adj=as.data.frame(prcpSimSeasonyear)%>%group_by(year)%>%summarise(rain=sum(rain),rain_adj=sum(rain_adj))

nonprcpSimSeasonyear=round(nonprcpSimSeasonyear,1)
data_sim=as.data.frame(cbind(prcpSimSeasonyear,nonprcpSimSeasonyear))
return(list(data_sim=data_sim,cm_seas_sim=cum_sim_year_normal_adj,parametre=para_gene_pluie))
stopCluster(cl)
}
