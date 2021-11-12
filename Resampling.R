Briggs_wilks_resampling=function(chemin.in,seas,prevu,n){

  #read station files
  kandi_val<- read_csv(chemin.in, col_types = cols(year = col_character()))
  Clim_real=kandi_val
  
  ####################Briggs and Wilks Methodology################## 
  #Select season
  saison_tampon=subset(Clim_real,mois%in%seas:(seas+2))
  saison_tampon=saison_tampon%>%filter(between(year,1981,2010))
  
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
                              rain>=tercile[2]&rain<=tercile[3],select = year)$year) 
  
  
  
  
  #Bootsrapping
  sample_below=sample(belowyear,n*prevu[1],replace=T)
  sample_normal=sample(normalyear,n*prevu[2],replace=T)
  sample_above=sample(aboveyear,n*prevu[3],replace = T)
  
  #Consitent climatology
  yearclimcond=c(sample_below,sample_normal,sample_above)
  return(yearclimcond)
}
