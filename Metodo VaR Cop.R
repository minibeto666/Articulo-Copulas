#-------------------------------------------------------------------------------
#              FUNCION PARA CALCULO DE VaR CON COPULA (2 ACTIVOS)
#-------------------------------------------------------------------------------

# CARGA LOS PAQUETES
  #loads packages
  
  library(quantmod)
  library(copula)
  library(VineCopula)
  library(rriskDistributions)  
  library(fitdistrplus)
  library(MASS)
  library(actuar) 
  library(univariateML)
  library(VGAM)
  library(lcopula)
  library(VC2copula)
  library(scatterplot3d)
  library(cvar)
  library(knitr)
  library(DescTools)
  library(fGarch)
  
  #carga el portafolio de 2 emisoras por rendimientos
  
  S1<-"BTC-USD"  #como aparecen en el url de yahoo finance
  S2<-"ETH-USD"
  
  portfolio_ticks<-c(S1,S2)
  
  #asigna fecha a partir de la cual se tomaran los precios 
  
  date<-'2021-01-01'  # yyyy-mm-dd
  
  #carga datos de open,high,low y cierre
  
  getSymbols(S1,src = "yahoo",from=date)
  n1<-nrow(get(S1)[,1])
  
  getSymbols(S2,src = "yahoo",from=date)
  n2<-nrow(get(S2)[,1])
  
  #carga los rendimientos diarios (ojo el primer dia será un NA si leading=False, si es True lo toma en cuenta)
  
  r1<- dailyReturn(get(S1),leading = TRUE) #sin get toma el - como operador
  r2<- dailyReturn(get(S2),leading = TRUE)
  
  #crea la base de datos de c/u con los rendimientos
  
  #data1<-as.data.frame(cbind(DB1,r1)) #rendimientos en la columna 
  #data2<-as.data.frame(cbind(DB2,r2))
  
  #cargamos vector de ultimos precios
  
  last_price<- as.vector(c( get(S1)[n1,6],get(S2)[n2,6] ))  #precios de cierre ajustados en columna 6 (ojo que puede ser NA)
  
    if (is.na(last_price[1]) | is.na(last_price[2])){
      
      last_price<- as.vector(c( get(S1)[n1-1,6],get(S2)[n2-1,6] ))  #toma los del cierre del dia anterior
      
    }  
  
  #ajustamos una distribución a los rendimientos usando criterio de AIC (investigar que metodo es mejor)
  
  bestr1<-model_select(r1,criterion="loglik") 
  
  bestr2<-model_select(r2,criterion="loglik")
  
  # OJO
  # una vez que se determina la mejor distribución, se va a usar su CDF, algunas necesitan de paquetes para usarla (No se puede automatizar :( )

  #este caso es una Skew Generalized Error
  
  # Calculamos los valores de u con la distribucion ajustada
  
  u<-psged(as.vector(r1), mean =  bestr1[1], sd = bestr1[2], nu = bestr1[3], xi = bestr1[4])

  #replicamos para v con r2
  
  v<-psged(as.vector(r2), mean =  bestr2[1], sd = bestr2[2], nu = bestr2[3], xi = bestr2[4])

  #grafico y prueba de correlacion
  
  plot(u,v)  
  
  cor.test(u,v, method="spearman")  
  
  #Ajustamos la mejor copula por AIC, BIC o loglik
  
  CopAdj<-BiCopSelect(u,v,familyset=NA,selectioncrit = "logLik") #cop rotada (el resultado es una lista)
  
  copula<-surBB1Copula(param = c(CopAdj[[2]],CopAdj[[3]]) )  #seleccionamos la copula correspondiente con los parametros obtenidos 
  
  
  #Creamos la distribución conjunta
  
  #en los parametros hay que ponerlosen lista dependiendo del nombre de cada parametro y el numero de parametros
  DistConj<-mvdc(copula,margins=c("sged","sged"),
                 paramMargins=list(
                   list(mean=coef(bestr1)[1],sd=coef(bestr1)[2],nu=coef(bestr1)[3],xi=coef(bestr1)[4]),
                   list(mean=coef(bestr2)[1],sd=coef(bestr2)[2],nu=coef(bestr2)[3],xi=coef(bestr2)[4])
                 )
  )
  
  # Simulamos puntos aleatorios de la copula 
  
  # Simulamos 1000 valores de X y Y 
  z <- rMvdc(n=1000,DistConj)
  
  # Creamos el grafico de la pdf
  pdf<-dMvdc(z,DistConj) 
  scatterplot3d(z[,1],z[,2],pdf,highlight.3d=T,main = "Densidad de la copula",xlab = "rendimientos 1",
                ylab = "rendimientos 2")
  
  # Creamos el grafico de la cdf
  cdf<-pMvdc(z,DistConj) 
  scatterplot3d(z[,1],z[,2],cdf,highlight.3d=T,main = "Distribucion de la copula",xlab = "rendimientos 1",
                ylab = "rendimientos 2")
  
  # gofCopula para pruebas de bondad de ajuste (no jala con todas las cop)
  
  valores=cbind(as.vector(r1),as.vector(r2))
  
  gofCopula(copula = copula,x=valores) #?
  
  
  #-----------------------------------------------------------------------------
  # Procedemos al calculo del VaR
  
  VaR_df<-data.frame()
  tVaR_df<-data.frame()
  
  for(k in 1:100){ #numero de veces que se calcularan los VaR y tVaR para sacarles promedio
    
    nsims<-1000 #numero de simulaciones sobre las cuales hace P&L, total de simulaciones 100x nsims
    
    sim_rend<-data.frame()
    
    sim_rend<- rMvdc(n=nsims,DistConj) #x y y sim con la copula   simula nsims puntos aleatorios de la distribución conjunta (rendimientos simulados)
    
    #calculamos reevaluacion
    reev_sim<-data.frame()
    
    for(j in 1:2){
      
      for( i in 1:nrow(sim_rend)){
        reev_sim[i,j]<-last_price[j]*(1+sim_rend[i,j])
      }#i
      
    }#j
    
    #calculamos P&L
    
    PL_sim<-data.frame()
    
    for (j in 1:2){
      
      for(i in 1:nrow(reev_sim)){
        
        PL_sim[i,j]<-last_price[j]-reev_sim[i,j]
      }#i
    }#j
    
    #añadimos PL del portafolio
    
    PL_port<-apply(PL_sim,1,sum) 
    
    PL_sim<-cbind(PL_sim,PL_port)
    
    #Calculamos VaR 95 y tVaR 95
    
    VaR95port<-apply(PL_sim,2,quantile,probs=.95)  #Aca le cambiamos para el 97,99,99.9 etc
    tVaR95port<-apply(PL_sim,2,ES,p_loss=.05)      #DE igual forma le tenemos que cambiar
    
    VaR_df[k,1]<-VaR95port[1]
    VaR_df[k,2]<-VaR95port[2]
    VaR_df[k,3]<-VaR95port[3]
    
    tVaR_df[k,1]<-tVaR95port[1]
    tVaR_df[k,2]<-tVaR95port[2]
    tVaR_df[k,3]<-tVaR95port[3]
    
  }#k
    
  # VaR y tVaR
  VaR_95<-colMeans(VaR_df)
  tVaR_95<-colMeans(tVaR_df) 
  
  copVaR_tabla<- kable( t(VaR_95), digits = 4, col.names = c("Activo1","Activo2","Portafolio"), caption = "Copula VaR 95%" )
  
  coptVaR_tabla<-kable( t(tVaR_95), digits = 4, col.names = c("Activo1","Activo2","Portafolio"), caption = "Copula tVaR 95% (Expected Shortfall)" )
  
  copVaR_tabla
  
  coptVaR_tabla
  
  
  
  
  
  
  
  
  
  
  
  
  