
##################
# Description outputs
##################

# obsFD: observed metrics of functional diversity

# nulldf: full dataframe with simmulated metrics

# nullsummary: summary of the null communities. mean, SE and 95%CI


##################
# 0. Preamble
##################

if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach,doParallel, vegan, SYNCSA, tidyr,tseries,strucchange, svMisc, dplyr)
'%ni%'<- Negate('%in%')

load('data_input/data1.rda')

# Load functions
source('functions/FD_df.R')
source('functions/multidim.R')

# Make directory for output
dir.create('data_output')

rep.nulls<-500 #number of null replicates
method.nulls<- 'r00_samp'  #for other methods see ?commsim


##################
# 1. Observed metrics
##################

df<-FD_df(data1, subcat = c(6,6,6))

.fdmetrics<- function(df){weight <- pivot_wider(df[order(as.numeric(df$time)),c('ecocode', 'time', 'abundance')],
                      names_from = 'ecocode', values_from = 'abundance',id_cols = 'time' )
weight<-data.frame(weight, check.names = FALSE)
rownames(weight)<-weight$time
weight<-weight[-1]

coord  <- subset(df,select = 1:(nchar(as.character(df$ecocode[1]))+1))
coord<- unique(coord)
rownames(coord)<- coord[,1]
coord<-coord[,-1]

for (n in 1:ncol(coord)) {
  coord[,n] = as.numeric(coord[,n])
}
coord<- scale(coord,center = T, scale = T )
obsFD<-as.data.frame(multidimFD(coord, as.matrix(weight)))
obsRao<-rao.diversity(traits=coord, weight)

obsFD$FRed <-obsRao$FunRedundancy
obsFD$raoQ <-obsRao$FunRao
obsFD$Simpson<- obsRao$Simpson

return(obsFD)
}

obsFD<- .fdmetrics(df)

#save ovserved functional diversity metrics as csv
write.csv(obsFD, file='data_output/obsFD.csv')


##################
# 2. Null models
##################

nm<-nullmodel(data1,method.nulls)
null<-simulate(nm, nsim =rep.nulls)

## WARNING: it may generate errors associated with convex hull. Those are removed from the final output.

.nulls<- function(null0, rep.nulls){
nullsdf<- data.frame()

foreach (p=1:rep.nulls) %do% {
  print(paste0('null iteration: ',p, ' out of ',rep.nulls ))
  #progress(p,max.value = rep.nulls )
  m1<- null0[,,p]
  #m1<-(t(m0)/rowSums(t(m0)))

  datam<-FD_df(as.data.frame(m1), subcat = c(6,6,6))
  df<-datam
  weight <- pivot_wider(df[order(as.numeric(df$time)),c('ecocode', 'time', 'abundance')],
                        names_from = 'ecocode', values_from = 'abundance',id_cols = 'time' )
  weight<-data.frame(weight, check.names = FALSE)
  rownames(weight)<-weight$time
  weight<-weight[-1]

  coord  <- subset(df,select = 1:(nchar(as.character(df$ecocode[1]))+1))
  coord<- unique(coord)
  rownames(coord)<- coord[,1]
  coord<-coord[,-1]

  for (n in 1:ncol(coord)) {
    coord[,n] = as.numeric(coord[,n])
  }
  coord<- scale(coord,center = T, scale = T )


  mdimDF<-as.data.frame(multidimFD(coord, as.matrix(weight),verb = TRUE))
  raoDF<-rao.diversity(traits=coord, weight)

  nullFunDiv<-data.frame(null=p, com= rownames(mdimDF),
                         Simpson=raoDF$Simpson,
                         FunRao= raoDF$FunRao,
                         FRed= raoDF$FunRedundancy,
                         FRic=mdimDF$FRic,
                         FEve=mdimDF$FEve,
                         FSpe= mdimDF$FSpe,
                         FDis= mdimDF$FDis,
                         FDiv= mdimDF$FDiv,
                         FOri= mdimDF$FOri,
                         Tax.rich=mdimDF$Nb_sp)


  nullsdf<-rbind(nullsdf, nullFunDiv)

}
return(nullsdf)
}

nulldf<- .nulls(null0 = null, rep.nulls)

.nullsummary<-function(truenull){
  nullmodelsummary<- truenull%>%
    group_by(com)%>%
    mutate(com=as.numeric(com))%>%
    summarise(mSimpson= mean(Simpson),
              mFunRao= mean(FunRao),
              mFRed = mean(FRed),
              mFRic = mean(FRic),
              mFEve= mean(FEve),
              mFSpe= mean(FSpe),
              mFDis = mean(FDis),
              mFDiv= mean(FDiv),
              mFOri = mean(FOri),
              mTaxRich = mean(Tax.rich),
              sd.Simpson= sd(Simpson, na.rm = TRUE),
              sd.FunRao= sd(FunRao, na.rm = TRUE),
              sd.FRed= sd(FRed, na.rm = TRUE),
              sd.FRic= sd(FRic, na.rm = TRUE),
              sd.FEve= sd(FEve, na.rm = TRUE),
              sd.FSpe= sd(FSpe, na.rm = TRUE),
              sd.FDis= sd(FDis, na.rm = TRUE),
              sd.FDiv= sd(FDiv, na.rm = TRUE),
              sd.FOri= sd(FOri, na.rm = TRUE),
              sd.TaxRich= sd(Tax.rich, na.rm = TRUE),
              n.FRic = n())%>%
    mutate(se.Simpson = sd.Simpson / sqrt(n.FRic),
           se.FunRao = sd.FunRao / sqrt(n.FRic),
           se.FRed = sd.FRed / sqrt(n.FRic),
           se.FRic = sd.FRic / sqrt(n.FRic),
           se.FEve = sd.FEve / sqrt(n.FRic),
           se.FSpe = sd.FSpe / sqrt(n.FRic),
           se.FDis = sd.FDis / sqrt(n.FRic),
           se.FDiv = sd.FDiv / sqrt(n.FRic),
           se.FOri = sd.FOri / sqrt(n.FRic),
           se.TaxRich = sd.TaxRich / sqrt(n.FRic),

           lower.ci.Simpson = mSimpson - qt(1 - (0.05 / 2), n.FRic - 1) * se.Simpson,
           upper.ci.Simpson = mSimpson + qt(1 - (0.05 / 2), n.FRic - 1) * se.Simpson,

           lower.ci.FunRao = mFunRao - qt(1 - (0.05 / 2), n.FRic - 1) * se.FunRao,
           upper.ci.FunRao = mFunRao + qt(1 - (0.05 / 2), n.FRic - 1) * se.FunRao,

           lower.ci.FRed = mFRed - qt(1 - (0.05 / 2), n.FRic - 1) * se.FRed,
           upper.ci.FRed = mFRed + qt(1 - (0.05 / 2), n.FRic - 1) * se.FRed,

           lower.ci.FRic = mFRic - qt(1 - (0.05 / 2), n.FRic - 1) * se.FRic,
           upper.ci.FRic = mFRic + qt(1 - (0.05 / 2), n.FRic - 1) * se.FRic,

           lower.ci.FEve = mFEve - qt(1 - (0.05 / 2), n.FRic - 1) * se.FEve,
           upper.ci.FEve = mFEve + qt(1 - (0.05 / 2), n.FRic - 1) * se.FEve,

           lower.ci.FSpe = mFSpe - qt(1 - (0.05 / 2), n.FRic - 1) * se.FSpe,
           upper.ci.FSpe = mFSpe + qt(1 - (0.05 / 2), n.FRic - 1) * se.FSpe,

           lower.ci.FDis = mFDis - qt(1 - (0.05 / 2), n.FRic - 1) * se.FDis,
           upper.ci.FDis = mFDis + qt(1 - (0.05 / 2), n.FRic - 1) * se.FDis,

           lower.ci.FDiv = mFDiv - qt(1 - (0.05 / 2), n.FRic - 1) * se.FDiv,
           upper.ci.FDiv = mFDiv + qt(1 - (0.05 / 2), n.FRic - 1) * se.FDiv,

           lower.ci.FOri = mFOri - qt(1 - (0.05 / 2), n.FRic - 1) * se.FOri,
           upper.ci.FOri = mFOri + qt(1 - (0.05 / 2), n.FRic - 1) * se.FOri,

           lower.ci.TaxRich = mTaxRich - qt(1 - (0.05 / 2), n.FRic - 1) * se.TaxRich,
           upper.ci.TaxRich = mTaxRich + qt(1 - (0.05 / 2), n.FRic - 1) * se.TaxRich)

  return(nullmodelsummary)
}

nullsummary<- .nullsummary(truenull=nulldf)

#save simulated Fd metrics for each null community
write.csv(nullsummary, file='data_output/nullsummary_NEE.csv')

#save summary of null models  as csv
write.csv(nulldf, file='data_output/nulls_NEE.csv')

