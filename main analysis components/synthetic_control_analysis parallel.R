# Clear workspace to get rid of old junk
#rm(list=ls(all=TRUE))

#These 2 packages need to be installed and loaded first
#install.packages('devtools')
require(devtools)
#install.packages('abind')
require(abind)

#read in function file--can read directly from github using devtools function 'source_url'
source_url("https://raw.githubusercontent.com/weinbergerlab/Brazil_state/master/functions%20glm%20aic%20mod%20ave.R")      
#source('functions glm aic mod ave.R')

packages <- c('RCurl','reshape2','RColorBrewer','matlib','lme4', 'dummy','knitr','plotly','MASS', 'parallel','splines', 'lubridate','devtools')
packageHandler(packages, update_packages = FALSE, install_packages = FALSE) #change to true if need to install any
sapply(packages, library, quietly = TRUE, character.only = TRUE)


#Set working directory: default to desktop--different path for windows vs Mac
if(.Platform$OS.type == "windows") {
  desktop<-file.path(Sys.getenv("USERPROFILE"),"Desktop")
} else {
  desktop<- "~/Desktop"
}

setwd(desktop)

#read in function file--can read directly from github
source_url("https://raw.githubusercontent.com/weinbergerlab/Brazil_state/master/functions%20glm%20aic%20mod%20ave.R")      
#source('functions glm aic mod ave.R')

###############USER SPECIFIED VALUES
country="Brazil"
all_cause_name <- 'ach_noj'
all_cause_pneu_name <- 'J12_18'
noj_name <- 'cJ20_J22'
date_name<-'date'
data_start_date <- as.Date('2004-01-01')   #When do you want the analysis to start? (yyyy-mm-01)
data_intervention_start<-as.Date('2010-01-01')   #When is the intervention introduced?  (yyyy-mm-01)
n_seasons=12 #12 for monthly, 4 for quarterly, 3 for triannually
N.sim=10000 #total number of random draws for the predictive distribution
bivariate=FALSE  #two control variables or 1 in regression?
season.control='dummy'  #can either by 'harmonic' or 'dummy'
#################################################
#IMPORT DATA 
#(sample data pulled directly from github)
git.data= getURL("https://raw.githubusercontent.com/weinbergerlab/Brazil_state/master/prelog_Brazil_state_processed_data.csv")      
d1 <- read.csv(text=git.data)
############################################  

output_directory<-paste0(getwd(),'/output/')
dir.create(output_directory, recursive=TRUE, showWarnings = FALSE)


d1.list<-split(d1,d1$age_group)
names(d1.list)

if(country %in% c("Brazil","Brazil_state")){
  agegrp02<-which(substr(names(d1.list),1,2)=="02")
  d1.list<-d1.list[-agegrp02]
}
##############################

#############################################
data.start<-1 #Determines which year to include; set to 1 to include full TS, to 61 to just 2008+      
reg.names<-names(d1.list)

#Initialize empty lists to store results
rr.post.q<- vector("list", length(reg.names)) 
rr.post.q.pca<- vector("list", length(reg.names)) 
rr.post.t<- vector("list", length(reg.names)) 
rr.post.t.pca<- vector("list", length(reg.names)) 
aic.df.save <- vector("list", length(reg.names)) 
quantiles<- vector("list", length(reg.names))
set.seed(123) #set random seed
###########LOOOP BY STATE
for(k in 1:length(reg.names)){
  in.data<-d1.list[[reg.names[k]]]
  ##################################################
  ##################################################
  #SECTION 1: IMPORTING AND FORMATTING TIME SERIES
  
  strata.label=unique(in.data$age_group)
  # output_directory<-paste('C:/Users/DMW63/Desktop/My documents h/Gates/CausalImpact code/Brazil_state/', strata.label,'/', sep='')
  #  dir.create(output_directory, recursive=TRUE, showWarnings = FALSE)
  ds1a<-in.data
  age_groups <- paste( unique(unlist(ds1a$age_group, use.names = FALSE)))
  
  #filter if  column mean<5?
  exclude<-c( which(colMeans(ds1a[,-c(1:3)])<5 )+3, which(is.na(colMeans(ds1a[,-c(1:3)])))+3)
  ds1a<-ds1a[,-exclude ] #Position of variables to filter out
  
  print(strata.label)
  potential_covars<-names(ds1a)[4:ncol(ds1a)]
  
  data_intervention_date <- data_intervention_start %m-% days(1)  #Day before intervention date
  
  data_end_date <- max(as.Date(ds1a[, date_name]))
  pre_period <- c(data_start_date, data_intervention_date)  #define training period
  post_period <- c(data_intervention_date + 1, data_end_date) #Define post-vaccine period
  #Define evaluation period
  
  eval_period <- c( (data_end_date %m-% months(24)), data_end_date) #take last 24m as evaluation period
  
  #Log transform covariates
  ds <- ds1a
  ds$date<-as.Date(ds$date)
  ds <- ds[, colSums(is.na(ds)) == 0] #Delete columns with NAs
  ds <- ds[match(data_start_date, ds$date):nrow(ds),]
  ds[ds == 0] <- 0.5
  ds[, 4:ncol(ds)] <- log(ds[, 4:ncol(ds)])
  
  data_start <- match(data_start_date, ds$date)
  time_points <- ds$date[data_start:nrow(ds)]
  post.start.index<-which(time_points==pre_period[2]+1)
  eval.start.index<-which(time_points==eval_period[1])
  
  if(ncol(ds)>=4){
    covars.raw <- as.data.frame(ds[data_start:nrow(ds), 4:ncol(ds)])
    names(covars.raw)<-names(ds)[4:ncol(ds)]
  }else{
    covars.raw<-as.data.frame(rep(1, nrow(ds)))
    names(covars.raw)<-'one'
  }
  if(ncol(ds)>4){
    covars.raw <- ds[data_start:nrow(ds), 4:ncol(ds)]
  }else{
    covars.raw <- as.data.frame(ds[data_start:nrow(ds), 4])
    names(covars.raw)<-names(ds)[4]
  }
  month_i <- as.factor(as.numeric(format(ds[, date_name][data_start:nrow(ds)], '%m')))
  if(country=="Brazil"){
    spline <- setNames(as.data.frame(bs(1:nrow(covars.raw), knots = 5, degree = 3)), c('bs1', 'bs2', 'bs3', 'bs4'))
    year_2008 <- numeric(nrow(covars.raw))
    year_2008[1:nrow(covars.raw) >= match(as.Date('2008-01-01'), ds[, date_name])] <- 1
    data <- cbind.data.frame(year_2008, spline, month_i)
    trend <- lapply(covars.raw, getTrend, data = data)
    covars.raw <- covars.raw - trend
  }
  
  pandemic <- ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0))
  variance.covars<-apply(as.data.frame(covars.raw[1:(post.start.index-1),]),2,var, na.rm=TRUE)
  with.variation<-which(variance.covars>0)
  no.variation<-which(variance.covars==0)
  names.keep<-names(covars.raw)[with.variation]
  if(length(no.variation)>0){covars.raw<-covin.data=ars.raw[,with.variation]} #eliminates covariates that have NO VARIABIITY IN PRE-PERIOD
  covars.raw<-as.data.frame(covars.raw)
  names(covars.raw)<-names.keep
  
  ####################################
  ####################################
  #SECTION 2: CREATING SMOOTHED VERSIONS OF CONTROL TIME SERIES AND APPENDING THEM ONTO ORIGINAL DATAFRAME OF CONTROLS
  #EXTRACT LONG TERM TREND WITH DIFFERENT LEVELS OF SMOOTHNESS USING STL
  # Set a list of parameters for STL
  t.windows <- c(5,25,59)
  covars.raw.compile<- vector("list", length(t.windows)) 
  s.windows <- "periodic"
  # STL
  for (value_t in 1:length(t.windows)) {
    for (value_s in 1:length(s.windows)) {
      t <- t.windows[value_t]
      s <- s.windows[value_s]
      stl.covars <- DoSTL_trend(covars.raw,t,s)
      covars.raw.compile[[value_t]] <-stl.covars
    }
  }
  covars.raw2<-do.call(cbind, covars.raw.compile)
  covars.raw2<-cbind(covars.raw,covars.raw2)
  covars<-covars.raw2 #COMBINE ALL VARIABLES WITH DIFFERENT SMOOTHING LEVELS, FROM RAW TO VERY SMOOTH
  ##########################################################################################
  ###########################################################################################SECTION 3: ADD MONTHLY DUMMY, TIME TREND VARIABLES TO DATASET
  covars$t<-1:nrow(covars)/120
  covars$nocovars<-rep(1, times=nrow(covars)) #constant
  data.sel<-covars
  covar<-data.sel
  #Create Monthly dummies#Monthly dummies
  if(n_seasons==4){x<-quarter(as.Date(time_points))}
  if(n_seasons==12){x<-month(as.Date(time_points))}
  if(n_seasons==3){
    x.m<-month(as.Date(time_points))
    x<-x.m
    x[x.m %in% c(1,2,3,4)]<-1
    x[x.m %in% c(5,6,7,8)]<-2
    x[x.m %in% c(9,10,11,12)]<-3
  }
  x.df<-as.data.frame(as.factor(x))
  names(x.df)<-'m'
  dummy<-dummy::dummy #resolves function name conflict with lme4
  season.dummies<-  as.numeric(as.character(as.matrix(dummy(x.df))))
  season.dummies<-matrix(season.dummies, nrow=length(time_points), ncol=n_seasons)
  season.dummies<-season.dummies[,-n_seasons]
  dimnames(season.dummies)[[2]]<-paste0('m',1:(n_seasons-1))
  predictors<-as.matrix(covars)
  #Altrnatively, use harmonics
  index=1:length(time_points)
  sint1<-sin(2*pi*index/n_seasons)
  cost1<-cos(2*pi*index/n_seasons)
  #predictors <-cbind(predictors, ifelse(time_points == '2009-08-01', 1, ifelse(time_points == '2009-09-01', 1, 0)))
  #dimnames(predictors)[[2]][ncol(predictors)]<-'pandemic'
  if(country=='Fiji'){
    outcome<-as.vector(ds[data_start:nrow(ds),'pneuvar']) 
  }else{
    outcome<-as.vector(ds[data_start:nrow(ds),all_cause_pneu_name]) 
  }
  outcome.pre<-outcome
  outcome.pre[post.start.index:length(outcome)]<-NA
  #COMBINE MONTHLY DUMMIES AND COVARIATES AND PANDEMIC INTO SINGLE DATAFRAME
  if(season.control=="dummy"){ covar.matrix<-cbind.data.frame(season.dummies,pandemic,predictors)
  }else{
    covar.matrix<-cbind.data.frame(sint1, cost1, pandemic,predictors)
  }
  covar.lab<-dimnames(covar.matrix)[[2]]
  covar.matrix<-apply(covar.matrix,2,scale)
  
  time<-1:nrow(data.sel)
  time_post<-(time[post.start.index:length(outcome)]-post.start.index+1)/100
  data.fit<-cbind.data.frame(outcome.pre, covar.matrix)
  data.fit$outcome.pre<-as.integer(data.fit$outcome.pre)
  one<-rep(1,times=nrow(data.fit))
  
  
  
  #INITIALIZE VARIOUS LISTS TO STORE RESULTS
  aic.test <- vector(mode="numeric", length=ncol(covars))
  V<- vector("list",  length=ncol(covars)) #combine models into a list
  pred.mean<- vector("list",  length=ncol(covars)) #combine models into a list
  pred.coefs<- vector("list",  length=ncol(covars)) #combine models into a list
  covar.test<- vector("list",  length=ncol(covars)) #combine models into a list
  coef1<- vector("list",  length=ncol(covars)) #combine models into a list
  preds.stage2<- vector("list",  length=ncol(covars)) #combine models into a list
  mod1<- vector("list",  length=ncol(covars)) #combine models into a list
  ds.fit<- vector("list",  length=ncol(covars)) #combine models into a list
  pred.mean<- matrix(NA, nrow=nrow(covars), ncol=ncol(covars)) #combine models into a list
  aic.test[]<-NA
  
  #Create a list that contains all of the dataset variants used for fitting
  #If using AIC-based selection, they need to all have exact same X variable with same N observations
  
  if(bivariate==TRUE){
    combos<-  t(combn(dimnames(predictors)[[2]],2))
    same<- substr(combos[,1],1,3)==substr(combos[,2],1,3)
    combos[-same,] #only allo combo of different covariates, not of same
    combos2<-paste0(combos[,1], ',', combos[,2] )
    combos3<-c(dimnames(predictors)[[2]],combos2 )
  }else{
    combos3<-dimnames(predictors)[[2]]
  }
  
  names(data.fit)[1]<-'y'
  for(p in 1:length(combos3)){
    if(grepl('nocovars',  combos3[p] )) {
      if (season.control=='dummy'){
        incl.names<-c('y','m1','m2','m3','m4','m5','m6', 'm7','m8','m9','m10','m11','pandemic' ) #NULL MODEL
      }else{
        incl.names<-c('y','sint1','cost1','pandemic' ) #NULL MODEL
      }
    }else{
      if(season.control=='dummy'){
        incl.names<-c('y','m1','m2','m3','m4','m5','m6', 'm7','m8','m9','m10','m11','pandemic', combos3[p] )
      }else{
        incl.names<-c('y','sint1','cost1','pandemic', combos3[p]  ) 
      }
    }
    keep.cols<-which(names(data.fit) %in% incl.names )
    ds.fit[[p]]<-data.fit[,keep.cols]
    comment(ds.fit[[p]])<-combos3[p]
  }
  ######################################################################
  ######################################################################
  ##SECTION 4: RUN INITIAL UNIVARIATE MODELS 
  n_cores <- detectCores()-1
  cl1 <- makeCluster(n_cores)
  clusterEvalQ(cl1, {library(lme4, quietly = TRUE)})
  clusterExport(cl1, c('ds.fit',  'glm.fun', 'time_points', 'n_seasons','post.start.index'), environment())
  
  glm.results<-parLapply(cl=cl1 ,     ds.fit, glm.fun )
  stopCluster(cl1)
  
  #glm.results<-lapply(   ds.fit, glm.fun )
  
  
  #glm.results<-tryCatch(lapply(ds.fit, glm.fun ), error=function(e) lapply(ds.fit, glm.fun.pois ) ) 
  
  #Extract AICs from list into dataframe
  aics<-unlist(lapply(glm.results, '[[', 'aic.test'))  # This returns a vector with AIC score
  vars<-unlist(lapply(glm.results, '[[', 'test.var'))  # This returns a vector with the variable names
  pred.mean<-lapply(glm.results, '[[', 'pred.mean') # This returns a vector with the variable names
  pred.mean<-do.call(cbind,pred.mean)
  pred.mean<-exp(pred.mean)
  aic.df<-cbind.data.frame(vars, aics)
  names(aic.df)<-c('covars','aic')
  aic.df$model.index<-1:nrow(aic.df)
  aic.df$grp<-as.numeric(as.factor(substr(aic.df$covars,1,3))) #for each smoothed or unsmoothed version of variable, assign it to a grouping 
  aic.df$delta.aic<-aic.df$aic-min(aic.df$aic)
  aic.df$w_aic<- exp(-0.5*aic.df$delta.aic)/sum( exp(-0.5*aic.df$delta.aic))
  aic.df<-aic.df[order(-aic.df$w_aic),]
  aic.df$cumsum<-cumsum(aic.df$w_aic) 
  aic.df$keep.high.weight<- aic.df$cumsum<=0.99  #only keep variables that contribute to 99% f weight
  aic.df$model.rank<- 1:nrow(aic.df)  #only keep variables that contribute to 99% f weight
  aic.df<-aic.df[order(aic.df$model.index),]
  aic.df$Nsamps<- rmultinom(n=1,size=N.sim, prob=aic.df$w_aic) #How many samples do we want from each model in total?
  for(i in 1:nrow(aic.df)){glm.results[[i]]$Nsamps<-aic.df$Nsamps[i,1]}
  aic.df<-aic.df[order(aic.df$grp, -aic.df$w_aic),]
  model.weight<-aic.df[,c('covars','w_aic', 'model.rank')]
  model.weight<-model.weight[order(-model.weight$w_aic),]
  model.weight$model.rank<-1:nrow(model.weight)
  
  #top.covar in each set
  aic.df2<-aic.df
  aic.df2<-aic.df2[order(-aic.df$w_aic),]
  aic.df2$deseason=0
  aic.df2$deseason[grep('trend',aic.df2$covars, fixed=TRUE)] <- 1
  aic.df2<-aic.df2[aic.df2$deseason==1,] #Only keep STL version of variable, not raw
  aic.df2$w_aic<-aic.df2$w_aic/sum(aic.df2$w_aic) #rescale weights
  top.covar.grp<-aic.df2[!duplicated(aic.df2$grp),]
  remove<-c('t','nocovars')
  top.covars<-as.character(top.covar.grp$covars[! top.covar.grp$covars %in% remove])
  
  
  ####Model averaging results
  preds.stage2<-lapply(glm.results,obs.uncertainty)
  all.preds<-do.call(cbind, preds.stage2)
  all.preds<-all.preds[ , colSums(is.na(all.preds)) == 0]
  preds.q<-t(apply(all.preds, 1,quantile, probs=c(0.025,0.5,0.975)))
  rr.post.t[[k]]<- outcome/preds.q
  #rr during eval period
  post.preds<-all.preds[eval.start.index:nrow(all.preds),]
  post.preds.sums<-apply(post.preds,2,sum)
  post.obs.sum<-sum(outcome[eval.start.index:nrow(all.preds)])
  post.rr<-post.obs.sum/post.preds.sums
  rr.post.q[[k]]<-quantile(post.rr,probs=c(0.025,0.5,0.975))            
  rr.mean.post.mod<-rr.post.q[[k]]
  #Pointwise CIs
  log_rr_full_t<-log((outcome+0.5)/(all.preds+0.5))
  log_rr_full_t_quantiles<-apply(log_rr_full_t,1, quantile, probs=c(0.025,0.5,0.975))
  log_rr_full_t_sd<-t(apply(log_rr_full_t, 1, sd, na.rm = TRUE))
  log_rr_full_t_samples.covar<-cov(t(log_rr_full_t))
  log_rr_full_t_samples.prec<-solve(log_rr_full_t_samples.covar)
  
  ###SETUP AND RUN MODELS WITH FIRST PC
  covars.keep.pca<-covars[,top.covars]
  pca <- prcomp(covars.keep.pca, scale=TRUE) # scale=TRUE should be added!!
  predictors2 <- as.data.frame(pca$x[,1]) # First "1" PC
  names(predictors2)<-'pca1'
  if(season.control=="dummy"){ covar.matrix.pca<-cbind.data.frame(outcome.pre,season.dummies,pandemic,predictors2)
  }else{
    covar.matrix.pca<-cbind.data.frame(outcome.pre,sint1, cost1, pandemic,predictors2)
  }
  covar.matrix.pca$obs<-as.factor(1:nrow(covar.matrix.pca))
  mod1<-glmer(outcome.pre~ m1+m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+pandemic+pca1 +(1|obs) ,data=covar.matrix.pca, family='poisson',control=glmerControl(optimizer="bobyqa",
                                                                                                                                                   optCtrl=list(maxfun=2e5)) )
  all.preds.pca<-simulate(mod1, nsim=10000, newdata=covar.matrix.pca, allow.new.levels=TRUE,re.form=NA)
  preds.q.pca<-t(apply(all.preds.pca, 1,quantile, probs=c(0.025,0.5,0.975)))
  rr.post.t.pca[[k]]<- outcome/preds.q.pca
  
  #rr during eval period from PCA model
  post.preds.pca<-all.preds.pca[eval.start.index:nrow(all.preds),]
  post.preds.sums.pca<-apply(post.preds.pca,2,sum)
  post.obs.sum<-sum(outcome[eval.start.index:nrow(all.preds)])
  post.rr.pca<-post.obs.sum/post.preds.sums.pca
  rr.post.q.pca[[k]]<-quantile(post.rr.pca,probs=c(0.025,0.5,0.975))
  rr.mean.post.mod.pca<-rr.post.q.pca[[k]]
  
  #point estimate for RR for each model
  log_rr_full_t.pca<-log((outcome+0.5)/(all.preds.pca+0.5))
  log_rr_full_t_quantiles.pca<-apply(log_rr_full_t.pca,1, quantile, probs=c(0.025,0.5,0.975))
  log_rr_full_t_sd.pca<-t(apply(log_rr_full_t.pca, 1, sd, na.rm = TRUE))
  log_rr_full_t_samples.covar.pca<-cov(t(log_rr_full_t.pca))
  log_rr_full_t_samples.prec.pca<-solve(log_rr_full_t_samples.covar.pca)
  quantiles[[k]]<-list(rr.mean.post.mod, outcome, log_rr_full_t_quantiles, log_rr_full_t_sd,log_rr_full_t_samples.prec,model.weight,
                       rr.mean.post.mod.pca, log_rr_full_t_quantiles.pca, log_rr_full_t_sd.pca,log_rr_full_t_samples.prec.pca)
  names(quantiles[[k]])<- c('rr.mean.post.mod', 'outcome', 'log_rr_full_t_quantiles', 'log_rr_full_t_sd','log_rr_full_t_samples.prec','model.weight',
                            'rr.mean.post.mod.pca', 'log_rr_full_t_quantiles.pca', 'log_rr_full_t_sd.pca','log_rr_full_t_samples.prec.pca')
  ##################################
  ##################################
  #STEP 5: PLOT RESULTS
  pdf(paste0(output_directory, 'plots_', country,strata.label, ".pdf"), width=7, height=5)
  
  #Plot observed vs expected +/- 95% CI
  matplot(preds.q.pca, type='l', lty=c(2,1,2), col='black', bty='l', ylim=c(0,max(preds.q)))
  points(outcome)
  abline(v=post.start.index,col='gray', lty=2)
  title(paste0(reg.names[k]," STL+PCA: Observed vs expected"))
  
  matplot(preds.q, type='l', lty=c(2,1,2), col='black', bty='l', ylim=c(0,max(preds.q)))
  points(outcome)
  abline(v=post.start.index,col='gray', lty=2)
  title(paste0(reg.names[k]," STL+BMA: Observed vs expected"))
  
  matplot(preds.q.pca, type='l', col='#af8dc3', lty=1, bty='l')
  points(outcome)
  points(preds.q[,'50%'], lty=1, col='#7fbf7b',type='l' )
  points(preds.q[,'2.5%'], lty=2, col='#7fbf7b',type='l' )
  points(preds.q[,'97.5%'], lty=2, col='#7fbf7b',type='l' )
  abline(v=post.start.index,col='gray', lty=2)
  title("Comparison of STL+BMA (green) with STL+PCA (purple)")
  
  matplot(t(log_rr_full_t_quantiles.pca),  type='l', lty=c(2,1,2), ylim=c(-1,1),col='black', bty='l', ylab="Log(Rate Ratio)")
  abline(h=0, col='gray', lty=2)
  abline(v=post.start.index,col='gray', lty=2) 
  title("STL+PCA:Log rate ratios")
  
  matplot(t(log_rr_full_t_quantiles),  type='l', lty=c(2,1,2), ylim=c(-1,1), col='black', bty='l', ylab="Log(Rate Ratio)")
  abline(h=0, col='gray', lty=2)
  abline(v=post.start.index,col='gray', lty=2) 
  title("STL+BMA:Log rate ratios")
  
  #
  #Plot expected from individual models vs observed
  col.plot<-rgb(0,0,0, alpha=aic.df$w_aic)
  #All preds
  matplot(pred.mean, type='l', col='gray', lty=1, bty='l')
  points(outcome)
  abline(v=post.start.index,col='gray', lty=2)
  title("Predictions from individual variables")
  
  #by weight
  matplot(pred.mean, type='l', col=col.plot, lty=1, bty='l')
  points(outcome)
  abline(v=post.start.index,col='gray', lty=2)
  title("Predictions from individual variables; weighted by fit")
  
  
  
  dev.off()
  
  
  
  ## MAKE THE INTERACTIVE PLOTLY PLOTS
  model_weights_2<-aic.df$w_aic
  model_weights_2[aic.df$w_aic<0.03]=0.03 #set a minimum threshold for weights so that even 0 weighted lines show up faintly
  mods.covars<-names(covars)
  format.mods.covars<-paste(round(rr.mean.post.mod,2) )
  
  col.plot<-rgb(0,0,0, alpha=model_weights_2)
  plot.order<-aic.df$model.rank
  n.plot<-sum(aic.df$keep.high.weight) #how many covars contriubuted to 99% of weight?
  n.plot<-min(n.plot,50) #limit to number of lines plotted
  outcome.list<-outcome
  plotlies <- plot_ly(x = ~time_points, y = ~outcome.list, type = 'scatter', mode = 'markers')%>%
    layout(
      yaxis=list(rangemode='tozero')
    )
  
  for(i in 1:n.plot ) {
    j<-which(plot.order==i)
    plotlies <- add_trace(plotlies,x=time_points,y = pred.mean[,j] ,name = mods.covars[j],  mode = "lines",
                          evaluate = TRUE,line = list(color = col.plot[j]),  hoverinfo = 'text',
                          text = paste0(mods.covars[j]," RR:",format.mods.covars[j]," Weight ",round(aic.df$w_aic,2)[j]))  
  }           
  htmlwidgets::saveWidget(plotlies, paste0(output_directory,'plottly.plots.', country,strata.label, '.html' ))
  
}

names(quantiles)<-reg.names
saveRDS(quantiles,'./output_STL/results.rds' )

rr.bma<-do.call(rbind, lapply(quantiles,function(xl) xl$'rr.mean.post.mod') )
rr.bma<-cbind.data.frame(reg.names,rr.bma)
rr.pca<-do.call(rbind, lapply(quantiles,function(xl) xl$'rr.mean.post.mod.pca') )
rr.pca<-cbind.data.frame(reg.names,rr.pca)
write.csv(rr.bma,'./output_STL/rr.bma.csv')
write.csv(rr.pca,'./output_STL/rr.pca.csv')

