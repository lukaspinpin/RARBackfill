##### packages
library(bcrm)
library(dfcrm)
library(extraDistr)
library(collapse)
library(mcp)
library(matrixStats)
library(Rlab)
library(rjags)
library(R2WinBUGS)

## functions
rm(list = ls(all.names = TRUE))

# SOME USEFUL FUNCTIONS

# logit function
logit_transform <- function(x) log(x/(1-x))

# "unlogit" function
unlogit_transform <- function(x) exp(x)/(1+exp(x))


### TS################################

# define a function to generate a beta distribution for a given arm
generate_beta <- function(alpha, beta) {
  rbeta(1, alpha, beta)
}

# define a function to select the arm with the highest expected value
select_arm <- function(num_arms, success_counts, failure_counts) {
  # generate a beta distribution for each arm
  beta_values <- sapply(1:num_arms, function(i) generate_beta(success_counts[i] + 1, failure_counts[i] + 1))
  # select the arm with the highest expected value
  return(which.max(beta_values))
}

# Function for Backfill
eff_mTS <- function(n=3, eff, doses){
  # num_trials <- n 
  levels <- sort(unique(doses))
  num_arms <-  length(levels)
  
  # initialize the success and failure counts for each arm
  success_counts <- rep(0, num_arms)
  failure_counts <- rep(0, num_arms)
  
  #Create Data Set
  data <- data.frame(eff, doses)
  
  #Count Successes and Failures
  success_counts <- sapply(1:num_arms, function(i) sum(data[data$doses==i,]$eff))
  failure_counts <- sapply(1:num_arms, function(i) length(data[data$doses==i,]$eff)) - success_counts
  
  #Choose Dose level 
  backfillSET <- c()
  for(i in 1:n){
    backfillSET <-  c(backfillSET, select_arm(num_arms, success_counts, failure_counts))
  }
  return(backfillSET)
}


################################################################################################################

### dose-finding plus randomisation

################################################################################################################

simulation <- function(n.cohorts=10, n.param ="1PM", n.sim=10, scenario="A", method = "TS"){
  
  TTL<-1/4 # Target Toxicity Level is 25%
  delta<-0.05 # Want DLT risk at MTD to be between TTL +- 5 pct points 
  n.doses<-7 # number of dose levels
  model.choice<- "empiric" # We will use the power (empiric) model
  
  #Scenario Specfic parameters
  if(scenario=="A"){
    expected.mtd<- 5 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- tox.skeleton
    eff.skeleton=c(0.05,0.15,0.25,0.25,0.25,0.25,0.25) # Scenario A
  }
  if(scenario=="E"){ 
    expected.mtd<- 5 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- tox.skeleton
    eff.skeleton=c(0.05,0.1,0.15,0.15,0.15,0.15,0.15) # Scenario B
  }
  if(scenario=="C"){ 
    expected.mtd<- 7 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- seq(0.05,0.35,0.05) # Scenario C,D 
    eff.skeleton=c(0.04,0.08,0.12,0.16,0.20,0.24,0.28) # Scenario C
  }
  if(scenario=="D"){ 
    expected.mtd<- 7 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- seq(0.05,0.35,0.05) # Scenario C,D 
    eff.skeleton=c(0.07,0.14,0.21,0.28,0.35,0.42,0.49) # Scenario D
  }
  if(scenario=="G"){ 
    expected.mtd<- 4 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- seq(0.025,0.475,0.075) # Scenario E,F
    eff.skeleton=c(0.05,0.1,0.15,0.2,0.2,0.2,0.2) # Scenario E 
  }
  if(scenario=="H"){ 
    expected.mtd<- 4 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- seq(0.025,0.475,0.075) # Scenario E,F
    eff.skeleton=c(0.04,0.08,0.12,0.16,0.2,0.24,0.24) # Scenario F 
  }
  if(scenario=="F"){ #F
    expected.mtd<- 7 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- seq(0.05,0.35,0.05) # Scenario G,H
    eff.skeleton=c(0.05,0.1,0.15,0.2,0.25,0.25,0.25) # Scenario G
  }
  if(scenario=="B"){ # B
    expected.mtd<- 7 # Prior belief is that dose level x is the MTD # A5, B5, C7, D7, E4, F4, G7, H7
    tox.skeleton<-getprior(target = TTL, halfwidth = delta, nlevel = n.doses, nu = expected.mtd, model = model.choice)
    p.tox0 <- seq(0.05,0.35,0.05) # Scenario G,H
    eff.skeleton=c(0.08,0.16,0.24,0.32,0.4,0.4,0.4) # Scenario H
  }
  
  round(p.tox0,2)
  
  truep <- tox.skeleton
  truep
  
  dose.labels <- seq(1:n.doses)
  
  ######################
  
  ## plot the resulting scenario
  
  ######################
  filename <- paste("Scenario",scenario, "plot.pdf",sep="_")
  pdf(file=filename)
 
  plot(eff.skeleton ~ dose.labels,type="b",col="blue", ylim=c(0,1),lwd=2, lty=2, xlab="dose level",ylab="probability of DLE / efficacy", axes=FALSE)
  axis(1, at=1:7)
  axis(2)
  lines(tox.skeleton ~ dose.labels,type="b", col="lightsalmon",lwd=2, lty=2)
  abline(h=0.25,lty=2)
  legend("topright", legend=c("true dose-DLE curve", "true dose-efficacy curve"), #,"target DLT level"
         col=c("lightsalmon", "blue", "black"), lty=c(2,2), lwd=c(2,2),
         box.lty=0)
  text(x=1.1,y=1,cex = 1.25,labels=paste("(",scenario,")",sep=""))
  
  dev.off()
  
  ######################
  
  ## define the first dose and set up the simulations
  
  ######################
  
  dose=1
  consecutive.cohort=1
  
  n.per.cohort=3
  
  if(method=="CRM"){
    n.per.backfill=0 
    n.cohorts=n.cohorts*2-1  
  } else{
    n.per.backfill=3 
  }
  
  doses <- c()
  DLTs <- c()
  responses <- c()
  sims <- c()
  IDs <- c()
  type.patients <- c()
  
  experimentations <- matrix(0,nrow=n.sim,ncol=n.doses)
  colnames(experimentations) <- dose.labels
  head(experimentations)
  
  experimentations.backfill <- matrix(0,nrow=n.sim,ncol=n.doses)
  colnames(experimentations.backfill) <- dose.labels
  head(experimentations.backfill)
  
  number.responses <- matrix(0,nrow=n.sim,ncol=n.doses)
  colnames(number.responses) <- dose.labels
  
  estimated.probas <- matrix(0,nrow=n.sim,ncol=n.doses)
  colnames(number.responses) <- dose.labels
   
  recommendations <- c()
  plateau.or.not <- c()
  start.plateau <- c()
  patients.per.dose <- data.frame(dose1=1,dose2=1,dose3=1,dose4=1,dose5=1,dose6=1,dose7=1)
  
  lowest.acceptable.dose.backfill.at.end <-c()
  
  ######################
  
  ## start of the loop 
  
  ######################
  
  # z=1
  # i=1
  max.dose <- 1
  
  for(z in 1:n.sim){
  
    u=1 # for the bayesian beta-binomial model comparing the lowest dose, with the combination of the other doses
    
    for(i in 1:n.cohorts){
    
    if(dose==1){
      doses.received.backfill <- c()
    }
    
    if(dose==2){
      doses.received.backfill <- rcat(n=n.per.backfill, prob = rep(1/length(doses.below),length(doses.below)) ) # Lukas this is where the allocation of backfill patients to the dose level is made
    }
    
    if(dose>2){
      if(method=="ER"){
        doses.received.backfill <- rcat(n=n.per.backfill, prob = rep(1/length(doses.below),length(doses.below)),labels=as.factor(doses.below) ) 
        doses.received.backfill <- as.numeric(levels(doses.received.backfill))[doses.received.backfill]
      }
      if(method=="TS"){
        dose.curr <- length(doses.below)
        eff.curr <- subset(data, dose<=dose.curr & sim==z)$efficacy
        doses.curr <- as.factor(subset(data, dose<=dose.curr & sim==z)$dose)
        doses.received.backfill <- eff_mTS(n=n.per.backfill, eff=eff.curr, doses=doses.curr)
      }
    }
    
    type.backfill=rep("backfill",length(doses.received.backfill))
    
    doses.received.regular=rep(dose,n.per.cohort)
    type.reg=rep("regular",length(doses.received.regular))
    
    type=c(type.backfill,type.reg)
    type.patients <- c(type.patients,type)
    
    if(dose>=2){
      doses <- c(doses,doses.received.backfill,doses.received.regular)
    }else{
      doses <- c(doses,doses.received.regular)
    }
    
    toxicities.regular <- rbinom(n.per.cohort,1,prob=truep[dose])
    toxicities.regular
    
    toxicities.backfill <- c()
    
    if(length(doses.received.backfill)>=1){
      for(t in 1:length(doses.received.backfill)){
        toxicities.backfill <- c(toxicities.backfill, rbinom(1,1,prob=truep[doses.received.backfill[t]]) )
      }
    }
    
    if(dose>=2){
      DLTs <- c(DLTs,toxicities.backfill,toxicities.regular)
    }else{
      DLTs <- c(DLTs,toxicities.regular)
    }
    
    eff.regular <- rbinom(n.per.cohort,1,prob=eff.skeleton[dose]) # Lukas this is where the efficacy results are generated
    eff.backfill <- c() 
    if(length(doses.received.backfill)>=1){
      for(t in 1:length(doses.received.backfill)){
        eff.backfill <- c(eff.backfill, rbinom(1,1,prob=eff.skeleton[doses.received.backfill[t]]) )
      }
    }
    
    if(dose>=2){
      responses <- c(responses,eff.backfill,eff.regular)
    }else{
      responses <- c(responses,eff.regular)
    }
    
    simulation=rep(z,length(c(eff.regular,eff.backfill)))
    sims=c(sims,simulation)
    
    data=data.frame(dose=doses,tox=DLTs,efficacy=responses, type=type.patients, 
    sim=sims 
    )

    for(j in 1:n.doses){
    assign(paste0("indicator_dose_greater",j), ifelse(data$dose>j,1,0) )
    data <- cbind(data, get(eval(paste0("indicator_dose_greater",j))) ) 
    last.col = dim(data)[2]
    colnames(data)[last.col] <- paste0("indic_dose_greater",j)
    }
    head(data)
    
    patient_all <- c()
    for(x in 1:z){
      patient_all <- c(patient_all, seq(1:length(subset(data,sim==x)[,1])))
    }
    
    data <- cbind(patient_all,data)
    
    data_crm <- subset(data,type=="regular")
    data_backfill <- subset(data,type=="backfill")
    head(data_crm)
    head(data_backfill)
    
    patient <- c()
    for(x in 1:z){
      patient <- c(patient, seq(1:length(subset(data_crm,sim==x)[,1])))
    }
    
    data_crm <- cbind(patient,data_crm)
    
    if(n.param == "1PM"){
      trial.output <- bcrm(stop = list(nmax = length(subset(data_crm,sim==z & type=="regular")[,1])), 
                           data = subset(data_crm,sim==z & type=="regular"), # constrain=FALSE,
                           p.tox0 = p.tox0,dose = dose.labels,  ff = "power", 
                           prior.alpha=list(1, 1, 1), target.tox = TTL)
    }
    if(n.param == "2PM"){
      ## Bivariate lognormal prior for two parameters
      mu <- c(2.15, 0.52)
      Sigma <- rbind(c(0.84^2, 0.134), c(0.134, 0.80^2))
      trial.output <- bcrm(stop = list(nmax = length(subset(data_crm,sim==z & type=="regular")[,1])),
                             data = subset(data_crm,sim==z & type=="regular"), # constrain=FALSE,
                             p.tox0 = p.tox0, dose = dose.labels,  ff = "logit2",
                             prior.alpha=list(4, mu, Sigma), target.tox = TTL, method="rjags")
    }
    
    dose = trial.output$ndose[[1]][[1]] #dose
    
    doses.below = seq(1,dose-1) # what are the doses below the current recommended dose
    #doses.below
    
    consecutive.cohort=i+1
    
    max.dose <- max(subset(data_crm,sim==z)$dose) #max.dose
    
    } # end of cohort loop 
    
    lowest.acceptable.dose.backfill.at.end <- c(lowest.acceptable.dose.backfill.at.end,u)
    
    ## fit efficacy models (one with change-point, one without)
    collapsed <- collap(subset(data,sim==z), efficacy + tox ~ dose, FUN = list(fsum))
    count.dose <- as.data.frame(table(subset(data,sim==z)$dose))[,2]
    collapsed.dose <- cbind(collapsed,count.dose)
    #collapsed.dose
    
    # model with change point
    model = list(
      efficacy | trials(count.dose) ~ 1+dose,  # constant rate
      ~ 0)
    fit = mcp(model, collapsed.dose, family = binomial())
    #plot(fit, q_fit = TRUE)
    
    # model without change point
    model_null = list(
    efficacy | trials(count.dose) ~ 1+dose
    )
    
    fit_null = mcp(model_null, collapsed.dose, family = binomial())
    
    fit$loo = loo(fit)
    fit_null$loo = loo(fit_null)
    
    loo.comparison <- loo::loo_compare(fit$loo, fit_null$loo)
    best.fit <- rownames(loo.comparison)[1]
    
    plateau.or.not <- c(plateau.or.not, best.fit)
    
    if(best.fit=="model1"){
    start <- mean(c(fit$mcmc_post[[1]][,1],fit$mcmc_post[[2]][,1],fit$mcmc_post[[3]][,1]))
    start.plateau <- c(start.plateau, ceiling(start))
    }else{
    start.plateau <- c(start.plateau, NA)
    }
    
    # extract tox data at end of trial  
    estimated.probas[z,] < - trial.output$ndose[[1]][[2]]
    
    final.recommendation <- trial.output$ndose[[1]][[1]]
    dose=1 # reset dose to level 1
    consecutive.cohort=1 # reset cohort order to first cohort
    
    recommendations <- c(recommendations,final.recommendation)
    
    
    for(j in 1:n.doses){
      assign(paste0("count_dose",j), length(which(subset(data,sim==z)$dose == j)) )
      experimentations[z,j] <- get(eval(paste0("count_dose",j)))
      assign(paste0("count_dose",j), length(which(subset(data_backfill,sim==z)$dose == j)) )
      experimentations.backfill[z,j] <- get(eval(paste0("count_dose",j)))
    }
    
    for(j in 1:n.doses){
      assign(paste0("count_eff",j), sum(subset((subset(data,sim==z)),dose==j)$efficacy) )
      number.responses[z,j] <- get(eval(paste0("count_eff",j)))
    }
  } ## end of simulations
  
  ######################
  
  ## process results of simulations and save in data fram
  result_rownames <- c("Dose-response curve", "Dose-DLT curve","Skeleton", "% of recommendations for dose", 
                       "% of patients receiving dose","% of backfill patients receiving dose","Average sample size")
  
  results <- data.frame(row.names=result_rownames)
  empty_cols <- c("Dose1", "Dose2", "Dose3", "Dose4", "Dose5", "Dose6", "Dose7")
  results[ , empty_cols] <- NA
  
  
  ######################
  
  results[1,] <- eff.skeleton
  results[2,] <- p.tox0
  results[3,] <- tox.skeleton
  
  ######
   
  head(number.responses)
  head(experimentations)
  
  proportions.responses <- number.responses/experimentations
  proportions.responses[is.nan(proportions.responses)] <- NA
  head(proportions.responses)
  
  percentage.experimentations <- matrix(0,nrow=n.sim,ncol=n.doses)
  colnames(percentage.experimentations) <- dose.labels
  for(t in 1:n.sim){
  percentage.experimentations[t,] <- experimentations[t,]/rowSums(experimentations)[t]
  }
  
  average.patients.per.dose <- colMeans(percentage.experimentations)
  results[5,] <- average.patients.per.dose
  round(average.patients.per.dose,4)
  sum(average.patients.per.dose) 
  
  
  percentage.experimentations.backfill <- matrix(0,nrow=n.sim,ncol=n.doses)
  colnames(percentage.experimentations.backfill) <- dose.labels
   for(t in 1:n.sim){
     percentage.experimentations.backfill[t,] <- experimentations.backfill[t,]/rowSums(experimentations.backfill)[t]
  }
   
   average.patients.per.dose.backfill <- colMeans(percentage.experimentations.backfill)
   results[6,] <- average.patients.per.dose.backfill
   round(average.patients.per.dose.backfill,4)
   sum(average.patients.per.dose.backfill) 
   
  # just reversing for better file names
  if(method=="CRM"){
    n.cohorts = (n.cohorts+1)/2
  }
  
  filename <- paste(method, n.param, scenario, n.cohorts,"Barplot", "AllPatients.pdf",sep="_")
  pdf(file=filename)
  barplot(average.patients.per.dose*100, ylim=c(0,40), xlab = "Doses", ylab = "% of patients assigned")
  dev.off()
  
  if(method!="CRM"){
  data <- average.patients.per.dose.backfill[1:6] * 100
  filename <- paste(method, n.param, scenario, n.cohorts,"Barplot", "BackfillPatients.pdf",sep="_")
  pdf(file=filename)
  barplot(data,  ylim=c(0,40), xlab = "Doses", ylab = "% of backfill patients assigned")
  
  # Add percentages above the bars
  for (i in 1:6) {
    text(x = i, y = data[i] + 2, labels = paste(round(data[i],1), "%"), pos = 3)
  }
  
  dev.off()
  }
  
  
  percentage.recommendations <- matrix(0,nrow=1,ncol=n.doses)
  colnames(percentage.recommendations) <- dose.labels
  
  for(j in 1:n.doses){
    assign(paste0("count_recommended_dose",j), length(which(recommendations == j)) )
    percentage.recommendations[1,j] <- get(eval(paste0("count_recommended_dose",j)))/n.sim
  } 
  percentage.recommendations
  recommendations 
  
  average.sample.size <- mean(rowSums(experimentations))
  average.sample.size
  results[7,1] <- average.sample.size
  
  start.plateau
  table(start.plateau)/n.sim
  
  matrix.recommendations <- data.frame(tox=recommendations, efficacy=start.plateau)
  
  matrix.recommendations <- transform(matrix.recommendations, min = pmin(tox, efficacy))
  
  matrix.recommendations$min <- ifelse(is.na(matrix.recommendations$min), matrix.recommendations$tox, matrix.recommendations$min)
  
  plateau.identification <- table(is.na(matrix.recommendations[,2]))
  plateau.identification
  
  sum(ifelse(matrix.recommendations$min<matrix.recommendations$tox,1,0))
  plateau.to.mtd<- sum(ifelse(matrix.recommendations$min<matrix.recommendations$tox,1,0))
  plateau.to.mtd
  
  table(matrix.recommendations$min)
  
  
  if(method=="CRM"){
    results[4,] <- percentage.recommendations
  } else {
    for(k in 1:n.doses){
      results[4,k] <- sum(matrix.recommendations$min == k)/n.sim
    }
  }
  
  filename <-  paste(method, n.param, scenario, n.cohorts,"Results.csv",sep="_")
  write.csv(round(results,4), filename, row.names=TRUE)
}


############### Execute Simulation #######################

##################  One Scenario #######################
simulation(n.cohorts=10, n.param="1PM", n.sim=10^1, scenario="A", method = "TS") #10^3


#################### All Scenarios ###################### 
########## To Run on a High-Performance-Cluster ##########
library(foreach)
library(doParallel)

scenarios <- c("A","B","C","D","E","F","G","H")
n.cohs <- 10
methods <- c("CRM","ER","TS")
#n.param <- c("1PM","2PM")
n.param <- c("1PM")


# Generate list of parameter combinations to simulate over
param_list <- expand.grid(n.cohs = n.cohs, n.param=n.param, scenario = scenarios, method = methods)

# Set up a parallel backend with 4 cores
cl <- makeCluster(48) 
registerDoParallel(cl)

# Run simulation study in parallel
foreach(param = iter(param_list, by='row'), 
        .packages=c("bcrm", "dfcrm", "extraDistr", "collapse",
                   "mcp", "matrixStats", "pseudorank","Rlab")) %dopar% { #dopar
  print(param)
  simulation(n.cohorts=param$n.cohs, n.param=param$n.param, n.sim=10^3, scenario=param$scenario, method = param$method) #10^3
}

# Clean up parallel backend
stopCluster(cl)
###############################################################################