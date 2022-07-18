#Importing Data file into Rstudio and creating working directory 28/5/22
dir.create("FulmerData")
#Needed first time, not to reload
fulmars <- read.table("C:/Users/Robert/Documents/FulmarLife/FulmarBird/FulmerData/2018CaptureMatrix.txt", header = TRUE, sep = "\t", stringsAsFactors = TRUE)

#Creating summary tables with data
summary(fulmars)
head(fulmars)
#change N/a to 0, changes 2s to 1s


#taken from rows 4-55 to only analyse capture recapture data from dataset
toydata <- fulmars[,4:55]
#visualising new dataset
View(toydata)

#finding which values are "2"
toydata[toydata == "2"] 
#replace values which are "2" into "1"
toydata[toydata == "2"] <- "1"
#replace values which are "na" into "0"

toydata[is.na(toydata)] <- 0
fulmars[,4:55]<-toydata
rm(toydata)

#change sex 1 to male and sex 2 to female
#creating datafile of sexes
sexdata <- fulmars[,2]
#Changing 1 to male 
sexdata[sexdata == "1"] <- "male"
#Changing 2 to female
sexdata[sexdata == "2"] <- "female"

#change known age 0 as adultring and known age as chickring
agedata <- fulmars[,3]
View(agedata)
agedata[agedata == "0"] <- "adultring"
agedata[agedata == "1"] <- "chickring"


#run model for survival only
naive.survival.model <- nimbleCode({
  # prior
  phi ~ dunif(0, 1)
  # likelihood
  y ~ dbinom(phi, n)
})
my.data <- list(n = 450, y = 55)
my.constants <- list(n = 450)
my.data <- list(y = 55)
initial.values <- function() list(phi = runif(1,0,1))
initial.values()
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
n.thin <- 1
parameters.to.save <- c("phi")
mcmc.output <- nimbleMCMC(code = naive.survival.model,
                          data = my.data,
                          constants = my.constants,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          thin = n.thin,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)
str(mcmc.output)
head(mcmc.output$chain1)
library(MCMCvis)
MCMCsummary(mcmc.output, round = 2)
MCMCtrace(mcmc.output,
          pdf = FALSE,
          ind = TRUE,
          Rhat = TRUE,
          n.eff = TRUE) 
#run model for recapture only
hmm.survival <- nimbleCode({
  phi ~ dunif(0, 1) # prior survival
  p ~ dunif(0, 1) # prior detection
  # parameters
  gamma[1,1] <- phi # Pr(alive t -> alive t+1)
  gamma[1,2] <- 1 - phi # Pr(alive t -> dead t+1)
  gamma[2,1] <- 0 # Pr(dead t -> alive t+1)
  gamma[2,2] <- 1 # Pr(dead t -> dead t+1)
  delta[1] <- 1 # Pr(alive t = 1) = 1
  delta[2] <- 0 # Pr(dead t = 1) = 0
  omega[1,1] <- 1 - p # Pr(alive t -> non-detected t)
  omega[1,2] <- p # Pr(alive t -> detected t)
  omega[2,1] <- 1 # Pr(dead t -> non-detected t)
  omega[2,2] <- 0 # Pr(dead t -> detected t)
  # likelihood
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2])
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})
first <- apply(y, 1, function(x) min(which(x !=0)))
my.constants <- list(N = nrow(y), T = 55, first = first)
my.constants
my.data <- list(y = y + 1)
zinits <- y + 1 # non-detection -> alive
zinits[zinits == 2] <- 1 # dead -> alive
initial.values <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
parameters.to.save <- c("phi", "p")
parameters.to.save
n.iter <- 5000
n.burnin <- 1000
n.chains <- 2
mcmc.output <- nimbleMCMC(code = hmm.survival,
                          constants = my.constants,
                          data = my.data,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains)




#run model by sex
#run model by known age

