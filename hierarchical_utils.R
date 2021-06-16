## This file includes functions to estimate the hierarchical Thurstone
## model, and perform diagnostics proposed by the paper
## Diagnostics for Pairwise Comparison Models,
## by Wu et al (2021).

suppressPackageStartupMessages(library(rstan))
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

fit.hierarchical.model <- function(stan.data,
                                   iter = 20000,
                                   thin = 10)
# Function to fit the hierarchical Thurstone model using STAN
# Inputs:
# stan.data: list that can be used to fit a STAN file
#            output of make.stan.data in wiki_utils.R
# iter: integer, number of iterations in the Markov Chain
# thin: integer, thinning parameter. Samples will be stored per 
#       thin iterations
# Output:
# stan.results: list that records the samples in the Markov Chain.
{
  stan.fit = stan(file = 'hierarchical_model.stan',
                  data = stan.data, 
                  iter = iter,
                  thin = thin)
  stan.results = extract(stan.fit)
  return(stan.results)
}

get.hierarchical.scores <- function(stan.results,obj.names)
# Function to get the object scores in the hierarchical Thurstone model
# Inputs:
# stan.results: the recorded samples in the Markov Chain
#               during STAN estimation
#               output of fit.hierarchical.model
# obj.names: names of the objects
# Output:
# hierarchical.scores: named vector, the scores of the objects
{
  hierarchical.scores = colMeans(stan.results$mu)
  names(hierarchical.scores) = obj.names[order(obj.names)]
  return(hierarchical.scores)
}

eval.hierarchical.deviance <- function(stan.data,theta.v)
# Function to calculate the deviance of the hierarchical Thurstone model
# (Mostly for internal use)
# Inputs:
# stan.data: list that can be used to fit a STAN file
#            output of make.stan.data in wiki_utils.R;
# theta.v: theta_v in the hierarchical model
# Output:
# deviance, the deviance of the model
{
  deviance = 0
  for(i in 1:stan.data$V)
  {
    p = pnorm(theta.v[stan.data$left.indices[i]] - 
                theta.v[stan.data$right.indices[i]])
    if(stan.data$y[i] == 0) p = 1 - p
    deviance = deviance - 2 * log(p)
  }
  return(deviance)
}

eval.hierarchical.gof <- function(stan.data,stan.results)
# Function to evaluate the goodness-of-fit of the hierarchical 
# Thurstone model
# Inputs:
# stan.data: list that can be used to fit a STAN file
#            output of make.stan.data in wiki_utils.R;
# stan.results: the recorded samples in the Markov Chain
#               during STAN estimation,
#               output of fit.hierarchical.model
# Output:
# gof, a list containing:
#      deviance, the deviance of the model;
#      eff.params: effective number of parameters;
#      DIC: the deviance information criterion (DIC);
{
  bar.D.theta = ave(stan.results$deviance)[1]
  theta.v.est = colMeans(stan.results$theta_v)
  D.bar.theta = eval.hierarchical.deviance(stan.data,theta.v)
  pD = bar.D.theta - D.bar.theta
  DIC = D.bar.theta + 2 * pD
  gof = list()
  gof$deviance = D.bar.theta
  gof$eff.params = pD
  gof$DIC = DIC
  return(gof)
}

ppc.hierarchical.model <- function(stan.results)
# Function to perform the posterior predictive check (PPC) of 
# the hierarchical Thurstone model
# Input:
# stan.results: the recorded samples in the Markov Chain
#               during STAN estimation,
#               output of fit.hierarchical.model
{
  deviances = stan.results$deviance
  deviance.rep = stan.results$deviance_rep
  ppc = list(deviances = deviances, deviance.rep = deviance.rep)
  ppc = as.data.frame(ppc)
  
  ggplot(ppc,aes(x = deviances,y=deviance.rep)) +
    geom_point() + geom_abline(slope = 1, intercept = 0) + 
    xlab(TeX('Deviance for real data, $D^{(k)}$')) +
    ylab(TeX('Deviance for simulated data, $D^{rep (k)}$'))
  
  p = sum(deviances > deviance.rep)/length(deviances)
  print(paste('The p-value of the posterior predictive test is',p))
}


