# This file contains utility functions for analyzing pairwise comparison
# data using the Bradley-Terry model (Bradley & Tery, 1952),
# presenting the results, and performing diagnostics
# according to the paper
# "Diagnostics for Pairwise Comparison Models",
# by Wu et al, 2021.

## Packages
suppressPackageStartupMessages(library(BradleyTerry2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(latex2exp))
suppressPackageStartupMessages(library(Rmisc))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(CompQuadForm))


sigmoid <- function(x)
# Sigmoid function: input x, output 1/(1+exp(-x))
# (Mostly for internal use)
{
  s = 1 / (1 + exp(-x))
  return(s)
}

qgenchisq <- function(p, lambda, max_iter = 1000, tol = 1e-4)
# Find the q-quantile of the generalized chi-square distribution
# (Depends on the CompQuadForm package)
# (Mostly for internal use)
{
  # Initialize
  min = 0
  cur_q = sum(lambda)
  max = 100 * cur_q
  cur_p = 1-davies(q = cur_q, lambda = lambda)$Qq
  for(num_iter in 1:max_iter)
  {
    #print(cur_q)
    if(abs(cur_p - p) <= tol)
    {
      print(paste('Quantile find after', as.character(num_iter),
                  'iterations', sep = ' '))
      break
    }
    if(cur_p - p > tol)
    {
      max = cur_q
      cur_q = (cur_q + min) / 2
      cur_p = 1-davies(q = cur_q, lambda = lambda)$Qq
    }
    if(cur_p - p < -tol)
    {
      min = cur_q
      cur_q = (cur_q + max) / 2
      cur_p = 1-davies(q = cur_q, lambda = lambda)$Qq
    }
  }
  if(num_iter == max_iter)
    print('Maximum iteration reached; result unreliable')
  return(cur_q)
}


fit.bt.model <- function(W)
## Fit the Bradley-Terry model under the sum constraint
# Input: W, named, n * n matrix, the contingency matrix
#        row and column names represents the objects
#        W[i,j]: number of times that object i beats object j
# Output: bt.scores, named vector with length n
# the fitted Bradley-Terry scores
{
    votes.to.btm = countsToBinomial(W)
    bt.model = BTm(outcome = cbind(win1,win2),
                   player1 = player1,
                   player2 = player2,
                   data = votes.to.btm)
    bt.scores = BTabilities(bt.model)
    bt.scores = as.data.frame(bt.scores)$ability
    bt.scores = bt.scores - ave(bt.scores)[1]
    names(bt.scores) = rownames(W)
    return(bt.scores)
}

bt.obj.diagnostics <- function(W, bt.scores, 
                               over.dispersion = TRUE,
                               r.boxplot = TRUE,
                               r.qqplot = TRUE,
                               r.vs.scores = TRUE,
                               r.vs.app = TRUE,
                               show.plots = TRUE)
# Perform object diagnostics for a Bradley-Terry model
# Inputs:
# W: named n*n matrix, the contingency matrix;
# bt.scores: named length n vector, the fitted Bradley-Terry scores
#            satisfying the sum constraint
# over.dispersion: boolean, whether to perform over-dispersion test
#                  and estimate a dispersion parameter
# r.boxplot: boolean, whether to construct the boxplot of object residuals
# r.qqplot: boolean, whether to construct the normal
#           Q-Q plot of the object residuals
# r.vs.scores: boolean, whether to construct the scatterplot
#              between object residuals and the Bradley-Terry scores
# r.vs.app: boolean, whether to construct the scatterplot
#           between object residuals and the number of appearances
# show.plots: boolean, whether to show the plots
# Output: diag.plots, a list, including
#         (1) r.boxplot, the boxplot of the object residuals;
#         (2) r.qqplot, the normal Q-Q plot of the object residuals;
#         (3) r.vs.scores, the scatterplot of the object residuals versus
#                          the Bradley-Terry scores;
#         (4) r.vs.app, the scatterplot of the object residuals versus
#                       the number of appearances
{
  # Calculate residuals
  n = nrow(W) # Number of objects
  V = W + t(W)
  A = matrix(data = 0, nrow = n, ncol = n)
  for(i in 1:n) A[,i] = bt.scores
  for(i in 1:n) A[i,] = A[i,] - bt.scores
  B = sigmoid(A)
  R = (W - V * B)/(sqrt((V+1e-8) * B * (1-B)))
  
  # Dispersion parameter
  phi = 1
  # Overdispersion test
  if(over.dispersion)
  {
    m = sum(V>0) / 2
    df = m - n + 1
    chi2 = sum(sum(R * R)) / 2
    p.value = 1 - pchisq(q = chi2, df = df)
    print(paste('Pearson X2 statistics is ', chi2, 
                ' with ', df, 'degrees of freedom.'))
    print(paste('p value for overdispersion test is ',
                p.value))
    if(p.value >= 0.05)
    {
      print('No significant evidence of overdispersion.')
    }
    if(p.value < 0.05)
    {
      phi = m / df
      print(paste('Siginificant evidence of overdispersion;',
                  ' overdispersion parameter estimate: ', phi))
    }
  }
  
  # Calculate object residuals
  r = R %*% bt.scores / sqrt(sum(bt.scores * bt.scores))
  r = r / phi
  
  obj.diag = matrix(data = 0, nrow = n, ncol = 3)
  colnames(obj.diag) = c('name','residual','score')
  obj.diag = as.data.frame(obj.diag)
  obj.diag$name = names(bt.scores)
  obj.diag$residual = r
  obj.diag$score = bt.scores
  obj.diag$num.app = rowSums(V)
  
  # Boxplot for object residuals
  p1 = ggplot(obj.diag,aes(x = factor(0), y=residual)) + 
    geom_boxplot() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ylab(TeX('residuals($r_i$)')) 
  
  # Normal Q-Q plot for object residuals
  p2 = ggplot(obj.diag,aes(sample = residual)) +
    stat_qq() + geom_abline(slope = 1, intercept = 0) + 
    xlab('Theoretical quantiles') + ylab(TeX('Residuals ($r_i$)')) 
  
  # Object residuals VS Bradley-Terry scores
  p3 = ggplot(obj.diag,
              aes(x = bt.scores, y = residual)) +
    geom_point() + 
    xlab(TeX('Bradley-Terry scores ($\\hat{\\beta}_i$)')) + 
    ylab(TeX('Residuals ($r_i$)')) + theme(legend.position = 'none')
  
  # Object resdiuals VS Number of appearances
  p4 = ggplot(obj.diag,
              aes(x = num.app, y = residual)) +
    geom_point() +
    xlab('Number of appearances') + ylab(TeX('Residuals ($r_i$)')) + 
    theme(legend.position = 'none')
  
  if(show.plots)
  {
    multiplot(p1,p3,p2,p4, cols = 2)
  }
  
  diag.plots = list()
  diag.plots$r.boxplot = NULL
  diag.plots$r.qqplot = NULL
  diag.plots$r.vs.scores = NULL
  diag.plots$r.vs.app = NULL
  if(r.boxplot) diag.plots$r.boxplot = p1
  if(r.qqplot) diag.plots$r.qqplot = p2
  if(r.vs.scores) diag.plots$r.vs.scores = p3
  if(r.vs.app) diag.plots$r.vs.app = p4
  return(diag.plots)
}

bt.sbj.diagnostics <- function(votes,
                               bt.scores,
                               diff.plot = TRUE,
                               infl.plot = TRUE,
                               rank.plot = TRUE,
                               show.plots = TRUE)
# Subject diagnostics of the Bradley-Terry model
# Inputs:
# votes: a data frame, each row representing a vote (comparison),
#        that includes at least the following columns:
#        Subject.ID: the subject that makes the comparison;
#        Winner.ID: the winner of the comparison;
#        Loser.ID: the loser of the comparison.
# bt.scores: named length n vector, the Bradley-Terry scores
# show.plots: boolean, whether to show the diagnostic plots
# Output: diag.plots, a list including
#         diff.plot, plot for "difference";
#         infl.plot, plot for "influence".
{
  total.votes = nrow(votes)
  sbj.list = unique(votes$Subject.ID)
  sbj.diag = matrix(nrow = length(sbj.list),ncol = 3)
  rownames(sbj.diag) = sbj.list
  colnames(sbj.diag) = c('num.votes','difference','influence')
  sbj.diag = as.data.frame(sbj.diag)
  
  obj.list = names(bt.scores)
  W = make.cont.matrix(votes,obj.list)
  n = nrow(W) # Number of objects
  V = W + t(W)
  A = matrix(data = 0, nrow = n, ncol = n)
  for(i in 1:21) A[,i] = bt.scores
  for(i in 1:21) A[i,] = A[i,] - bt.scores
  B = sigmoid(A)
  H.beta = V * B * (1-B)
  tmp = rowSums(H.beta)
  diag(H.beta) = -tmp
  #H.beta is the Hessian matrix of \beta
  H.gamma = H.beta[2:n,2:n]
  #H.gamma is the Hessian matrix of \gamma
  Var.gamma = solve(-H.gamma)
  #Var.gamma is the variance of \gamma
  P = matrix(data = 0, nrow = n, ncol = n-1) #\beta = P\gamma
  for(j in 1:(n-1)) P[j+1,j] = 1
  P = P - 1 / n
  Var.beta = P %*% Var.gamma %*% t(P)
  
  for(sbj in sbj.list)
  {
    delta.votes = votes[votes$Subject.ID == sbj,]
    delta.W = make.cont.matrix(delta.votes,obj.list)
    delta.V = delta.W + t(delta.W)
    delta.w = rowSums(delta.W)
    E.delta.w = rowSums(delta.V * B)
    diff = delta.w - E.delta.w
    infl = Var.beta %*% diff
    sbj.diag[sbj,'num.votes'] = nrow(delta.votes)
    sbj.diag[sbj,'difference'] = sum(diff**2)
    sbj.diag[sbj,'influence'] = sum(infl**2)
  }
  
  sbj.diag = sbj.diag[order(sbj.diag$num.votes,decreasing = T),]
  sbj.diag$rank = seq(1,nrow(sbj.diag))
  rank.plot = ggplot(sbj.diag,aes(x = rank, y = num.votes)) +
    geom_point() +
    xlab('Voters ranked by number of votes') + 
    ylab('Number of votes cast per voter') + 
    theme(axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))
  
  # Reference for difference
  lambdas = eigen(-H.beta,symmetric = TRUE, 
                  only.values = TRUE)$values[1:n-1]
  upper.dif = qgenchisq(p = 0.975, lambda = lambdas, max_iter = 100)
  lower.dif = qgenchisq(p = 0.025, lambda = lambdas, max_iter = 100)
  middle.dif = qgenchisq(p = 0.5, lambda = lambdas, max_iter = 100)
  
  # Plot for difference
  diff.plot = 
    ggplot(data = sbj.diag, aes(x = num.votes, y = difference))+
    geom_point() +
    geom_abline(intercept = 0, slope = lower.dif / total.votes,
                color = 'green') +
    geom_abline(intercept = 0, slope = middle.dif / total.votes,
                color = 'blue') +
    geom_abline(intercept = 0, slope = upper.dif / total.votes,
                color = 'red') +
    xlab(TeX('Number of votes $V^{(s)}$')) + 
    ylab(TeX('Difference $D^{(s)}$'))  +
    scale_color_discrete('')
  
  # Reference for influence
  mus = 1 / lambdas
  upper.inf = qgenchisq(p = 0.975, lambda = mus, max_iter = 100)
  lower.inf = qgenchisq(p = 0.025, lambda = mus, max_iter = 100)
  middle.inf = qgenchisq(p = 0.5, lambda = mus, max_iter = 100)
  
  # Plot for influence
  infl.plot = 
    ggplot(data = sbj.diag, aes(x = num.votes, y = influence))+
    geom_point() +
    geom_abline(intercept = 0, slope = lower.inf / total.votes,
                color = 'green') +
    geom_abline(intercept = 0, slope = middle.inf / total.votes,
                color = 'blue') +
    geom_abline(intercept = 0, slope = upper.inf / total.votes,
                color = 'red') +
    xlab(TeX('Number of votes $V^{(s)}$')) + 
    ylab(TeX('Influence $I^{(s)}$')) +
    scale_color_discrete('')
  
  diag.plots = list()
  diag.plots$diff.plot = NULL
  diag.plots$infl.plot = NULL
  diag.plots$rank.plot = NULL
  if(diff.plot) diag.plots$diff.plot = diff.plot
  if(infl.plot) diag.plots$infl.plot = infl.plot
  if(rank.plot) diag.plots$rank.plot = rank.plot
  if(show.plots)
  {
    multiplot(diff.plot, NULL, infl.plot, NULL, cols = 2)
  }
  return(diag.plots)
}

bt.visualize <- function(W, bt.scores, 
                         top = 10, 
                         obj.names = NULL,
                         rank.graph = TRUE,
                         prob.heatmap = TRUE,
                         zscore.heatmap = TRUE)
# Visualize the results of the Bradley-Terry model
# Inputs:
# W: named n*n matrix, the contingency matrix
# bt.scores: named length n vector, the Bradley-Terry scores
# top: integer, number of objects included in the visualization
# obj.names: list of characters, names of the objects to be shown
# rank.graph: boolean, whether to construct the ranking graph
# prob.heatmap: boolean, whether to construct the probability heatmap
# zscore.heatmap: boolean, whether to construct the z-score heatmap
{
  n = nrow(W) # Number of objects
  V = W + t(W)
  A = matrix(data = 0, nrow = n, ncol = n)
  for(i in 1:21) A[,i] = bt.scores
  for(i in 1:21) A[i,] = A[i,] - bt.scores
  B = sigmoid(A)
  H.beta = V * B * (1-B)
  tmp = rowSums(H.beta)
  diag(H.beta) = -tmp
  #H.beta is the Hessian matrix of \beta
  H.gamma = H.beta[2:n,2:n]
  #H.gamma is the Hessian matrix of \gamma
  Var.gamma = solve(-H.gamma)
  #Var.gamma is the variance of \gamma
  P = matrix(data = 0, nrow = n, ncol = n-1) #\beta = P\gamma
  for(j in 1:(n-1)) P[j+1,j] = 1
  P = P - 1 / n
  Var.beta = P %*% Var.gamma %*% t(P)
  
  obj.list = rownames(W)
  bt.results = matrix(nrow = length(obj.list),ncol = 4)
  rownames(bt.results) = obj.list
  colnames(bt.results) = c('est','lower','upper','names')
  bt.results = as.data.frame(bt.results)
  bt.results$est = bt.scores
  bt.results$lower = bt.scores - 2 * diag(Var.beta)
  bt.results$upper = bt.scores + 2 * diag(Var.beta)
  bt.results = bt.results[order(bt.results$est,
                                decreasing = T),]
  bt.results = bt.results[1:top,]
  bt.results$Summary = obj.names
  
  bt.visual = list()
  bt.visual$rank.graph = NULL
  bt.visual$prob.heatmap = NULL
  bt.visual$zscore.heatmap = NULL
  
  if(rank.graph)
  {
    p = ggplot(bt.results,
               aes(x = Summary, y = est)) 
    p = p + geom_pointrange(aes(ymin = lower, ymax = upper)) + 
      labs(x = 'Object', y = 'Score')
    p = p + coord_flip() + 
      scale_x_discrete(limits = 
                         bt.results$Summary[order(bt.results$est,
                                                  decreasing = F)])
    p = p + theme(axis.title = element_text(size = 10),
                  axis.text = element_text(size = 10),
                  plot.margin = unit(c(3,0.3,3,0.3),'cm')) 
    bt.visual$rank.graph = p
  }
  
  if(prob.heatmap)
  {
    B = B[order(bt.scores,decreasing = TRUE), order(bt.scores,decreasing = TRUE)]
    x = as.character(seq(1,top))
    y = paste(x,obj.names,sep = '.')
    heat.data = expand.grid(X=x,Y=y)
    heat.data$Probability = as.vector(t(B[1:top,1:top]))
    p = ggplot(heat.data, aes(X, Y, fill= Probability)) + 
      geom_tile() +
      scale_fill_gradient2(high="red", mid= 'white',low = 'navy') +
      theme(axis.text.x = element_text(angle=0, hjust=1, vjust=.5)) + 
      scale_x_discrete(limits = x[order(bt.results$est,decreasing = T)],
                       position = 'top') +
      scale_y_discrete(limits = 
                         y[order(bt.results$est,decreasing = F)]) +
      xlab('') + ylab('') + theme(plot.margin = unit(c(2,0,2,0),'cm'))
    bt.visual$prob.heatmap = p
  }
  
  if(zscore.heatmap)
  {
    V.beta = matrix(data = 0, nrow = n, ncol = n)
    for(i in 1:n) V.beta[,i] = V.beta[,i] + diag(Var.beta)
    for(i in 1:n) V.beta[i,] = V.beta[i,] + diag(Var.beta)
    V.beta = V.beta - 2 * Var.beta
    heat.matrix = A / sqrt(V.beta + 1e-8)
    
    zscore.matrix = heat.matrix[order(bt.scores,decreasing = T),
                                order(bt.scores,decreasing = T)]
    
    zscore.matrix = zscore.matrix[1:top,1:top]
    
    x = as.character(seq(1,top))
    y = paste(x,obj.names,sep = '.')
    heat.data = expand.grid(X=x,Y=y)
    heat.data$Z = as.vector(t(zscore.matrix))
    p = ggplot(heat.data, aes(X, Y, fill= Z)) + 
      geom_tile() +
      scale_fill_gradient2(high="red", mid= 'white',low = 'navy') +
      theme(axis.text.x = element_text(angle=0, hjust=1, vjust=.5)) + 
      scale_x_discrete(limits = x[order(bt.results$est,decreasing = T)],
                       position = 'top') +
      scale_y_discrete(limits = 
                         y[order(bt.results$est,decreasing = F)]) +
      xlab('') + ylab('') + theme(plot.margin = unit(c(2,0,2,0),'cm'))
    bt.visual$zscore.heatmap = p           
  }
  return(bt.visual)
}

eval.bt.gof <- function(W, bt.scores)
# Evaluate the goodness-of-fit of the Bradley-Terry model
# Inputs:
# W: named n*n matrix, the contingency matrix
# bt.scores: named length n vector, the Bradley-Terry scores
# Output:
# gof, a list showing the goodness-of-fit of the model, including
#      deviance, the deviance of the model;
#      num.params, the number of parameters;
#      AIC, the Akaike Information Criterion of the model.
{
  n = nrow(W)
  V = W + t(W)
  A = matrix(nrow = n, ncol = n)
  for(i in 1:n) A[,i] = bt.scores
  for(i in 1:n) A[i,] = A[i,] - bt.scores
  B = log(sigmoid(A))
  D = -2 * sum(sum(W * B))
  deviance = D - sum(sum(log(choose(V,W)))) # deviance
  p = n - 1
  AIC = deviance + 2 * (n-1)
  gof = list()
  gof$deviance = deviance
  gof$num.params = p
  gof$AIC = AIC
  return(gof)
}


