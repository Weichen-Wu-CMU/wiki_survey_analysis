# This file contains code to analyze the PlaNYC wiki-survey
# using the framework proposed by Wu et al (2021).

# Set working directory to files pane location
# setwd("~/Desktop/ADA/code")

# Import utility files
source("wiki_utils.R")
source("bt_utils.R")
source("hierarchical_utils.R")

prefix <- 'planyc'

if (prefix == "planyc") {
  question.text <- "Which do you think is a better idea for creating a greener, greater New York City?"
  cex.main.val <- 0.8
  question.text.2line <- "Which do you think is a better idea for creating\na greener, greater New York City?"
  wiki.survey.num <- 608
  votes.filename <- paste("data/wikisurvey_608_votes","_cleaned.csv", sep="")
  ideas.filename <- paste("data/wikisurvey_608_ideas", "_cleaned.csv", sep="")
  nonvotes.filename <- paste("data_precleaned/wikisurvey_608_nonvotes", "_cleaned.csv", sep="")
  start.date <- as.Date("2010-10-07")
  stop.date <- as.Date("2011-01-30")
}

## Load data
votes <- read.csv(votes.filename, header=TRUE, sep = ",", dec=".")
ideas <- read.csv(ideas.filename, header=TRUE, sep=",", dec=".")

## Preprocessing
ideas = preprocess.ideas(ideas,thres = 0)
votes = preprocess.votes(votes, ideas)

# Distribution of the number of times each idea
# appeared in comparisons
# (Fig.5 on page 22 of the paper)
show.num.app(ideas)

#Bradley-Terry model fit
W = make.cont.matrix(votes,ideas)
bt.scores = fit.bt.model(W)
print('Model fitted.')

# Summary of the top-10 ideas
Summary = c('Ban fracking',
            'Invest in transportation',
            'Plug ships into electricity grids',
            'Enhance bike lane network',
            'Year-long greenmarkets',
            'Local food in public schools',
            'Protected bike paths',
            'Transit service outside Manhattan',
            'Support community gardens',
            'Upgrade building energy efficiency')

# Visualize the results of the Bradley-Terry model
bt.visual = bt.visualize(W,bt.scores, obj.names = Summary)
# Ranking graph of the top-10 ideas
# (Fig.6 on page 23 of the paper)
bt.visual$rank.graph

# Z-score heatmap of the top-10 ideas
# (Fig.7 on page 24 of the paper)
bt.visual$zscore.heatmap

# Object(idea) diagnostics
# (Fig.8 on page 25 of the paper)
obj.diag.plots = bt.obj.diagnostics(W, bt.scores)

# Subject (voter) diagnostics
# (Fig.9 on page 27 of the paper)
sbj.diag.plots = bt.sbj.diagnostics(votes, bt.scores)

# Distribution of the number of votes cast by each voter
# (Fig.4 on page 22 of the paper)
sbj.diag.plots$rank.plot

# Evaluate the goodness-of-fit of the Bradley-Terry model
# (Line 1 of Table 1 on page 31 of the paper)
bt.gof = eval.bt.gof(W,bt.scores)
print(paste("Deviance:",bt.gof$deviance))
print(paste("Effective number of parameters:",bt.gof$num.params))
print(paste("AIC (DIC):",bt.gof$AIC))

# Prepare data for the hierarchical model
stan.data = make.stan.data(votes)
# Fit the hierarchical Thurstone model
stan.results = fit.hierarchical.model(stan.data)
# Evaluate the goodness-of-fit of the hierarchical Thurstone model
# (Line 3 of Table 1 on page 31 of the paper)
gof = eval.hierarchical.gof(stan.data,stan.results)
print(paste("Deviance:",gof$deviance))
print(paste("Effective number of parameters:",gof$eff.params))
print(paste("DIC:",gof$DIC))

# Posterior predictive check of the hierarchical Thurstone model
# (Fig.13 on page 31 of the paper)
ppc.hierarchical.model(stan.results)






