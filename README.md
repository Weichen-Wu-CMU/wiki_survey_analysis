# README

This repository contains utility files to analyze wiki-survey data with pairwise comparison models like the Bradley-Terry model and the hierarchical Thurstone model, and perform model diagnostics and visualize the results. It also contains code that can be used to reproduce the results and figures in the paper *Diagnostics for Pairwise Comparison Models*, by Wu et al (2021).



## Utility Files

The following files contain utility functions that can be used to process wiki-survey data and analyze pairwise comparison data using the Bradley-Terry model and the hierarchical Thurstone model. 



### wiki_utils.R

This file contains functions to process the wiki-survey data and prepare for analysis. This file depends on the R package `ggplot2`. It has the following functions:



<u>get.active.ideas</u>: get the active ideas (mostly for internal use);

<u>rename.ideas</u>: add "o." at the beginning of the IDs of the ideas (mostly for internal use);

<u>change.type.ideas</u>: change the type of the text of the ideas to characters (mostly for internal use);

<u>threshold.ideas</u>: get the ideas that win and lose in the comparisons for at least once (or another threshold that can be decided by the users, mostly for internal use) ;

<u>preprocess.ideas</u>: preprocess the ideas by running the four functions above (for external use);

<u>rename.votes</u>: add "s." at the beginning of the Session IDs and "o." at the beginning of the idea IDs of each vote (mostly for internal use);

<u>get.valid.votes</u>: get the valid votes that were cast between two active ideas (mostly for internal use);

<u>give.vote.results</u>: add a boolean attribute called "results" to the data frame representing the votes. Results are set to 1 if and only if the idea shown on the left won the comparison (mostly for internal use);

<u>preprocess.votes</u>: preprocess the votes by running the three function above (for external use);

<u>make.cont.matrix</u>: make the contingency matrix $W$ for the Bradley-Terry model. $W[i,j]$ represents the number of times that idea $i$ beats idea $j$ in the comparisons (for external use);

<u>get.theta.v.indexes</u>: provide indexes for $\boldsymbol{\theta}_v$ in the hierarchical Thurstone model. The vector $\boldsymbol{\theta}_v$ represents all the elements in the opinion matrix that are estimated by votes. For the detailed description of this concept, refer to *Wiki surveys: Open and quantifiable social data collection*, by Salganik and Levy, 2012 (mostly for internal use);

<u>make.stan.data</u>: the wrapper function that provides the data for STAN estimation of the hierarchical Thurstone model (for external use).



### bt_utils.R

This file contains the functions to analyze pairwise comparison data using the Bradley-Terry model and perform model diagnostics. This file depends on the following R packages: `ggplot2, BradleyTerry2, latex2exp, Rmisc, ComQuadForm`. It has the following functions:



<u>sigmoid</u>: the sigmoid function, $\sigma(x) = 1/(1+\exp(-x))$ (mostly for internal use);

<u>qgenchisq</u>: find the quantile of a generalized $\chi^2$ distribution (mostly for internal use);

<u>fit.bt.model</u>: fit the Bradley-Terry model under the sum constraint (for external use);

<u>bt.obj.diagnostics</u>: generate (and show) the object diagnostic plots of the Bradley-Terry model, including the box plot and normal Q-Q plot of the object residuals, the scatterplot between object residuals and fitted Bradley-Terry scores, and the scatterplot between object residuals and the number of appearances (for details, please refer to Wu et al, 2021). The users can choose which plots to generate and whether to show them (for external use);

<u>bt.sbj.diagnostics</u>: generate (and show) the subject diagnostic plots of the Bradley-Terry model, including the difference plot, the influence plot and the ranking graph of the number of coparisons made by each subject (for details, please refer to Wu et al, 2021). The users can choose which plots to generate and whether to show them (for external use);

<u>bt.visualize</u>: visualize the results of the Bradley-Terry model. The users can choose to generate and show the ranking graph, probability heatmap and z-score heatmap of the top objects (for details, please refer to Wu et al, 2021; for external use);

<u>bt.dic</u>: calculate the deviance information criterion (DIC) of the Bradley-Terry model (for external use).



### hierarchical_utils.R

This file contains the functions to analyze multi-subject pairwise comparison data with the hierarchical Thurstone model proposed by Salganik & Levy (2012). It also includes functions to perform model diagnostics. This file depends on the R package `rstan` , and the installation of STAN. It has the following functions:



<u>fit.hierarchical.model</u>: fit the hierarchical Thurstone model with STAN (for external use);

<u>get.hierarchical.scores</u>: get the object scores in the hierarchical Thurstone model (for external use);

<u>eval.hierarchical.deviance</u>: calculate the deviance of the hierarchical Thurstone model (mostly for internal use);

<u>eval.hierarchical.gof</u>: evaluate the goodness-of-fit of the hierarchical Thurstone model (for external use).



### hierarchical_model.stan

This is a STAN file that specifies the hierarchical Thurstone model proposed by Salganik & Levy (2012). This file is called by the function <u>fit.hierarchical.model</u> in the file **hierarchical_utils.R**, and mostly for internal use.



## Files to replicate the results in Wu et al (2021)

The following files contain code to replicate the results in Wu et al (2021). 



### simulation.R

This file contains code for experiments 1 and 2 in the paper. Specifically, it can be used to produce:



<u>Figure 1 on Page 13 of the paper</u>, which shows the normal Q-Q plot of the object residuals in Example 1;

<u>FIgure 2 on Page 13 of the paper</u>, which shows the relationship between object residuals and the fitted Bradley-Terry scores in Example 2;

<u>Figure 3 on Page 14 of the paper</u>, which gives an explanation to the relationship shown in Figure 2.



### planyc_analysis.R

This file contains the code to analyze the PlaNYC wiki-survey data using the framework discribed by Wu et al (2021). Specifically, it can be used to preprocess the data, fit the Bradley-Terry and hierarchical Thurstone models, perform model diagnostics accordingly, and replicate the following results:

<u>Figure 4 on page 22 of the paper</u>, which illustrates the distribution of number of votes cast by each voter;

<u>Figure 5 on page 22 of the paper</u>, which illustrates the distribution of the number of times each idea appeared in comparisons;

<u>Figure 6 on page 23 of the paper</u>, a ranking graph of the top-10 ideas in the wiki-survey found by the Bradley-Terry model;

<u>Figure 7 on page 24 of the paper</u>, a z-score heatmap of the difference in scores among the top-10 ideas in the wiki-survey found by the Bradley-Terry model;

<u>Figure 8 on page 25 of the paper</u>, the idea (object) diagnostic plots for the Bradley-Terry model fitted on the PlaNYC wiki-survey;

<u>Figure 9 on page 27 of the paper</u>, the voter (subject) diagnostic plots for the Bradley-Terry model fitted on the PlaNYC wiki-survey;

<u>Table 1 on page 31 of the paper</u>, which compares the goodness-of-fit among the Bradley-Terry model, the Bayesian Bradley-Terry model, and the Hierarchical Thurstone model;

<u>Figure 13 on page 31 of the paper</u>, which shows the relationship between the deviance for real and simulated data in a posterior predictive check of the hierarchical Thurstone model.



### appendix.R

This file contains code to produce <u>Figures 14 and 15 on page 35 in the appendix of the paper</u>, which illustrates the relationship among the number of appearances ($N_i$) , the Bradley-Terry scores $\hat{\beta}_i$, and their estimated variances $\widehat{\text{Var}}(\hat{\beta}_i)$. 



### bootstrap.R

This file contains code to perform the parameter bootstrap experiments on Section 5.2 of the paper.  **This file needs further cleaning at this point.** Specifically, it can be used to reproduce:



<u>Figure 10 on page 28 of the paper</u>, the voter (subject) diagnostic plots for Experiment 1;

<u>Figure 11 on page 29 of the paper</u>, the voter (subject) diagnostic plots for Experiment 2;

<u>Figure 12 on page 30 of the paper</u>, the voter (subject) diagnostic plots for Experiment 3.