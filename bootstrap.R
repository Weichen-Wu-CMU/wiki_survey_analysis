rownames(B) = act.ideas.list
colnames(B) = act.ideas.list
total.votes = nrow(votes)
# Every voter is faced with the same choices
sim.voter.diag = matrix(nrow = length(voters.list),
                        ncol = 3)
rownames(sim.voter.diag) = voters.list
colnames(sim.voter.diag) = c('num.votes','difference','influence')
sim.voter.diag = as.data.frame(sim.voter.diag)
for(voter in voters.list)
{
  delta.votes = votes[votes$Session.ID == voter,]
  delta.W = matrix(data = 0, nrow = n, ncol = n)
  rownames(delta.W) = act.ideas.list
  colnames(delta.W) = act.ideas.list
  for(i in 1:nrow(delta.votes))
  {
      left = delta.votes$Left.Choice.ID[i]
      right = delta.votes$Right.Choice.ID[i]
      x = runif(1)
      if(x < B[left,right]) # left wins
        delta.W[left,right] = delta.W[left,right] + 1
      else
        delta.W[right,left] = delta.W[right,left] + 1
  }
  delta.V = delta.W + t(delta.W)
  delta.w = rowSums(delta.W)
  E.delta.w = rowSums(delta.V * B)
  diff = delta.w - E.delta.w
  infl = Var.beta %*% diff
  sim.voter.diag[voter,'num.votes'] = nrow(delta.votes)
  sim.voter.diag[voter,'difference'] = sum(diff**2)
  sim.voter.diag[voter,'influence'] = sum(infl**2)
}

diff_bootstrap3 = 
  ggplot(data = sim.voter.diag, aes(x = num.votes, y = difference))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.dif / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.dif / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.dif / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Difference $D^{(s)}$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
diff_bootstrap3



inf_bootstrap3 = 
  ggplot(data = sim.voter.diag, aes(x = num.votes, y = influence))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.inf / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.inf / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.inf / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Influence $I^{(s)}$'))+ 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
inf_bootstrap3

# Every voter faces new choices
sim.voter.diag2 = matrix(nrow = length(voters.list), ncol = 3)
colnames(sim.voter.diag2) = c('num.votes','difference','influence')
rownames(sim.voter.diag2) = voters.list
sim.voter.diag2 = as.data.frame(sim.voter.diag2)

for(voter in voters.list)
{
  delta.votes = votes[votes$Session.ID == voter,]
  num.votes = nrow(delta.votes)
  sim.voter.diag2[voter,'num.votes'] = num.votes
  indices = sample(x = total.votes, size = num.votes)
  delta.votes = votes[indices,]
  delta.W = matrix(data = 0, nrow = n, ncol = n)
  rownames(delta.W) = act.ideas.list
  colnames(delta.W) = act.ideas.list
  for(i in 1:nrow(delta.votes))
  {
    left = delta.votes$Left.Choice.ID[i]
    right = delta.votes$Right.Choice.ID[i]
    x = runif(1)
    if(x < B[left,right]) # left wins
      delta.W[left,right] = delta.W[left,right] + 1
    else
      delta.W[right,left] = delta.W[right,left] + 1
  }
  delta.V = delta.W + t(delta.W)
  delta.w = rowSums(delta.W)
  E.delta.w = rowSums(delta.V * B)
  diff = delta.w - E.delta.w
  infl = Var.beta %*% diff
  
  sim.voter.diag2[voter,'difference'] = sum(diff**2)
  sim.voter.diag2[voter,'influence'] = sum(infl**2)
}

diff_bootstrap2 = 
  ggplot(data = sim.voter.diag2, aes(x = num.votes, y = difference))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.dif / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.dif / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.dif / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Difference $D^{(s)}$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
diff_bootstrap2



inf_bootstrap2 = 
  ggplot(data = sim.voter.diag2, aes(x = num.votes, y = influence))+
  geom_point() +
  geom_abline(intercept = 0, slope = lower.inf / total.votes,
              color = 'green') +
  geom_abline(intercept = 0, slope = middle.inf / total.votes,
              color = 'blue') +
  geom_abline(intercept = 0, slope = upper.inf / total.votes,
              color = 'red') +
  xlab(TeX('Number of votes $V^{(s)}$')) + 
  ylab(TeX('Influence $I^{(s)}$')) + 
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  scale_color_discrete('')
inf_bootstrap2
