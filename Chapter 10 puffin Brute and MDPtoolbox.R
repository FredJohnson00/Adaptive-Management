setwd('C:\\Users\\fjohn\\OneDrive\\Documents\\PROJECTS\\AM Primer Johnson-Nichols\\Chapter 10 solve MDPs')
rm(list=ls())

library(MDPtoolbox)  # package for MDPtoolbox for R by Chadès et al. 2014. Ecography 37:916-920.
library(matrixcalc)  # package for matrix calculations
library(popbio)      # package for eigen analysis

# Brute force method for finding optimal policy ####################################################
#...define transition matrices and reward vectors
p = 0.3 # probability of recovery from a collapsed state in the absence of harvest
#p=0.6   # alternative probability of recovery from a collapsed state in the absence of harvest
(A1 = matrix(c(1-p,0,0,p,0.5,0.4,0,0.5,0.6),byrow=TRUE,nrow=3))  # no-harvest transitions
(R1 = matrix(c(0,0,0),ncol=1))                                   # no-harvest rewards
(A2 = matrix(c(1,0.6,0,0,0.4,0.7,0,0,0.3),byrow=TRUE,nrow=3))    # harvest transitions
(R2 = matrix(c(-99,1,1),ncol=1))                                 # harvest rewards

#...compile transitions for all possible 2^3 = 8 policies by pulling appropriate columns from matrices A1 and A2
pndx = expand.grid(1:2,1:2,1:2)
(pndx = as.matrix(pndx[,3:1]))
A = array(c(A1,A2),dim=c(3,3,2))
pols = array(dim=c(3,3,8))
for (i in 1:8) pols[,,i] = matrix(c(A[,1,pndx[i,1]],A[,2,pndx[i,2]],A[,3,pndx[i,3]]),nrow=3)
pols

#...compute longterm distributions of states for each potential policy
LTD = array(dim=c(3,1,8))
for (i in 1:8) LTD[,,i] = matrix.power(pols[,,i],200)[,1]
LTD[,1,1]
eigen.analysis((pols[,,1]))$stable.stage # gives equivalent results

#...compute immediate rewards for each potential policy by pulling appropriate rows of vectors R1 and R2
R = array(c(R1,R2),dim=c(3,1,2))
returns = array(dim=c(3,1,8))
for (i in 1:8) returns[,,i] = c(R[1,,pndx[i,1]],R[2,,pndx[i,2]],R[3,,pndx[i,3]])
returns

#...compute long-term expected rewards for each potential policy
Ereturns = NULL
for (i in 1:8) Ereturns[i] = LTD[,,i]%*%returns[,,i]

#...bundle results
results = matrix(nrow=8,ncol=4)
for (i in 1:8) results[i,] = c(LTD[,,i],Ereturns[i])
results = round(results,3)
(results = cbind(pndx,results))

#...which strategy is optimal?
results[which.max(results[,7]),] 

# Solve for optimal policy using dynamic programming
#...with backward iteration algorithm
#...time horizon = 10 and discount factor = 1
P = array(c(t(A1),t(A2)),dim=c(3,3,2))               # ROW-stochastic transitions
R = cbind(R1,R2)                                     # combine action-specific returns
mdp_check(P,R)                                       # check P and R
mdp = mdp_finite_horizon(P, R, 1, 10)                # solve the MDP
mdp$policy                                           

#...with average value algorithm
#...time horizon = infinite and discount factor = 1
(mdp = mdp_relative_value_iteration(P,R))

# Passive adaptive management ##############################################
w = seq(0,1,0.02) # probabilities of the null model                 
p = 0.3 # null model
A1 = matrix(c(1-p,0,0,p,0.5,0.4,0,0.5,0.6),byrow=TRUE,nrow=3)  # col stochastic
NullMod = t(A1) # row stochastic
p = 0.6 # alternative model
A1 = matrix(c(1-p,0,0,p,0.5,0.4,0,0.5,0.6),byrow=TRUE,nrow=3)  # col stochastic
AltMod = t(A1) # row stochastic

polAM = VAM = policyAM = list()
for (i in 1:length(w))
{ 
  P = array(data=NA,dim=c(3,3,2))
  P[,,1] = (NullMod*w[i]+AltMod*(1-w[i]))   # calculate average model for each belief weight, w
  P[,,2] = t(A2)
  mdp_check(P,R)
  results = mdp_finite_horizon(P, R, 1, 200)
  VAM[[i]] <- results$V[,1]
  policyAM[[i]] <- results$policy[,1]
}
#policyAM

ppi=400
dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
x = 1:3
y = unique(w)
z = matrix(nrow=length(y),ncol=length(x))
for (i in 1:length(w)) z[i,] = policyAM[[i]]
image(x=x,y=y,z=t(z),xlab="State",ylab='P(NULL Model)',main='',col=gray(c(0.9,0.1)),
      cex.lab=1.25,cex.axis=1.25,
      cex.main=2,xaxt='n',las=1);box()
axis(1,at=1:3,labels=c('Collapsed','Vulnerable','Robust'))
abline(h=seq(0,1,0.1),v=1:3,col='lightgray',lty=3)
legend('topleft',legend=c('No harvest','Harvest'),fill=gray(c(0.9,0.1)),bty='n',
       border='white',text.col='black',cex=1.2)

# EVPI
max = c(VAM[[1]][2],VAM[[51]][2])
for (i in 1:length(w)) avg.max = (1-w)*max[1] + w*max[2]
max.avg = NULL
for (i in 1:length(w)) max.avg[i] = VAM[[i]][2]
EVPI = avg.max - max.avg
#... search for maximum EVPI
maxEVPI = max(EVPI)
maxEVPI.ndx = which(EVPI==max(EVPI))
x=cbind(w[maxEVPI.ndx], max.avg[maxEVPI.ndx],EVPI[maxEVPI.ndx])
(maxEVPI.improvement = x[,3]/x[,2])

dev.new(height=6*ppi,width=9*ppi,noRStudioGD = TRUE)
plot(w,avg.max,las=1,type='n',ylab='Cumulative return',xlab='P(NULL model)',lwd=2,ylim=c(max[2],max[1]))
grid()
lines(w,avg.max,lwd=2,lty=2)
lines(w,max.avg,lwd=2)
text(0.75,90,'avg(max)')
text(0.3,90,'max(avg)')
lines(c(w[maxEVPI.ndx],w[maxEVPI.ndx]),c(max.avg[maxEVPI.ndx],avg.max[maxEVPI.ndx]),lwd=3,lty=2)
text(0.60,85,paste('max(EVPI) = ',round(maxEVPI,1),sep=''),pos=4)

# Simulations of passive adaptive management ####################################
#... # Bayes theorem for 2 models
bayes = function(prior,l1,l2) (prior*l1)/(prior*l1+(1-prior)*l2) 
action = c('N','H') # action labels

# select 'true' model
model = 1 # null model
M = array(dim=c(3,3,2))
M[,,1] = if(model==0) NullMod else AltMod
M[,,2] = t(A2)
P0 = array(dim=c(3,3,2))
P0[,,1] = NullMod
P0[,,2] = t(A2)
P1 = array(dim=c(3,3,2))
P1[,,1] = AltMod
P1[,,2] = t(A2)

T = 15                      # timeframe for simulations
sim = matrix(nrow=3,ncol=T) # container for simulation results
sim[1,1] = 3                # initialize system state
sim[2,1] = 0.5              # initialize weight

for (j in 2:T) {            # determine AM policy, choose action, update system state and belief weights
  P[,,1] = (NullMod*sim[2,j-1]+AltMod*(1-sim[2,j-1]))                   # calculate average model 
  results = mdp_finite_horizon(P, R, 1, 2000)               # calculate AM policy
  policyAM = results$policy[,1]                             # store AM policy
  sim[3,j-1] = policyAM[sim[1,j-1]]                         # implement optimal action
  sim[1,j] = sample(1:3,1,prob=M[sim[1,j-1],,sim[3,j-1]])   # update state
  sim[2,j] = round(bayes(sim[2,j-1],P0[sim[1,j-1],sim[1,j], # update belief weight
                                       sim[3,j-1]],P1[sim[1,j-1],sim[1,j],sim[3,j-1]]),4)
}

#... display results: row 1 = states, row 2 = weight; row 3 = actions
plot(1:T,sim[2,],xlab='Time',ylab='P(Null model)',main='Truth = alternative model',
     las=1,type='l',
     lwd=3,ylim=c(0,1));grid()
text(1:T,0,action[sim[3,]])
text(1:T,1,sim[1,])

ppi = 600
tiff(file='alt model sim1.tif',res=ppi,width = 9*ppi, height = 6*ppi)
plot(1:T,sim[2,],xlab='Time',ylab='P(Null model)',main='Truth = null model',
     las=1,type='l',
     lwd=3,ylim=c(0,1));grid()
text(1:T,0,action[sim[3,]])
text(1:T,1,sim[1,])
dev.off()
