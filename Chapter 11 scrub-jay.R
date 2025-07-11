# Scrub-Jay.R

#... Scrub-jay habitat management, with state-action transitions given by:
#... Williams, B. K., and F. A. Johnson. 2018. Value of sample information in dynamic, 
#... structurally uncertain resource systems. PLOS ONE 13:e0199326.

library(MDPtoolbox)
library(matrixcalc)

# Function for normalizing returns to [0,1]
norm01 <- function(x){(x-min(x))/(max(x)-min(x))}

#-------------------------------------------------------------
#                       SECTION A                             
#-------------------------------------------------------------
# State & actions
scrub=c('short-open','short-closed','optimal-open','optimal-closed','tall')
action=c('nothing','routine','intensive')

#-------------------------------------------------------------
# State-action returns
# Demographic performance for each scrub class is characterized by the expected
#...number of yearlings produced per breeding pair minus the number of breeders 
#...per pair that die during the annual cycle. Relative to the do-nothing alternative 
#...(action 1), routine and intensive burns (actions 2 and 3, respectively) were 
#...arbitrarily assigned a 10% and 40% reduction in demographic performance, 
#...respectively, to account for the costs of burning. 
(return = matrix(c(-0.310,-0.310,0.490,0.150,-0.240,
                   -0.341,-0.341,0.441,0.135,-0.264,
                   -0.434,-0.434,0.294,0.090,-0.336),nrow=5))
(R = norm01(return))  # normalize returns

#-------------------------------------------------------------
# State-action transitions (row stochastic)

#... do-nothing
(n = t(matrix(c(0.276,0.310,0.146,0.259,0.009,0.000,0.484,0.043,0.473,0.000,0.000,0.028,0.448,0.476,0.048,
                0.000,0.011,0.000,0.822,0.167,0.000,0.003,0.002,0.002,0.993),nrow=5)))

#... routine burn
(r = t(matrix(c(0.346,0.346,0.192,0.116,0.000,0.203,0.322,0.051,0.407,0.017,0.156,0.000,0.558,0.169,0.117,
                0.150,0.010,0.163,0.461,0.216,0.045,0.002,0.012,0.006,0.935),nrow=5)))

#... intensive burn
(i = t(matrix(c(0.438,0.156,0.187,0.219,0.000,0.438,0.156,0.187,0.219,0.000,0.273,0.000,0.568,0.091,0.068,
                0.272,0.006,0.296,0.260,0.166, 0.098,0.005,0.025,0.010,0.862),nrow=5)))

#-------------------------------------------------------------
# Solve the MDP for Model 1
P1 = array(data=NA,dim=c(5,5,3))
P1[,,1] = n
P1[,,2] = r
P1[,,3] = i

mdp_check(P1,R)
results1 = mdp_finite_horizon(P1, R, 1, 2000)
V1 <- results1$V[,1]
policy1 <- results1$policy[,1]
pol1 = data.frame(cbind(scrub,action[policy1],round(V1,2)))
colnames(pol1)=c('shrub','action','value')
pol1
results1$policy[,1:20]

#-------------------------------------------------------------
# Solve the MDP for a null Model 0, where intensive burn only as effective as routine burn
P0 = array(data=NA,dim=c(5,5,3))
P0[,,1] = n
P0[,,2] = r
P0[,,3] = r 

mdp_check(P0,R)
results0 = mdp_finite_horizon(P0, R, 1, 2000)
V0 <- results0$V[,1]
policy0 <- results0$policy[,1]
pol0 = data.frame(cbind(scrub,action[policy0],round(V0,2)))
colnames(pol0)=c('shrub','action','value')
pol0

#--------------------------------------------------------------
# Long-run characteristics of doing nothing and of the optimal, model-specific mgmt policies
(LRstates.nomgmt = matrix.power(n,100)[1,])
(LRdemo.nomgmt = sum(LRstates.nomgmt*return[,1]))

# policy0
P0opt = matrix(nrow=5,ncol=5)
for (i in 1:5) P0opt[i,] = P0[i,,policy0[i]]
(LRstates.P0opt = matrix.power(P0opt,100)[1,])
(LRdemo.P0opt = sum(LRstates.P0opt*return[,1]))

# policy1
P1opt = matrix(nrow=5,ncol=5)
for (i in 1:5) P1opt[i,] = P1[i,,policy1[i]]
(LRstates.P1opt = matrix.power(P1opt,100)[1,])
(LRdemo.P1opt = sum(LRstates.P1opt*return[,1]))

#-------------------------------------------------------------
#                       SECTION B                             
#-------------------------------------------------------------
# Passive adaptive management
w = seq(0,1,0.05) # belief weights for Model 0                    
# w = 0.5

polAM = VAM = policyAM = list()
for (i in 1:length(w))
{ 
P = array(data=NA,dim=c(5,5,3))
P[,,1] = n
P[,,2] = r
P[,,3] = (P0[,,3]*w[i]+P1[,,3]*(1-w[i]))   # calculate average model for each belief weight, w

mdp_check(P,R)
results = mdp_finite_horizon(P, R, 1, 2000)
VAM[[i]] <- results$V[,1]
policyAM[[i]] <- results$policy[,1]
polAM[[i]] = data.frame(cbind(scrub,action[policyAM[[i]]],round(VAM[[i]],2)))
colnames(polAM[[i]])=c('scrub','action','value')
}
polAM

#...graph AM policy
scrub.lab = c('SO','SC','OO','OC','T')
x = 1:5
y = unique(w)
z = matrix(nrow=length(y),ncol=length(x))
for (i in 1:length(w)) z[i,] = policyAM[[i]]

ppi=400
# jpeg(file='C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\final figures\\Johnson_figure11.1.jpg',
#      res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,4,2))
image(x=x,y=y,z=t(z),xlab="State",ylab='P(Model 0)',main='',col=gray(c(0.9,0.5,0.1)),cex.lab=1.5,cex.axis=1.25,
      cex.main=2,xaxt='n');box()
axis(1,at=1:5,labels=scrub.lab)
abline(h=seq(0,1,0.1),v=1:5,col='lightgray',lty=3)
legend('bottomleft',legend=c('do nothing','routine burn','intensive burn'),fill=gray(c(0.9,0.5,0.1)),bty='n',
       border='white',text.col='black',cex=1.2)
#dev.off()

#-------------------------------------------------------------
#                       SECTION C                             
#-------------------------------------------------------------
# EVPI
# avgmax = matrix(nrow=length(w),ncol=5)
# for (i in 1:length(w)) avgmax[i,] = VAM[[1]]*(1-w[i])+VAM[[length(w)]]*w[i]  # expected performance under certainty
# 
# maxavg = matrix(nrow=length(w),ncol=5)
# for (i in 1:length(w)) maxavg[i,] = VAM[[i]]                                 # expected performance under uncertainty
# 
# EVPI = cbind(w,avgmax-maxavg)
# colnames(EVPI) = c('w',scrub.lab)
# EVPI
# 
# plot(EVPI[,1],EVPI[,6],type='n',las=1,xlab='P(Model 0)',ylab='EVPI(Tall)')
# abline(v=seq(0,1,.1),h=seq(0,40,5),col='lightgray',lty=3)
# lines(EVPI[,1],EVPI[,6],lwd=3)

#-------------------------------------------------------------
#                       SECTION D                             
#-------------------------------------------------------------
# Simulations of adaptive management
bayes = function(prior,l1,l2) (prior*l1)/(prior*l1+(1-prior)*l2) # Bayes theorem for 2 models
action = c('N','R','I')

model = 0                   # select 'true' model
M = if(model==0) P0 else P1
T = 25                      # timeframe for simulations
sim = matrix(nrow=3,ncol=T) # container for simulation results
sim[1,1] = 1                # initialize system state
sim[2,1] = 0.5              # initialize weight

#(seed = sample(1:100,1))
seed = 69
set.seed(seed)
for (j in 2:T) {            # determine AM policy, choose action, update system state and belief weights
  P[,,3] = (P0[,,3]*sim[2,j-1]+P1[,,3]*(1-sim[2,j-1]))      # calculate average model
  results = mdp_finite_horizon(P, R, 1, 2000)               # calculate AM policy
  policyAM = results$policy[,1]                             # store AM policy
  sim[3,j-1] = policyAM[sim[1,j-1]]                         # implement optimal action
  sim[1,j] = sample(1:5,1,prob=M[sim[1,j-1],,sim[3,j-1]])   # update state
  sim[2,j] = round(bayes(sim[2,j-1],P0[sim[1,j-1],sim[1,j], # update belief weight
            sim[3,j-1]],P1[sim[1,j-1],sim[1,j],sim[3,j-1]]),4)
}
# sim                                                     # display results: row 1 = states, row 2 = weight; row 3 = actions
(sum.demo.perf = sum(return[sim[1,]]))                    # calculate cumulative real demographic performance

ppi=400
# jpeg(file='C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\final figures\\Johnson_figure11.2.jpg',
#      res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,4,2))
plot(1:T,sim[2,],xlab='Time',ylab='P(model 0)',main=paste("Truth = Model",model),type='l',
     sub=paste('Demo =',sum.demo.perf),lwd=3,ylim=c(0,1));grid()
text(1:T,0,action[sim[3,]])
text(1:T,1,sim[1,])
#dev.off()

###########################
model = 1                   # select 'true' model
M = if(model==0) P0 else P1
T = 25                      # timeframe for simulations
sim = matrix(nrow=3,ncol=T) # container for simulation results
sim[1,1] = 1                # initialize system state
sim[2,1] = 0.5              # initialize weight

seed = 77
set.seed(seed)
for (j in 2:T) {            # determine AM policy, choose action, update system state and belief weights
  P[,,3] = (P0[,,3]*sim[2,j-1]+P1[,,3]*(1-sim[2,j-1]))      # calculate average model
  results = mdp_finite_horizon(P, R, 1, 2000)               # calculate AM policy
  policyAM = results$policy[,1]                             # store AM policy
  sim[3,j-1] = policyAM[sim[1,j-1]]                         # implement optimal action
  sim[1,j] = sample(1:5,1,prob=M[sim[1,j-1],,sim[3,j-1]])   # update state
  sim[2,j] = round(bayes(sim[2,j-1],P0[sim[1,j-1],sim[1,j], # update belief weight
                                       sim[3,j-1]],P1[sim[1,j-1],sim[1,j],sim[3,j-1]]),4)
}
# sim                                                     # display results: row 1 = states, row 2 = weight; row 3 = actions
(sum.demo.perf = sum(return[sim[1,]]))                    # calculate cumulative demographic performance

ppi=400
 # jpeg(file='C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\final figures\\Johnson_figure11.3.jpg',
 #      res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,4,2))
plot(1:T,sim[2,],xlab='Time',ylab='P(model 0)',main=paste("Truth = Model",model),type='l',
     sub=paste('Demo =',sum.demo.perf),lwd=3,ylim=c(0,1));grid()
text(1:T,0,action[sim[3,]])
text(1:T,1,sim[1,])
#dev.off()

#####################################
seed = 69
set.seed(seed)
for (j in 2:T) {            # determine AM policy, choose action, update system state and belief weights
  P[,,3] = (P0[,,3]*sim[2,j-1]+P1[,,3]*(1-sim[2,j-1]))      # calculate average model
  results = mdp_finite_horizon(P, R, 1, 2000)               # calculate AM policy
  policyAM = results$policy[,1]                             # store AM policy
  sim[3,j-1] = policyAM[sim[1,j-1]]                         # implement optimal action
  sim[1,j] = sample(1:5,1,prob=M[sim[1,j-1],,sim[3,j-1]])   # update state
  sim[2,j] = round(bayes(sim[2,j-1],P0[sim[1,j-1],sim[1,j], # update belief weight
                                       sim[3,j-1]],P1[sim[1,j-1],sim[1,j],sim[3,j-1]]),4)
}
# sim                                                     # display results: row 1 = states, row 2 = weight; row 3 = actions
(sum.demo.perf = sum(return[sim[1,]]))                    # calculate cumulative demographic performance

ppi=400
# jpeg(file='C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\final figures\\Johnson_figure11.4.jpg',
#      res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,4,2))
plot(1:T,sim[2,],xlab='Time',ylab='P(model 0)',main=paste("Truth = Model",model),type='l',
     sub=paste('Demo =',sum.demo.perf),lwd=3,ylim=c(0,1));grid()
text(1:T,0,action[sim[3,]])
text(1:T,1,sim[1,])
#dev.off()

