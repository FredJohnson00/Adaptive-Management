setwd('C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\Chapter 13 learning')

# Single-loop learning in coin tosses

p1=0.25 # candidate model 1 (probability of heads)
p2=0.50 # candidate model 2 (probability of heads)
p3=0.75 # candidate model 3 (probability of heads)
p = p2  # true model

T=25 # timeframe (# of coin flips)
n=1  # number of tosses each time step


b=matrix(nrow=3,ncol=T+1) # container for posterior probs
b[,1]=rep(1/3,3)          # priors for the three models

set.seed(1) # to ensure repeatability
for (t in 1:T)
 {
  toss=rbinom(n,1,p)               # coin tosses
  x=sum(toss)                      # number of heads
  l1=p1^x*(1-p1)^(n-x)             # likelihood under model 1
  l2=p2^x*(1-p2)^(n-x)             # likelihood under model 2
  l3=p3^x*(1-p3)^(n-x)             # likelihood under model 3
  norm.constant = (l1*b[1,t] + l2*b[2,t] + l3*b[3,t])
  b[1,t+1]=l1*b[1,t]/norm.constant # posterior for model 1
  b[2,t+1]=l2*b[2,t]/norm.constant # posterior for model 2
  b[3,t+1]=l3*b[3,t]/norm.constant # posterior for model 2
 }

# plot posterior model probabilities
# ppi=400
# tiff(file='cointoss3models1toss.tif',res=ppi,width = 9*ppi, height = 6*ppi)
# par(mar=c(5,4,2,2))
ppi = 400
# jpeg(file='C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\final figures\\Johnson_figure13.3.jpg',
#      res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,2,2) )
plot(1:(T+1),b[1,],type='n',ylab='Information state (model probability)',xlab='Time',las=1,ylim=c(0,1),lwd=3)
grid()
lines(1:(T+1),b[1,],lwd=3,lty=1)
lines(1:(T+1),b[2,],lwd=3,lty=2)
lines(1:(T+1),b[3,],lwd=3,lty=3)
legend('topleft',lwd=3,lty=1:3,bty='n',legend=c('Model 1','Model 2','Model 3'))
#dev.off()


# Single-loop learning in 2 coin tosses

p1=0.25 # candidate model 1 (probability of heads)
p2=0.50 # candidate model 2 (probability of heads)
p3=0.75 # candidate model 3 (probability of heads)
p = p2  # true model

T=25 # timeframe (# of coin flips)
n=2  # number of tosses each time step


b=matrix(nrow=3,ncol=T+1) # container for posterior probs
b[,1]=rep(1/3,3)          # priors for the three models

set.seed(1) # to ensure repeatability
for (t in 1:T)
{
  toss=rbinom(n,1,p)               # coin tosses
  x=sum(toss)                      # number of heads
  l1=p1^x*(1-p1)^(n-x)             # likelihood under model 1
  l2=p2^x*(1-p2)^(n-x)             # likelihood under model 2
  l3=p3^x*(1-p3)^(n-x)             # likelihood under model 3
  norm.constant = (l1*b[1,t] + l2*b[2,t] + l3*b[3,t])
  b[1,t+1]=l1*b[1,t]/norm.constant # posterior for model 1
  b[2,t+1]=l2*b[2,t]/norm.constant # posterior for model 2
  b[3,t+1]=l3*b[3,t]/norm.constant # posterior for model 2
}

# plot posterior model probabilities
# ppi=400
# tiff(file='cointoss3models1toss.tif',res=ppi,width = 9*ppi, height = 6*ppi)
# par(mar=c(5,4,2,2))
ppi = 400
#tiff(file='Fig12.3.tif',res=ppi,width = 9*ppi, height = 6*ppi)
# jpeg(file='C:\\Users\\fjohn\\OneDrive - Aarhus universitet\\U-Drive\\Documents\\AM Primer Johnson-Nichols\\final figures\\Johnson_figure13.4.jpg',
#      res=ppi,width = 9*ppi, height = 6*ppi)
par(mar=c(5,5,2,2) )
plot(1:(T+1),b[1,],type='n',ylab='Information state (model probability)',xlab='Time',las=1,ylim=c(0,1),lwd=3)
grid()
lines(1:(T+1),b[1,],lwd=3,lty=1)
lines(1:(T+1),b[2,],lwd=3,lty=2)
lines(1:(T+1),b[3,],lwd=3,lty=3)
legend('topleft',lwd=3,lty=1:3,bty='n',legend=c('Model 1','Model 2','Model 3'))
#dev.off()



