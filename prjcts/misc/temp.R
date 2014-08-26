require(rcolgem)
parms_truth <- list(gamma0 = 1
 , gamma1 = 1/7
 , gamma2 = 1/2
 , mu = 1/30
 , b=.036
 , beta0 = 12./10
 , beta1=3./100
 , beta2=9./100
 , S_0=3000
 , I0_0=1, I1_0=0.01, I2_0=0.01
 , m=3, mm=1) 

INFECTEDNAMES <- c('I0', 'I1', 'I2')

births <- rbind( 
  c('parms$beta0 * S * I0 / (S + I0 + I1 + I2)', '0', '0'), 
  c('parms$beta1 * S * I1 / (S + I0 + I1 + I2)', '0', '0'), 
  c('parms$beta2 * S * I2 / (S + I0 + I1 + I2)', '0', '0')
)
rownames(births)=colnames(births)<- INFECTEDNAMES
migrations <- rbind(
  c('0', 'parms$gamma0 * I0', '0'),
  c('0', '0', 'parms$gamma1 * I1'),
  c('0', '0', '0')
)
rownames(migrations)=colnames(migrations) <- INFECTEDNAMES
deaths <- c(
 'parms$mu*I0'
 , 'parms$mu*I1'
 , 'parms$mu*I2 + parms$gamma2 * I2'
)
names(deaths) <- INFECTEDNAMES
nonDemeDynamics <- paste(sep='',
'-parms$mu*S + parms$mu*(S + I0 + I1 + I2)'
 ,'- S * (parms$beta0*I0+parms$beta1*I1+parms$beta2*I2) / (S + I0 + I1 + I2)'
)
names(nonDemeDynamics) <- 'S'
n <- 100
sampleTimes <- seq(40, 50, length.out=n) 
sampleStates <- c( rep(1, 10) , rep(2, 70), rep(3, 20))
sampleStates <- sample(sampleStates, n , replace=FALSE) #shuffle elements
names(sampleTimes) <- 1:n
print(system.time( bdt <- simulate.binary.dated.tree(births, deaths, nonDemeDynamics
  ,  t0=0
  , x0=c(I0=1, I1=.01, I2=.01, S = parms_truth$S_0)
  , sampleTimes=sampleTimes
  , sampleStates=sampleStates
  ,  migrations=migrations
  ,  parms=parms_truth
  , fgyResolution = 2000
  , integrationMethod = 'rk4')
))
