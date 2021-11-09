set.seed(85) 
library(readxl)
library(dplyr) 

state.data <- function(dir_input,state_input){
  print(paste(dir_input,state_input,"_daily_cases.csv",sep=""))
  state_daily_cases <- read.table(paste(dir_input,state_input,"_daily_cases.csv",sep=""), sep = ",", header = T)
  state_daily_cases <- state_daily_cases[,1:3]
  state_cum_cases <- read.table(paste(dir_input,state_input,"_cum_cases.csv",sep=""), sep = ",", header = T)
  state_dates <- rev(seq(as.Date("2020/01/23"), by="day", length.out=nrow(state_cum_cases)))
  statedf <- data.frame(state_daily_cases, state_cum_cases, state_dates)
  statedf <- statedf %>% select("date"= state_dates, 
                                "cases" = New.Cases, 
                                "cumulative_cases" = Total.Cases) %>%
    arrange(date)
  statedf$day <- as.Date(statedf$date)-as.Date("2020-02-01")
  (statedf <- statedf[statedf$day >= 0,])
}

drw = 0 
f = function(p){ ## This is my function to compute the rms error for a given parameter value. 
    k <<- k+1
    cat(k," ") 
    if(min(p)<.0000001) return(9e+20) 
    if(sum(p[38:52]) > 1) return(9e+20) 
    expectedn = rep(0,days)  
    mu = p[1] 
    K = p[2:37]
    g = c(p[38:52],1-sum(p[38:52])) 
    for(i in 1:days){
	if(i>16) tmp = t3[i-c(1:16)]
	if(i<17) tmp = c(t3[(i-1):1],rep(0,(17-i))) 
	if(i<2) tmp = rep(0,16) 
	ind = ceiling(i/16) 
	expectedn[i] = mu + K[ind] * sum(g[1:16]*tmp) 
	}
    if(drw>1) lines(1:16,g) 
    ans = sqrt(mean((expectedn - t3[1:days])^2)) 
    cat(ans,"\n") 
    ans 
} 

hawkesN.run=function(dir,state)
{
  statedf.hawkes <- state.data(dir, state)
  t1 <<- statedf.hawkes
  t2 <<- t1[[3]] 
  t2 = as.numeric(gsub(",","",t2)) ## need this for Alabama and California 
  t3 <<- as.numeric(t2[(length(t2)-1):1]) 
  days <<- 36*16 ## I'll just use this part of the data since I want 16 day segments. 
  k <<- 1 
  theta1 = c(rep(1,37),rep(1/16,15)) 
  drw = 0
  b1 = optim(theta1,f,method="Nelder-Mead",control=list(maxit=100000)) 
  b2 = optim(b1$par,f,method="Nelder-Mead",control=list(maxit=100000))   
  return(b2$par)
}
results=hawkesN.run("/Users/sashafarzin-nia/Desktop/UCLA/Year 2/Thesis/coviddata/", "ca")
results

mu=results[1]
kappa=results[2:37]
g=results[38:52]
lc.g=1-sum(g)
g=c(g,lc.g)
g <- data.frame(g)
nl.result=nls(g~beta*exp(-beta*(1:16)), data = g, start=list(beta=0.01))
beta=summary(nl.result)$par[[1]]

### Forecasting
hawkesN.covid.forecasting=function(background,beta,K,susceptible,sim_num,T_end)
{
  N=susceptible
  result.matrix <- matrix(NA, nrow = sim_num, ncol = T_end)
  total=c()
  for (p in 1:sim_num)
  {
    mu=background
    num.background = rpois(1, mu*T_end)
    bgpoints = sort(runif(num.background)*T_end)
    acc_pt=bgpoints
    sj_I=acc_pt
    length(acc_pt)
    n=1
    while (n<2e5)
    {
      n=length(acc_pt)
      cand_pt=c()
      for (k in 1:length(sj_I))
      {
        r_0=K
        gamma=beta
        if (length(sj_I)>10)
        {
          nsup = rpois(1, r_0)
        }
        if (length(sj_I)<=10)
        {
          nsup = max(rpois(1, r_0),1)
        }
        if (nsup>0)
        {
          newpt=rexp(nsup,rate=mu)+sj_I[k]
          cand_pt=c(cand_pt,newpt)
        }
        if (nsup==0)
        {
          cand_pt=cand_pt
        }
      }
      sim_pt=sort(c(acc_pt,cand_pt))
      acc_pt=sort(acc_pt)
      cand_pt=sort(cand_pt)
      n2=length(sim_pt)
      sl_E=c()
      if (length(cand_pt)==0)
      {
        break
      }
      if (length(cand_pt)>0)
      {
        acc_counter=0
        for (i in 1:length(cand_pt))
        {
          N=susceptible
          j1 = max((1:n)[acc_pt<cand_pt[i]])
          lambda=(1-(j1/N))*sum(r_0*gamma*exp(-gamma*(cand_pt[i]-acc_pt[1:j1])))
          j2 = max((1:n2)[sim_pt<cand_pt[i]])
          nu=sum(r_0*gamma*exp(-gamma*(cand_pt[i]-sim_pt[1:j2])))
          a=lambda/nu
          keep=runif(1)<a
          if (keep==1 & cand_pt[i]<T_end)
          {
            acc_pt=sort(c(acc_pt, cand_pt[i]))
            n=n+1
            sl_E=sort(c(sl_E, cand_pt[i]))
            if (n%%1000==0)
            {
              acc_pt2=acc_pt[acc_pt>0]
              acc_trunc=trunc(acc_pt2,0)
              acc_trunc_ct=as.numeric(table(acc_trunc))
              print(n)
            }
          }
        }
        sj_I=sort(rexp(length(sl_E),rate=gamma)+sl_E)
      }
    }
    acc_pt=acc_pt[acc_pt>0]
    acc_trunc=trunc(acc_pt,0)
    acc_trunc_ct=as.numeric(table(acc_trunc))
    print(acc_trunc_ct)
    result.matrix[p,] <- acc_trunc_ct
  }
  med.vec <- apply(result.matrix, 2, median)
  return(med.vec)
}

forecast_results=matrix(0,length(kappa),30)
for (i in 1:2)
{
  forecast_results[i,]=hawkesN.covid.forecasting(background=mu*5, beta=beta, K=kappa[i], susceptible=1000000, sim_num=3,T_end=30)
}