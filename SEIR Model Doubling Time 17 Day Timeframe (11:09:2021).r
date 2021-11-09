# SEIR Model

#Packages
library(rvest)
library(dplyr)
library(readr)
library(optimx)


#State Population Dataset
popdfpage <- read_html("https://en.wikipedia.org/wiki/List_of_U.S._states_and_territories_by_population")
popdf <- html_node(popdfpage, ".wikitable")
popdf <- html_table(popdf, fill = T) 
popdf <- popdf[-1,]
popdf <- popdf[,3:4]
popdf <- popdf %>% select(c("state" = "State or territory","population"="Census population[7][a]"))
popdf <- popdf[-4,]
popdf$population=parse_number(popdf$population)
initials=c("CA","TX","FL","PA","IL","OH","GA","NC", "MI", "NJ", "VA", "WA", "AZ", "MA", "TN", "IN", "MD", "MO", "WI", "CO", "MN", "SC", "AL", "LA", "KY", "OR", "OK", "CT", "PR", "UT", "IA", "NV", "AR", "MS", "KS", "NM", "NE", "ID", "WV", "HI", "NH", "ME", "RI", "MT", "DE", "SD", "ND", "AK", "DC", "VT", "WY", "GU", "VI", "MP", "AS")
initials=tolower(initials)
popdf=cbind(initials,popdf[1:55,2])
popdf %>% add_row(initials="nyc",population=8804190) %>%
  add_row(initials="ny",population=11397059) %>%
  add_row(initials="mh", population=59194) %>%
  add_row(initials="pw", population=18092)


#State Data Table with Population, Cases and Deaths
state_data=function(dir_input,state_input,state_population)
{
  library(dplyr)
  print(paste(dir_input,state_input,"_daily_cases.csv",sep=""))
  state_daily_cases <- read.table(paste(dir_input,state_input,"_daily_cases.csv",sep=""), sep = ",", header = T)
  state_daily_cases <- state_daily_cases[,1:3]
  state_daily_deaths <- read.table(paste(dir_input,state_input,"_daily_deaths.csv",sep=""), sep = ",", header = T)
  state_daily_deaths <- state_daily_deaths[,1:3]
  state_cum_cases <- read.table(paste(dir_input,state_input,"_cum_cases.csv",sep=""), sep = ",", header = T)
  state_cum_deaths <- read.table(paste(dir_input,state_input,"_cum_deaths.csv",sep=""), sep = ",", header = T)
  state_dates <- rev(seq(as.Date("2020/01/23"), by="day", length.out=nrow(state_cum_cases)))
  state_population <- rep(state_population, length(state_dates))
  statedf <- data.frame(state_daily_cases, state_daily_deaths, state_cum_cases, state_cum_deaths, state_dates, state_population)
  statedf <- statedf %>% select("date"= state_dates, 
                        "population" = state_population, 
                        "cases" = New.Cases, 
                        "cumulative_cases" = Total.Cases, 
                        "deaths" = New.Deaths, 
                        "cumulative_deaths" = Total.Deaths) %>%
  arrange(date)
  statedf$day <- as.Date(statedf$date)-as.Date("2020-02-01")
  (statedf <- statedf[statedf$day >= 0,])
}


## Fitting
SEIR_fit <- function(x,data,begin) {
  end=as.Date(begin)+16
  data.sub=data[(which(as.Date(data$date)==as.Date(begin)):which(as.Date(data$date)==as.Date(end))),]
  a=x[1]
  beta=x[2]
  delta=x[3]
  epsilon=x[4]
  gamma=x[5]
  rho=x[6]
  alpha=x[7]
  U=x[8]
  q1=x[9]
  q2=x[10]
  q3=x[11]
  q4=x[12]
  Q=0
  I=0
  D=0
  E=0
  R=0
  S=1-Q-U-I-D-E-R
  time=begin
  I_cum=0
  D_cum=0
  counter=1
  N=40000000
  ls1=0
  ls2=0
  ls3=0
  if (S>1|Q>1|U>1|I>1|D>1|E>1|R>1)
  {
    ls1=1e20
    ls2=1e20
  }
  if (S+Q+U+I+D+E+R>1)
  {
    ls1=1e20
    ls2=1e20
  }
  if (S<0|Q<0|U<0|I<0|D<0|E<0|R<0)
  {
    ls1=1e20
    ls2=1e20
  }
  if (a<=0.98|beta<0|delta<0|epsilon<0|gamma<0|rho<0|alpha<0|q1<0|q2>0|q3<0|q4>0)
  {
    ls1=1e20
    ls2=1e20
  }
  if (a>1.02|beta>50|delta>1|epsilon>1|gamma>1|rho>1|alpha>1|q1>1|abs(q2)>1|q3>1|abs(q4)>1)
  {
    ls1=1e20
    ls2=1e20
  }
  c=matrix(nrow=(as.Date(end)-as.Date(begin)+1),ncol=5)
  comp1=c()
  comp2=c()
  while (time<=end)
  {
    if (time==begin)
    {
      U=x[8]
    }
    if (time==as.Date("2020-03-19"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- (-1*(beta*S*(U^a)))-(q1*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q1+epsilon+delta))
      dQ <- q1*(U+S)
      Q=max(Q+dQ,0)
      I=max(I+dI,0)
      U=max(U+dU,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-(Q+U+I+D+E+R)
    }
    if (time==as.Date("2020-05-18"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- (-1*(beta*S*(U^a)))-(q2*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q2+epsilon+delta))
      dQ <- q2*(U+S)
      Q=max(Q+dQ,0)
      I=max(I+dI,0)
      U=max(U+dU,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-(Q+U+I+D+E+R)
    }
    if (time==as.Date("2020-11-19"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- (-1*(beta*S*(U^a)))-(q3*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q3+epsilon+delta))
      dQ <- q3*(U+S)
      Q=max(Q+dQ,0)
      I=max(I+dI,0)
      U=max(U+dU,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-(Q+U+I+D+E+R)
    }
    if (time==as.Date("2021-01-25"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- (-1*(beta*S*(U^a)))-(q4*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q4+epsilon+delta))
      dQ <- q4*(U+S)
      Q=max(Q+dQ,0)
      I=max(I+dI,0)
      U=max(U+dU,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-(Q+U+I+D+E+R)
    }
    else
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- (-1*(beta*S*(U^a)))+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(epsilon+delta))
      dQ <- 0
      Q=max(Q+dQ,0)
      I=max(I+dI,0)
      U=max(U+dU,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-(Q+U+I+D+E+R)
    }
    I_cum=I_cum+I/14
    D_cum=D_cum+dD
    compare=which(as.Date(data.sub$date)==time)
    c[counter,1]=format(time, format="%Y-%m-%d")
    c[counter,2]=I_cum*N
    c[counter,3]=data.sub$cumulative_cases[compare]
    c[counter,4]=D_cum*N
    c[counter,5]=data.sub$cumulative_deaths[compare]
    #print(list(time,I_cum*N))
    #print(c[counter,3])
    comp1[counter]=(as.numeric(c[counter,2])-as.numeric(c[counter,3]))^2
    comp2[counter]=(as.numeric(c[counter,4])-as.numeric(c[counter,5]))^2
    if (time==as.Date(end))
    {
      S_final=S
      Q_final=Q
      U_final=U
      I_final=I
      D_final=D
      E_final=E
      R_final=R
      I_pt=I/14
    }
    time=as.Date(time)+1
    counter=counter+1
  }
  ls1=ls1+sum(comp1)
  ls2=ls2+sum(comp2)
  ld=data.sub$cumulative_deaths[length(data.sub$cumulative_deaths)]
  lc=data.sub$cumulative_cases[length(data.sub$cumulative_cases)]
  ls=ls1+((lc+ld)/ld)*ls2
  #print(x)
  #print(ls)
  return(list(ls,I_cum*N,D_cum*N,S_final,Q_final,U_final,I_final,D_final,E_final,R_final,I_pt*N))
}

## Forecasting
SEIR_forecast <- function(x,data,begin,end,I_cum_input,D_cum_input,Input_List) {
  data.sub=data[(which(as.Date(data$date)==as.Date(begin)):which(as.Date(data$date)==as.Date(end))),]
  a=x[1]
  beta=x[2]
  delta=x[3]
  epsilon=x[4]
  gamma=x[5]
  rho=x[6]
  alpha=x[7]
  q1=x[9]
  q2=x[10]
  q3=x[11]
  q4=x[12]
  S=Input_List[[1]]
  Q=Input_List[[2]]
  U=Input_List[[3]]
  I=Input_List[[4]]
  D=Input_List[[5]]
  E=Input_List[[6]]
  R=Input_List[[7]]
  N=40000000
  I_cum=I_cum_input/N
  D_cum=D_cum_input/N
  time=begin
  counter=1
  c=matrix(nrow=(as.Date(end)-as.Date(begin)+1),ncol=4)
  while (time<=end)
  {
    if (time==as.Date("2020-03-19"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- -1*(beta*S*(U^a))-(q1*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q1+epsilon+delta))
      dQ <- q1*(U+S)
      Q=max(Q+dQ,0)
      U=max(U+dU,0)
      I=max(I+dI,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-Q-U-I-D-E-R
    }
    if (time==as.Date("2020-05-18"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- -1*(beta*S*(U^a))-(q2*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q2+epsilon+delta))
      dQ <- q2*(U+S)
      Q=max(Q+dQ,0)
      U=max(U+dU,0)
      I=max(I+dI,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-Q-U-I-D-E-R
    }
    if (time==as.Date("2020-11-19"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- -1*(beta*S*(U^a))-(q3*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q3+epsilon+delta))
      dQ <- q3*(U+S)
      Q=max(Q+dQ,0)
      U=max(U+dU,0)
      I=max(I+dI,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-Q-U-I-D-E-R
    }
    if (time==as.Date("2021-01-25"))
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- -1*(beta*S*(U^a))-(q4*S)+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(q4+epsilon+delta))
      dQ <- q4*(U+S)
      Q=max(Q+dQ,0)
      U=max(U+dU,0)
      I=max(I+dI,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-Q-U-I-D-E-R
    }
    else
    {
      dI <- delta*U-(I*(gamma+alpha))
      dR <- alpha*I-rho*R
      dD <- gamma*I
      dE <- epsilon*U-rho*E 
      dS <- -1*(beta*S*(U^a))+(rho*(E+R))
      dU <- beta*S*(U^a)-(U*(epsilon+delta))
      dQ <- 0
      Q=max(Q+dQ,0)
      U=max(U+dU,0)
      I=max(I+dI,0)
      D=max(D+dD,0)
      E=max(E+dE,0)
      R=max(R+dR,0)
      S=1-Q-U-I-D-E-R
    }
    
    I_cum=I_cum+I/14
    D_cum=D_cum+dD
    compare=which(as.Date(data.sub$date)==time)
    c[counter,1]=format(time, format="%Y-%m-%d")
    c[counter,2]=I_cum*N
    c[counter,3]=D_cum*N
    c[counter,4]=(I/14)*N
    time=as.Date(time)+1
    counter=counter+1
  }
  return(c)
}

## Fit 10 different 17 day periods
SEIR_results=function(dataset,begin_input,i_max_input,end_forecast)
{
  i_max=i_max_input
  fitted_param=matrix(NA,nrow=i_max,ncol=14)
  forecasted=list()
  for (i in 1:i_max)
  {
    fit=optimx(par=c(1,0.7,0.5,0.1,0.02,0.01,0.3,1e-5,0.1,-0.5,0.1,-0.5),
             function(x,data,begin) SEIR_fit(x,dataset,as.Date(begin_input)+(17*(i-1)))[[1]], 
             method = c("Nelder-Mead"), control = list(maxit = 20000))
  
  pars=c(fit$p1,fit$p2,fit$p3,fit$p4,fit$p5,fit$p6,fit$p7,fit$p8,fit$p9,fit$p10,fit$p11,fit$p12)
  
  I_cum=SEIR_fit(pars,dataset,as.Date(begin_input)+(17*(i-1)))[[2]]
  
  D_cum=SEIR_fit(pars,dataset,as.Date(begin_input)+(17*(i-1)))[[3]]
  
  Input_List=SEIR_fit(pars,dataset,as.Date(begin_input)+(17*(i-1)))[4:10]
  
  I_pt=SEIR_fit(pars,dataset,as.Date(begin_input)+(17*(i-1)))[[11]]
  
  fitted_param[i,]=c(as.numeric(pars),as.numeric(I_cum),as.numeric(D_cum)) 

  forecasted[[i]] = SEIR_forecast(x=fitted_param[i,1:12], dataset, as.Date(begin_input)+(17*(i))[[1]], 
                                  as.Date(end_forecast), I_cum, D_cum, Input_List)
  print(i)
  }
  return(list(fitted_param,forecasted))
}

s=state_data("/Users/sashafarzin-nia/Desktop/UCLA/Year 2/Thesis/coviddata/","al",popdf[popdf$initials=="al",2])
results=SEIR_results(s,"2021-02-01",3, "2021-05-01")
results

#Fitted parameters include in order: a, beta, delta, epsilon, gamma, rho, alpha, U, q1, q2, q3, q4.
#q1 and q2 are spring 2020 quarantine start and end whereas q3 and q4 are fall 2020 quarantine start and end

#Output includes forecast date, total infections, total deaths, daily infection rate

