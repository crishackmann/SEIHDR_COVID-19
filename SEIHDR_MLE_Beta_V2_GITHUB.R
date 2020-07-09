# EDO SEIHDR com beta estimado por MLE 

require(deSolve)
require(ggplot2)
require(reshape2)
library(tidyverse)
library(lubridate)
options(stringsAsFactors = FALSE)

number_reported = c(52,75,86,85,85,97,88,93,94,100,104,97,102,110,102,105,108,108,110,109,125,129,133,127,135,136,139,165,158,173,179,175,187,184,206,212,215,214,241,247,274,275,287,311,313,330,367,335,348,370,398,428)
time = c(19: (18 + length(number_reported)))
reported_data = data.frame(time, number_reported)

#### Dates
date_ini = dmy("25/Abr/2020") 
num_days = 200 # mais o primeiro dia
series_days = date_ini + days(0:num_days)
times  = seq(0, num_days, by = 1)

num_days_hospi = length(number_reported)
date_ini_hospi = dmy("14/May/2020") 
###
# INPUT
sigma  = 1/5.1

gammah = 1/5.72 
gammad = 0 
gammar = 1/5.6

etad   = 1/15.56 
etar   = 1/10.93 

Theta  = 0.01 
Lambda = 0.154 

N = 1014009
E = 2 * 1318
I = 1318
H = 52
S = N - I - H - E
D = 0
R = 0
initial_state_values = c(S = S, E = E, I = I, H = H, D = D, R = R)

times = seq(from = 0, to = num_days, by = 1)

# SEIDR 
seidr_model = function(times, state, parameters) {  
  
  with(as.list(c(state, parameters)), {
    
    
    dS = - (beta * I) * (S / N)
    
    dE =   (beta * I) * (S / N) - (sigma * E)
    
    dI =   (sigma * E) - (Theta * gammah * I) - 
      ((1 - Theta) * (1 - Lambda) * gammar * I) - 
      ((1 - Theta) * Lambda *  gammad * I)
    
    dH =   (Theta * gammah * I) - (Lambda * etad * H) - 
      ((1 - Lambda) * etar * H)
    
    dD =   ((1 - Theta) * Lambda *  gammad * I) + 
      (Lambda * etad * H) 
    
    dR =   ((1 - Theta) * (1 - Lambda) * gammar * I) + 
      ((1 - Lambda) * etar * H)
    
    return(list(c(dS, dE, dI, dH, dD, dR))) 
  })
}

# Distance Function
loglikelihood_fun = function(parameters, dat) { 
  
  beta = parameters[1]
  
  output = as.data.frame(ode(y = initial_state_values, 
                             times = times, 
                             func = seidr_model,
                             parms = c(beta = beta)))  
  
  # Calculate log-likelihood 
  LL = -sum(dpois(x = dat$number_reported, lambda = output$H[output$time %in% dat$time], log = TRUE))
  return(LL) 
}

# Optimization
param_otimim = optim(par = c(0.4),
                     fn = loglikelihood_fun,
                     dat = reported_data, 
                     control = list(fnscale = 1),
                     method = 'Brent',
                     lower = 0,
                     upper = 2) 

beta_estimated = round(param_otimim$par[1], 4)
######################
parameters = c(beta = beta_estimated, 
               sigma  = sigma,
               gammah = gammah, 
               gammad = gammad, 
               gammar = gammar,
               etad   = etad,
               etar   = etar, 
               Theta  = Theta,
               Lambda = Lambda)

output = as.data.frame(ode( y     = initial_state_values, 
                            times = times, 
                            func  = seidr_model,
                            parms = parameters))

### Dates
# Change time index
output$time = series_days
time = date_ini_hospi + days(1:num_days_hospi - 1)
reported_data = data.frame(time, number_reported)

last_day_hosp = tail(time, n = 1)
d = day(last_day_hosp)
m = months(last_day_hosp, abbreviate = T)
y = year(last_day_hosp)
day_month_year_last_day_hosp = paste(d,"/",m,"/",y)
###
##### Plot Resíduos
out_res = output[output$time %in% reported_data$time, output_column(c('time','H'))]
res = reported_data$number_reported - out_res$H
sum_res = sum(res)/num_days_hospi
plot(res, main = sum_res)
abline(h = 0)
hist(res)
#####

# PLot the model fit
base_1 = ggplot() +
  geom_line(data = output, aes(x = time, y = H), colour = '#119299', size = 2) +
  geom_point(data = reported_data, aes(x = time, y = number_reported, colour = "Hospitalizados reais"), size = 2) + 
  xlab("Data") +                                              
  ylab("Número de Hospitalizados") +                                 
  labs(title = paste("Modelo Calibrado até",day_month_year_last_day_hosp), colour = "") +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(angle = 30, vjust = 0.5, size = 13)) +
  theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(angle = 0, vjust = 0.5, size = 13)) +
  theme(legend.position = "bottom") +
  theme(text = element_text(size = 18))

base_1 + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")

####################
output_long = melt(as.data.frame(output), id = "time") 

base_2 = ggplot() +
  geom_line(data = output_long, aes(x = time, y = value, colour = variable, group = variable)) +                              
  geom_point(data = reported_data, aes(x = time, y = number_reported, colour = "Casos")) + 
  xlab("Tempo (dias)")+                                              
  ylab("Número de Pessoas") +                                 
  labs(title = paste("Beta Calibrado = ",beta_estimated, colour = ""))+
  labs(colour = "")

base_2 + scale_x_date(date_breaks = '2 week', date_labels = "%d %b")
#############
# R0 (next-generation matrix) - Ini

F1 = quote(beta * S * I / N)
F2 = 0
F3 = 0

###################################################
Vm1 = quote(sigma * E)
Vm2 = quote(Theta * gammah * I + (1 - Theta) * 
              (1 - Lambda) * gammar * I + (1 - Theta) * 
              Lambda * gammad * I)
Vm3 = quote(Lambda * etad * H + (1 - Lambda) * etar * H)

###################################################
Vp1 = 0
Vp2 = quote(sigma * E)
Vp3 = quote(Theta * gammah * I)

###################################################
V1 = substitute(a - b, list(a = Vm1, b = Vp1))
V2 = substitute(a - b, list(a = Vm2, b = Vp2))
V3 = substitute(a - b, list(a = Vm3, b = Vp3))

###################################################
f11 = D(F1, "E"); f12 = D(F1, "I"); f13 = D(F1, "H")
f21 = D(F2, "E"); f22 = D(F2, "I"); f23 = D(F2, "H")
f31 = D(F3, "E"); f32 = D(F3, "I"); f33 = D(F3, "H")

v11 = D(V1, "E"); v12 = D(V1, "I"); v13 = D(V1, "H")
v21 = D(V2, "E"); v22 = D(V2, "I"); v23 = D(V2, "H")
v31 = D(V3, "E"); v32 = D(V3, "I"); v33 = D(V3, "H")

###################################################

paras = list(N = N, S = N, E = 0, I = 0, 
             H = 0, D = 0, R = 0, 
             beta = beta_estimated,
             sigma = sigma, Theta = Theta, 
             Lambda = Lambda, gammah = gammah,
             gammad = gammad, gammar = gammar, 
             etad = etad, etar = etar)

f = with(paras, 
         matrix(c(eval(f11), eval(f12), eval(f13),
                  eval(f21), eval(f22), eval(f23),
                  eval(f31), eval(f32), eval(f33)),
                nrow = 3, byrow = T))

v = with(paras, 
         matrix(c(eval(v11), eval(v12), eval(v13),
                  eval(v21), eval(v22), eval(v23),
                  eval(v31), eval(v32), eval(v33)),
                nrow = 3, byrow = T))
###################################################
R0 = max(eigen(f %*% solve(v))$values)
# R0 (next-generation matrix) - End
###################################################

# R effective
Re = R0 * (output$S[output$time == last_day_hosp]/N)