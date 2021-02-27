#Zombie Display
# zombie apocalypse modeling
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import integrate
from scipy.optimize import fmin
#=====================================================


#The Model
#=======================================================
def eq(par,initial_cond,start_t,end_t,incr):
     #-time-grid-----------------------------------
     t  = np.linspace(start_t, end_t,incr)
     #differential-eq-system----------------------
     def funct(y,t):
        Si=y[0]
        Zi=y[1]
        Ri=y[2]
        P,d,B,G,A=par
        # the model equations (see Munz et al. 2009)
        f0 = P - B*Si*Zi - d*Si
        f1 = B*Si*Zi + G*Ri - A*Si*Zi
        f2 = d*Si + A*Si*Zi - G*Ri
        return [f0, f1, f2]
     #integrate------------------------------------
     ds = integrate.odeint(funct,initial_cond,t)
     return (ds[:,0],ds[:,1],ds[:,2],t)
#1.Get Data
#====================================================
Td=np.array([0.5,1,1.5,2,2.2,3,3.5,4,4.5,5])#time
Zd=np.array([0,2,2,5,2,10,15,50,250,400])#zombie pop
#====================================================

#2.Set up Info for Model System
#===================================================
# model parameters
#----------------------------------------------------
P = 0       # birth rate
d = 0.0001  # natural death percent (per day)
B = 0.0095  # transmission percent  (per day)
G = 0.0001  # resurect percent (per day)
A = 0.0001  # destroy perecent (per day)
rates=(P,d,B,G,A)

# model initial conditions
#---------------------------------------------------
S0 = 500.               # initial population
Z0 = 0                  # initial zombie population
R0 = 0                  # initial death population
y0 = [S0, Z0, R0]      # initial condition vector

# model steps
#---------------------------------------------------
start_time = 0.0
end_time = 5.0
intervals = 1000
mt=np.linspace(start_time,end_time,intervals)

# model index to compare to data
#----------------------------------------------------
findindex = lambda x:np.where(mt>=x)[0][0]
mindex = map(findindex,Td)
#=======================================================


#3.Score Fit of System
#=========================================================
def score(parms):
    #a.Get Solution to system
    F0,F1,F2,T=eq(parms,y0,start_time,end_time,intervals)
    #b.Pick of Model Points to Compare
    Zm=F1[mindex]
    #c.Score Difference between model and data points
    ss=lambda data,model:((data-model)**2).sum()
    return ss(Zd,Zm)
#========================================================


#4.Optimize Fit
#=======================================================
fit_score=score(rates)
answ=fmin(score,(rates),full_output=1,maxiter=1000000)
bestrates=answ[0]
bestscore=answ[1]
P,d,B,G,A=answ[0]
newrates=(P,d,B,G,A)
#=======================================================

#5.Generate Solution to System
#=======================================================
F0,F1,F2,T=eq(newrates,y0,start_time,end_time,intervals)
Zm=F1[mindex]
Tm=T[mindex]
#======================================================




#6. Plot Solution to System
#=========================================================
plt.figure()
plt.plot(T,F1,'b-',Tm,Zm,'ro',Td,Zd,'go')
plt.legend(('Zombies','Model Points','Data Points'),
           'upper center')
plt.xlabel('days')
plt.ylabel('population')
title='Zombie Apocalypse  Fit Score: '+str(bestscore)
plt.title(title)