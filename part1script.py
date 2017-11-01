#Exercise9, 10/27/17

import numpy
import pandas
from scipy.optimize import minimize 
from scipy.stats import norm
from plotnine import *

data = pandas.read_csv('ponzr1.csv', header=0,sep=',')

#subset data into three different data frames 
sub1=data.loc[data.mutation.isin(['WT','M124K']),:]
sub2=data.loc[data.mutation.isin(['WT','V456D']),:]
sub3=data.loc[data.mutation.isin(['WT','I213N']),:]


#Make new data frame with 'group' column (your x=0 or x=1)
#var2=pandas.DataFrame({'y':var1.col2name, 'x':})
sub1frame=pandas.DataFrame({'y':sub1.ponzr1Counts,'x':0})
sub2frame=pandas.DataFrame({'y':sub2.ponzr1Counts,'x':0})
sub3frame=pandas.DataFrame({'y':sub3.ponzr1Counts,'x':0})
#Designate 'treatment' group as x=1
#var2.loc[var1.col1name=='name of treatment group', 'x']=1
sub1frame.loc[sub1.mutation=='M124K','x']=1
sub2frame.loc[sub2.mutation=='V456D','x']=1
sub3frame.loc[sub3.mutation=='I213N','x']=1
#print(sub3frame)

#### Define null function
def nllikeNull(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

# estimate parameters by minimizing the NLL for sub1frame
initialGuess=numpy.array([1,1])
fit=minimize(nllikeNull,initialGuess,method="Nelder-Mead",options={'disp': True},args=sub1frame)
nullM124K= fit.fun #M124K
#gives NLL value for affect of mutation 1 in null model
#print(nullM124K)

initialGuess=numpy.array([1,1])
fit=minimize(nllikeNull,initialGuess,method="Nelder-Mead",options={'disp': True},args=sub2frame)
nullV456D= fit.fun #V456D

initialGuess=numpy.array([1,1])
fit=minimize(nllikeNull,initialGuess,method="Nelder-Mead",options={'disp': True},args=sub3frame)
nullI213N = fit.fun #I213N

#### Define function y=B0+B1*treat+error
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    
    expectedAlt=B0+B1*obs.x
    nll=-1*norm(expectedAlt,sigma).logpdf(obs.y).sum()
    return nll
   
#estimate parameters by minimizing the NLL
initialGuess=numpy.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=sub1frame)
altM124K = fit.fun #M124K

initialGuess=numpy.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=sub2frame)
altV456D = fit.fun #V456D

initialGuess=numpy.array([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=sub3frame)
altI213N = fit.fun #I213N




