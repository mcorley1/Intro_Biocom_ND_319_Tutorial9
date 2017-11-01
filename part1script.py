#Exercise9, 10/27/17

import numpy
import pandas
from scipy.optimize import minimize 
from scipy.stats import norm
from plotnine import *


data = pandas.read_csv('ponzr1.csv', header=0,sep=',')

sub1=data.loc[data.mutation.isin(['WT','M124K']),:]
sub2=data.loc[data.mutation.isin(['WT','V456D']),:]
sub3=data.loc[data.mutation.isin(['WT','I213N']),:]

def nllike(p,obs):
    B0=p[0]
    sigma=p[1]
    
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll
    
initialGuess=numpyarray([1,1,1])
fit=minimize(nllike,initialGuess,method="Nelder-Mead",options={'disp': True},args=df)
print(fit.x)