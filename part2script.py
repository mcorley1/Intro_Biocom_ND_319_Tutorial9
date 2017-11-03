###Exercise 9 Part 2
import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *
data=pandas.read_csv("MmarinumGrowth.csv")
def Grate(p,obs):
    um=p[0]
    K=p[1]
    sigma=p[2]
    ue=um*(obs.S/(obs.S+K))
    nll=-1*norm(ue,sigma).logpdf(obs.u).sum()
    return nll
probswrong=numpy.array([1,1,1])
best=minimize(Grate,probswrong,method="Nelder-Mead",options={"disp":True},args=data)
print (best.x)
