import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *

leaf=pd.read_csv("leafDecomp.csv")


ggplot(leaf, aes(x="Ms", y="decomp")) + geom_point()


def nllike_null (p,obs):
    a=p[0]
    b=p[1]
    c=p[2]
    sigma=p[3]
    
    expected=a
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return(nll)
    
def nllike_linear(p,obs):
    a=p[0]
    b=p[1]
    c=p[2]
    sigma=p[3]
    
    expected=a+b*obs.Ms
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return(nll)    
    
def nllike_hump(p,obs):
    a=p[0]
    b=p[1]
    c=p[2]
    sigma=p[3]
    
    expected=a+b*obs.Ms+c*obs.Ms**2
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return(nll)


guess_hump=np.array([1,100,500,50])
fit_null=minimize(nllike_null,guess_hump, method="Nelder-Mead", options={'disp':True},args=leaf)
fit_linear=minimize(nllike_linear,guess_hump, method="Nelder-Mead", options={'disp':True},args=leaf)
fit_hump=minimize(nllike_hump,guess_hump, method="Nelder-Mead", options={'disp':True},args=leaf)

chi2_linear=2*(fit_null.fun-fit_linear.fun)
1-sp.stats.chi2.cdf(x=chi2_linear,df=1)

chi2_hump=2*(fit_null.fun-fit_hump.fun)
1-sp.stats.chi2.cdf(x=chi2_hump,df=2)

#based on the chi2 test, the hump-shaped model is the best fit for the leaf decomposition data.
