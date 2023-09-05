import numpy as np
from Constitutrel import *

def objectfunc(d):
    global S, Breq, nmemb, L, A, funcname, member, deltasupp, ns, multiplier
    
    delta, fdiff = Constitutrel(L, A, funcname, member, nmemb, d[0:nmemb])
    errvec = np.dot(S, d) - multiplier * delta - deltasupp
    errval = np.sqrt(np.dot(errvec.T, errvec))
    graderrval = (np.concatenate([-np.diag(fdiff * multiplier) @ errvec, np.dot(Breq, errvec), np.zeros(ns)]) / errval)
    
    return errval, graderrval