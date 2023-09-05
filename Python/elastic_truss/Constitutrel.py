import numpy as np
from model import *

def Constitutrel(L, A, funcname, member, nmemb, forcememb):
    """
    Constitutrel: Evaluates the value of strain for assumed stress values
    of each truss member.

    Parameters:
    L: array of member lengths
    A: array of cross-sectional areas
    funcname: list of function names
    member: array of member properties
    nmemb: number of members
    forcememb: array of member forces

    Returns:
    delta: array of elongations
    fdiff: array of derivatives of elongations with respect to stress
    """

    fdiff = np.zeros((nmemb, 1))
    strain = np.zeros((nmemb, 1))
    delta = np.zeros((nmemb, 1))

    for i in range(nmemb):
        sigma = (forcememb[i]*1) / A[i]
        K1 = sigma
        K2 = sigma**2
        namestr = funcname[i + 1, 0].strip()
        beta, dbeta = globals()[namestr](K1, K2, member[i, 4:])

        strain[i] = beta[0] + beta[1] * sigma + beta[2] * sigma**2
        delta[i] = strain[i] * L[i]
        fdiff[i] = (dbeta[0, 0] + dbeta[0, 1] * 2 * sigma + dbeta[1, 0] * sigma +
                   dbeta[1, 1] * 2 * sigma**2 + beta[1] + dbeta[2, 0] * sigma**2 +
                   2 * dbeta[2, 1] * sigma**3 + 2 * beta[2] * sigma) * L[i] / A[i]

    return delta, fdiff
