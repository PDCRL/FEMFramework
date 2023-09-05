import numpy as np
def model(K1, K2, matpar):
    # Model to find out the invariants of stress tensor
    # Input arguments:
    # K1, K2, and matpar - Material Parameters
    # Output arguments:
    # beta - values of beta
    # dbeta - derivative of beta
    alpha1 = matpar[0]
    alpha2 = matpar[1]
    alpha3 = matpar[2]
    alpha4 = matpar[3]
    
    beta = np.zeros(3)
    beta[0] = alpha1 * K1
    beta[1] = alpha2 * K2**alpha3 + alpha4
    beta[2] = 0
    
    dbeta = np.zeros((3, 3))
    dbeta[0, 0] = alpha1
    dbeta[1, 1] = alpha2 * alpha3 * K2**(alpha3 - 1)
    dbeta[2, 2] = 0
    
    return beta, dbeta