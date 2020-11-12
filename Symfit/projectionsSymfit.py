import symfit
import sympy
from sympy.core.numbers import Infinity as inf
import json
import numpy as np

from GeneralPDFSymfit import  *
import Signals_Symfit as Signals
from Backgrounds_Symfit import *

#Jun 2020
# Here are the proyections of the first implementation of the Complete Model:

# - Independent Variables

# - Signal Side:
#   - Angular:
#     + Theory*Efficiency
#   - Mass:
#     + Gaussian

# - Bckground Side:
#   - Angular:
#     + Chebyshev+Gaussian
#   - Mass:
#     + Exponential
def angularProy(cos, frac, AFB, FH, coefsEff, fraction, mu, sigma, coefs):
    """
    NonExtended complete PDF on the angular varibale (Integrating the mass variable)
    """
    signal = Signals.angularSignal(cos, AFB, FH, coefsEff)
    background = angularBackground(cos, fraction, mu, sigma, coefs)
    return frac*signal + (1-frac)*background
    
    
def massProy(mass, mu, sigma, lambda_, frac, limits):
    """
    NonExtended complete PDF on the mass varibale (Integrating the angular variable)
    """
    signal = gaussian(mass, mu, sigma, limits)
    background = exponential(mass, lambda_, limits)
    return frac*signal + (1-frac)*background







# Nov 2020
# Here are the proyections of the second implementation of the Complete Model:

# - Independent Variables

# - Signal Side:
#   - Angular:
#     + Theory*Efficiency
#   - Mass:
#     + Gaussian+Gaussian+CrystallBall

# - Background Side:
#   - Angular:
#     + LeftSB(Chebyshev+Gaussian) + RightSB(Chebyshev+Gaussian)
#   - Mass:
#     + Exponential+Gaussian
def angularProy_v1(cos, frac, AFB, FH, coefsEff, fraction_Left, Left, Right, **kwargs):
    """
    NonExtended complete PDF on the angular varibale (Integrating the mass variable)
    """
    signal = Signals.angularSignal(cos, AFB, FH, coefsEff)
    background = angularBackground_SideBands(cos, fraction_Left, Left, Right, **kwargs)
    return frac*signal + (1-frac)*background

def massProy_v1(mass, frac, mu, sigma, lambda_, fraction_exp,
                muCB, sigmaCB, alphaCB, nCB,  
                muGauss,sigmaGauss1, sigmaGauss2,
                fracCB, fracGauss1, fracGauss2=False, 
                limits=[5, 5.7], **kwargs):
    """
    NonExtended complete PDF on the mass varibale (Integrating the angular variable)
    """
    print("algo")
    signal = Signals.CrystalBall_2Gauss(mass, limits, muCB, sigmaCB, alphaCB, nCB, 
                       muGauss,sigmaGauss1, sigmaGauss2,
                       fracCB, fracGauss1, fracGauss2, **kwargs)
    background = Gauss_Exp(mass,  mu, sigma, lambda_, fraction_exp, limits, **kwargs)
    
    return frac*signal + (1-frac)*background
