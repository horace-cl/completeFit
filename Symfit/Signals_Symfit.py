import symfit
import sympy
from sympy.core.numbers import Infinity as inf
import json
import numpy as np
from GeneralPDFSymfit import  *
from Signals_Symfit import *
from Backgrounds_Symfit import *
from projectionsSymfit import *


###Theory
def decayWidth(cos, AFB, FH, **kwargs):
    """Angular Model taken from AN2014-244 [TODO ADD REFS]"""
    pdf = (3/4)*(1-FH)*(1-cos**2)
    pdf+= (1/2)*FH
    pdf+= AFB*cos
    return pdf

###AngularPart of the signal side
def angularSignal(cos, AFB, FH, coefsEff, **kwargs):
    """Angular part of the signal from the complete PDF\nIt is the PRODUCT of the Decay Width times the Efficiency"""
    efficiency = chebyPol(cos, coefsEff)
    decay = decayWidth(cos, AFB, FH)
    return productPDF(cos, [-1,1], [decay,efficiency])



#### Mass Parts of the signal side
def CrystalBall_2Gauss(mass, limits, muCB, sigmaCB, alphaCB, nCB, 
                       muGauss,sigmaGauss1, sigmaGauss2,
                       fracCB, fracGauss1, fracGauss2=False, **kwargs):

    CB = crystal_ball(mass, muCB, sigmaCB, alphaCB, nCB, limits)
    G1 = gaussian(mass, muGauss, sigmaGauss1, limits)
    G2 = gaussian(mass, muGauss, sigmaGauss2, limits)
    if fracGauss2:
        tot = fracCB + fracGauss1 + fracGauss2
        fracCB, fracGauss1, fracGauss2 = fracCB/tot, fracGauss1/tot, fracGauss2/tot
        
    mass = fracCB*CB + fracGauss1*G1 + (1-fracCB-fracGauss1)*G2
    return mass



### COMPLETE Signal  ### COMPLETE Signal### ### COMPLETE Signal###
### COMPLETE Signal  ### COMPLETE Signal### ### COMPLETE Signal###
### COMPLETE Signal  ### COMPLETE Signal### ### COMPLETE Signal###


#Adding mass variable (Gausian model deprecated)
def completeSignal(mass, cos, AFB, FH, coefsEff, mu, sigma, limits):
    """Signal part of the complete PDF"""
    angular = angularSignal(cos, AFB, FH, coefsEff)
    mass = gaussian(mass, mu, sigma, limits)
    return angular*mass

#Adding mass variable (2 Gaussians + CrystalBall)
def completeSignalCrystalBall_2Gauss(mass, cos, AFB, FH, coefsEff, 
                                     muCB, sigmaCB, alphaCB, nCB,  
                                    muGauss,sigmaGauss1, sigmaGauss2,
                                    fracCB, fracGauss1, fracGauss2=False, 
                                    limits=[5, 5.7], **kwargs):
    """Signal part of the complete PDF"""

    
    mass = CrystalBall_2Gauss(mass, limits, muCB, sigmaCB, alphaCB, nCB, 
                       muGauss,sigmaGauss1, sigmaGauss2,
                       fracCB, fracGauss1, fracGauss2, **kwargs)
    angular = angularSignal(cos, AFB, FH, coefsEff)
    
    return angular*mass
