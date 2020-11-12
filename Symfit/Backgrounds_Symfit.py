import symfit
import sympy
from sympy.core.numbers import Infinity as inf
import json
import numpy as np
from GeneralPDFSymfit import  *




### Angular Components START ### Angular Components START ### Angular Components START ### Angular Components START 
### Angular Components START ### Angular Components START ### Angular Components START ### Angular Components START 
### Angular Components START ### Angular Components START ### Angular Components START ### Angular Components START 
### Angular Components START ### Angular Components START ### Angular Components START ### Angular Components START 
### Angular Components START ### Angular Components START ### Angular Components START ### Angular Components START 
### Angular Components START ### Angular Components START ### Angular Components START ### Angular Components START 

#Initial Guess
def angularBackground(cos, fraction, mu, sigma, coefs):
    """
    Angular part of the background of the complete PDF.
    It is a SUM of a gaussian plus a Chebysev polynomial
    """
    chebAngular = chebyPol(cos, coefs)
    gaussAngular = gaussian(cos, mu, sigma, [-1,1])
    return fraction*gaussAngular + (1-fraction)*chebAngular


def chebyPlusGauss(cos, mu, sigma, coefs, fraction_Cheby, **kwargs):
    cheb = chebyPol(cos, coefs)
    if fraction_Cheby==1:
        return cheb
    gauss = gaussian(cos, mu, sigma, [-1,1])
    return fraction_Cheby*cheb + (1-fraction_Cheby)*gauss   
    

#Heribertos IDEA after looking how bad projections looked
def angularBackground_SideBands(cos, fraction_Left, Left, Right, **kwargs):
    """In this model we have a pdf for each sideband composed by a chebyshev polynomial and possibly a gaussian the free parameter is the fraction between them"""
    #template = {
    #                'coefs':[0.1,0.1],   'fraction_Cheby':0.5,
    #                'mu':0.75,           'sigma':0.1
    #            }
    

    #error_str = 'If only chebyshev, fill gaussian slots with anything, but set fraction_Cheby=0\n'
    #error_str+= 'If only gaussian, not implemented yet!!!'
         #   raise NameError(error_str)
    #checkParams(template, Left, error_str)
    #checkParams(template, Right, error_str)
    
    
    LeftModel = chebyPlusGauss(cos, **Left)
    RightModel = chebyPlusGauss(cos, **Right)
        
    return fraction_Left*LeftModel + (1-fraction_Left)*RightModel

### Angular Components END ### Angular Components END ### Angular Components END ### Angular Components END 
### Angular Components END ### Angular Components END ### Angular Components END ### Angular Components END 
### Angular Components END ### Angular Components END ### Angular Components END ### Angular Components END 
### Angular Components END ### Angular Components END ### Angular Components END ### Angular Components END 
### Angular Components END ### Angular Components END ### Angular Components END ### Angular Components END 
### Angular Components END ### Angular Components END ### Angular Components END ### Angular Components END 







### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 

# Gaussian plus Exponential
def Gauss_Exp(mass,  mu, sigma, lambda_, fraction_exp, limits, **kwargs):
    #error_str= 'Fill all slots, and use fraction to set only gauss or only exponential'
    template = {
        'gauss':{
            'mu':5.1,
            'sigma':0.04},
        'lambda' : -0.2,
        'fraction_Exp' : 0.5
    } 
    #checkParams(template, **kwargs)
    if fraction_exp>1 or fraction_exp<0:
        raise NotImplementedError(f'\nFractions must be in [0, 1]\nWhat to do when fraction_exp = {fraction_exp}?')
    mass_exp = exponential(mass, lambda_, limits)
    mass_gauss = gaussian(mass, mu, sigma, limits)
    mass = fraction_exp*mass_exp + (1-fraction_exp)*mass_gauss
    return mass

### Mass Components ### Mass Components ### Mass Components ### Mass Components
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 
### Mass Components ### Mass Components ### Mass Components ### Mass Components 












### COMPLETE Background ### COMPLETE Background ### ### COMPLETE Background ###
### COMPLETE Background ### COMPLETE Background ### ### COMPLETE Background ###
### COMPLETE Background ### COMPLETE Background ### ### COMPLETE Background ###
    
    
    
def completeBackground(mass, cos, fraction, mu, sigma, coefs, lambda_, limits):
    angular = angularBackground(cos, fraction, mu, sigma, coefs)
    mass = exponential(mass, lambda_, limits)
    return mass*angular




def completeBackground_SB_GausExp(mass, cos, fraction_Left, Left, Right, 
                                  mu, sigma, lambda_, fraction_exp, **kwargs):
    """In this version of the model we have background angular pdf for each SB, and for the mass variable a gaussian+exponential"""
    angular = angularBackground_SideBands(cos, fraction_Left, Left, Right, **kwargs)
    masspdf = Gauss_Exp(mass,  mu, sigma, lambda_, fraction_exp, [-1,1], **kwargs)
    
    return masspdf*angular