import symfit
import sympy
import scipy.integrate as integrate
from sympy.core.numbers import Infinity as inf
import json
import numpy as np
import pdb


def chebyPol(variable:"A symfit.Variable instance", 
             coefs, 
             norm:"If False, return an unnormalized polynomial"=True):
    """
    Create a chebyshev polynomial given a list of coeficients starting from the 0th coeficient.
    """
    _pols = [coefs[n]*symfit.chebyshevt(n,variable) for n in range(len(coefs))]
    unnPDF =sum(_pols)
    if norm:
        norm = sympy.integrate(unnPDF, (variable, -1, 1))
        return unnPDF/norm
    return unnPDF


def productPDF(variable:"symfit.Variable", 
               limits:"To renormalize", 
               pdfs:"As a python iterable"):
    """The product of a PDF is not necesarilly a PDF if they depend on the same variable; at least, one must normlalize"""
    pdf = pdfs[0]
    for i in range(1, len(pdfs)):
        pdf*=pdfs[i]
    norm = sympy.integrate(pdf, (variable, limits[0], limits[1]))
    return pdf/norm


def gaussian(var:"symfit.Variable", mu, sigma, limits=None, norm=True):
    """A sympy gaussian to be used by symfit"""
    unnPDF_ = symfit.exp(-(var-mu)**2/(2*sigma**2))
    if not norm:
        return unnPDF_
    if limits:
        #norm=sympy.integrate(unnPDF_, (var, limits[0], limits[1])).doit()
        ini, fin = symfit.sympify(limits[0]), symfit.sympify(limits[1]) 
        norm = integrate.quad(unnPDF_, ini, fin)[0]
        #norm = symfit.sympify(norm)
    else:
        norm =(sigma*sympy.sqrt(2*sympy.pi))
    return unnPDF_/norm


def tail_power(var:"symfit.Variable", mu, sigma, alpha, n):
    """To be used by the Crystal Ball implementation on symfit"""
    A = ((n/symfit.Abs(alpha))**n) *  symfit.exp(-0.5*symfit.Abs(alpha)**2)
    B = n/symfit.Abs(alpha) - symfit.Abs(alpha)
    A, B = symfit.sympify(A), symfit.sympify(B)

    return A*(B-(var-mu)/sigma)**(-n)

def power_law(var:"symfit.Variable", a, k):
    """To be used by the Crystal Ball implementation on symfit"""
    return a*var**k


def crystal_ball(var:"symfit.Variable", mu, sigma, alpha, n, limits=None, norm=True):
    """Crystal Ball implementation on symfit *Beware the case with n=1 since the normalization could be faulty*
    \n--https://en.wikipedia.org/wiki/Crystal_Ball_function """
    mu = symfit.sympify(mu)
    alpha = symfit.Abs(alpha)
    sigma = symfit.sympify(sigma)
    n = symfit.sympify(round(n))
    print(alpha, limits)

    gauss =  gaussian(var, mu, sigma, norm=False)
    tail  = tail_power(var, mu, sigma, alpha, n)
    
    t = (var-mu)/sigma

    conditionG = symfit.simplify(t>-alpha)
    conditionT = symfit.simplify(t<=-alpha)
    #condition = sympy.simplify((t>-symfit.Abs(alpha))
    unnPDF_ = symfit.Piecewise((tail, t<=-alpha), (gauss,t>-alpha), evaluate=True)
    #unnPDF_ = symfit.Piecewise((gauss,condition>-alpha), (tail, True))
    
    if not norm:
        return unnPDF_
    
    if limits:
        print(limits, var)
        try:
            norm= (sympy.integrate(unnPDF_, (var, limits[0], limits[1]))).evalf()
        except Exception as err:
            print(err)
            print('falling back to numeric integration')
            norm = integrate.quad(unnPDF_, limits[0], limits[1])[0]
            print(norm)
            
    else:
        if n==1:
            norm=sympy.integrate(unnPDF_, (var, mu-(100*sigma), mu+(100*sigma)) )
        else:
            C = (n/symfit.Abs(alpha)) * (1/(n-1)) * symfit.exp(-symfit.Abs(alpha)**2/2)
            D = sympy.sqrt(symfit.pi/2) * (1 + symfit.erf(symfit.Abs(alpha)/symfit.sqrt(2)) )
            norm = (sigma*(C+D))        
    return unnPDF_/norm


def crystal_ball_z(var:"symfit.Variable", mu, sigma, alpha, n, limits=None, norm=True):
    """Crystal Ball implementation on symfit *Beware the case with n=1 since the normalization could be faulty*
    \n--https://en.wikipedia.org/wiki/Crystal_Ball_function """
    mu = symfit.sympify(mu)
    abs_alpha = symfit.Abs(symfit.sympify(alpha))
    sigma = symfit.sympify(sigma)
    n = symfit.sympify(n)
    if alpha<0:
        t = (-var+mu)/sigma
    else:
        t = (var-mu)/sigma
        
    A = ((n/abs_alpha)**n) *  symfit.exp(-0.5*symfit.Abs(alpha)**2)
    B = (n/abs_alpha) - abs_alpha
    
    gauss =  exponential(t**2, 0.5, norm=False)
    tail  = power_law(var, A, -n)
    
    #condition = ((var-mu)/sigma)
    condition = sympy.simplify(t>-alpha)
    #unnPDF_ = symfit.Piecewise((gauss,condition>-alpha), (tail, condition<= -alpha))
    unnPDF_ = symfit.Piecewise((gauss,condition), (tail, True))
    
    if not norm:
        return unnPDF_
    
    if limits:
        try:
            norm=1/(sympy.integrate(unnPDF_, (var, limits[0], limits[1]))).evalf()
        except Exception as err:
            print(err)
            print('falling back to numeric integration')
            norm = 1/(integrate.quad(unnPDF_, limits[0], limits[1]))[0]
            print(norm)
            
    else:
        if n==1:
            norm=1/(sympy.integrate(unnPDF_, (var, mu-(100*sigma), mu+(100*sigma)) ))
        else:
            C = (n/symfit.Abs(alpha)) * (1/(n-1)) * symfit.exp(-symfit.Abs(alpha)**2/2)
            D = sympy.sqrt(symfit.pi/2) * (1 + symfit.erf(symfit.Abs(alpha)/symfit.sqrt(2)) )
            norm = 1/(sigma*(C+D))        
    return unnPDF_*norm

                               
def exponential(var, lambda_, limits = [0,1], norm=True):
    """ EXPONENTIAL PDF """
    unnPDF_ = symfit.exp(var*lambda_)
    if not norm:
        return unnPDF_
    norm = 1/(sympy.integrate(unnPDF_, (var, limits[0], limits[1])))
    return norm*unnPDF_




############################################################
################## TOOLS TOOLS TOOLS TOOLS #################
############################################################

def checkParams(template, dicti, text=''):
    
    error_str = '\nPARAM DICTIONARY STRUCTURE\n\n{\n'
    for k,v in template.items():
        print(k, type(v))
        if type(v)==dict:
            error_str+='\t' +k+' : {\n'
            for l in v: 
                error_str+=f'\t\t{l} : {type(v[l])}\n'
            error_str+= '\t}\n'
        else:
            error_str+=f'\t{k} : {type(v)}\n'
    error_str+= '}\n'

    try:
        for k,v in template.items():
            _ = dicti[k]
            if type(v)==dict:
                for l in v:
                    _ = dicti[k][l]
            elif type(v)==list and not type(_) == list:
                    raise KeyError
    except KeyError:
        raise NameError(error_str+'\n\n'+text)
        