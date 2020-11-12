import symfit
import sympy
from sympy.core.numbers import Infinity as inf
import json
from functions import polyToChebCoefs
import numpy as np

with open("/home/horacio/Documents/hcl/toysParams.json") as f:
    toysParams = json.load(f)


    
    
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
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
    """The product of a PDF is not a PDF if they depend on the same variable; at least one must renormalize"""
    pdf = pdfs[0]
    for i in range(1, len(pdfs)):
        pdf*=pdfs[i]
    norm = sympy.integrate(pdf, (variable, limits[0], limits[1]))
    return pdf/norm


def gaussian(var:"symfit.Variable", mu, sigma, limits=None, norm=True):
    """A sympy gaussian to be used by symfit"""
    unnPDF_ = symfit.exp(-(var-mu)**2/(2*sigma**2))
    if limits:
        norm=1/(sympy.integrate(unnPDF_, (var, limits[0], limits[1])))
    else:
        norm = 1/(sigma*sympy.sqrt(2*sympy.pi))
    return norm*unnPDF_

def tail_power(var:"symfit.Variable", mu, sigma, alpha, n):
    """To be used by the Crystal Ball implementation on symfit"""
    A = ((n/symfit.Abs(alpha))**n) *  symfit.exp(-symfit.Abs(alpha)**2/2)
    B = n/symfit.Abs(alpha) - symfit.Abs(alpha)
    
    return A*(B-(var-mu)/sigma)**(-n)


def crystal_ball(var:"symfit.Variable", mu, sigma, alpha, n, limits=None):
    """Crystal Ball implementation on symfit *Beware the case with n=1 since the normalization could be faulty*
    \n--https://en.wikipedia.org/wiki/Crystal_Ball_function """
    gauss =  gaussian(var, mu, sigma, norm=False)
    tail  = tail_power(var, mu, sigma, alpha, n)
    
    condition = (var-mu)/sigma 
    unnPDF_ = symfit.Piecewise((gauss,condition>-alpha), (tail, condition<= -alpha))
    
    if limits:
        norm=1/(sympy.integrate(unnPDF_, (var, limits[0], limits[1])))
    else:
        if n==1:
            norm=1/(sympy.integrate(unnPDF_, (var, mu-(100*sigma), mu+(100*sigma)) ))
        else:
            C = (n/symfit.Abs(alpha)) * (1/(n-1)) * symfit.exp(-symfit.Abs(alpha)**2/2)
            D = sympy.sqrt(symfit.pi/2) * (1 + symfit.erf(symfit.Abs(alpha)/symfit.sqrt(2)) )
            norm = 1/(sigma*(C+D))        
    return unnPDF_*norm


def exponential(var, lambda_, limits = [0,1]):
    """ EXPONENTIAL PDF """
    unnPDF_ = symfit.exp(var*lambda_)
    norm = 1/(sympy.integrate(unnPDF_, (var, limits[0], limits[1])))
    return norm*unnPDF_
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS
#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS#TOOLS




###SIGNAL
###SIGNAL
def decayWidth(cos, AFB, FH):
    """Angular Model taken from AN2014-244 [TODO ADD REFS]"""
    pdf = (3/4)*(1-FH)*(1-cos**2)
    pdf+= (1/2)*FH
    pdf+= AFB*cos
    return pdf


def angularSignal(cos, AFB, FH, coefsEff):
    """Angular part of the signal from the complete PDF\nIt is the PRODUCT of the Decay Width and the Efficiency"""
    efficiency = chebyPol(cos, coefsEff)
    decay = decayWidth(cos, AFB, FH)
    return productPDF(cos, [-1,1], [decay,efficiency])


def completeSignal(mass, cos, AFB, FH, coefsEff, mu, sigma, limits):
    """Signal part of the complete PDF"""
    angular = angularSignal(cos, AFB, FH, coefsEff)
    mass = gaussian(mass, mu, sigma, limits)
    return angular*mass

def completeSignalCrystalBall_2Gauss(mass, cos, AFB, FH, coefsEff, limits, **kwargs):
    """Signal part of the complete PDF"""
    needed_params = ['CB_mu', 'CB_sigma', 'CB_alpha', 'CB_n', 'CB_frac', 
                     'G_mu', 
                     'G1_sigma', 'G1_frac', 
                     'G2_sigma', 'G2_frac']
    
    for k in needed_params:
        if k not in kwargs: raise KeyError(f'Key: {k} Not found in **kwargs')
    
    angular = angularSignal(cos, AFB, FH, coefsEff)
    CB = crystal_ball(mass, kwargs['CB_mu'], kwargs['CB_sigma'], 
                      kwargs['CB_alpha'], kwargs['CB_n'], limits)
    G1 = gaussian(mass, kwargs['G_mu'], kwargs['G1_sigma'], limits)
    G2 = gaussian(mass, kwargs['G_mu'], kwargs['G2_sigma'], limits)
    
    mass = (CB*kwargs['CB_frac'] + G1*kwargs['G1_frac'] + G2*kwargs['G2_frac'])
    mass/= (kwargs['CB_frac']+kwargs['G1_frac']+kwargs['G2_frac'])
    
    return angular*mass
###SIGNAL
###SIGNAL




###BACKGROUND
###BACKGROUND
#JUN 2020
def angularBackground(cos, fraction, mu, sigma, coefs):
    """
    Angular part of the background of the complete PDF.
    It is a SUM of a gaussian plus a Chebysev polynomial
    """
    chebAngular = chebyPol(cos, coefs)
    gaussAngular = gaussian(cos, mu, sigma, [-1,1])
    return fraction*gaussAngular + (1-fraction)*chebAngular
    
def completeBackground(mass, cos, fraction, mu, sigma, coefs, lambda_, limits):
    angular = angularBackground(cos, fraction, mu, sigma, coefs)
    mass = exponential(mass, lambda_, limits)
    return mass*angular




#NOV 2020
def angularBackground_SideBands(cos, fraction, paramsLeft, paramsRight):
    """In this model we have a pdf for each sideband composed by a chebyshev polynomial and possibly a gaussian the free parameter is the fraction between them"""
    chebLeft = chebyPol(cos, paramsLeft['cheby'])
    if 'gauss' in paramsLeft:
        gaussLeft = gaussian(cos, paramsLeft['gauss']['mu'], paramsLeft['gauss']['sigma'], [-1,1])
        angularLeft = chebLeft*gaussL
    else:
        angularLeft = chebLeft
        
    chebRight = chebyPol(cos, paramsRight['cheby'])
    if 'gauss' in paramsRight:
        gaussRight = gaussian(cos, paramsRight['gauss']['mu'], paramsRight['gauss']['sigma'], [-1,1])
        angularRight = chebRight*gaussRight
    else:
        angularRight = chebRight
        
    return fraction*angularLeft + (1-fraction)*angularRight

def completeBackground_SB_GausExp(mass, cos, frac_mass, mu, sigma, lambda_, frac_ang, paramsLeft, paramsRight, limits):
    """In this version of the model we have background angular pdf for each SB, and for the mass variable a gaussian+exponential"""
    angular = angularBackground_SideBands(cos, frac_ang, paramsLeft, paramsRight)
    
    mass_exp = exponential(mass, lambda_, limits)
    mass_gauss = gaussian(mass, mu, sigma, limits)
    mass = frac_mass * mass_gauss + (1-frac_mass)*mass_exp
    
    return mass*angular
###BACKGROUND
###BACKGROUND
    
    
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
#Jun 2020
def angularProy(cos, frac, AFB, FH, coefsEff, fraction, mu, sigma, coefs):
    """
    NonExtended complete PDF on the angular varibale (Integrating the mass variable)
    """
    signal = angularSignal(cos, AFB, FH, coefsEff)
    background = angularBackground(cos, fraction, mu, sigma, coefs)
    return frac*signal + (1-frac)*background
    
def massProy(mass, mu, sigma, lambda_, frac, limits):
    """
    NonExtended complete PDF on the mass varibale (Integrating the angular variable)
    """
    signal = gaussian(mass, mu, sigma, limits)
    background = exponential(mass, lambda_, limits)
    return frac*signal + (1-frac)*background



#Nov 2020
def angularProy_SB(cos, frac, AFB, FH, coefsEff, fraction, paramsLeft, paramsRight):
    """
    NonExtended complete PDF on the angular varibale (Integrating the mass variable)
    """
    signal = angularSignal(cos, AFB, FH, coefsEff)
    background = angularBackground_SideBands(cos, fraction, paramsLeft, paramsRight)
    return frac*signal + (1-frac)*background

def massProy_CB2G(mass, mu, sigma, lambda_, frac, limits):
    """
    NonExtended complete PDF on the mass varibale (Integrating the angular variable)
    """
    CB = crystal_ball(mass, kwargs['CB_mu'], kwargs['CB_sigma'], 
                      kwargs['CB_alpha'], kwargs['CB_n'], limits)
    G1 = gaussian(mass, kwargs['G_mu'], kwargs['G1_sigma'], limits)
    G2 = gaussian(mass, kwargs['G_mu'], kwargs['G2_sigma'], limits)
    
    signal = (CB*kwargs['CB_frac'] + G1*kwargs['G1_frac'] + G2*kwargs['G2_frac'])
    signal/= (kwargs['CB_frac']+kwargs['G1_frac']+kwargs['G2_frac'])
    
    background_exp = exponential(mass, lambda_, limits)
    background_gauss = gaussian(mass, mu, sigma, limits)
    background = frac_mass * background_gauss + (1-frac_mass)*background_exp
    
    return frac*signal + (1-frac)*background
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###
###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###PROYECCIONS###





def Complete(mass, cos, AFB, FH, frac, coefsEff, mu, sigma,fraction, muB, sigmaB, coefs, lambda_, limitsAng, limitsMass):
    

    signal = completeSignal(mass, cos, AFB, FH, coefsEff, mu, sigma, limitsMass)
    back = completeBackground(mass, cos, fraction, muB, sigmaB, coefs, lambda_, limitsMass)
    complete = frac*signal+(1-frac)*back
    return complete






class complete(object):
    """
    TODO : We need a method to evaluate the model only giving cos and mass
    """
    def __init__(self, nBin, AFB, FH, 
                limits = [5,5.7], 
                limitsAFB=None, 
                limitsFH=None, 
                fixed_params = [],
                params = None,
                efficiency='OriginalEffyCoefs'):
        
        if not limitsAFB:
            d1 = [0.5*FH-AFB, AFB+0.5*FH]
            limitsAFB = [AFB-min(d1)*0.9, AFB+min(d1)*0.9]
        
        if not limitsFH:
            d2 = [FH, 3-FH]
            limitsFH = [FH-min(d2), FH+min(d2)]
            
        self.fixed_params = [f.upper() for f in fixed_params]
        
        if params:
            self.parameters = params[str(nBin)]
        else:
            self.parameters = toysParams[str(nBin)]  
        
        try:
            self.efficiencyCoefs = self.parameters['fixed'][efficiency]
            print('Using :', efficiency)
        except Exception as e:
            self.efficiencyCoefs = self.parameters['fixed']["OriginalEffyCoefs"]
            print(e)
            print(efficiency + ' Not found in json, falling back to "OriginalEffyCoefs"')
        
        S, B = self.parameters['yield'].values()
        if S>B : S,B = B,S
        self._frac = S/(S+B)
        self._AFB = AFB
        self._FH = FH
        
        free = self.parameters['free']
        self._lambda = free['LAMBDA']
        self._mu = free['MU']
        self._sigma = free['SIGMA']
        
        
        self.__init_Vars()
        self.__init_Params(limitsAFB, limitsFH)
                              
        self.create_pdf()
           
            
            

    def __init_Params(self, limitsAFB, limitsFH):
        
        self.frac  = symfit.Parameter('FRAC', self._frac, 0, 1)
        
        if 'AFB' in self.fixed_params:
            self.AFB = self._AFB
        else:
            self.AFB = symfit.Parameter('AFB', self._AFB, limitsAFB[0], limitsAFB[1])        
        
        
        if 'FH' in self.fixed_params:
            self.FH = self._FH                                              
        else:
            self.FH = symfit.Parameter('FH', self._FH, limitsFH[0], limitsFH[1])                                                          
        
        
        lambLims = [-4,-0.1]
        if self._lambda>-0.1:
            lambLims[1]=0
        if self._lambda<-4:
            lambLims[0]=self._lambda-1
            
        if 'LAMBDA' in self.fixed_params:
            self.lambdA = self._lambda
        else:
            self.lambdA = symfit.Parameter('LAMBDA', self._lambda, lambLims[0], lambLims[1])                             
            
        if 'MU' in self.fixed_params:
            self.mu = self._mu
        else:
            self.mu = symfit.Parameter('MU', self._mu, 5.2, 5.4 )                           
        
            
        if 'SIGMA' in self.fixed_params:
            self.sigma = self._sigma
        else:
            self.sigma = symfit.Parameter('SIGMA', self._sigma, 0.01, 0.1)
            
        
#         self._Signal()
#         self._Back()
        
        
    def _check_valuesInside(self):
        """I would like a way to guess the best range of the initial values. Only for MC"""
        pass
        

    def __init_Vars(self):
        self.cos = symfit.Variable('cos')
        self.mass = symfit.Variable('mass')
    
    
#     def _Signal(self):
#         """TODO! INCOMPLETE"""
#         fix = self.parameters['fixed']
#         coefsEff = polyToChebCoefs(fix['OriginalEffyCoefs'])
#         self.angular_S = angularSignal(self.cos, self.AFB, self.FH, coefsEff)
#         self.mass_S = gaussian(self.mass, self.mu, self.sigma, [5.0,5.7])
#         self.pdf_S = self.angular_S*self.mass_S
    
#     def _Back(self):
#         """TODO!"""
#         fix = self.parameters['fixed']
#         self.angular_B = angularBackground(self.cos, fix['frac'], fix['mu'], fix['sigma'], fix['sidebandsCoefs'])
#         self.mass_B = exponential(self.mass, self.lambdA, [5.0,5.7])
#         self.pdf_B = self.angular_B*self.mass_B

    def create_pdf(self):
        fix = self.parameters['fixed']
        #if 'HoracioEffyCoefs' in self.parameters:
         #   print('Using HoracioEffyCoefs')
          #  coefsEff = polyToChebCoefs(fix['HoracioEffyCoefs'])
        #else:
         #   print('Using OriginalEffyCoefs')
          #  coefsEff = polyToChebCoefs(fix['OriginalEffyCoefs'])
        complete = Complete(self.mass, self.cos, self.AFB, self.FH, 
                        self.frac, self.efficiencyCoefs, mu=self.mu, 
                        sigma=self.sigma, fraction=fix['frac'],muB=fix['mu'],
                        sigmaB=fix['sigma'], coefs=fix['sidebandsCoefs'],
                        lambda_=self.lambdA, limitsMass = [5.0, 5.7], limitsAng=[-1,1])
        
        self.projection_mass= massProy(self.mass, self.mu, self.sigma, 
                                       self.lambdA, self.frac, [5.0, 5.7])
        self.projection_angular= angularProy(self.cos, self.frac, self.AFB, self.FH, 
                                    self.efficiencyCoefs, fix['frac'], fix['mu'],
                                    fix['sigma'], fix['sidebandsCoefs'])
        
        self.signal_angular = self.frac*angularSignal(self.cos, 
                                self.AFB, self.FH, self.efficiencyCoefs)
        self.signal_mass = self.frac*gaussian(self.mass, self.mu, self.sigma, [5.0,5.7])
        
        self.background_angular = (1-self.frac)*angularBackground(self.cos,
                                            fix['frac'], fix['mu'], 
                                            fix['sigma'], fix['sidebandsCoefs'])
        self.background_mass = (1-self.frac)*exponential(self.mass, 
                                                    self.lambdA, [5.0,5.7])
        
        self.only_signal = self.frac*decayWidth(self.cos, self.AFB, self.FH)
    
        self.pdf=complete
        #self.signal_pdf = signal
        #self.background_pdf = back
        #self.angularProyection = angular
        #self.massProyection = mass
        

        
    def symParams(self):
        
        return {'AFB': self.AFB,
               'FH': self.FH,
               'LAMBDA' : self.lambdA,
                'MU' : self.mu ,
                'SIGMA' : self.sigma,
                'FRAC' : self.frac,
               }

                              
        
        
                              
                              
class angularModel(object):
    
    def __init__(self, nBin, AFB, FH, limitsAFB=[-0.1,0.1], limitsFH=[0,3]):
        
        
        self.parameters = toysParams[str(nBin)]['fixed']
        
        self.cos = symfit.Variable('cos')
            
        S, B = toysParams[str(nBin)]['yield'].values()
        if S>B : S,B = B,S
        self._frac = S/(S+B)
        self.frac = symfit.Parameter('f',self._frac, 0, 1)
        
        self._AFB = AFB
        self.AFB = symfit.Parameter('AFB', self._AFB, limitsAFB[0], limitsAFB[1])
        
        self._FH = FH
        self.FH = symfit.Parameter('FH', self._FH, limitsFH[0], limitsFH[1])
        
        self.pdf = self._pdf()
        
        
    
    def _pdf(self):
        if 'HoracioEffyCoefs' in self.parameters:
            print('Using HoracioEffyCoefs')
            coefsEff = polyToChebCoefs(self.parameters['HoracioEffyCoefs'])
        else:
            print('Using OriginalEffyCoefs')
            coefsEff = polyToChebCoefs(self.parameters['OriginalEffyCoefs'])
        return angularComplete(self.cos, self.frac, self.AFB, self.FH, coefsEff, self.parameters['frac'], self.parameters['mu'], self.parameters['sigma'], self.parameters['sidebandsCoefs'])
    
    def __call__(self, **kwargs):
        cos = kwargs.get('cos',np.linspace(-1,1,100))
        AFB = kwargs.get('AFB', self._AFB)
        FH = kwargs.get('FH', self._FH)
        f = kwargs.get('f', self._frac)
        return self.pdf(cos=cos, AFB=AFB, FH=FH, f=f)
    
    
    
    def sample(self, nEvents = 10000):
        pass
        
        
        
        
# class model(object):
    
#     def __init__(self, params, variables, pdf):
#         self.params = params
#         self.variables = variables
#         self.pdf = pdf
    
#     def __call__(self, **params):
        
        
if __name__=='__main__':
    completePDF = complete(nBin=10, AFB=0.5, FH=1.5, limitsAFB=[-0.1,0.1], limitsFH=[0,3])
    
        
        
        
