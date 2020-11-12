import symfit
import sympy
from sympy.core.numbers import Infinity as inf
import json
from copy import deepcopy
import numpy as np

from GeneralPDFSymfit import  *
import Signals_Symfit as Signals
import Backgrounds_Symfit as Backgrounds
import projectionsSymfit as projections


# Inital version of the complete pdf (deprecated)
# ~Jun 2020
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
def Complete(mass, cos, AFB, FH, frac, coefsEff, mu, sigma,fraction, muB, sigmaB, coefs, lambda_, limitsAng, limitsMass):
    

    signal = Signals.completeSignal(mass, cos, AFB, FH, coefsEff, mu, sigma, limitsMass)
    back = completeBackground(mass, cos, fraction, muB, sigmaB, coefs, lambda_, limitsMass)
    complete = frac*signal+(1-frac)*back
    return complete









### THIS CLASSES TRY TO MIMIC (A LITTLE BIT THO) THE BEHAVIOUR OF ZFIT (AND THUS ROOFIT?)
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
    
    
    def create_pdf(self):
        fix = self.parameters['fixed']
    
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
        

        
    def symParams(self):
        
        return {'AFB': self.AFB,
               'FH': self.FH,
               'LAMBDA' : self.lambdA,
                'MU' : self.mu ,
                'SIGMA' : self.sigma,
                'FRAC' : self.frac,
               }

                              
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
# Nov 2020
# Updated Verison of the complete PDF
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

def Complete_V1(mass, cos, frac, masslimits, dictionary_params):

    signalParams = dictionary_params['Signal']
    backParams = dictionary_params['Background']
        
    signal = Signals.completeSignalCrystalBall_2Gauss(mass,cos,
                                            **signalParams['angle'],
                                            **signalParams['mass']['2Gaussian+CrystalBall'],
                                             limits = masslimits)
    back = Backgrounds.completeBackground_SB_GausExp(mass, cos,
                                     **backParams['mass']['Exponential+Gauss'],
                                     **backParams['angle']['2SideBands'])
    
    complete = frac*signal+(1-frac)*back
    return complete


def massProj_warpper(mass, cos, masslimits, dictionary_params):
        pass
        
        


# V1 OF THE COMPLETE PDF!
### THIS CLASSES TRY TO MIMIC (A LITTLE BIT THOU) THE BEHAVIOUR OF ZFIT (AND THUS ROOFIT?)
class complete_v1(object):
    """
    TODO : We need a method to evaluate the model only giving cos and mass
    """
    def __init__(self, 
                nBin:"Bin number (0-10) Bins 3 and 5 are resonances, so, restricted", 
                params:"Dictionary with all params needed to define the complete PDF",
                limitsAFB=None, 
                limitsFH=None, 
                fraction_Floating = True,
                free_Params:"Which Params will be set to float" = {'Signal:angle:AFB':[0,-3,3], 'Signal:angle:FH':[0.2,0,3]}):
        
        

        ## OBTAIN A DICTIONARY WITH FIXED PARAMETERS REPLACING FREE PARAMS WITH syfit.Parameter free_Params dict
        correctParams = deepcopy(params)
        for name in free_Params:
            path_ = name.split(':')
            if len(path_) == 5:
                correctParams[path_[0]][path_[1]][path_[2]][path_[3]][path_[4]] = symfit.Parameter(path_[4], *free_Params[name]) 
            elif len(path_) == 4:
                correctParams[path_[0]][path_[1]][path_[2]][path_[3]] = symfit.Parameter(path_[3], *free_Params[name]) 
            elif len(path_) == 3:
                correctParams[path_[0]][path_[1]][path_[2]] = symfit.Parameter(path_[2], *free_Params[name]) 
            elif len(path_) == 2:
                print('WARNING!!!  you want to remove :' , name)
                correctParams[path_[0]][path_[1]] = symfit.Parameter(path_[1], *free_Params[name])
            else:
                raise NotImplementedError('Please double check your json or edit the script')
        
        self.Params = correctParams
        self.FreeParams = [k for k in free_Params]
        
        S, B = params['Signal']['yield'], params['Background']['yield']
        if fraction_Floating:
            self.fraction = symfit.Parameter('Frac', 0.5, 0, 1)
            self.FreeParams+=['Frac']
        else:
            self.fraction = S/(S+B)
        
        #self.limits = params['limits']
        self.cos = symfit.Variable('cos')
        self.mass = symfit.Variable('mass')
    
        self.create_pdf()
            
            
    
    
    def create_pdf(self):        
        complete = Complete_V1(self.mass, self.cos, self.fraction, self.Params['limits'], self.Params)
        self.pdf=complete
        

        projection_mass= projections.massProy_v1(self.mass, self.fraction, 
                                          limits = self.Params['limits'],
                                          **self.Params['Background']['mass']['Exponential+Gauss'], 
                                          **self.Params['Signal']['mass']['2Gaussian+CrystalBall'])

        self.projection_mass = projection_mass
        self.mass_signal = self.fraction*Signals.CrystalBall_2Gauss(self.mass,
                                        **self.Params['Signal']['mass']['2Gaussian+CrystalBall'], 
                                        limits = self.Params['limits'])
        self.mass_background = (1-self.fraction)*Backgrounds.Gauss_Exp(self.mass,
                                        **self.Params['Background']['mass']['Exponential+Gauss'], 
                                        limits = self.Params['limits'])
               
        
        
        
        self.projection_angular = projections.angularProy_v1(self.cos, self.fraction,  
                                                            **self.Params['Signal']['angle'],  
                                                            **self.Params['Background']['angle']['2SideBands'])
        self.angular_signal = self.fraction*Signals.angularSignal(self.cos, **self.Params['Signal']['angle'])
        self.angular_back = (1-self.fraction)*Backgrounds.angularBackground_SideBands(self.cos, **self.Params['Background']['angle']['2SideBands'])
        
        
        
    
    
    def symParams(self, param=None):
        
        symP  = dict()
        
        for k in self.FreeParams:
            path_ = k.split(':')
                        
            if len(path_) == 5:
                param_ = self.Params[path_[0]][path_[1]][path_[2]][path_[3]][path_[4]]
                symP[param_.name] = param_
            elif len(path_) == 4:
                param_ = self.Params[path_[0]][path_[1]][path_[2]][path_[3]]
                symP[param_.name] = param_
            elif len(path_) == 3:
                param_ = self.Params[path_[0]][path_[1]][path_[2]]
                symP[param_.name] = param_
            elif len(path_) == 2:
                param_ = self.Params[path_[0]][path_[1]]
                symP[param_.name] = param_
            elif len(path_)==1 and path_[0]=='Frac':
                symP['Frac'] = self.fraction
            else:
                raise NotImplementedError('Please double check your json or edit the script')

        
        if param:
            return symP[param]
        return symP
        
        
