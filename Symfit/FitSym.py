import argparse


parser = argparse.ArgumentParser(prog='Measurement with symfit',
            description='Fit a given bin with data obtained from Params.json and data', 
            epilog='Hope this works :V')

parser.add_argument('--bin',
                    action='store',
                   type=int,
                    required=True,
                   help='Bin number!')

parser.add_argument('--events',
                  type=int,
                  help="Number of times a new fit is produced",
                  default=1
                  )

parser.add_argument('--iterations',
                   type=int, 
                   help='Number of iterations of the fit default to 2',
                   default=1
                  )

parser.add_argument('--debug',
                   type=bool,
                   default=False
                  )

args = parser.parse_args()


import numpy as np
import matplotlib.pyplot as plt

import scipy.integrate as integrate
from time import time
import pandas as pd
import json
import pdb

import mplhep as hep
plt.style.use(hep.styles.CMS)
from matplotlib.backends.backend_pdf import PdfPages

from symfit import Fit
from symfit.core.objectives import LogLikelihood
from symfit.core.minimizers import  BasinHopping, SLSQP
import sympy

import sys
sys.path.append('/home/horacio/Documents/hcl/')
import completeSymfit as Complete

import os
from copy import deepcopy
from random import random


data_dir = '/home/horacio/Documents/hcl/DataSelection/Data/AntiRad_Resonance/'
dataDF = pd.read_pickle(data_dir+f'df_Bin{args.bin}OnlyResonanceRej.pkl')
antiRad_JP_1 = (np.abs((dataDF.BMass - 5.27934)- (dataDF.DiMuMass-3.0969)) > 0.130) | (dataDF.DiMuMass>3.0969)
antiRad_JP_2 = ((np.abs((dataDF.BMass - 5.27934)- (dataDF.DiMuMass-3.0969)) > 0.127) | (dataDF.DiMuMass > 3.35))  | (dataDF.DiMuMass<3.0969)
antiRad_PP_1 = (np.abs((dataDF.BMass - 5.27934)- (dataDF.DiMuMass-3.6861)) > 0.085) | (dataDF.DiMuMass>3.6861)
antiRad_PP_2 = ((np.abs((dataDF.BMass - 5.27934)- (dataDF.DiMuMass-3.6861)) > 0.054) | (dataDF.DiMuMass > 3.92))  | (dataDF.DiMuMass<3.6861)
dataDF = dataDF[antiRad_JP_1 & antiRad_JP_2 &antiRad_PP_1 &antiRad_PP_2]
dataDF = dataDF[(dataDF.R_XGB>0.96) & (dataDF.L_XGB>0.98)]
dataDF = dataDF[(dataDF.BMass>5) & (dataDF.BMass<5.7)]


parameter_path = f'/home/horacio/Documents/hcl/Experiments_Fitting/CompleteFit/Symfit/Params/Bin{args.bin}.json'
with open(parameter_path) as file:
    parameters = json.load(file)

if args.debug:
    pdb.set_trace()

#RENAME THE KEY FOR THE EFFICIENCY COEFFICIENTS
parameters['Signal']['angle']['coefsEff'] = parameters['Signal']['angle']['Chi2EffyCoefs'] 
    
def create_models(afb, fh, bin_data) :   

    free = {'Signal:angle:AFB':[afb,-1.5,1.5], 'Signal:angle:FH':[fh,0,3],
        'Background:angle:2SideBands:fraction_Left':[0.5, 0, 1]}
    print(free)
    model = Complete.complete_v1(args.bin, params = parameters, free_Params = free)
    
    AFB, FH = model.symParams('AFB'), model.symParams('FH')
    _constraints = [
                sympy.LessThan(FH, 3),
                sympy.GreaterThan(FH, 0),
                sympy.LessThan(sympy.Abs(AFB), 0.5*FH)
                ]


    fitConstraintSLSQP = Fit(model.pdf, 
                                cos=np.array(bin_data.cosThetaKMu), 
                                mass=np.array(bin_data.BMass), 
                                objective=LogLikelihood, 
                                constraints=_constraints,
                                minimizer=SLSQP)


    fitConstraintBH_SLSQP = Fit(model.pdf, 
                                cos=np.array(bin_data.cosThetaKMu), 
                                mass=np.array(bin_data.BMass), 
                                objective=LogLikelihood, 
                                constraints=_constraints,
                                minimizer=[BasinHopping, SLSQP])
    
    return fitConstraintSLSQP, fitConstraintBH_SLSQP, model



def update_to_save(result, diction, model, j):
    for k in result.params:
        
        if not k+str(j) in diction: diction[k+str(j)]=list()
        if not 'err'+k+str(j) in diction: diction['err'+k+str(j)]=list()
            
        diction[k+str(j)].append(result.params[k])
        diction['err'+k+str(j)].append(result.stdev(model.symParams(k)))
    
    if not 'Status'+str(j) in diction:   diction['Status'+str(j)]=list()
    diction['Status'+str(j)].append(result.status_message)
        
    if not 'NLL'+str(j) in diction:   diction['NLL'+str(j)]=list()
    diction['NLL'+str(j)].append(result.minimizer_output['fun'])

    
# Used to produce the pandas dataframe
dictSLSQP = dict()
path = sys.argv[0].split('/')[:-1]
print(path)
path = '/'.join(path)
path += 'Results/'
print(path)


for i in range(args.events):
    
    fh_ = [round(random()*1.5, 3)]
    afb_ = [round(random()*fh_[0] - fh_[0]/2,4)]
    
    for j in range(args.iterations):
        
        fitConstraintSLSQP, fitConstraintBH_SLSQP, model = create_models(afb_[-1], fh_[-1], dataDF)
        resultSLSQP = fitConstraintSLSQP.execute()

        #Append the result of fit to use it as input for the next iteration
        if np.abs(resultSLSQP.params['AFB'])<0.5*resultSLSQP.params['FH']:
            afb_.append(resultSLSQP.params['AFB'])
        else:
            # When the fitted value result on the boundary, 
            # the precision make it look like is "just" outside 
            # the boundary
            # Here we try to correct this issue
            sign = resultSLSQP.params['AFB']/np.abs(resultSLSQP.params['AFB'])
            afb_.append(sign*resultSLSQP.params['FH']/2)
        fh_.append(resultSLSQP.params['FH'])
        
        update_to_save(resultSLSQP, dictSLSQP, model, j )

    if 'afb_ini' not in dictSLSQP: dictSLSQP['afb_ini']=list()
    dictSLSQP['afb_ini'].append(afb_[0])

    if 'fh_ini' not in dictSLSQP: dictSLSQP['fh_ini']=list()
    dictSLSQP['fh_ini'].append(fh_[0])
    
    print(dictSLSQP)
    
    for method in ['SLSQP', 'BH+SLSQP']:
        if method == 'SLSQP': ress = dictSLSQP
        else: continue #ress = dictBH_SLSQP
        with PdfPages(path+method+'.pdf') as pdf:
            for k in resultSLSQP.params:
                plt.figure()
                plt.hist(ress[k+str(j)])
                plt.xlabel(k)
                plt.title('Last value : '+ method+'_'+sys.argv[1])
                #plt.savefig(path+k+'_'+method+'.png')
                pdf.savefig()
                plt.close()
    
    plt.figure()
    x_ = np.linspace(-1.5,1.5, 3000)
    plt.plot(x_, np.abs(x_)*2, ls='--', color='grey', linewidth=3)
    plt.scatter(dictSLSQP['afb_ini'], dictSLSQP['fh_ini'])
    plt.savefig(path+'values.png')
    plt.close()
    
    

    pd.to_pickle(pd.DataFrame.from_dict(dictSLSQP), path+'/SLSQP.pkl')
    
    
