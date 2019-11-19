import os
import stochpy
import random
import numpy as numpy
from scipy import stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--suffix", help="suffix for files")
args = parser.parse_args()
suffix = args.suffix

workingdir = os.getcwd()

# Simulation parameters
start_time = 0.0
end_time = 8760
n_runs = 1

# Run is a single run of the model that returns the number of incident cases
def SensRun(filename,k,parameter,pdict):
    model = stochpy.SSA()
    model.Model(model_file=filename, dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('rho_N',pdict['rho_n']*3.973)
    model.ChangeParameter('rho_D',pdict['rho_d']*0.181)
    model.ChangeParameter('sigma',pdict['sigma']*0.054)
    model.ChangeParameter('psi',pdict['psi']*0.0815306)
    model.ChangeParameter('theta',pdict['theta']*0.00949)
    model.ChangeParameter('nu',pdict['nu']*0.0779)
    model.ChangeParameter('iota_N',pdict['iota_n']*6.404)
    model.ChangeParameter('iota_D',pdict['iota_d']*1.748)
    model.ChangeParameter('tau_N',pdict['tau_n']*2.728)
    model.ChangeParameter('tau_D',pdict['tau_d']*0.744)
    model.ChangeParameter('mu',pdict['mu']*0.002083)
    model.ChangeParameter(parameter,k)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    return cases

def Batch(filename,k,parameter,pdict):
    batch_cases = numpy.empty([n_runs,1])
    for i in range(0,n_runs):
        batch_cases[i,0] = SensRun(filename,k,parameter,pdict)
        batch_avg = numpy.mean(batch_cases)
    return batch_avg

def AcceptReject(objective,result,tolerance):
    rate = result/(365*18)
    lnrate = numpy.log(rate)
    lnobjective = numpy.log(objective/(365*18))
    distance = abs(lnobjective - lnrate)
    if distance <= tolerance:
        accept = 1
    else:
        accept = 0
    return accept

def ABCFit(config,objective,tol,iterations,parameter,priorhi,priorlow):
    results = numpy.zeros([iterations,3])
    
    pdict = {'rho_n':random.uniform(0.5,1.5),
            'rho_d':random.uniform(0.5,1.5),
             'sigma':random.uniform(0.5,1.5),
             'psi':random.uniform(0.5,1.5),
             'theta':random.uniform(0.5,1.5),
             'nu':random.uniform(0.5,1.5),
             'iota_n':random.uniform(0.5,1.5),
             'iota_d':random.uniform(0.5,1.5),
             'tau_n':random.uniform(0.5,1.5),
             'tau_d':random.uniform(0.5,1.5),
             'mu':random.uniform(0.5,1.5)
            }   
    
    for i in range(0,iterations):
        draw = random.uniform(priorhi,priorlow)
        results[i,0] = draw
        sim_avg = Batch(config,draw,parameter,pdict)
        results[i,1] = sim_avg
        decision = AcceptReject(objective=objective,result=sim_avg,tolerance=tol)
        results[i,2] = decision
    if len(results[:,0][results[:,2]==1])==0:
        median_val = 0
    else: 
        median_val = numpy.median(results[:,0][results[:,2]==1])
    return median_val,pdict

def GlobalSens(iterations):
    parameters = numpy.zeros([iterations,11])
    for k in range(0,iterations):
        fit = ABCFit(config='MRSA_Decol_Sens.psc',objective=22.36691,tol=0.05,iterations=200,parameter='delta',priorhi=1.0,priorlow=0.0)
        parameters[k,0] = fit[0]
        parameters[k,1] = fit[1]['rho_n']
        parameters[k,2] = fit[1]['rho_d']
        parameters[k,3] = fit[1]['psi']
        parameters[k,4] = fit[1]['theta']
        parameters[k,5] = fit[1]['nu']
        parameters[k,6] = fit[1]['iota_n']
        parameters[k,7] = fit[1]['iota_d']
        parameters[k,8] = fit[1]['tau_n']
        parameters[k,9] = fit[1]['tau_d']
        parameters[k,10] = fit[1]['mu']
        print("*** Iteration %i of %i ***" % (k+1,iterations))
    return parameters

# Fit delta
global_sweep = GlobalSens(n_runs)

numpy.savetxt(''.join(['global_sweep_',suffix,'.csv']),global_sweep,delimiter=',',comments=',')
