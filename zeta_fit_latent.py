import os
import stochpy
import random
import numpy as numpy
from scipy import stats
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--suffix", help="suffix for files")
parser.add_argument("--delta", help="value for delta")
args = parser.parse_args()
suffix = args.suffix
fitted_delta = float(args.delta)

workingdir = os.getcwd()

# Simulation parameters
start_time = 0.0
end_time = 8760
n_runs = 1

# Run is a single run of the model that returns the number of incident cases
def Run(filename,k,parameter):
    model = stochpy.SSA()
    model.Model(model_file=filename, dir=workingdir)
    model.Endtime(end_time)
        
    latent_duration = random.randrange(24,96)
    gamma_val = 1/latent_duration
    model.ChangeParameter('gamma',gamma_val)
    
    model.ChangeParameter('delta',fitted_delta)
    model.ChangeParameter(parameter,k)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[17][0][-1]
    return cases

def Batch(filename,k,parameter):
    batch_cases = numpy.empty([n_runs,1])
    for i in range(0,n_runs):
        batch_cases[i,0] = Run(filename,k,parameter)
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
    for i in range(0,iterations):
        draw = random.uniform(priorhi,priorlow)
        results[i,0] = draw
        sim_avg = Batch(config,draw,parameter)
        results[i,1] = sim_avg
        decision = AcceptReject(objective=objective,result=sim_avg,tolerance=tol)
        results[i,2] = decision
        print("*** Iteration %i of %i ***" % (i+1,iterations))
    return results

# Fit zeta
zeta_fit = ABCFit(config='MRSA_Decolonization_Latent.psc',objective=22.36691,tol=0.05,iterations=22728,parameter='zeta',priorhi=1.0,priorlow=0.0)

zeta_estimate = zeta_fit[:,0][zeta_fit[:,2]==1]

numpy.savetxt(''.join(['zetafit_latent_',suffix,'.csv']),zeta_estimate,delimiter=',',comments=',')
