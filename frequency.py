import os
import stochpy
import random
import numpy as numpy
from scipy import stats
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--delta", help="value for delta")
parser.add_argument("--zeta", help="value for zeta")
args = parser.parse_args()
fitted_delta = float(args.delta)
fitted_zeta = float(args.zeta)

workingdir = os.getcwd()

# Simulation parameters
start_time = 0.0
end_time = 8760
n_runs = 1000
scenarios = 6

# Model output storage arrays
results = numpy.empty([n_runs*scenarios,2])

def Control(iteration,scen):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Decolonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('eta',0)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    results[iteration+(n_runs*scen),0] = 0
    results[iteration+(n_runs*scen),1] = (cases/(365*18*5))*1000
    
def OneDay(iteration,scen):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Decolonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('delta',fitted_delta)
    model.ChangeParameter('zeta',fitted_zeta)
    model.ChangeParameter('eta',1/24)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    results[iteration+(n_runs*scen),0] = 24
    results[iteration+(n_runs*scen),1] = (cases/(365*18*5))*1000

def TwoDay(iteration,scen):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Decolonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('delta',fitted_delta)
    model.ChangeParameter('zeta',fitted_zeta)
    model.ChangeParameter('eta',1/48)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    results[iteration+(n_runs*scen),0] = 48
    results[iteration+(n_runs*scen),1] = (cases/(365*18*5))*1000
    
def ThreeDay(iteration,scen):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Decolonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('delta',fitted_delta)
    model.ChangeParameter('zeta',fitted_zeta)
    model.ChangeParameter('eta',1/72)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    results[iteration+(n_runs*scen),0] = 72
    results[iteration+(n_runs*scen),1] = (cases/(365*18*5))*1000

def FourDay(iteration,scen):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Decolonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('delta',fitted_delta)
    model.ChangeParameter('zeta',fitted_zeta)
    model.ChangeParameter('eta',1/96)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    results[iteration+(n_runs*scen),0] = 96
    results[iteration+(n_runs*scen),1] = (cases/(365*18*5))*1000

def FiveDay(iteration,scen):
    model = stochpy.SSA()
    model.Model(model_file='MRSA_Decolonization.psc', dir=workingdir)
    model.Endtime(end_time)
    model.ChangeParameter('delta',fitted_delta)
    model.ChangeParameter('zeta',fitted_zeta)
    model.ChangeParameter('eta',1/120)
    model.DoStochSim()
    model.GetRegularGrid(n_samples=end_time)
    outcomes = model.data_stochsim_grid.species
    cases = outcomes[16][0][-1]
    results[iteration+(n_runs*scen),0] = 120
    results[iteration+(n_runs*scen),1] = (cases/(365*18*5))*1000
    
for i in range(0,n_runs):
	print("*** Iteration %i of %i ***" % (i+1,n_runs))
	Control(i,0)
	OneDay(i,1)
	TwoDay(i,2)
	ThreeDay(i,3)
	FourDay(i,4)
	FiveDay(i,5)

numpy.savetxt('frequency.csv',results,delimiter=','
,header="Frequency,Cases",comments='')



