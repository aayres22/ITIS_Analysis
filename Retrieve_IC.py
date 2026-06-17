import numpy as np
from ITIS_Model_funcs import get_nominal_param, ITIS
import pickle
from scipy.integrate import odeint

param_log,IC = get_nominal_param()

##Define time spaces
tstart = 0
t_end =24*8
dt = .05
tspace = np.arange(tstart,t_end,dt)
tfinal = np.arange(0,24+dt,dt)

##Solve the model at nominal values
ode_options = {'rtol': 1e-6}
ODE_sol = odeint(ITIS, IC, tspace, args = (param_log,), **ode_options)

# MJC - this is where we can specify that the infection begins during this final time point
IC_2    = ODE_sol[-1,:]
#print(IC_2)
D = '3'
if D == '1':
    output_ids = [6,7]
    param_ids = [33, 12, 36, 32, 43, 39, 16, 7, 31, 37, 8, 15, 9, 28, 14]
elif D == '2':
    output_ids = [3,6,7]
    param_ids = [12, 33, 16, 7, 36, 43, 15, 39, 32, 8, 31, 14, 9, 1, 34]
elif D == '3':
    output_ids = [3,4,6,7]
    param_ids =[12, 33, 7, 36, 16, 43, 15, 39, 32, 8, 9, 31, 21, 14, 10, 34]

with open('paramRankings\\rankingDesign' + D + '.pkl', 'rb') as f:
    results = pickle.load(f)

rank_value = results['rank_value']      # sorted ranking values
param_sorted_all = results['param_sorted']  #accordingly sorted params as strings
rank_cut = np.array(rank_value)[np.array(rank_value) > .25 * np.array(rank_value)[0]]
param_sorted = np.array(param_sorted_all)[:len(rank_cut)]
print(param_sorted)

param_titles = ['d1',
'k1','k2','h1','h2','h3','d2',
'k3','k4','h4','d3',
'h5','h6','k5','k6','h7','d4',
'b1','k7','h8','k8','h9','d5','h10',
'b2','k9','k10','k11','d6',
'k12','k13','k14','h11','d7',
'k15','k16','d8',
'alpha', 'k', 'beta', 'L', 'eps', 'delta', 'T', 'Nc']

print([param_titles.index(param_sorted[i]) for i in range(len(param_sorted))])