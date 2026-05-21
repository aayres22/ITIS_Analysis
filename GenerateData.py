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
IC_2[0] = 2.0
# Resolve model only over 24 hours at subsampled time points. This will be used as our data
true_sol = odeint(ITIS, IC_2, tfinal, args = (param_log,), **ode_options)


#t_data = np.linspace(0,24,24+1) # time points every hour
t_data = np.linspace(0,24,12+1) # less data, every 30 minutes


#Generate data
#output_ids = [6, 7] # Neuroendocrine = ACTH + Cort
#output_ids = [3, 6, 7] # TNF + ACTH + Cort
output_ids = [3,4,6,7] #TNF + IL10 + ACTH + Cort

# Define data
nt = len(t_data)
y_data = np.zeros((len(t_data),len(output_ids)))

# Here, we can add Gaussian, iid noise to the observations [REMOVE NOISE?]
for i, time in enumerate(t_data):
    index = np.where(tfinal==time)
    for j,out_id in enumerate(output_ids):
        # This was causing an issue - it gave back a vector of 1s
        y_data[i,j] = true_sol[index,out_id][0][0]

all_results = {
     'y_data': y_data,
     't_data': t_data,
     'tfinal': tfinal,
  }

with open('syntheticData\\TRUE_SOL32.pkl', 'wb') as f:
      pickle.dump(all_results, f)