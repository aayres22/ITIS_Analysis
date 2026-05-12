import pickle
import numpy as np
import matplotlib.pyplot as plt
from ITIS_Model_funcs import get_nominal_param,OLS_res,ITIS

param_log,IC = get_nominal_param()

pickles = ['11', '12', '21', '22', '31', '32']

with open('OLS_Results\\ResAnalysis11.pkl', 'rb') as f:
    results = pickle.load(f)

        ##Reference
        # all_results = {
        #     'output_ids': output_ids,
        #     'param_ids' : param_ids,
        #     'param_in' : param_in,
        #     'twonorms' : twonorms,
        #     't_data' : t_data,
        #     'y_data' : y_data,
        #     'opt_model' : opt_model,
        #     'true_sol' : true_sol,
        #     'param_opt' : param_opt,
        #     'optimal_solution' :np.exp(least_sq_sol.x),
        #     'objvalue' : least_sq_sol.cost
        # }


output_ids = results['output_ids']
t_data = results['t_data']
param_opt = results['param_opt']
y_data = results['y_data']
param_ids = results['param_ids']

##Sensitvity analysis of residuals Starting from optimized param values?

##SENSITIVITY ANALYSIS
h = 1e-6  #amount to perturb parameters
n_param = len(param_opt)
n_states = len(output_ids)

Sensitivity_Mat = np.zeros((n_param, len(t_data) * n_states))  ##Initialize shape of sensitivity matrix.


##Reference
# for i in range(n_param):  #calculate the relative sensitivity to each 45 parameters
#
#         param_in = param_log[i]
#         param_delta = param_in + h ##    [HOW EXACTLY SHOULD I BE PERTURBING MY PARAM?]
#         ##THIS THING NEXT DO THIS!!!
#         Sensitivity_Mat[i, :] = ((1 / h) * (OLS_res(param_delta, y_data, t_data, i, output_ids, param_all, IC)
#                                             - OLS_res(param_in, y_data, t_data, i, output_ids, param_all, IC)))
for i in range(n_param):  #calculate the relative sensitivity to each 45 parameters

        param_in = param_opt[i]
        param_delta = param_in + h ##    [HOW EXACTLY SHOULD I BE PERTURBING MY PARAM?]
        ##THIS THING NEXT DO THIS!!!
        Sensitivity_Mat[i, :] = ((1 / h) * (OLS_res(param_delta, y_data, t_data, i, output_ids, param_log, IC)
                                            - OLS_res(param_in, y_data, t_data, i, output_ids, param_log, IC)))


print(Sensitivity_Mat)
#
# twonorms = np.zeros(n_param)
#
# for i in range(n_param):
#         twonorms[i] = np.power(sum(np.power(Sensitivity_Mat[i, :], 2)), 1 / 2)
#
