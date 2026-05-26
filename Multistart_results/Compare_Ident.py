import pickle
import numpy as np
import matplotlib.pyplot as plt

# all_results = {
#     'output_ids': output_ids,
#     'param_ids' : param_ids,
#     'param_in' : param_in,
#     't_data' : t_data,
#     'y_data' : y_data,
#     'opt_model' : opt_model,
#     'true_sol' : true_sol,
#     'param_opt' : param_opt,
#     'opt_cost' : opt_cost,
#     'opt_solutions' : opt_solutions,
#     'S_opt' : S_opt,
#     'var' : var,
#     'rel_err' : rel_err,
#
# }
paramtitles = ['d1',
'k1','k2','h1','h2','h3','d2',
'k3','k4','h4','d3',
'h5','h6','k5','k6','h7','d4',
'b1','k7','h8','k8','h9','d5','h10',
'b2','k9','k10','k11','d6',
'k12','k13','k14','h11','d7',
'k15','k16','d8',
'alpha', 'k', 'beta', 'L', 'eps', 'delta', 'T', 'Nc']

D = '1'
with open('OLSMultistart_nonoise_' + D + '1' + '.pkl', 'rb') as f:
    results = pickle.load(f)

rel_err1 = results['rel_err']
param_ids1 = results['param_ids']
opt_cost1 = results['opt_cost']

with open('Multistart_results\\OLSMultistart_nonoise_' + D + '2' + '.pkl', 'rb') as f:
    results = pickle.load(f)

rel_err2 = results['rel_err']
param_ids2 = results['param_ids']
opt_cost2 = results['opt_cost']

for i in range(len(param_ids1)):
    plt.subplot(len(param_ids1),2, i+1)
    plt.scatter(rel_err1[:,i],opt_cost1,marker='o')
    plt.axvline(0,linestyle='--')
    plt.title(paramtitles[param_ids1[i]])
    plt.xlim([-0.4,0.4])

    plt.subplot(len(param_ids1),2, i + 1 + len(param_ids1))
    plt.scatter(rel_err2[:, i], opt_cost2, marker='o')
    plt.axvline(0, linestyle='--')
    plt.title(paramtitles[param_ids1[i]])
    plt.xlim([-0.4, 0.4])

plt.show()
