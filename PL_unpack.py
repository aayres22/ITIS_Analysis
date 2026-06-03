import pickle
import matplotlib.pyplot as plt

paramtitles = ['d1',
'k1','k2','h1','h2','h3','d2',
'k3','k4','h4','d3',
'h5','h6','k5','k6','h7','d4',
'b1','k7','h8','k8','h9','d5','h10',
'b2','k9','k10','k11','d6',
'k12','k13','k14','h11','d7',
'k15','k16','d8',
'alpha', 'k', 'beta', 'L', 'eps', 'delta', 'T', 'Nc']

from PL_theiaprep import param_global_all

D = '1'         # '1' for ACTH + Cort,
                # '2' for ACTH + Cort + TNF-a
                # '3' for ACTH + Cort + TNF-a + IL10
dpoints = '1'   # '1' for 25, '2' for 13

if D == '1':
    output_ids = [6,7]
elif D == '2':
    output_ids = [3,6,7]
else:
    output_ids = [3,4,6,7]

with open('PL_' + D + dpoints + '_nonoise.pkl', 'rb') as f:
    results = pickle.load(f)

# all_results = {
#     'param_ids' : param_ids,
#     't_data' : t_data,
#     'y_data' : y_data,
#     'true_sol' : true_sol,
#     'J_save': J_save,
#     'par_save': par_save,
#     'param_global_all': param_global_all
# }

param_ids = results['param_ids']
t_data = results['t_data']
y_data = results['y_data']
true_sol = results['true_sol']
J_save = results['J_save']
par_save = results['par_save']

for i in range(len(param_ids)):

    plt.figure()
    plt.title(paramtitles[param_ids[i]])
    plt.plot(par_save[i,i,:],J_save[i,:])
    plt.show()