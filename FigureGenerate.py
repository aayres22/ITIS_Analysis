import pickle
import numpy as np
import matplotlib.pyplot as plt
from ITIS_Model_funcs import get_nominal_param

param_log,IC = get_nominal_param()

pickles = ['11', '12', '21', '22', '31', '32']

for pic in pickles:
    with open('OLS_Results\\ResAnalysis' + pic + '.pkl', 'rb') as f:
        results = pickle.load(f)


    ##Retrieve inferred parameters
    paramtitles = ['d1',
                   'k1', 'k2', 'h1', 'h2', 'h3', 'd2',
                   'k3', 'k4', 'h4', 'd3',
                   'h5', 'h6', 'k5', 'k6', 'h7', 'd4',
                   'b1', 'k7', 'h8', 'k8', 'h9', 'd5', 'h10',
                   'b2', 'k9', 'k10', 'k11', 'd6',
                   'k12', 'k13', 'k14', 'h11', 'd7',
                   'k15', 'k16', 'd8',
                   'alpha', 'k', 'beta', 'L', 'eps', 'delta', 'T', 'Nc']

    param_ids = results['param_ids']
    inferred = [paramtitles[i] for i in param_ids]
    print(inferred)


    # print(results)

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

    ##Retrieve residuals

    # optimal_sol = results['optimal_solution']
    # true_sol = results['true_sol']
    # objvalue = results['objvalue']
    # #param_in = results['param_in'] ##this recorded param_in from sensitivity analysis, useless :/
    # param_ids = results['param_ids']
    #
    # print("optimal solution",optimal_sol)
    # print("true solution",np.exp(param_log[param_ids]))
    # print("objective function value",objvalue)


    # dt = .05
    # tfinal = np.arange(0,24+dt,dt)
    #
    ##Unpack
    param_opt = np.exp(results['param_opt'][results['param_ids']])
    param_true = np.exp(param_log[results['param_ids']])

    print(pic + " Residuals = ", abs(param_true - param_opt))

    t_data = results['t_data']
    y_data = results['y_data']
    output_ids = results['output_ids']
    opt_model = results['opt_model']
    true_sol = results['true_sol']
    twonorms = results['twonorms']

    ## ROUTINE FOR LOCAL SENSITIVITY ANALYSIS
    fig, ax = plt.subplots(figsize = (12,5))
    data = sorted(zip(paramtitles, twonorms), key=lambda v: v[1], reverse = True)
    for (i, c) in data:  # unpack and plot each tuple in sorted order
        bars = plt.bar(i, c, align='center', color = '#1f77b4')
    plt.xlim(-0.9, len(paramtitles) - 1 + 0.9)

    #plt.savefig("C:\\Users\\ayres\\Desktop\\RESEARCH, GRANTS, LIFE!!!!\\ITIS Analysis 25-26\\Figures\\LSorted" + pic + ".png")
    fig.tight_layout()
    plt.show()



### ROUTINE FOR OPTIMAL VS TRUE MODEL FITS
#
# fig, ax = plt.subplots(4,1, figsize = (6,11))
# fig.supylabel("Concentration")
#
# ax[0].scatter(t_data,y_data[:,0],marker='x',color='k', label='data')
# ax[0].plot(tfinal,opt_model[:,output_ids[0]],label='opt fit')
# ax[0].plot(tfinal,true_sol[:,output_ids[0]],label='true model')
# ax[0].set_title('TNF-\u03B1')
# #ax[0].legend()
#
# ax[1].scatter(t_data,y_data[:,1],marker='x',color='k', label='data')
# ax[1].plot(tfinal,opt_model[:,output_ids[1]],label='opt fit')
# ax[1].plot(tfinal,true_sol[:,output_ids[1]],label='true model')
# ax[1].set_title('IL10')
#
# ax[2].scatter(t_data,y_data[:,2],marker='x',color='k', label='data')
# ax[2].plot(tfinal,opt_model[:,output_ids[2]],label='opt fit')
# ax[2].plot(tfinal,true_sol[:,output_ids[2]],label='true model')
# ax[2].set_title('ACTH')
#
# ax[3].scatter(t_data,y_data[:,3],marker='x',color='k', label='data')
# ax[3].plot(tfinal,opt_model[:,output_ids[3]],label='opt fit')
# ax[3].plot(tfinal,true_sol[:,output_ids[3]],label='true model')
# ax[3].set_title('Cortisol')
#
#
# ax[3].set_xlabel("Time (hrs)", size = 12)
#
# fig.tight_layout()
# plt.savefig("C:\\Users\\ayres\\Desktop\\RESEARCH, GRANTS, LIFE!!!!\\ITIS Analysis 25-26\\Figures\\32.png")
# plt.show()

