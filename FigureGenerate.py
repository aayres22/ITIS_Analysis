import pickle
import numpy as np
import matplotlib.pyplot as plt

def get_nominal_param():
    d1 = 1.35e-7 # min N-unit inverse
    # Parameters for phagocytes
    k1 = 4.9956e7 #N-unit Kg / min pg
    k2 = 12.94907 # ND
    h1 = 1693.9509 #pg/mL
    h2 = 0.07212   # pg/mL
    h3 = 147.68    # pg/mL
    d2 = 0.1439   # 1/min

    # Parameters for TGF beta
    k3 = 1.546*1e-9 # mL/pg N-unit min
    k4 = 0.5     # mL/pg min
    h4 = 500#*1e6/100     # pg / mL
    d3 = 0.031777 # 1/min

    # Parameters for TNF alpha
    h5 = 550e4   # N-unit
    h6 = 0.1589  # pg/mL
    k5 = 25.5194 # pg/mL min
    k6 = 3.5514e4 # pg/mL min
    h7 = 1.5495e3 # pg/mL
    d4 = 0.0293 #0.0307/1   #mL/pg min

    # Parameters for IL-10
    b1 = 1187.2 # pg/mL min
    k7 = 267480 #  pg/mL min
    h8 = 8.0506e7 # N-unit
    k8 = 43875 #  pg/mL min
    h9 = 0.38  # pg/mL
    d5 = 98.932 # 1/min
    h10 = 791.27 #pg/mL

    # Parameters for CRH
    b2 = 0.001 # pg/mL min
    k9 = 6.84e9 # pg/ml min
    k10 = 1.7558e9#/(1e12))*(100^2) #(mL/pg^2)
    k11 = 0.0667 # 1/min
    d6 = 0.032 #1/min

    # Parameters for ACTH
    k12 = 2.3688e4 #1/min
    k13 = 1.7778e5#/1e6)*100 #mL/pg
    k14 = 112 #pg/mL min
    h11 = 80 #pg/mL
    d7  = 0.016 # 1/min

    # Parameters for Cortisol
    k15 = 4.2541e-04#*1e6/100 #pg mL / mL min pg
    k16 = 12 # pg/mL
    d8  = 0.0266 #1/min

    # Parameters for the circadian rhythm
    alpha = 300 #min
    k     = 5 #ND
    beta = 950 #min
    L = 6 #ND
    eps = 0.01 #ND
    delta = 70#76.37 #min
    T = 1440 # min
    Nc = 1.9168


    param = [d1,            #0
    k1,k2,h1,h2,h3,d2,      #1-6
    k3,k4,h4,d3,            #7-10
    h5,h6,k5,k6,h7,d4,      #11-16
    b1,k7,h8,k8,h9,d5,h10,  #17-23
    b2,k9,k10,k11,d6,       #24-28
    k12,k13,k14,h11,d7,     #29-33
    k15,k16,d8,             #34 - 36
    alpha, k, beta, L, eps, delta, T, Nc #37-43
    ]
    param_log = np.log(param)

    # Initial Conditions
    endo_0 = 0
    phag_0 = 0
    TGFB_0 = 0.004 / d3
    TNF_0 = 0
    IL10_0 = b1 * h10 / (d5 * h10 - b1)
    CRH_0 = 22
    ACTH_0 = 150
    cort_0 = 1.9

    # MJC - start this without the infection so we reach a steady state, then we can implement the infection
    IC = [0.0, 0., 0.1595967, 0., 14.68266827,
          41.49891056, 39.97751646, 11.29827208]

    return param_log, IC

param_log,IC = get_nominal_param()

pickles = ['11', '12', '21', '22', '31', '32']

for pic in pickles:
    with open('ResAnalysis' + pic + '.pkl', 'rb') as f:
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

    plt.savefig("C:\\Users\\ayres\\Desktop\\RESEARCH, GRANTS, LIFE!!!!\\ITIS Analysis 25-26\\Figures\\LSorted" + pic + ".png")
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

