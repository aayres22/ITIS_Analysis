#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pickle
from scipy.optimize import least_squares

D = '1'         # '1' for ACTH + Cort,
                # '2' for ACTH + Cort + TNF-a
                # '3' for ACTH + Cort + TNF-a + IL10
dpoints = '1'   # '1' for 25, '2' for 13

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
    beta  = 950 #min
    L     = 6 #ND
    eps   = 0.01 #ND
    delta = 70#76.37 #min
    T     = 1440 # min
    Nc    = 1.9168


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
    #IC = [ 0.,          0.,          0.14759607,  0.,         13.75484703,  3.13344442, 15.73028492,  1.38649737] ##Pre-inj VALUES AT t = 24*8
    IC =  [endo_0, phag_0, TGFB_0, TNF_0, IL10_0, CRH_0, ACTH_0, cort_0]
    return param_log, IC

def ITIS(y,t,param):
 endo,phag,TGFB,TNF,IL10,CRH,ACTH,cort = y

 (d1,
k1,k2,h1,h2,h3,d2,
k3,k4,h4,d3,
h5,h6,k5,k6,h7,d4,
b1,k7,h8,k8,h9,d5,h10,
b2,k9,k10,k11,d6,
k12,k13,k14,h11,d7,
k15,k16,d8,
alpha, k, beta, L, eps, delta, T, Nc) = np.exp(param)

#Circadian Rhythm
 tm = (t * 60 - delta) % T
 R_1 = tm**k / (tm**k + alpha**k)
 R_2 = ((T - tm)**L) / ((T - tm)**L + beta**L)
 R_func = Nc*(R_1*R_2 + eps)

#Other things
 phag1 = 1 + k2 * (TNF / (h1 + TNF))
 phag2 = (h2 / (h2+TGFB)) * (h3 / (h3 + IL10))

 tnf1 = (phag / (h5 + phag)) * (h6**4 / (h6**4 + TGFB**4))
 tnf2 = k5 + k6 * (TNF / (h7 + TNF))

 il10_1 = k7 * (phag**3 / (h8**3 + phag**3)) + k8 * (TGFB**6 / (h9**6 + TGFB**6))

 crh1 = R_func * k9 * (CRH / (1+k10*cort**2)) ##[RESOLVE]

#Differential equations
 dendo = -d1 * endo * phag
 dphag = k1 * (phag1 * phag2) * endo - d2 * phag
 dTGFB = k3 * phag + k4 * cort / (h4 + cort) - d3 * TGFB
 dTNF = tnf1 * tnf2 - d4 * (TNF**2)
 dIL10 = b1 + il10_1 - d5 * (h10 / (h10 + IL10)) * IL10
 dCRH = b2 * 60 + crh1 * 60 + k11 * TNF - d6 * CRH * 60
 dACTH = k12 * (CRH / (1 + k13 * cort)) * 60 + k14 * (TNF**2 / (h11**2 + TNF**2)) - d7 * ACTH * 60
 dcort = k15 * (ACTH **2 / (1 + k16 * TGFB)) * 60 - d8 * cort * 60

 dy = [dendo,dphag,dTGFB,dTNF,dIL10,dCRH,dACTH,dcort]
 return dy

#%%
def OLS_res(param_in,y_data,t_data,param_ids,output_ids,param_fix,IC):
    # First, we may only want to update certain parameters. This will be defined in param_ids, which should include
    # the indices for which parameters are updated. Param_fix will have the "fixed" values
    param_eval = param_fix.copy()
    param_eval[param_ids] = param_in.copy()
    ##Define time spaces
    # tstart = 0
    # t_end =24*8
    dt = .05
    #tspace = np.arange(tstart,t_end,dt) # The time frame to solve the model to reach steady state (before we grab any output)
    tfinal = np.arange(0,24+dt,dt)

    ##Solve the model at nominal values
    ode_options = {'rtol': 1e-6} # sets the numerical accuracy of the ODE solver (i.e., how fine of a time step do we need)
    #ODE_sol = odeint(ITIS, IC, tspace, args = (param_eval,), **ode_options)

    # MJC - this is where we can specify that the infection begins during this final time point
    #IC_2    = ODE_sol[-1,:]
    IC[0] = 2.0
    # Resolve model only over 24 hours at subsampled time points. This will be used as our data
    ODE_sol = odeint(ITIS, IC, tfinal, args = (param_eval,), **ode_options)

    # In a similar way, the data may be only a subset of the observations. We will assume that
    # we have the same number of time points for now, but can modify this as well
    model_output = np.zeros((len(t_data),len(output_ids)))

    for i, time in enumerate(t_data):
        t_index = np.where(tfinal==time) #Grabs the index correpsonding to where we have data
        for j, mod_out in enumerate(output_ids): #Loop over which outputs we want to compare to data
            model_output[i,j] = ODE_sol[t_index,mod_out][0][0] ##THIS IS BEHAVING WEIRDLY

    # Now, we have model outputs at the time points and recorded states. We can construct the residual and then the sum of squares
    # in the future, we can look at applying weights to the data that correspond to measurement noise
    # Here we just scale by the data
    # MJC - Since we pass in the data and then take the mean here, we don't need to do it at the beginning
    residual = (y_data - model_output)/np.mean(y_data,axis=0)

    return residual.flatten()

import numpy as np
from scipy.optimize import least_squares


def calc_PL_fixpar(
    f,
    params,
    bounds,
    par_fix,
    N=21,
    profile_width=0.5,
    max_nfev=300,
    finite_diff_step=1e-4,
    verbose=False,
):
    """
    Profile likelihood calculation with one globally fixed parameter.

    Parameters
    ----------
    f : callable
        residual function. Should accept a 1D parameter vector and return model output.
        Example: yhat = f(params)-data or yhat = f(params)-data / data



    params : array_like
        Initial parameter guess, shape (num_params,).

    bounds : array_like
        Parameter bounds, shape (num_params, 2).
        First column is lower bound, second column is upper bound.


    par_fix : int
        Index of the parameter to exclude from profiling.
        Uses Python 0-based indexing.

    N : int, optional
        Number of profile points per profiled parameter. Default is 201.

    profile_width : float, optional
        Profiles each parameter over:
            [(1 - profile_width) * p, (1 + profile_width) * p]
        clipped to bounds.
        Default 0.5 gives the MATLAB range [0.5*p, 1.5*p].

    max_nfev : int, optional
        Maximum function evaluations for each least-squares solve.

    finite_diff_step : float, optional
        Finite difference step size.

    verbose : bool, optional
        Whether to print progress.

    Returns
    -------
    J_save : ndarray
        Sum of squared scaled residuals for each profile.
        Shape: (num_profiled_params, N)

    par_save : ndarray
        Saved parameter values for each profile.
        Shape: (num_profiled_params, num_profiled_params, N)

        Important: this stores only the profiled subset, excluding `par_fix`,
        just like the MATLAB code after `params = xall(par_ids)`.

    param_global_all : ndarray
        Globally optimized full parameter vector.
    """

    params = np.asarray(params, dtype=float).ravel()
    bounds = np.asarray(bounds, dtype=float)

    if bounds.shape != (params.size, 2):
        raise ValueError("bounds must have shape (num_params, 2).")

    num_par = params.size

    # if not (0 <= par_fix < num_par):
    #     raise ValueError("par_fix must be a valid 0-based parameter index.")

    if N < 3 or N % 2 == 0:
        raise ValueError("N should be an odd integer >= 3 so there is a midpoint.")

    lower_all = bounds[:, 0]
    upper_all = bounds[:, 1]

    par_ids = np.arange(num_par)
    par_ids = np.delete(par_ids, par_fix)

    def eval_PL(qest, qfix=None, par_i=None, par_ids=None, param_global_all=None):
        """
        Residual function for least_squares.

        If par_i is None:
            qest is the full parameter vector.

        Otherwise:
            qest contains the free parameters within the profiled subset,
            qfix is the fixed profile value,
            par_i is the index inside the profiled subset,
            par_ids maps profiled-subset indices back to full-parameter indices.
        """

        if par_i is None:
            q_eval = np.asarray(qest, dtype=float)
        else:
            q_eval = np.array(param_global_all, dtype=float, copy=True)

            profiled_full_id = par_ids[par_i]
            free_profile_ids = np.delete(np.arange(len(par_ids)), par_i)
            free_full_ids = par_ids[free_profile_ids]

            q_eval[free_full_ids] = qest
            q_eval[profiled_full_id] = qfix

        residual = np.asarray(f(q_eval), dtype=float).ravel()

        return residual

    # ------------------------------------------------------------------
    # 1. Global optimization over all parameters
    # ------------------------------------------------------------------
    global_res = least_squares(
        fun=lambda q: eval_PL(q),
        x0=params,
        bounds=(lower_all, upper_all),
        method="trf",  # supports bounds; SciPy's lm does not support bounds
        diff_step=finite_diff_step,
        max_nfev=max_nfev,
    )

    param_global_all = global_res.x.copy()

    # Profile only the non-fixed parameters.
    profile_params = param_global_all[par_ids]
    lower = lower_all[par_ids]
    upper = upper_all[par_ids]

    num_profile = len(profile_params)
    Nhalf = N // 2

    # ------------------------------------------------------------------
    # 2. Build profile ranges
    # ------------------------------------------------------------------
    temp_max = upper#np.minimum((1.0 + profile_width) * profile_params, upper)
    temp_min = lower#np.maximum((1.0 - profile_width) * profile_params, lower)

    PL_max = np.maximum(temp_max, temp_min)
    PL_min = np.minimum(temp_max, temp_min)

    J_save = np.full((num_profile, N), np.nan)
    par_save = np.full((num_profile, num_profile, N), np.nan)

    # ------------------------------------------------------------------
    # 3. Profile each non-fixed parameter
    # ------------------------------------------------------------------
    for par_i in range(num_profile):
        if verbose:
            print(f"Profiling parameter {par_i + 1} of {num_profile}")

        par_space = np.linspace(PL_min[par_i], PL_max[par_i], N)

        free_mask = np.ones(num_profile, dtype=bool)
        free_mask[par_i] = False

        qest0 = profile_params[free_mask]
        lower_i = lower[free_mask]
        upper_i = upper[free_mask]

        def run_profile_fit(qfix, q_start):
            res = least_squares(
                fun=lambda q: eval_PL(
                    qest=q,
                    qfix=qfix,
                    par_i=par_i,
                    par_ids=par_ids,
                    param_global_all=param_global_all,
                ),
                x0=q_start,
                bounds=(lower_i, upper_i),
                method="trf",
                diff_step=finite_diff_step,
                max_nfev=max_nfev,
            )

            J = np.sum(res.fun ** 2)
            return res.x, J

        # --------------------------------------------------------------
        # Midpoint optimization
        # --------------------------------------------------------------
        qfix_mid = par_space[Nhalf]

        try:
            x_mid, J_mid = run_profile_fit(qfix_mid, qest0)

            par_save[par_i, free_mask, Nhalf] = x_mid
            par_save[par_i, par_i, Nhalf] = qfix_mid
            J_save[par_i, Nhalf] = J_mid

            qest_half = x_mid.copy()

        except Exception as exc:
            if verbose:
                print(f"Midpoint failed for parameter {par_i}: {exc}")

            # Fall back to the initial guess so the rest of the profile can still run.
            qest_half = qest0.copy()

        # --------------------------------------------------------------
        # Move right from midpoint
        # --------------------------------------------------------------
        qest = qest_half.copy()

        for N_i in range(Nhalf + 1, N):
            qfix = par_space[N_i]

            if verbose:
                print(f"  par_i={par_i}, N_i={N_i}, qfix={qfix}")

            try:
                x, J = run_profile_fit(qfix, qest)

                par_save[par_i, free_mask, N_i] = x
                par_save[par_i, par_i, N_i] = qfix
                J_save[par_i, N_i] = J

                qest = x.copy()

            except Exception as exc:
                if verbose:
                    print(f"  Right profile failed at N_i={N_i}: {exc}")

        # --------------------------------------------------------------
        # Move left from midpoint
        # --------------------------------------------------------------
        qest = qest_half.copy()

        for N_i in range(Nhalf - 1, -1, -1):
            qfix = par_space[N_i]

            if verbose:
                print(f"  par_i={par_i}, N_i={N_i}, qfix={qfix}")

            try:
                x, J = run_profile_fit(qfix, qest)

                par_save[par_i, free_mask, N_i] = x
                par_save[par_i, par_i, N_i] = qfix
                J_save[par_i, N_i] = J

                qest = x.copy()

            except Exception as exc:
                if verbose:
                    print(f"  Left profile failed at N_i={N_i}: {exc}")

    return J_save, par_save, param_global_all
#%%
paramtitles = ['d1',
'k1','k2','h1','h2','h3','d2',
'k3','k4','h4','d3',
'h5','h6','k5','k6','h7','d4',
'b1','k7','h8','k8','h9','d5','h10',
'b2','k9','k10','k11','d6',
'k12','k13','k14','h11','d7',
'k15','k16','d8',
'alpha', 'k', 'beta', 'L', 'eps', 'delta', 'T', 'Nc']

stateTitles = ['Endotoxin', 'Phagocytes', 'TGF-\u03b2','TNF-\u03B1','IL10','CRH','ACTH','Cortisol']
#%%
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


# t_data = np.linspace(0,24,24+1) # time points every hour
# t_data = np.linspace(0,24,12+1) # less data, every 30 minutes
#%%
#Generate data
if D == '1':
    output_ids = [6,7]
    param_ids = [33, 12, 36, 32, 43, 39, 16, 7, 31, 37, 8, 15, 9, 28, 14]
elif D == '2':
    output_ids = [3,6,7]
    param_ids = [12, 33, 16, 7, 36, 43, 15, 39, 32, 8, 31, 14, 9, 1, 34]
elif D == '3':
    output_ids = [3,4,6,7]
    param_ids =[12, 33, 7, 36, 16, 43, 15, 39, 32, 8, 9, 31, 21, 14, 10, 34]


##use generated sets instead
with open('/work/aayres/syntheticData/TRUE_SOL' + D + dpoints + '.pkl', 'rb') as f:
    results = pickle.load(f)
y_data = results['y_data'].copy()
t_data = results['t_data'].copy()

# plt.figure()
# plt.subplot(4,1,1)
# plt.scatter(t_data,y_data[:,0],marker='x',color='k',label='obs')
# plt.plot(tfinal,true_sol[:,output_ids[0]],label='true')
# plt.title('True model + observations')
# plt.legend()
# #
# plt.subplot(4,1,2)
# plt.scatter(t_data,y_data[:,1],marker='x',color='k')
# plt.plot(tfinal,true_sol[:,output_ids[1]])
#
# # plt.subplot(4,1,3)
# # plt.scatter(t_data,y_data[:,2],marker='x',color='k')
# # plt.plot(tfinal,true_sol[:,output_ids[2]])
#
# # plt.subplot(4,1,4)
# # plt.scatter(t_data,y_data[:,3],marker='x',color='k')
# # plt.plot(tfinal,true_sol[:,output_ids[3]])
#
# plt.show()

#%%
# MJC - make these bounds based on the parameter values
UB = np.log(np.exp(param_log[param_ids])*10)
LB =  np.log(np.exp(param_log[param_ids])*0.001)
bounds = np.array([LB,UB]).T

#%%
def ols_call(x):
    r = OLS_res(x,y_data,t_data,param_ids,output_ids,param_log,IC_2)
    return(r)


par_fix = np.round(np.linspace(0,int(len(param_log)-1),int(len(param_log))))
par_fix = np.delete(par_fix,param_ids)
print(par_fix)
print(param_ids)
#%%
J_save, par_save, param_global_all = calc_PL_fixpar(ols_call,
    params=param_log[param_ids],
    bounds=bounds,
    par_fix=[],
    N=201,
    profile_width=0.5,
    max_nfev=300,
    finite_diff_step=1e-4,
    verbose=True
    #n_jobs=12
            )
#%%
# for i in range(len(param_ids)):
#
#     plt.figure()
#     plt.title(paramtitles[param_ids[i]])
#     plt.plot(par_save[i,i,:],J_save[i,:])
#     plt.show()

# plt.figure()
# plt.plot(par_save[1,1,:],J_save[1,:])
# plt.show()
#
# plt.figure()
# plt.plot(par_save[2,2,:],J_save[2,:])
# plt.show()
#%%
all_results = {
    'param_ids' : param_ids,
    't_data' : t_data,
    'y_data' : y_data,
    'true_sol' : true_sol,
    'J_save': np.asarray(J_save),
    'par_save': np.asarray(par_save),
    'param_global_all' : np.asarray(param_global_all)
}

with open('PL_' + D + dpoints + 'nonoise.pkl', 'wb') as f:
    pickle.dump(all_results, f)
