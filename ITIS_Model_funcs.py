import numpy as np
from scipy.integrate import odeint

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
    #IC = [0.0, 0., 0.1595967, 0., 14.68266827, 41.49891056, 39.97751646, 11.29827208] (UNSURE WHAT TIMEPOINT THIS WAS PULLED FROM)
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
    param_eval = param_fix
    param_eval[param_ids] = param_in
    ##Define time spaces
    tstart = 0
    t_end =24*8
    dt = .05
    tspace = np.arange(tstart,t_end,dt) # The time frame to solve the model to reach steady state (before we grab any output)
    tfinal = np.arange(0,24+dt,dt)

    ##Solve the model at nominal values
    ode_options = {'rtol': 1e-6} # sets the numerical accuracy of the ODE solver (i.e., how fine of a time step do we need)
    ODE_sol = odeint(ITIS, IC, tspace, args = (param_eval,), **ode_options)

    # MJC - this is where we can specify that the infection begins during this final time point
    IC_2    = ODE_sol[-1,:]
    IC_2[0] = 2.0
    # Resolve model only over 24 hours at subsampled time points. This will be used as our data
    ODE_sol = odeint(ITIS, IC_2, tfinal, args = (param_eval,), **ode_options)

    # In a similar way, the data may be only a subset of the observations. We will assume that
    # we have the same number of time points for now, but can modify this as well
    model_output = np.zeros((len(t_data),len(output_ids)))

    for i, time in enumerate(t_data):
        t_index = np.where(tspace==time) #Grabs the index correpsonding to where we have data
        for j, mod_out in enumerate(output_ids): #Loop over which outputs we want to compare to data
            model_output[i,j] = ODE_sol[t_index,mod_out][0][0] ##THIS IS BEHAVING WEIRDLY

    # Now, we have model outputs at the time points and recorded states. We can construct the residual and then the sum of squares
    # in the future, we can look at applying weights to the data that correspond to measurement noise
    # Here we just scale by the data
    # MJC - Since we pass in the data and then take the mean here, we don't need to do it at the beginning
    residual = (y_data - model_output)/np.mean(y_data,axis=0)

    return residual.flatten()
