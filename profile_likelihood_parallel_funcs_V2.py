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