from scipy.optimize import curve_fit
import numpy as np
import numpy.typing as npt
from typing import Union, List, Tuple


def biexponential_model(t: float, A: float, tau1: float, tau2: float) -> float:
    """Fit data to a biexponential function (sum of two exponential decays). A and B should sum to 1.0
    C(t) = A*exp(-t/tau1) + B*exp(-t/tau2)

    Args:
        t (float): time (independent variable)
        A (float): Preexponential factor
        tau1 (float): time constant
        B (float): Preexponential factor
        tau2 (float): time constant

    Returns:
        float: Result of the expression, used for fitting
    """
    return A * np.exp(-t / tau1) + (1 - A) * np.exp(-t / tau2)


def fit_biexponential(
    tau_timeseries: Union[List[float], npt.NDArray],
    tcf_timeseries: Union[List[float], npt.NDArray],
    initial_guess: Union[List[float], npt.NDArray] = [0.5, 1, 2],
) -> Tuple[Tuple[float, float, float], npt.NDArray, npt.NDArray, float]:
    """Fit a biexponential function of the form (A*exp(-t/tau1) + B*exp(-t/tau2)) to the tau values and time autocorrelation function values (using the continuous bond definition).
    Using the relation from Gowers et al. (2015).

    Args:
        tau_timeseries (Union[List[float], npt.NDArray]): Lag times for the TCF
        tcf_timeseries (Union[List[float], npt.NDArray]): Time autocorrelation function values (this fit seems to work for the continuous bond definition).
        Once a bond breaks, it is considered broken even if reformed.
        initial_guess (Union[List[float], npt.NDArray], optional): Initial guess for parameters A, tau1 B, and tau2. Defaults to [0.5, 1, 2].

    Returns:
        Tuple[Tuple[float, float, float], npt.NDArray, npt.NDArray, float]: Returns a tuple of fitted parameters (A, tau1 and tau2), times, TCF values and the lifetime.
    """

    # Initial guess for the parameters [A, tau1, tau2]

    # Increase the maxfev parameter to allow more iterations
    params, params_covariance = curve_fit(
        biexponential_model,
        tau_timeseries,
        tcf_timeseries,
        p0=initial_guess,
        maxfev=10000,
    )

    fit_t = np.linspace(tau_timeseries[0], tau_timeseries[-1], 1000)
    fit_ac = biexponential_model(fit_t, *params)

    # Extract the parameters
    A, tau1, tau2 = params

    # Compute the lifetime
    lifetime = A * tau1 + (1 - A) * tau2

    return params, fit_t, fit_ac, lifetime
