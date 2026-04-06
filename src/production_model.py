"""
New Keynesian production economy for a Python replication of:

    Eggertsson & Mehrotra (2014) "A Model of Secular Stagnation"
    sections 4–6

extended with rational bubbles following:

    Lu (2015) "Bubbles in a Secular Stagnation Economy"
    section 3

The steady state is characterised by two curves in (Y, Pi) space:
aggregate supply (AS) and aggregate demand (AD), each with two regimes.
Equilibrium is found by numerically equating AS and AD.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


# ---------------------------------------------------------------------------
# Potential output
# ---------------------------------------------------------------------------

def potential_output(L_bar: float, alpha: float) -> float:
    """Return full-employment (potential) output.

    Parameters
    ----------
    L_bar : float
        Labour endowment of the middle-aged (e.g. calibration.L_BAR = 1.0).
    alpha : float
        Labour share in production (e.g. calibration.ALPHA = 0.7).

    Returns
    -------
    float
        Potential output Y^f = L_bar^alpha.

    Notes
    -----
    Eggertsson & Mehrotra (2014), equation (24):

        Y^f = L_bar^alpha
    """
    return L_bar ** alpha


# ---------------------------------------------------------------------------
# Aggregate supply
# ---------------------------------------------------------------------------

def aggregate_supply(
    Pi: float | np.ndarray,
    Y_f: float,
    gamma: float,
    alpha: float = 0.7,
) -> float | np.ndarray:
    """Return output on the aggregate supply (AS) curve.

    Two regimes depending on gross inflation Pi:

    - **Normal regime** (Pi >= 1): wages are flexible upward, labour market
      clears at full employment.
    - **ZLB / deflation regime** (Pi < 1): downward nominal wage rigidity
      prevents real wages from falling fast enough; employment and output
      fall below potential.

    Parameters
    ----------
    Pi : float or ndarray
        Gross inflation rate (P_t / P_{t-1}).
    Y_f : float
        Potential output, as returned by :func:`potential_output`.
    gamma : float
        Degree of downward nominal wage rigidity (calibration.GAMMA).
        0 = fully flexible, 1 = fully rigid.
    alpha : float, optional
        Labour share in production (default 0.7, calibration.ALPHA).
        Determines the curvature of the AS curve below Pi = 1.

    Returns
    -------
    float or ndarray
        Output Y on the AS curve.  Shape matches ``Pi``.

    Notes
    -----
    Eggertsson & Mehrotra (2014), equations (24)–(25):

    Normal regime (Pi >= 1):

        Y = Y^f                                                     (24)

    Deflation regime (Pi < 1): wage rigidity binds, so

        gamma / Pi = 1 - (1 - gamma) * (Y / Y^f)^{(1-alpha)/alpha}  (25)

    Solving (25) for Y:

        Y = Y^f * [(1 - gamma/Pi) / (1 - gamma)]^{alpha/(1-alpha)}
    """
    scalar_input = np.ndim(Pi) == 0
    Pi = np.atleast_1d(np.asarray(Pi, dtype=float))

    exponent = alpha / (1.0 - alpha)

    # Deflation branch: solve eq (25) for Y
    ratio = (1.0 - gamma / Pi) / (1.0 - gamma)
    # Clamp to avoid negative base when Pi is very small
    ratio = np.maximum(ratio, 0.0)
    Y_deflation = Y_f * ratio ** exponent

    Y = np.where(Pi >= 1.0, Y_f, Y_deflation)

    return float(Y[0]) if scalar_input else Y


# ---------------------------------------------------------------------------
# ZLB kink point
# ---------------------------------------------------------------------------

def pi_kink(i_star: float, phi_pi: float, Pi_star: float) -> float:
    """Return the gross inflation rate at which the ZLB becomes binding.

    The Taylor rule sets the nominal rate to

        1 + i = (1 + i_star) * (Pi / Pi_star)^{phi_pi}

    The ZLB binds when this rule would prescribe i < 0, i.e. when Pi falls
    below Pi_kink, where 1 + i = 1 exactly.

    Parameters
    ----------
    i_star : float
        Target nominal interest rate in the Taylor rule (calibration.I_STAR).
    phi_pi : float
        Taylor-rule coefficient on inflation (calibration.PHI_PI).
    Pi_star : float
        Gross inflation target (calibration.PI_STAR).

    Returns
    -------
    float
        Pi_kink = Pi_star * (1 / (1 + i_star))^{1 / phi_pi}.

    Notes
    -----
    Eggertsson & Mehrotra (2014), equation (28):

        Pi_kink = (1 / (1 + i_star))^{1 / phi_pi} * Pi_star
    """
    return Pi_star * (1.0 / (1.0 + i_star)) ** (1.0 / phi_pi)


# ---------------------------------------------------------------------------
# Aggregate demand
# ---------------------------------------------------------------------------

def aggregate_demand(
    Pi: float | np.ndarray,
    D: float,
    g: float,
    beta: float,
    phi_pi: float,
    Pi_star: float,
    i_star: float,
    A_bubble: float = 0.0,
) -> float | np.ndarray:
    """Return output on the aggregate demand (AD) curve.

    Two regimes depending on whether the ZLB is binding:

    - **Normal regime** (Pi >= Pi_kink): Taylor rule is active, the nominal
      rate i > 0 responds to inflation.
    - **ZLB regime** (Pi < Pi_kink): nominal rate is stuck at zero; the real
      rate rises with deflation, further depressing demand.

    Parameters
    ----------
    Pi : float or ndarray
        Gross inflation rate.
    D : float
        Collateral / borrowing limit (calibration.D or calibration.D_L).
    g : float
        Population growth rate per OLG period (calibration.G).
    beta : float
        Household discount factor (calibration.BETA).
    phi_pi : float
        Taylor-rule coefficient on inflation (calibration.PHI_PI).
    Pi_star : float
        Gross inflation target (calibration.PI_STAR).
    i_star : float
        Target nominal interest rate (calibration.I_STAR).
    A_bubble : float, optional
        Bubble size per capita (default 0.0, fundamental economy).

    Returns
    -------
    float or ndarray
        Output Y on the AD curve.  Shape matches ``Pi``.

    Notes
    -----
    Eggertsson & Mehrotra (2014), equations (26)–(27), extended by
    Lu (2015), equations (29)–(30):

    Let  Gamma_star = Pi_star^{phi_pi} / (1 + i_star).

    Normal regime (Pi >= Pi_kink):

        Y = D + (1+beta)/beta * (1+g)*D * Gamma_star / Pi^{phi_pi-1}
              + (1+beta)/beta * A_bubble                               (26/29)

    ZLB regime (Pi < Pi_kink):

        Y = D + (1+beta)/beta * (1+g)*D * Pi
              + (1+beta)/beta * A_bubble                               (27/30)

    The bubble term is additive in both regimes because it enters the
    underlying aggregate demand equation before any policy rule is
    substituted.
    """
    scalar_input = np.ndim(Pi) == 0
    Pi = np.atleast_1d(np.asarray(Pi, dtype=float))

    Pi_k = pi_kink(i_star, phi_pi, Pi_star)
    Gamma_star = Pi_star ** phi_pi / (1.0 + i_star)
    coeff = (1.0 + beta) / beta * (1.0 + g) * D
    bubble_term = (1.0 + beta) / beta * A_bubble

    Y_normal = D + coeff * Gamma_star / Pi ** (phi_pi - 1.0) + bubble_term
    Y_zlb    = D + coeff * Pi + bubble_term

    Y = np.where(Pi >= Pi_k, Y_normal, Y_zlb)

    return float(Y[0]) if scalar_input else Y


# ---------------------------------------------------------------------------
# Equilibrium solver
# ---------------------------------------------------------------------------

def find_equilibrium(
    D: float,
    g: float,
    beta: float,
    phi_pi: float,
    Pi_star: float,
    i_star: float,
    gamma: float,
    alpha: float,
    L_bar: float,
    A_bubble: float = 0.0,
) -> dict:
    """Find the steady-state equilibrium (Y, Pi) numerically.

    Solves aggregate_supply(Pi) == aggregate_demand(Pi) for Pi using
    Brent's method, then evaluates Y and classifies the regime.

    Parameters
    ----------
    D : float
        Collateral / borrowing limit.
    g : float
        Population growth rate per OLG period.
    beta : float
        Household discount factor.
    phi_pi : float
        Taylor-rule coefficient on inflation.
    Pi_star : float
        Gross inflation target.
    i_star : float
        Target nominal interest rate.
    gamma : float
        Degree of downward nominal wage rigidity.
    alpha : float
        Labour share in production.
    L_bar : float
        Labour endowment of the middle-aged.
    A_bubble : float, optional
        Bubble size per capita (default 0.0).

    Returns
    -------
    dict
        Keys:

        - ``Y``         : equilibrium output
        - ``Pi``        : equilibrium gross inflation
        - ``regime``    : ``'normal'`` or ``'zlb'``
        - ``r_natural`` : natural real interest rate (net), R^n - 1
        - ``i_nominal`` : equilibrium nominal interest rate (net)

    Notes
    -----
    The root is bracketed on (gamma, 1.5) because:

    - Pi cannot fall below gamma (AS is not defined for Pi < gamma
      when (1 - gamma/Pi) < 0).
    - Pi = 1.5 is safely above any reasonable equilibrium inflation.

    Equilibrium condition: AS(Pi) - AD(Pi) = 0.
    """
    Y_f = potential_output(L_bar, alpha)

    def residual(Pi: float) -> float:
        return (aggregate_supply(Pi, Y_f, gamma, alpha)
                - aggregate_demand(Pi, D, g, beta, phi_pi, Pi_star,
                                   i_star, A_bubble))

    try:
        Pi_eq = brentq(residual, gamma + 1e-8, 2.0, xtol=1e-10, rtol=1e-10)
    except ValueError:
        raise ValueError(
            "No equilibrium found in search range — try adjusting parameters."
        )
    Y_eq  = aggregate_supply(Pi_eq, Y_f, gamma, alpha)

    Pi_k  = pi_kink(i_star, phi_pi, Pi_star)
    regime = 'normal' if Pi_eq >= Pi_k else 'zlb'

    # Natural rate (gross) from the endowment-economy formula at Y_m = Y_f
    R_n = natural_rate_production(D, g, Y_f, beta, A_bubble)
    r_natural = R_n - 1.0

    # Nominal rate: Taylor rule if normal, zero if ZLB
    if regime == 'normal':
        i_nominal = (1.0 + i_star) * (Pi_eq / Pi_star) ** phi_pi - 1.0
    else:
        i_nominal = 0.0

    return {
        'Y':         Y_eq,
        'Pi':        Pi_eq,
        'regime':    regime,
        'r_natural': r_natural,
        'i_nominal': i_nominal,
    }


# ---------------------------------------------------------------------------
# Natural rate in the production economy
# ---------------------------------------------------------------------------

def natural_rate_production(
    D: float,
    g: float,
    Y_f: float,
    beta: float,
    A_bubble: float = 0.0,
) -> float:
    """Return the natural (gross) real interest rate in the production economy.

    Evaluated at full employment (Y_m = Y_f) and steady state (D_prev = D),
    using the loan-market clearing condition from the endowment economy.

    Parameters
    ----------
    D : float
        Collateral / borrowing limit (current and previous period equal in
        steady state).
    g : float
        Population growth rate per OLG period.
    Y_f : float
        Potential output (middle-aged income at full employment).
    beta : float
        Household discount factor.
    A_bubble : float, optional
        Bubble size per capita (default 0.0).

    Returns
    -------
    float
        Gross natural rate R^n = (1 + r^n).

    Notes
    -----
    Lu (2015), Proposition 2.  At full employment the endowment-economy
    loan market clears at:

        R^n = (1 + g) * D / [beta/(1+beta) * (Y_f - D) - A_bubble]

    A bubble (A_bubble > 0) reduces the denominator and raises R^n,
    potentially lifting the economy out of secular stagnation.
    """
    numerator   = (1.0 + g) * D
    denominator = beta / (1.0 + beta) * (Y_f - D) - A_bubble
    return numerator / denominator
