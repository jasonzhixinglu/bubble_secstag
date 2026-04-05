"""
Endowment economy for a Python replication of:

    Eggertsson & Mehrotra (2014) "A Model of Secular Stagnation"

extended with a rational bubble asset following:

    Lu (2015) "Bubbles in a Secular Stagnation Economy"

The economy is a 3-period OLG model (young, middle-aged, old).
Only the middle-aged receive endowment income Y_m; the old receive Y_o
(zero in the EM baseline).  The young are borrowing-constrained at D.

All functions accept explicit arguments and are usable standalone.
Calibration constants from src.calibration are noted in the signatures
for reference but are not imported as defaults here.
"""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Loan market
# ---------------------------------------------------------------------------

def loan_demand(D: float, g: float, r: float) -> float:
    """Return the young generation's loan demand.

    The young borrow up to their binding debt limit D, discounted by the
    population growth factor (1+g) and the gross real interest rate (1+r).

    Parameters
    ----------
    D : float
        Collateral / debt limit on the young (e.g. calibration.D).
    g : float
        Population growth rate per OLG period (e.g. calibration.G).
    r : float
        Net real interest rate (not gross).

    Returns
    -------
    float
        Aggregate loan demand L^d (per middle-aged household).

    Notes
    -----
    Eggertsson & Mehrotra (2014), equation (9):

        L^d = (1+g)/(1+r) * D
    """
    return (1.0 + g) / (1.0 + r) * D


def loan_supply(
    Y_m: float,
    D_prev: float,
    Y_o: float,
    r: float,
    beta: float,
    A_bubble: float = 0.0,
) -> float:
    """Return the middle-aged generation's loan supply.

    The middle-aged save out of their endowment Y_m after repaying last
    period's debt D_prev, discounting old-age income Y_o, and net of any
    funds absorbed by the bubble asset.

    Parameters
    ----------
    Y_m : float
        Middle-aged endowment income (e.g. 1.0, normalised).
    D_prev : float
        Debt repaid this period (equals D from the previous period, per
        the binding constraint; e.g. calibration.D or calibration.D_L).
    Y_o : float
        Old-generation endowment income (e.g. calibration.Y_O = 0.0).
    r : float
        Net real interest rate (not gross).
    beta : float
        Household discount factor (e.g. calibration.BETA).
    A_bubble : float, optional
        Bubble asset size per capita (default 0.0, i.e. fundamental
        equilibrium).  Lu (2015) extension.

    Returns
    -------
    float
        Aggregate loan supply L^s (per middle-aged household).

    Notes
    -----
    Eggertsson & Mehrotra (2014), equation (10):

        L^s = beta/(1+beta) * (Y_m - D_prev)
            - 1/(1+beta) * Y_o/(1+r)

    Extended with the bubble term from Lu (2015), equation (15):

        L^s = beta/(1+beta) * (Y_m - D_prev)
            - 1/(1+beta) * Y_o/(1+r)
            - A_bubble
    """
    saving_motive = beta / (1.0 + beta) * (Y_m - D_prev)
    old_income_pv = 1.0 / (1.0 + beta) * Y_o / (1.0 + r)
    return saving_motive - old_income_pv - A_bubble


# ---------------------------------------------------------------------------
# Equilibrium natural rate
# ---------------------------------------------------------------------------

def natural_rate(
    D: float,
    D_prev: float,
    g: float,
    Y_m: float,
    Y_o: float,
    beta: float,
    A_bubble: float = 0.0,
) -> float:
    """Solve for the natural (gross) real interest rate analytically.

    Sets loan_demand == loan_supply and solves for (1+r).  The solution
    is closed-form because both sides are linear in 1/(1+r).

    Parameters
    ----------
    D : float
        Current-period collateral limit (borrowing by today's young).
    D_prev : float
        Previous-period collateral limit (debt repaid by today's middle-aged).
    g : float
        Population growth rate per OLG period.
    Y_m : float
        Middle-aged endowment income.
    Y_o : float
        Old-generation endowment income.
    beta : float
        Household discount factor.
    A_bubble : float, optional
        Bubble asset size per capita (default 0.0).

    Returns
    -------
    float
        Gross real interest rate R = (1+r) in equilibrium.

    Notes
    -----
    Eggertsson & Mehrotra (2014), equation (11):

        (1+g)/(1+r) * D  =  beta/(1+beta) * (Y_m - D_prev)
                          -  1/(1+beta) * Y_o/(1+r)

    Rearranging to isolate 1/(1+r):

        [(1+g)*D + Y_o/(1+beta)] / (1+r)
            = beta/(1+beta) * (Y_m - D_prev) - A_bubble

        1/(1+r) = [beta/(1+beta) * (Y_m - D_prev) - A_bubble]
                  / [(1+g)*D + Y_o/(1+beta)]

    Hence:

        R = (1+r) = [(1+g)*D + Y_o/(1+beta)]
                    / [beta/(1+beta) * (Y_m - D_prev) - A_bubble]
    """
    numerator = (1.0 + g) * D + Y_o / (1.0 + beta)
    denominator = beta / (1.0 + beta) * (Y_m - D_prev) - A_bubble
    return numerator / denominator


# ---------------------------------------------------------------------------
# Bubble sustainability
# ---------------------------------------------------------------------------

def bubble_upper_bound(Y_m: float, D: float, beta: float) -> float:
    """Return the maximum sustainable bubble size per capita.

    A bubble can exist in equilibrium only if it does not violate the
    middle-aged household's budget constraint.  The upper bound is
    derived from the condition that loan supply remains non-negative
    at r = g (the golden-rule interest rate).

    Parameters
    ----------
    Y_m : float
        Middle-aged endowment income.
    D : float
        Collateral limit (same period, i.e. D_prev = D in steady state).
    beta : float
        Household discount factor.

    Returns
    -------
    float
        Maximum bubble size A_max >= 0.

    Notes
    -----
    Lu (2015), Proposition 1:

        A_max = (beta * Y_m - (2 + beta) * D) / (1 + beta)
    """
    return (beta * Y_m - (2.0 + beta) * D) / (1.0 + beta)


def is_bubble_sustainable(
    A_bubble: float,
    Y_m: float,
    D: float,
    beta: float,
) -> bool:
    """Return True if the given bubble size is sustainable in equilibrium.

    A bubble is sustainable when it is non-negative and does not exceed
    the upper bound derived in Lu (2015) Proposition 1.

    Parameters
    ----------
    A_bubble : float
        Candidate bubble size per capita.
    Y_m : float
        Middle-aged endowment income.
    D : float
        Collateral limit (steady-state value, D_prev = D).
    beta : float
        Household discount factor.

    Returns
    -------
    bool
        True if 0 <= A_bubble <= bubble_upper_bound(Y_m, D, beta).
    """
    return 0.0 <= A_bubble <= bubble_upper_bound(Y_m, D, beta)


def fundamental_rate(
    D: float,
    g: float,
    Y_m: float,
    Y_o: float,
    beta: float,
) -> float:
    """Return the natural rate in the fundamental economy (no bubble).

    Convenience wrapper around natural_rate with A_bubble=0 and
    D_prev=D (steady state).  This is r^f in Lu (2015) Proposition 1,
    the threshold below which bubbles can exist in equilibrium.

    Parameters
    ----------
    D : float
        Steady-state collateral limit (D_prev = D = D in steady state).
    g : float
        Population growth rate per OLG period.
    Y_m : float
        Middle-aged endowment income.
    Y_o : float
        Old-generation endowment income.
    beta : float
        Household discount factor.

    Returns
    -------
    float
        Gross real interest rate R^f = (1 + r^f) in the fundamental
        steady state.  A bubble can exist when R^f <= 1 + g (i.e. when
        the economy is dynamically inefficient).
    """
    return natural_rate(D, D, g, Y_m, Y_o, beta, A_bubble=0.0)
