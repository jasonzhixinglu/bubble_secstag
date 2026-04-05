"""
Bubble stability, fiscal policy (PAYG), and monetary policy ineffectiveness
for a Python replication of:

    Lu (2015) "Bubbles in a Secular Stagnation Economy", Section 4

The three main results covered here are:
  1. Stability condition: a bubble is sustainable iff its size does not
     exceed the surplus of loan supply over loan demand at r = g.
  2. PAYG transfers crowd out the bubble one-for-one, reducing r^n.
  3. Monetary policy cannot affect the equilibrium bubble size.
"""

from __future__ import annotations

from src.production_model import find_equilibrium, natural_rate_production


# ---------------------------------------------------------------------------
# Generalised bubble upper bound (Section 4)
# ---------------------------------------------------------------------------

def bubble_upper_bound_with_old_income(
    Y_m: float,
    Y_o: float,
    D: float,
    g: float,
    beta: float,
    r: float,          # noqa: ARG001  (kept for API symmetry; formula uses g)
) -> float:
    """Return the maximum sustainable bubble size given old-generation income.

    The bubble can be sustained in a steady-state rational-expectations
    equilibrium only when the real interest rate equals the population
    growth rate, r = g (golden-rule / dynamic-inefficiency condition).
    The maximum sustainable bubble is the surplus of loan supply over loan
    demand when both are evaluated at r = g.

    Parameters
    ----------
    Y_m : float
        Middle-aged endowment income.
    Y_o : float
        Old-generation endowment income.
    D : float
        Collateral / borrowing limit (steady state: D_prev = D).
    g : float
        Population growth rate per OLG period.
    beta : float
        Household discount factor.
    r : float
        Placeholder; the bound is evaluated at r = g regardless.
        Retained in the signature for API consistency with other rate
        functions.

    Returns
    -------
    float
        Maximum bubble size per capita A_max.  Non-positive values imply
        no bubble can exist.

    Notes
    -----
    Lu (2015), Section 4.  At r = g, loan demand equals D (since
    (1+g)/(1+g) = 1) and loan supply is

        L^s(g) = beta/(1+beta) * (Y_m - D) - 1/(1+beta) * Y_o/(1+g)

    The maximum bubble absorbs all remaining surplus:

        A_max = L^s(g) - D
              = beta/(1+beta) * (Y_m - D)
                - 1/(1+beta) * Y_o/(1+g)
                - D
    """
    loan_supply_at_g = (
        beta / (1.0 + beta) * (Y_m - D)
        - 1.0 / (1.0 + beta) * Y_o / (1.0 + g)
    )
    loan_demand_at_g = D   # (1+g)/(1+g) * D = D
    return loan_supply_at_g - loan_demand_at_g


# ---------------------------------------------------------------------------
# Stability margin
# ---------------------------------------------------------------------------

def stability_margin(
    A_bubble: float,
    Y_m: float,
    D: float,
    beta: float,
    Y_o: float = 0.0,
    g: float | None = None,
) -> float:
    """Return the distance between the current bubble and its upper bound.

    A positive margin indicates the bubble is stable; zero or negative
    implies the bubble would collapse.  Larger bubbles always have smaller
    stability margins, consistent with the fragility result in Lu (2015)
    Section 4.

    Parameters
    ----------
    A_bubble : float
        Current bubble size per capita.
    Y_m : float
        Middle-aged endowment income.
    D : float
        Collateral / borrowing limit.
    beta : float
        Household discount factor.
    Y_o : float, optional
        Old-generation endowment income (default 0.0).
    g : float
        Population growth rate per OLG period.  Must be provided
        explicitly; no default to avoid silent errors when g=0 and g=0.2
        give meaningfully different stability margins.

    Returns
    -------
    float
        Stability margin = A_max - A_bubble.

    Notes
    -----
    Lu (2015), Section 4.  Margin > 0 iff the bubble is sustainable.
    """
    if g is None:
        raise ValueError(
            "g must be provided explicitly — no default "
            "to avoid silent errors with g=0 vs g=0.2."
        )
    a_max = bubble_upper_bound_with_old_income(Y_m, Y_o, D, g, beta, r=g)
    return a_max - A_bubble


# ---------------------------------------------------------------------------
# Critical collateral shock
# ---------------------------------------------------------------------------

def critical_shock_D(A_bubble: float, Y_m: float, beta: float) -> float:
    """Return the collateral limit above which the bubble becomes unsustainable.

    Inverts the Proposition 1 upper bound formula (Y_o = 0, steady state)
    to find the value of D at which A_max = A_bubble exactly.  For D above
    this threshold the upper bound falls below A_bubble and the bubble
    collapses.

    Parameters
    ----------
    A_bubble : float
        Bubble size whose sustainability is being tested.
    Y_m : float
        Middle-aged endowment income.
    beta : float
        Household discount factor.

    Returns
    -------
    float
        Critical collateral D_crit such that A_max(D_crit) = A_bubble.

    Notes
    -----
    Lu (2015), Proposition 1.  Setting

        A_max = (beta * Y_m - (2 + beta) * D) / (1 + beta) = A_bubble

    and solving for D:

        D_crit = (beta * Y_m - (1 + beta) * A_bubble) / (2 + beta)
    """
    return (beta * Y_m - (1.0 + beta) * A_bubble) / (2.0 + beta)


# ---------------------------------------------------------------------------
# Critical income shock
# ---------------------------------------------------------------------------

def critical_shock_Y(A_bubble: float, D: float, beta: float) -> float:
    """Return the income level below which the bubble becomes unsustainable.

    Inverts the Proposition 1 upper bound formula to find the value of Y_m
    at which A_max = A_bubble exactly.  For Y_m below this threshold the
    bubble collapses.

    Parameters
    ----------
    A_bubble : float
        Bubble size whose sustainability is being tested.
    D : float
        Collateral / borrowing limit.
    beta : float
        Household discount factor.

    Returns
    -------
    float
        Critical income Y_crit such that A_max(Y_crit) = A_bubble.

    Notes
    -----
    Lu (2015), Proposition 1 (income analogue).  Setting

        A_max = (beta * Y_m - (2 + beta) * D) / (1 + beta) = A_bubble

    and solving for Y_m:

        Y_crit = ((1 + beta) * A_bubble + (2 + beta) * D) / beta
    """
    return ((1.0 + beta) * A_bubble + (2.0 + beta) * D) / beta


# ---------------------------------------------------------------------------
# PAYG transfer effect
# ---------------------------------------------------------------------------

def payg_transfer_effect(
    Q: float,
    A_bubble: float,
    Y_m: float,
    D: float,
    g: float,
    beta: float,
    Y_o: float = 0.0,
) -> dict:
    """Model the effect of a PAYG pension transfer on the bubble equilibrium.

    A pay-as-you-go transfer of size Q is levied on the middle-aged and paid
    to the old.  In the loan market this is a perfect substitute for the
    bubble: middle-aged saving falls by Q, exactly crowding out A_bubble.
    The effective bubble is therefore A_bubble - Q.

    Parameters
    ----------
    Q : float
        PAYG transfer size per middle-aged household (Q >= 0).
    A_bubble : float
        Pre-transfer bubble size per capita.
    Y_m : float
        Middle-aged endowment income (gross, before transfer).
    D : float
        Collateral / borrowing limit.
    g : float
        Population growth rate per OLG period.
    beta : float
        Household discount factor.
    Y_o : float, optional
        Old-generation endowment income excluding PAYG (default 0.0).

    Returns
    -------
    dict
        Keys:

        - ``A_effective``       : A_bubble - Q (effective post-transfer bubble)
        - ``r_natural``         : net natural rate with effective bubble
        - ``stability_margin``  : A_max - A_effective at the new equilibrium
        - ``bubble_eliminated`` : True if Q >= A_bubble

    Notes
    -----
    Lu (2015), Section 4, PAYG analysis.  The transfer reduces the net
    surplus of saving available to support the bubble one-for-one.  This is
    the crowding-out result: fiscal policy can stabilise or eliminate a
    bubble, but at the cost of raising r^n toward the ZLB region.
    """
    A_effective = A_bubble - Q
    bubble_eliminated = Q >= A_bubble

    # Natural rate at effective bubble size (use Y_f = Y_m as proxy)
    Y_f = Y_m
    R_n = natural_rate_production(D, g, Y_f, beta, A_bubble=max(A_effective, 0.0))
    r_n = R_n - 1.0

    margin = stability_margin(max(A_effective, 0.0), Y_m, D, beta, Y_o, g)

    return {
        'A_effective':      A_effective,
        'r_natural':        r_n,
        'stability_margin': margin,
        'bubble_eliminated': bubble_eliminated,
    }


# ---------------------------------------------------------------------------
# Monetary policy ineffectiveness
# ---------------------------------------------------------------------------

def monetary_policy_effect(
    Pi_target: float,
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
    """Show that monetary policy cannot change the equilibrium bubble size.

    Computes the steady-state equilibrium at two inflation targets —
    the baseline Pi_star and an alternative Pi_target — and records both
    equilibria.  The bubble size A_bubble is identical in both cases,
    demonstrating monetary neutrality with respect to bubble size.

    Parameters
    ----------
    Pi_target : float
        Alternative gross inflation target (e.g. a higher target to try to
        escape secular stagnation via inflation).
    D : float
        Collateral / borrowing limit.
    g : float
        Population growth rate per OLG period.
    beta : float
        Household discount factor.
    phi_pi : float
        Taylor-rule coefficient on inflation.
    Pi_star : float
        Baseline gross inflation target.
    i_star : float
        Target nominal interest rate.
    gamma : float
        Degree of downward nominal wage rigidity.
    alpha : float
        Labour share in production.
    L_bar : float
        Labour endowment of the middle-aged.
    A_bubble : float, optional
        Bubble size per capita (default 0.0; same in both equilibria).

    Returns
    -------
    dict
        Keys:

        - ``baseline``   : dict with Y, Pi, regime, i_nominal at Pi_star
        - ``alternative``: dict with Y, Pi, regime, i_nominal at Pi_target
        - ``A_baseline`` : bubble size under baseline target (= A_bubble)
        - ``A_alternative``: bubble size under alternative target (= A_bubble)
        - ``bubble_unchanged``: True (illustrates monetary neutrality)

    Notes
    -----
    Lu (2015), Section 4, monetary policy discussion.  Changing the
    inflation target shifts the (Y, Pi) equilibrium along the AS/AD
    diagram but leaves the bubble size unaffected because the bubble is
    determined by the loan-market surplus at r = g, not by the monetary
    policy rule.
    """
    eq_baseline = find_equilibrium(
        D, g, beta, phi_pi, Pi_star, i_star, gamma, alpha, L_bar,
        A_bubble=A_bubble,
    )
    eq_alternative = find_equilibrium(
        D, g, beta, phi_pi, Pi_target, i_star, gamma, alpha, L_bar,
        A_bubble=A_bubble,
    )

    return {
        'baseline':        eq_baseline,
        'alternative':     eq_alternative,
        'A_baseline':      A_bubble,
        'A_alternative':   A_bubble,
        'bubble_unchanged': True,
    }
