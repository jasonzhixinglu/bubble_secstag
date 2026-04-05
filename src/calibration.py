"""
Baseline parameter calibration for a Python replication of:

    Eggertsson & Mehrotra (2014) "A Model of Secular Stagnation"

extended with rational asset-price bubbles following:

    Lu (2015) "Bubbles in a Secular Stagnation Economy"

Parameter values are taken from Table 1 of Eggertsson & Mehrotra (2014)
unless otherwise noted.
"""

# ---------------------------------------------------------------------------
# Household
# ---------------------------------------------------------------------------
BETA = 0.77         # Discount factor (Table 1)

# ---------------------------------------------------------------------------
# Endowment / collateral
# ---------------------------------------------------------------------------
D = 0.28            # Collateral constraint, high steady state D_H (Table 1)
D_L = 0.259         # Low collateral constraint post-deleveraging shock
                    #   (7.5% tightening of D_H, as in EM Figure 3)
G = 0.2             # Population growth rate per OLG period (Table 1)
Y_O = 0.0           # Old-generation endowment income (zero in baseline)

# ---------------------------------------------------------------------------
# Production
# ---------------------------------------------------------------------------
ALPHA = 0.7         # Labour share in production (Table 1)
L_BAR = 1.0         # Labour endowment of middle-aged, normalised (Table 1)
GAMMA = 0.3         # Degree of downward nominal wage rigidity (Table 1)
                    #   (0 = flexible, 1 = fully rigid)

# ---------------------------------------------------------------------------
# Taylor rule
# ---------------------------------------------------------------------------
PHI_PI = 2.0        # Taylor-rule coefficient on inflation (Table 1)
PI_STAR = 1.01      # Gross inflation target (Table 1)
I_STAR = 0.0        # Targeted nominal interest rate in Taylor rule; EM assume
                    #   i_star = r^f so the kink shifts with the natural rate.
                    #   Set to 0.0 as baseline (consistent with PI_STAR = 1.01
                    #   and a near-zero natural rate).

# ---------------------------------------------------------------------------
# Bubble (Lu 2015 extension)
# ---------------------------------------------------------------------------
A_BUBBLE = 0.0      # Bubble size per capita (0.0 => fundamental equilibrium)
