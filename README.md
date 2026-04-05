# Bubbles in a Secular Stagnation Economy — Python Replication

Replication of Lu (2015), "Bubbles in a Secular Stagnation Economy," building on
Eggertsson & Mehrotra (2014), "A Model of Secular Stagnation," NBER Working Paper 20574.

## Overview

This project replicates the theoretical and quantitative results of Lu (2015) in Python.
The paper extends the Eggertsson-Mehrotra (EM) three-period OLG framework to allow for
trade in a rational bubble asset (Tirole, 1985). The central finding is that bubbles increase
the natural rate of interest by absorbing excess savings, thereby alleviating secular
stagnation — but larger bubbles are less stable, and their collapse can trigger it.

## Model structure

The replication proceeds in three layers, matching the paper's structure:

**Layer 1 — Endowment economy (Section 2)**
A three-period OLG model (young, middle-aged, old) with a bond market and a bubbly
asset. The young borrow subject to a binding debt limit D; the middle-aged save via bonds
and by purchasing the bubble from the old. The equilibrium real interest rate is derived
analytically and shown to be increasing in bubble size A. Proposition 1 characterises the
set of sustainable bubble equilibria and the upper bound on A.

**Layer 2 — Production economy with secular stagnation (Section 3)**
New Keynesian ingredients are added: a production function, downward nominal wage
rigidity, a Taylor rule, and the ZLB. Steady-state aggregate supply and aggregate demand
curves are derived in (Y, Pi) space. Proposition 2 characterises how the natural rate of
interest maps to output, and how bubbles shift aggregate demand to restore full employment.

**Layer 3 — Bubble stability and policy (Section 4)**
Comparative statics on the upper bound of A with respect to shocks to Y, D, and g.
Analysis of a PAYG fiscal transfer as a bubble substitute. Monetary policy ineffectiveness
at the ZLB.

## Repository structure

```
bubble_secstag/
├── src/
│   ├── __init__.py
│   ├── calibration.py          # EM (2014) Table 1 parameter values
│   ├── endowment_model.py      # Layer 1: loan market, natural rate, bubble bound
│   ├── production_model.py     # Layer 2: AS/AD curves, ZLB, equilibrium solver
│   └── stability.py            # Layer 3: bubble stability, PAYG, monetary policy
├── notebooks/
│   ├── 01_endowment_economy.ipynb   # Loan market diagrams, bubble effects (EM Fig 1)
│   ├── 02_production_economy.ipynb  # AS/AD diagrams (EM Fig 2–3, Lu Fig 3)
│   └── 03_stability_and_policy.ipynb # Stability margin, critical shocks, PAYG
├── figures/                    # Output charts (PDF and PNG)
├── requirements.txt
└── README.md
```

## Installation and usage

```bash
git clone https://github.com/jasonzhixinglu/bubble_secstag.git
cd bubble_secstag
pip install -r requirements.txt
jupyter notebook
```

Then open any notebook under `notebooks/` and run all cells (Kernel → Restart & Run All).
The notebooks add the repo root to `sys.path` automatically, so no installation of the
`src` package is required.

A `.devcontainer/devcontainer.json` is provided for one-click setup in GitHub Codespaces
or VS Code Dev Containers.

## Key results (EM 2014 Table 1 calibration)

| Quantity | Value | Notes |
|----------|-------|-------|
| β (discount factor) | 0.77 | EM Table 1 |
| D_H (collateral limit) | 0.28 | Pre-shock |
| D_L (post-deleveraging) | 0.259 | 7.5% tightening |
| g (population growth) | 0.20 | Per OLG period |
| A_max at D_L | 0.0634 | Max sustainable bubble |
| r^f at D_H (fundamental) | +7.3% | Normal regime |
| r^f at D_L (fundamental, no bubble) | −3.6% | Secular stagnation |
| Equilibrium Y at D_H | 1.000 (Y^f) | Full employment |
| Equilibrium Pi at D_H | 1.094 | Above target |
| Equilibrium Y at D_L, A=0 | 0.918 | Below potential |
| Equilibrium Pi at D_L, A=0 | 0.923 | Deflation |
| Bubble restoration (A = A_max/2) | Y=1.000, Pi=1.010 | Full employment restored |

## Key equations

Loan market equilibrium (endowment economy, with bubble):

    1 + r = (1 + g) * D / [ beta/(1+beta) * (Y_m - D) - A ]

Bubble upper bound (Proposition 1):

    A_max = (beta * Y_m - (1 + 2*beta) * D) / (1 + beta)

Aggregate demand (ZLB binding, i = 0):

    Y = D + (1+beta)/beta * (1+g)*D * Pi + (1+beta)/beta * A

Aggregate supply (wage norm binding, Pi < 1):

    gamma/Pi = 1 - (1-gamma) * (Y/Y^f)^((1-alpha)/alpha)

Natural rate of interest:

    r^n = (1+g)*D / [ beta/(1+beta) * (Y^f - D) - A ] - 1

## Known deviations from the paper

**Bubble upper bound formula** (`src/endowment_model.py`):
The formula given in Lu (2015) Proposition 1 is sometimes cited as
`(beta*Y_m - (2+beta)*D) / (1+beta)`. This is incorrect. The correct
derivation — setting loan supply equal to loan demand at r = g — gives:

    A_max = beta/(1+beta) * (Y_m - D) - D
           = (beta * Y_m - (1 + 2*beta) * D) / (1 + beta)

The coefficient on D is `(1 + 2*beta)`, not `(2 + beta)`. The two
coincide only when beta = 1. This replication uses the corrected formula.

## References

- Lu, J. (2015). Bubbles in a Secular Stagnation Economy. University of Cambridge.
- Eggertsson, G. B. and Mehrotra, N. R. (2014). A Model of Secular Stagnation. NBER WP 20574.
- Tirole, J. (1985). Asset Bubbles and Overlapping Generations. Econometrica, 53(5), 1071–1100.
- Galí, J. (2014). Monetary Policy and Rational Asset Price Bubbles. AER, 104(3), 721–752.
