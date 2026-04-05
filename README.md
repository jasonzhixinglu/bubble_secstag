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
secular-stagnation-bubbles/
├── src/
│   ├── __init__.py
│   ├── endowment_model.py      # Layer 1: loan market equilibrium, Proposition 1
│   ├── production_model.py     # Layer 2: AS/AD, Proposition 2, ZLB dynamics
│   ├── stability.py            # Layer 3: bubble stability, policy analysis
│   ├── calibration.py          # Parameter values
│   └── plots.py                # All figures
├── notebooks/
│   ├── 01_endowment_economy.ipynb
│   ├── 02_production_economy.ipynb
│   └── 03_stability_and_policy.ipynb
├── figures/                    # Output charts
├── requirements.txt
└── README.md
```

## Key equations

Loan market equilibrium (endowment economy, with bubble):

    1 + r_t = (1 + g) * D / [ beta/(1+beta) * (Y - D) - A ]

Aggregate demand (ZLB binding, i = 0):

    Y = D + (1+beta)/beta * (1+g)*D * Pi + (1+beta)/beta * A

Aggregate supply (wage norm binding, Pi < 1):

    gamma/Pi = 1 - (1-gamma) * (Y/Y^n)^((1-alpha)/alpha)

Natural rate of interest:

    r^n = (1+g)*D / [ beta/(1+beta) * (Y^n - D) - A ] - 1

## References

- Lu, J. (2015). Bubbles in a Secular Stagnation Economy. University of Cambridge. (First version August 2015; this replication uses the March 2017 draft.)
- Eggertsson, G. B. and Mehrotra, N. R. (2014). A Model of Secular Stagnation. NBER WP 20574.
- Tirole, J. (1985). Asset Bubbles and Overlapping Generations. Econometrica, 53(5), 1071-1100.
- Gali, J. (2014). Monetary Policy and Rational Asset Price Bubbles. AER, 104(3), 721-752.