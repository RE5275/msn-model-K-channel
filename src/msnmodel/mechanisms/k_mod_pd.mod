COMMENT
This NMODL file defines a potassium current (k_mod_pd) adapted for Parkinson's disease modeling.

Enhancements:
- Simulates dopamine depletion via DA_level
- A2A receptor activation increases cAMP production based on DA level
- cAMP dynamically modulates potassium conductance (g_mod)
- Optional cAMP kinetics (production/degradation)
- MODIFIED: A2A activation now *INCREASES* potassium conductance (facilitatory effect on current)
- NEW: Includes a 'basal_cAMP_prod' to maintain a stable baseline cAMP level in control conditions.
- TEMPORARY DEBUGGING: 'n' gate is set to 1 to isolate g_mod modulation.
ENDCOMMENT

NEURON {
    SUFFIX k_mod_pd
    USEION k READ ek WRITE ik VALENCE 1
    RANGE gkbar, ik, g_mod, n_inf, tau_n
    RANGE A2A_on, DA_level, base_cAMP, pd_cAMP_gain
    GLOBAL nq, prod_rate, deg_rate, basal_cAMP_prod
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
}

PARAMETER {
    gkbar = 0.036 (S/cm2)       : Max potassium conductance (this is now the BASELINE for modulation)
    A2A_on = 0                  : A2A receptor activation toggle (0 = off, 1 = on)
    DA_level = 1                : Dopamine level (1 = normal, 0 = PD)
    pd_cAMP_gain = 0.5          : Sensitivity of g_mod to cAMP in PD (now represents FACILITATION strength)
    base_cAMP = 1               : Baseline cAMP level
    prod_rate = 0.5 (/ms)       : A2A-driven cAMP production rate (when A2A_on and DA loss)
    deg_rate = 0.1 (/ms)        : cAMP degradation rate
    basal_cAMP_prod = 0.1 (/ms) : Basal cAMP production rate (always on)
    nq = 4                      : Gating exponent
}

ASSIGNED {
    v (mV)
    ek (mV)
    ik (mA/cm2)
    g_mod (S/cm2)
    n_inf
    tau_n (ms)
}

STATE {
    n                               : Activation gate
    cAMP                            : Intracellular modulator (dynamic)
}

INITIAL {
    rates(v)
    n = 1.0  : TEMPORARY CHANGE: Set n to always be open for debugging
    cAMP = base_cAMP
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    rates(v)

    : Clamp cAMP to safe physiological range
    if (cAMP < 0) {
        cAMP = 0
    } else if (cAMP > 3) {
        cAMP = 3
    }

    : Modulate conductance based on cAMP (with A2A receptor activity)
    : A2A activation (via cAMP increase) should *INCREASE* potassium conductance
    if (A2A_on == 1) {
        : g_mod increases as cAMP increases above base_cAMP
        g_mod = gkbar * (1 + pd_cAMP_gain * (cAMP - base_cAMP))
        : Optional: Add a check for maximum conductance if needed,
        : e.g., if you don't want it to increase indefinitely beyond a certain point.
        : Example: if (g_mod > gkbar * 2) { g_mod = gkbar * 2 }
    } else {
        : In control, A2A_on is 0, so g_mod is just gkbar (no A2A modulation)
        g_mod = gkbar
    }

    ik = g_mod * n^nq * (v - ek)
}

DERIVATIVE states {
    : cAMP kinetics: A2A-driven production (when A2A_on and DA loss) + basal production - degradation
    cAMP' = (A2A_on * (1 - DA_level) * prod_rate) + basal_cAMP_prod - (deg_rate * cAMP)
    n' = (n_inf - n) / tau_n  : Keep n' equation, but n is forced to 1 in INITIAL
}

PROCEDURE rates(v (mV)) {
    LOCAL alpha_n, beta_n

    alpha_n = 0.01 * vtrap(10 - v, 10)
    beta_n = 0.125 * exp(-v / 80)

    tau_n = 1 / (alpha_n + beta_n)
    n_inf = alpha_n / (alpha_n + beta_n)
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
    if (fabs(x / y) < 1e-6) {
        vtrap = y * (1 - x / y / 2)
    } else {
        vtrap = x / (exp(x / y) - 1)
    }
}