## Variational QC / Gross-error penalty block (GNSS RSPD)

This is the **nonlinear “variational QC” (gross-error check)** part of the cost function for a GNSS RSPD observation in GSI. It computes:

- A **penalty term** `valqc` that replaces the simple quadratic penalty when gross-error probability is enabled.
- A corresponding **weight** `wgt` (roughly the derivative/linearized weight applied in the minimization), which downweights large innovations.

Conceptually it uses a **two-component mixture model**:

- “Good obs” likelihood: Gaussian `~ exp(-0.5 * (normalized innovation)^2)`
- “Gross obs” component: a broader/flat component controlled by `cvar_pg` (probability of gross error) and `cvar_b` (width/shape parameter)

This is standard in GSI: small departures behave like normal least-squares; big departures get **robustly downweighted**.

---

## Line-by-line meaning

### 1) Form the Gaussian likelihood for the (scaled) innovation

```fortran
val      = error*ddiff
val2     = val*val
exp_arg  = -half*val2
...
arg  = exp(exp_arg)
```

- `ddiff` is (typically) the innovation (obs−guess) in some normalized form.
- `error` is the inverse/scale factor used in GSI (often effectively `1/σ` or similar depending on convention).
- So `val = error*ddiff` is the **normalized departure**.
- `exp_arg = -0.5 * val^2`
- `arg = exp(exp_arg) = exp(-0.5 * val^2)` is proportional to the **Gaussian likelihood** of the observation given the background.

If you had no gross-error model, the cost term would just be `exp_arg` (i.e., quadratic penalty).

---

### 2) Bring in the gross-error mixture parameters

```fortran
rat_err2 = ratio_errors**2
if (cvar_pg(ikx) > tiny_r_kind .and. error > tiny_r_kind) then
   arg  = exp(exp_arg)
   wnotgross= one-cvar_pg(ikx)
   cg_gnssrspd=cvar_b(ikx)
   wgross = cg_term*cvar_pg(ikx)/(cg_gnssrspd*wnotgross)
```

- `cvar_pg(ikx)` is the **a priori probability of gross error** for this obs type (`ikx`).
- `wnotgross = 1 - pg` is probability of “not gross”.
- `cvar_b(ikx)` (here stored in `cg_gnssrspd`) is the **gross-error model scale/shape** parameter for this obs type.
- `wgross` becomes a **relative weight** of the gross-error component vs the “good” Gaussian component.

---

### 3) Compute the robustified cost contribution

```fortran
term = log((arg+wgross)/(one+wgross))
...
valqc = -two*rat_err2*term
```

Key identity:

- If `wgross = 0`, then:
  - `term = log(arg) = exp_arg`
  - so you recover the normal quadratic form.

When `wgross > 0`:

- `term = log( (arg + wgross) / (1 + wgross) )`
- This is the **log of the mixture likelihood ratio** (mixture vs normalization), which yields a **bounded / softened penalty** for large departures.

Then `valqc = -2 * (ratio_errors^2) * term` is the final penalty contribution (scaled by `ratio_errors^2`; in GSI this usually accounts for prescribed error inflation or ratioing).

---

### 4) Compute the effective weight (downweighting outliers)

```fortran
wgt  = one-wgross/(arg+wgross)
wgt  = wgt/wgtlim
```

Algebraically:

- `wgt = 1 - wgross/(arg+wgross) = arg/(arg+wgross)`

So:

- If the innovation is small: `arg ≈ 1`, then `wgt ≈ 1/(1+wgross)` (near 1 if `wgross` small).
- If the innovation is huge: `arg → 0`, then `wgt → 0` (obs gets essentially ignored).

So `wgt` is basically the **posterior probability of “not gross”** (up to constants), used as the **linearized weight** in the minimization.

`wgtlim` is a limiter/normalization factor used elsewhere in the routine (often to keep weights within a controlled range).

---

### 5) Fallback when variational QC is off or error is degenerate

```fortran
else
   term = exp_arg
   wgt  = wgtlim
   rwgt = wgt/wgtlim
endif
```

If `pg` is ~0 (QC off) or `error` is ~0, it reverts to:

- pure Gaussian penalty (`term = -0.5*val^2`)
- default weight (`wgt = wgtlim`)

---

## Quick interpretation summary

- **`term`**: log-likelihood-based robust penalty term (equals `-0.5*val^2` when gross-error model is disabled)
- **`wgross`**: how “available” the gross-error alternative is (based on `pg` and `b`)
- **`wgt = arg/(arg+wgross)`**: the downweighting factor; goes to **0 for outliers**
- **`valqc`**: the resulting QC penalty added into the cost
