# PESolver Algorithm Reference

## Purpose and Scope

This document explains how `PESolver` computes the diffraction efficiency of a grating for a single (wavelength, incidence angle) pair. Or to be more acurate, it shows the physical validation how and why you can compute something likes this and also how it is (numerical) implemented. It does not covers the geometry aspect: it consumes this input from `PEG`. It also does not include the MPI orchestration in `mainMPI`, or the material refractive-index database lookup — those are documented separately.

The method implemented here is a **differential (coupled-wave) method**: the structure is treated as a stack of thin horizontal layers; within each layer, a coupled ordinary differential equation (ODE) system is integrated numerically along the vertical coordinate `y`; the per-layer results are combined into a numerically stable overall response via **S-matrix recursion**. This is closely related to the family of methods used in Rigorous Coupled-Wave Analysis (RCWA) and multilayer optics/ellipsometry.

---

## Table of Contents

0. [Physical Setup](#00-physical-setup)
    1. [Maxwell to Scalar Wavefunction](#01-maxwell-to-scalar-wavefunction)
    2. [Fourier & Floquet Expansion](#02-fourier-and-floquet-expansion)
    3. [Coupled Mode Equation Derivation](#03-coupled-mode-equation-derivation)
    4. [Physical Interpretation](#04-physical-interpretation)

---

## 0.0 Physical Setup

Overview of the Theoretical Pipeline:

$$\text{Maxwell's Equations} \longrightarrow \text{Helmholtz Equation} \longrightarrow \text{Scalar Wave Equation (TE)}$$
$$\downarrow$$
$$\text{Periodic Fourier/Floquet Expansion} \longrightarrow \text{Coupled ODEs in } y \longrightarrow \text{Local Layer Response}$$
$$\downarrow$$
$$\text{Transfer Matrix} \longrightarrow \text{Scattering Matrix} \longrightarrow \text{Reflection Amplitudes} \longrightarrow \text{Diffraction Efficiencies}$$


Parameter Categorization:

- Geometry, Grating, and Material Structure: Defined by the spatial permittivity distribution $\varepsilon_r(x,y)$ and non-magnetic permeability ($\mu = \mu_0$). This encompasses the grating period $d$, layer/coating thicknesses along $y$ and the physical profile/surface topology of the grating.
- Incident Field Parameters: Free-space wavelength $\lambda$ and angle of incidence $\theta$


### 0.1 Maxwell to Scalar Wavefunction

Faraday's law of induction and Ampère's circuital law:
$$\vec{\nabla} \times \vec{E} = -\frac{\partial \vec{B}}{\partial t}, \quad \vec{\nabla} \times \vec{H} = -\vec{i} + \frac{\partial \vec{D}}{\partial t}$$

with $\vec{D} = \varepsilon_0 \varepsilon_r \vec{E}$ and $\vec{B} = \mu \vec{H}$ (assuming non-magnetic media where $\mu = \mu_0$ and $\vec i = \vec 0$), which can be formulated as:
$$\vec{\nabla} \times \vec{E} = -\mu_0 \frac{\partial \vec{H}}{\partial t}, \quad \vec{\nabla} \times \vec{H} = \varepsilon_0 \varepsilon_r \frac{\partial \vec{E}}{\partial t}$$

Separation of time and space assuming time-harmonic fields $\vec{E} e^{-i\omega t}$ and $\vec{H} e^{-i\omega t}$:
$$\vec{\nabla} \times \vec{E} = i\omega\mu_0\vec{H}, \quad \vec{\nabla} \times \vec{H} = -i \omega \varepsilon_0\varepsilon_r\vec{E}$$

To eliminate $\vec{H}$, take the curl of Faraday's law:
$$\vec{\nabla} \times (\vec{\nabla} \times \vec{E}) = i\omega\mu_0 (\vec{\nabla} \times \vec{H})$$

Substituting Ampère's law into the right-hand side yields:
$$\vec{\nabla} \times (\vec{\nabla} \times \vec{E}) = i\omega\mu_0 \left( -i\omega\varepsilon_0\varepsilon_r\vec{E} \right) = \omega^2 \mu_0 \varepsilon_0 \varepsilon_r \vec{E}$$

Using the identity $\vec{\nabla} \times (\vec{\nabla} \times \vec{E}) = \vec{\nabla}(\vec{\nabla} \vec{E}) - \Delta \vec{E}$ and setting $k_0^2 = \omega^2 \mu_0 \varepsilon_0$, we obtain the 3D vector wave equation:
$$\vec{\nabla}(\vec{\nabla} \vec{E}) - \Delta \vec{E} = k_0^2 \varepsilon_r \vec{E}$$

For Transverse Electric (TE) polarization in a 1D grating invariant along the $z$-axis ($\frac{\partial}{\partial z} = 0$), the electric field simplifies to $u(x,y) \equiv E_z(x,y)$ with $\vec{\nabla} \vec{E} = 0$. This reduces the system to the scalar wave equation:
$$\left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + k_0^2 \varepsilon_r(x,y) \right) u(x,y) = 0, \quad k_0=\frac{2\pi}{\lambda}$$


### 0.2 Fourier and Floquet Expansion

Due to the spatial periodicity of the material, $\varepsilon_r(x+d,y) = \varepsilon_r(x,y)$, the problem is invariant under translation $x \rightarrow x+d$. According to Bloch-Floquet theory ($\dot {x} =A(t)x$), all scattered fields must acquire the same phase shift upon translation by one period $d$. Given an incident tangential wavevector component $\alpha_{\text{inc}} = k_{\text{top}}\sin\theta$, the allowed tangential wavenumbers are discretized as:
$$\alpha_n = k_{\text{top}}\sin\theta + n\frac{2\pi}{d}, \quad n \in \mathbb{Z}$$

For the TE-Polarisation, the field $u(x,y)$ is expanded by separating the known $x$-dependence from the unknown transverse profiles $u_n(y)$:
$$u(x,y) = \sum_{n=-N}^{N} u_n(y) e^{i\alpha_n x}$$

Here, $u_n(y)$ represents the complex amplitude of the $n$-th Floquet harmonic at depth $y$.

### 0.3 Coupled Mode Equation Derivation

Substituting this series into the 2D scalar Helmholtz equation:

$$\left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + k_0^2\varepsilon_r(x,y) \right) \left( \sum_{n=-N}^{N} u_n(y) e^{i\alpha_n x} \right) = 0$$

- **Second derivative with respect to $x$:**
  $$\frac{\partial^2}{\partial x^2} u(x,y) = -\sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x}$$

- **Second derivative with respect to $y$:**
  $$\frac{\partial^2}{\partial y^2} u(x,y) = \sum_m u_m''(y) e^{i\alpha_m x}$$

Combining these gives:

$$\sum_m u_m''(y) e^{i\alpha_m x} - \sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x} + k^2(x,y) u(x,y) = 0$$

where $k^2(x,y) \equiv k_0^2 \varepsilon_r(x,y)$ is the periodic wavenumber profile. Due to spatial periodicity, $k^2(x,y)$ can also be expanded into a Fourier series:

$$k^2(x,y) = \sum_p k_p^2(y) e^{i p K x}, \quad K = \frac{2\pi}{d}$$

Inserting this expansion into the material product term gives:

$$\sum_m u_m''(y) e^{i\alpha_m x} - \sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x} + \sum_p k_p^2(y) e^{i p K x} \sum_{m} u_m(y) e^{i\alpha_m x} = 0$$

$$\sum_m u_m''(y) e^{i\alpha_m x} - \sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x} + \sum_p \sum_m k_p^2(y) u_m(y) e^{i (p K + \alpha_m) x} = 0$$

$$\sum_m \left( u_m''(y) - \alpha_m^2 u_m(y) \right) e^{i\alpha_m x} + \sum_p \sum_m k_p^2(y) u_m(y) e^{i (p K + \alpha_m) x} = 0$$

Using the relation for the allowed tangential wavenumbers $\alpha_m = k\sin\theta + m K$ ([Fourier](#02-fourier)), we observe that $p K + \alpha_m = \alpha_{m+p}$. Performing an index substitution by setting $n = m + p$ (or equivalently $p = n - m$), the second summation becomes:

$$\sum_p \sum_m k_p^2(y) u_m(y) e^{i (p K + \alpha_m) x} = \sum_n \left( \sum_m k_{n-m}^2(y) u_m(y) \right) e^{i \alpha_n x}$$

Substituting this back into the full differential equation gives:

$$\sum_m \left( u_m''(y) - \alpha_m^2 u_m(y) \right) e^{i\alpha_m x} + \sum_n \left( \sum_m k_{n-m}^2(y) u_m(y) \right) e^{i \alpha_n x} = 0$$

Renaming $m \to n$ in the first sum allows factoring out $e^{i\alpha_n x}$:

$$\sum_n \left( u_n''(y) - \alpha_n^2 u_n(y) + \sum_m k_{n-m}^2(y) u_m(y) \right) e^{i \alpha_n x} = 0$$

Since the set of spatial harmonics $\{ e^{i\alpha_n x} \}$ forms an orthogonal basis, this equality must hold independently for each harmonic mode $n$:

$$u_n''(y) = \alpha_n^2 u_n(y) - \sum_m k_{n-m}^2(y) u_m(y)$$

By introducing the Kronecker delta:

$$\delta_{nm} = \begin{cases} 1 & \text{if } m = n \\ 0 & \text{if } m \neq n \end{cases}$$

we can rewrite $\alpha_n^2 u_n(y)$ as a sum $\sum_m \alpha_n^2 \delta_{nm} u_m(y)$ and combine both terms:

$$u_n''(y) = \sum_m \left( \alpha_n^2 \delta_{nm} u_m(y) \right) - \sum_m \left( k_{n-m}^2(y) u_m(y) \right)$$

Factoring out $u_m(y)$ yields the final system of coupled differential equations:

$$u_n''(y) = \sum_m \left[ \alpha_n^2 \delta_{nm} - k^2_{n-m}(y) \right] u_m(y)$$

By defining the matrix operator $M_{nm}(y) \equiv \alpha_n^2 \delta_{nm} - k^2_{n-m}(y)$, this system can be expressed compactly in full vector-matrix notation, this simplifies to the linear system of second-order ordinary differential equations:

$$\mathbf{u}''(y) = \mathbf{M}(y) \mathbf{u}(y)$$

### 0.4 Physical Interpretation

The value of $\alpha_n^2 \delta_{nm}$ is only related to the tangential wave number, due $k^2_{n-m}(y)$ only the material structure.






-------


### 1.1 The grating equation — `computeAlphaAndBeta()`

A grating with period `d` illuminated at incidence angle `theta` can only scatter light into a discrete set of directions, indexed by an integer diffraction order `n` (from `-N` to `N`). This is the grating equation.

Every wave in the problem — incident, or any diffracted order — has a wavevector with a **tangential** component (parallel to the surface, along `x`) and a **normal** component (perpendicular, along `y`). The code:

```cpp
double alpha = k_2 * sin(theta_2) + 2 * M_PI * n / d;
alpha_[i] = alpha;
```

Physically: the incident wave contributes a fixed tangential wavevector `k·sin(theta)`; the grating's periodicity contributes an additional `2π·n/d` per order `n` (the "grating momentum kick"). This is the diffraction equation in wavevector form:

$$\alpha_n = k \sin\theta + n \cdot \frac{2\pi}{d}$$

`alpha_n` is the same in *every* layer of the structure — it depends only on geometry and the incident wave, never on `y` or on the local material. That is why it is computed once, before the per-layer loop even starts.

The normal component `beta_n` follows from the dispersion relation `k² = alpha_n² + beta_n²` in whichever homogeneous medium is being asked about:

```cpp
// beta2_: rayleigh expansion above grating.
double k22minusAn2 = k_2*k_2 - alpha*alpha;
if(k22minusAn2 >= 0)
    betaM_[i] = gsl_complex_rect(sqrt(k22minusAn2), 0);
else
    betaM_[i] = gsl_complex_rect(0, sqrt(-k22minusAn2));
```

$$\beta_n^{(M)} = \begin{cases}\sqrt{k^2 - \alpha_n^2} & \text{if } k^2 \ge \alpha_n^2 \quad\text{(propagating order)} \\[4pt] i\sqrt{\alpha_n^2 - k^2} & \text{if } k^2 < \alpha_n^2 \quad\text{(evanescent order)}\end{cases}$$

There are **two** beta arrays, computed with two different values of `k`, because they serve as boundary conditions at the two ends of the structure:

- `betaM_` — computed with the vacuum wavenumber. This is the boundary condition at the *top* (where the outgoing diffracted orders live — "M" for the top medium in the layer numbering, see §1.2).
- `beta1_` — computed with the substrate refractive index `v_1_`. This is the boundary condition at the *bottom*.

```cpp
// beta1_: rayleigh expansion inside grating
gsl_complex k12minusAn2 = gsl_complex_sub_real(gsl_complex_mul(k_1, k_1), alpha*alpha);
beta1_[i] = complex_sqrt_upperComplexPlane(k12minusAn2);
```

`complex_sqrt_upperComplexPlane()` (§6) picks the branch of the complex square root with `Im(w) ≥ 0` — physically, this enforces that waves *decay* (rather than grow) going deeper into an absorbing substrate.

Note that **no coating index (`v_c_`) appears anywhere in this function.** Coating only matters for the *local* material distribution inside the structure (§2.1, §3), not for these two asymptotic boundary conditions above and below everything.

### 1.2 Layering for numerical stability — `computeLayers()`

```cpp
// How many layers do we need?
double a = g_.totalHeight();

double magicNumber = 3;
// should be ln(1e15). However, emperically this is not enough to maintain
// stability (ex: REIXS LEG).  7 = ln(1e3) seems stable for all tests so far.

// How many layers to use? In order to keep size of exp(i betaM_{±N}) < 1e15 to avoid losing
// precision in double values compared with unity-size numbers.
numLayers_ = std::max( gsl_complex_abs(betaM_[0])*a/magicNumber, gsl_complex_abs(betaM_[2*N_])*a/magicNumber );
if(numLayers_ < 1)
    numLayers_ = 1;
    // we need at least one layer, in addition to the substrate.
```

Evanescent orders grow/decay as `exp(beta_n · y)` while propagating through the structure. A `double` carries about 15–16 significant decimal digits; if `beta_n · (total height)` gets too large, `exp(beta_n · height)` exceeds that dynamic range, and subsequent matrix operations lose essentially all precision comparing it against unity-sized numbers.

The fix: don't integrate the whole height `a` in one shot. Split it into `numLayers_` thinner slices, chosen so that even the *most* evanescent order (largest `|beta_n|`, which occurs at the highest or lowest computed diffraction order — hence checking both `betaM_[0]` and `betaM_[2*N_]`) stays within a safe exponential range over a single layer:

$$\text{numLayers} = \max\left(\frac{|\beta_0^{(M)}| \cdot a}{\text{magicNumber}}, \; \frac{|\beta_{2N}^{(M)}| \cdot a}{\text{magicNumber}}\right)$$

The comment here is worth calling out explicitly: the theoretically "correct" bound would be `magicNumber ≈ ln(10^15) ≈ 34.5`, but the code uses `3` — a much more conservative choice — because `34.5`, and even the less strict `ln(1000) ≈ 7`, turned out to be *insufficient* in a real test case (referenced in the comment as "REIXS LEG"). This heuristic is not a cosmetic detail: it is the mechanism that keeps the whole algorithm numerically stable regardless of how many "real" material regions (coating, profile, substrate) sit inside the structure. A thin coating barely changes `a`, so it barely changes `numLayers_` — the scheme adapts automatically rather than becoming unstable.

```cpp
y_ = new double[M_];
// To be consistent with the text, we number the lowest layer as 1.  y_[0] is
// therefore unused, so that we can use y_m = y_[m], with lowest m=1, highest m=M-1.

for(int m=1; m<M_; ++m) {
    y_[m] = double(m-1)/numLayers_*a;
}
```
`y_` holds the layer boundaries, evenly spaced from `y=0` (substrate interface) to `y=a` (top of the structure, into vacuum). `M_ = numLayers_ + 2`, matching the convention that layer boundary `m` runs from 1 to `M_-1`.

---

## 2. Inside a Layer: the Coupled-Wave ODE

### 2.1 Field ansatz and the differential equation

Inside the grating region, the local permittivity (via `k² = (2π/λ)² · v(x,y)²`) is not uniform — it varies periodically in `x`, and (for non-rectangular or coated profiles) may also vary with `y`. A single plane wave cannot satisfy the wave equation there. Instead, the field is expanded as a Fourier series in `x`, with **height-dependent coefficients**:

$$u(x, y) = \sum_{n=-N}^{N} u_n(y)\, e^{i \alpha_n x}$$

Substituting into the (scalar) wave equation and matching Fourier components gives a set of coupled second-order ODEs in `y` for the `u_n(y)`:

$$u_n''(y) = \sum_{m=-N}^{N} \Big[ \alpha_n^2\, \delta_{nm} \;-\; k^2_{n-m}(y) \Big]\, u_m(y)$$

where `k²_{n-m}(y)` is the `(n-m)`-th Fourier coefficient of the local permittivity profile `k²(x,y)` at height `y` (computed by `computeGratingExpansion()`, §3). This is exactly what `odeFunction()` implements:

```cpp
// w contains the last values of u_n{re, im} and u'_n{re, im}, in that order.
// need to compute f = dw/dy = u'_n{re, im} followed by u''_n{re, im}
...
gsl_complex M_nm = gsl_complex_rect(0,0);
if(n-m >= -N_ && n-m <= N_)
    M_nm = gsl_complex_mul_real(localK2[n-m + N_], -1.0);	// -k^2_{n-m}
if(n == m)
    M_nm = gsl_complex_add_real(M_nm, alpha2);

upp_n = gsl_complex_add(upp_n, gsl_complex_mul(M_nm, u_m));
```

**This is the physical heart of the whole method.** `k²_{n-m}(y)` is exactly where the grating's cross-sectional shape — including whether and how it's coated — enters the calculation. Change the shape, and this coupling matrix changes, and so does the solution. Everything else described in this document (T-matrix, S-matrix, alpha/beta) is generic machinery built around this one physical ingredient.

### 2.2 State vector layout

GSL's ODE solver works with a flat array of `double`, not complex vectors, so `u_n` and `u_n'` are packed into one array `w` of size `8N+4`:

```
w = [ u_{-N}.re, u_{-N}.im, ..., u_N.re, u_N.im,   u'_{-N}.re, u'_{-N}.im, ..., u'_N.re, u'_N.im ]
    \_____________ fourNp2_ values _____________/ \_______________ fourNp2_ values _______________/
```

`odeFunction()`'s job ("compute dw/dy") then splits naturally in two:
```cpp
// fourNp2 divides the top and bottom of the arrays w, dwdy.  Top of w is u; bottom of w is u' = v.
// Top of dwdy is u';  bottom of dwdy is u''.
// total size is (2N+1)x2x2, ie: 8N+4.
for(int i=0; i<fourNp2_; ++i) {
    // working on computing u'_n [top of dwdy array]. Just copy from u' = [bottom half of w array]
    dwdy[i] = w[i + fourNp2_];
}
```
The top half of `dwdy` is simply a copy of the bottom half of `w` — by definition, `d(u_n)/dy = u_n'`. The bottom half of `dwdy` is the coupling sum from §2.1.

### 2.3 The Jacobian, and its known limitation

Implicit ODE methods (like the Adams-Moulton corrector used here — §2.4) need the Jacobian `∂f/∂w` to iterate efficiently. `odeJacobian()` builds it from exactly the same `M_nm` coupling matrix, laid out as a block structure — an identity block (since `∂u'_n/∂u'_n = 1`) plus the coupling block (`∂u''_n/∂u_m = M_nm`):

```cpp
// go through top rows of jac [i=0,fourNp2]. Set ident. matrix in upper right-hand block.
for(int i=0; i<fourNp2_; ++i) {
    dfdw[i*eightNp4_+fourNp2_+i] = 1.0;		// set at index dfdw(i, fourNp2+i).
}
...
// set 4x4 matrix here: [M_re, -M_im; M_im, M_re] at (i,j), (i, j+1); (i+1, j), (i+1, j+1)
dfdw[i*eightNp4_ + j] = GSL_REAL(M_nm);
dfdw[i*eightNp4_ + j + 1] = -GSL_IMAG(M_nm);
dfdw[(i+1)*eightNp4_ + j] = GSL_IMAG(M_nm);
dfdw[(i+1)*eightNp4_ + j + 1] = GSL_REAL(M_nm);
```

The one place this is **not** exact: `dfdy` — the *explicit* `y`-dependence of the right-hand side, holding `u` fixed — is left at zero, with the code's own flag on it:

```cpp
/// \todo IMPORTANT! Leaving dfdy = 0 for now. This is only true in case of rectangular grating...
```

This is exact when `k²(x,y)` doesn't vary continuously with `y` within a layer — true for **rectangular** profiles (vertical walls: `xIntersection1()`/`xIntersection2()` don't depend on `y`), including a rectangular profile with an interpenetrating coating (only the *value* of `k²` jumps at specific `y`, never the *x*-positions). It is **not** exact for blazed, sinusoidal, trapezoidal, or custom profiles, where the wall position genuinely depends on `y`. In practice, an approximate Jacobian mainly affects the convergence speed/robustness of the implicit solver's inner iteration rather than the correctness of a converged result — but it remains a real, documented gap for non-rectangular shapes.

### 2.4 Solving one layer — `integrateTrialSolutionAlongY()`

```cpp
// define ode solving system, with our function to evaluate dw/dy, the Jacobian, and 8*N_+4 components.
gsl_odeiv2_system odeSys = {odeFunctionCB, odeJacobianCB, eightNp4_, this};

// initial starting step in y: choose grating height / 200.
double hStart = (yEnd - yStart)/200;

// setup driver
gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_standard_new(&odeSys, gsl_odeiv2_step_msadams,
    hStart, integrationTolerance_, integrationTolerance_, 0.5, 0.5);
// Variable-coefficient linear multistep Adams method in Nordsieck form. Uses explicit
// Adams-Bashforth (predictor) and implicit Adams-Moulton (corrector) methods in P(EC)^m
// functional iteration mode.
```

GSL's variable-coefficient linear multistep Adams method (predictor-corrector), starting with an initial step of `(layer height)/200`. `integrationTolerance_` (from `PEMathOptions`) controls both the absolute and relative error targets. If the solver can't converge to that tolerance, or the geometry function returns an error mid-integration (invalid profile, `y` above the profile height), the caller gets back a `ConvergenceFailure`:

```cpp
if (status != GSL_SUCCESS) {
    if(status == GSL_EBADFUNC)
        std::cout << "ODE: Integration failure: Invalid Geometry. Check your grating geometry specification." << std::endl;
    else if(status == GSL_FAILURE)
        std::cout << "ODE: Integration failure: Can't achieve step tolerance required. Try changing the integration tolerance." << std::endl;
    ...
    return PEResult::ConvergenceFailure;
}
```

### 2.5 Trial solutions and the T-matrix — `computeTMatrixBelowLayer()`

A single layer, by itself, is a linear boundary-value problem. For `2N+1` orders, you need `2·(2N+1)` independent solutions to span the full solution space — one "trial solution" per basis vector `u_n(y_start) = δ_{n,p}` (§2.6 covers why: this becomes the columns of a transfer matrix). `setIntegrationStartingValues()` sets up exactly this delta-function initial condition, once per trial index `p`:

```cpp
// set u[j] = 1.  Multiplication by 2 is due to {re,im}.
w[2*j] = 1.0;
gsl_complex uprime = gsl_complex_mul_imag(m == 1 ? beta1_[j] : betaM_[j], secondRound ? 1 : -1);
```
(using `beta1_` instead of `betaM_` specifically for the bottom-most layer, `m == 1` — the substrate boundary condition from §1.1).

Each trial solution is integrated **independently** — embarrassingly parallel, which is exactly what the surrounding `#pragma omp parallel for` parallelizes over:

```cpp
// We now need 2*(2N+1) trial solutions.  j will be the loop index over p, but ranging from [0,4*N+1].
#pragma omp parallel for num_threads(numThreads_) schedule(dynamic) reduction(||:integrationFailureOccurred)
for(int j=0; j<fourNp2_; ++j) {
    double* w = wVectorForP(j);
    setIntegrationStartingValues(w, j, m-1);
    PEResult::Code status = integrateTrialSolutionAlongY(w, y_[m-1], y_[m]);
    ...
}
```

The results are assembled column-by-column into four `(2N+1)×(2N+1)` blocks `T11_, T12_, T21_, T22_` — together forming the transfer matrix relating field values at the bottom of the layer to field values at the top:

```cpp
// T12_ij = 0.5(u_ij - u'_ij / (i*betaM_n) )
// T22_ij = 0.5(u_ij + u'_ij / (i*betaM_n) )
gsl_matrix_complex_set(T12_, i, jj, gsl_complex_mul_real(gsl_complex_sub(*u_ij, temp), 0.5));
gsl_matrix_complex_set(T22_, i, jj, gsl_complex_mul_real(gsl_complex_add(*u_ij, temp), 0.5));
```

### 2.6 Why not just multiply T-matrices? — the S-matrix recursion

The straightforward way to combine many layers' T-matrices is to multiply them together. This is exactly the numerically unstable approach multilayer optics (and ellipsometry) is well known for: evanescent solutions grow exponentially in one direction and decay in the other, so naive multiplication mixes enormous and tiny numbers in the same matrix, and rounding error swamps the result.

The code avoids this by recursing on the **S-matrix** (the scattering relationship between incoming and outgoing amplitudes) instead of the T-matrix, which stays numerically well-conditioned regardless of how many layers there are:

```cpp
// For the first layer, we have Zinv_ = T11_.
// S12_ = T21_ Zinv_^{-1}
// S22_ = Zinv_^{-1}
...
// (for subsequent layers m = 3 .. M_-1:)
// Compute Zinv_ = T11_ + T12_ S12_.
gsl_matrix_complex_memcpy(Zinv_, T11_);
gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, one, T12_, S12_, one, Zinv_);
// Invert Zinv_...
// S12 = (T21 + T22 S12) Z
// S22 = S22 Z
```
(`Zinv_`/`Z_` naming note, from the header: *"we don't actually compute any inverses [of Z directly]; Zinv_ is directly calculated from Zinv^{q+1} = T11^{q+1} + T12^{q+1} S12^{1}, and then we use LU decomp and multiplication to avoid loss of precision in computing Z = Zinv_^{-1}."*)

Combined with the adaptive layering from §1.2, this S-matrix recursion is what lets the algorithm handle an arbitrary number of thin layers — however many a coating requires — without the result blowing up numerically.

---

## 3. Grating Geometry Input — `computeGratingExpansion()`

This is the bridge between the pure-geometry world of `PEG.cpp` and the pure-numerics world of §2.

### 3.1 Where the coating (and everything else about the shape) actually enters

```cpp
// wave number in the coating layer:
gsl_complex k_c = gsl_complex_mul(v_c_, k_M);
...
// Compute multistep function from grating:
int numSteps = g_.computeK2StepsAtY(y, k2_M, k2_1, k2_c, stepsX, stepsK2);
```

`PEGrating::computeK2StepsAtY()` (in `PEG.cpp`) answers one question: *"at this specific height `y`, sweeping across one period in `x`, which materials do I cross, and where?"* It returns a **multi-step function**: a list of `x`-positions and the `k²` value immediately to the right of each. Depending on the coating thickness relative to the profile height, this dispatches to `computeK2StepsAtY_noCoating()` (2 steps: vacuum, substrate), `_interpenetratingCoating()` (2 or 4 steps depending on whether `y` is below/within/above the coating-thickness range), or `_thickCoating()` (1, 2, or 2 steps across its three `y`-regions). This dispatch, and therefore everything downstream of it, is where a change in coating thickness or material actually changes the physics.

### 3.2 Fourier-decomposing the step function

The ODE (§2.1) needs `k²_n(y)`, the Fourier coefficients of that step function — not the step function itself. For a piecewise-constant function with jumps `σ_p = k²_{p+1} - k²_p` at positions `x_p` (with wraparound, `σ_{last} = k²_0 - k²_{last}`), the Fourier coefficients have a closed form:

$$k^2_0 = \frac{1}{d}\left(k^2_0 \cdot d \;-\; \sum_p \sigma_p\, x_p\right), \qquad\qquad k^2_n = \frac{-1}{2\pi n}\sum_p \sigma_p\, e^{i n K x_p} \quad (n \ne 0)$$

with `K = 2π/d`. This is exactly:

```cpp
// compute sigma values at crossings:
// sigma[p] = stepsK2[p+1] - stepsK2[p] for p<numSteps-1; sigma[numSteps-1]=sigma[0]-sigma[numSteps-1]
...
if(n == 0) {
    gsl_complex f0 = gsl_complex_mul_real(stepsK2[0], d);
    for(int p=0; p<numSteps; ++p)
        f0 = gsl_complex_sub(f0, gsl_complex_mul_real(sigma[p], stepsX[p]));
    k2[i] = gsl_complex_div_real(f0, d);
}
else {
    gsl_complex fn = gsl_complex_rect(0,0);
    for(int p=0; p<numSteps; ++p) {
        double nKx = n*K*stepsX[p];
        fn = gsl_complex_add(fn, gsl_complex_mul(sigma[p], gsl_complex_rect(sin(nKx), cos(nKx))));
    }
    k2[i] = gsl_complex_div_real(fn, -2*M_PI*n);
}
```

There's a special-cased shortcut for `numSteps == 1` — a fully homogeneous layer, e.g. deep inside a thick coating: all `n ≠ 0` coefficients are exactly zero, and `k²_0` is just that one material's value directly, no sum needed:
```cpp
// Optimization for numSteps = 1: f_n = 0 (n!=0).   f_0 = stepsK2[0].
if(numSteps == 1) {
    for(int i=0; i<twoNp1_; ++i)
        k2[i] = gsl_complex_rect(0,0);
    k2[N_] = stepsK2[0];
    return;
}
```

`\warning assumes numSteps is in [1, PEG_MAX_PROFILE_CROSSINGS]` — this and the sibling stack-array buffers (`stepsX`, `stepsK2`, `sigma`, all sized `PEG_MAX_PROFILE_CROSSINGS = 60`) are the reason coating geometry currently supports at most a handful of material crossings per height slice; see the note on raw fixed-size buffers in the companion modernization notes for this file.

---

## 4. Assembling the Result — `getEff()`

`getEff()` is the orchestrator; here's its flow, matching the numbered comments already in the code:

1. **Setup** — look up `v_1_` (always) and `v_c_` (only if `g_.coatingThickness() > 0`); fail with `PEResult::MissingRefractiveDataFailure` if either lookup comes back empty (`(0,0)`, the sentinel `PEGrating::refractiveIndex()` uses for "not found").
2. **`computeAlphaAndBeta()` then `computeLayers()`** — §1.1, §1.2. (Comment in the code: *"Calculates how many vertical layers we need, and the division into slices at y_."*)
3. **Recursive S-matrix build** — *"Recursive computation of S-matrix below each layer."* The first layer (`m=2`) is a special case (comment: *"Handle first layer separately, as a special case."* — `Zinv_ = T11_` directly, since there is no previous S-matrix to fold in yet); layers `m = 3 .. M_-1` use the general recursion from §2.6 (comment: *"First layer done. Handle subsequent layers"*).
4. **`computeBMFromSMatrix()`** — *"Calculate B_n^M from center column of S matrix * exp(...)."* Extracts the outgoing reflected amplitudes `BM_n` from the finished S-matrix:
```cpp
BM_[i] = gsl_complex_mul(
            gsl_matrix_complex_get(S12_, i, N_),
            gsl_complex_exp(gsl_complex_mul_imag(gsl_complex_add(betaM_[i], betaM_[N_]), -a)));
```
(the `exp(...)` factor re-references the phase to a common origin at `y = 0`).

Then, *"Now we have BM_. Compute efficiency and put into result structure"*:

```cpp
result.eff[i] =  gsl_complex_abs2(BM_[i])*GSL_REAL(betaM_[i])/GSL_REAL(betaM_[N_]);
// is this a non-propagating order?  Then the real part of beta2_n will be exactly 0, so the
// efficiency will come out as 0.
```

`|BM_n|²` is the *intensity* of order `n`. But different orders leave the structure at different angles, so intensity alone is not power — you need the "obliquity" projection factor `cos(angle_n)/cos(angle_0) = Re(β_n)/Re(β_0)`:

$$\eta_n = |B_n|^2 \cdot \frac{\mathrm{Re}(\beta_n^{(M)})}{\mathrm{Re}(\beta_N^{(M)})}$$

(`betaM_[N_]` is index `N_` in the zero-based array, i.e. order `n=0` — the specular/zeroth order — used as the reference angle.) For evanescent orders, `Re(β_n) = 0`, so `η_n` comes out to exactly zero automatically, with no special-casing needed — exactly the property the code's own comment above points out.

If `rmsRoughnessNm > 0`, every order is finally scaled by `PEGrating::roughnessFactor()` — a Névot–Croce-style Debye–Waller-type correction, using the coating index if present, the substrate index otherwise:
```cpp
double roughnessFactor = g_.roughnessFactor(rmsRoughnessNm/1000., wl, g_.coatingThickness() > 0 ? v_c_ : v_1_, incidenceDeg);
for(int i=0; i<twoNp1_; ++i)
    result.eff[i] = roughnessFactor*result.eff.at(i);
```
The formula itself (`PEGrating::roughnessFactor()`) is documented in `PEG.h`, not repeated here.

When `measureTiming_` is set, `getEff()` reports where the time actually went, in the same seven buckets it accumulates into throughout the function (`timing_[0]` allocation, `[1]` refractive-index lookup, `[2]` alpha/beta/layers, `[3]` numerical integration, `[4]` matrix operations, `[5]` computing `BM_`, `[6]` packaging efficiencies) — useful for judging where a slow calculation (e.g. many thin layers from a demanding coating) is actually spending its time.

---

## 5. Threading Model

- `numThreads_` OpenMP threads are used in three places: `computeAlphaAndBeta()`, the trial-solution loop inside `computeTMatrixBelowLayer()`, and `computeBMFromSMatrix()`.
- Each thread needs its *own* `k2[]` scratch array, since `computeGratingExpansion()` gets called from inside the parallel trial-solution loop at a different `y` per solution. That's what `k2_` (one array per thread, allocated in the constructor: *"we need one k2_ array for each thread, since they will be used simultaneously"*) and `k2ForCurrentThread()` (*"Return based on OpenMP current thread"*) are for.
- Similarly, each of the `4N+2` trial solutions gets its own non-overlapping slice of the shared `wVectors_` buffer via `wVectorForP(j) = wVectors_ + eightNp4_*j` — "Returns the wVector for a trial solution `p` at index `j`."

---

## 6. Small Helpers

- **`complex_sqrt_upperComplexPlane()`** — complex square root with the branch cut chosen so `Im(w) ≥ 0`, the physically correct choice for a wave that *decays* (rather than grows) going into an absorbing medium. This differs from the "principal" square-root branch, which instead guarantees `Re(w) ≥ 0`.
- **`linalg_LU_complex_solve()`** — solves `AX = B` for a whole matrix `X` at once, given a pre-computed LU decomposition of `A`, by looping GSL's single-column solve (`gsl_linalg_complex_LU_solve`) over each column of `B`/`X`.
- **`conditionNumber()`** — **incomplete**, per its own doc comment (`INCOMPLETE!`). It computes an LU decomposition and inverse of the input matrix but never actually computes a norm from them, and unconditionally returns `1` (with a leftover comment: *"how to get the norm of a complex matrix?"*). It does not appear to be called anywhere in the solver's own code path.

---

## 7. Open Items (relevant to future work, including TM)

- **`dfdy = 0` in `odeJacobian()`** is only exact for rectangular profiles (§2.3) — a documented, not-yet-resolved approximation for blazed/sinusoidal/trapezoidal/custom shapes.
- **The entire derivation in §2.1 is the scalar wave equation — this is TE polarization specifically.** TM polarization needs the analogous equation formulated with a spatially-varying `1/ε(x,y)` rather than `k²(x,y)` (since it is `H`, not `E`, whose boundary conditions differ across a material interface) — a genuinely different coupled ODE, not a parameter tweak to the one derived here. This is the natural entry point for a TM implementation: a parallel `odeFunction`/`odeJacobian` pair, sharing the geometry/Fourier-expansion machinery in §3.
- **`conditionNumber()`** is unused and incomplete (§6) — candidate for removal or completion, but not on the critical path of any calculation today.
- /// \warning assumes numSteps is in [1, PEG_MAX_PROFILE_CROSSINGS]
- /// \todo IMPORTANT! Leaving dfdy = 0 for now. This is only true in case of rectangular grating...
