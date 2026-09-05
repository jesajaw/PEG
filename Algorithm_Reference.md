# TESolver Algorithm Reference

## Purpose and Scope

This document explains how `TESolver` computes the diffraction efficiency of a grating for a single (wavelength, incidence angle) pair. Or to be more acurate, it shows the physical validation how and why you can compute something likes this and also how it is (numerical) implemented. It does not covers the geometry aspect: it consumes this input from `PEG`. It also does not include the MPI orchestration in `mainMPI`, or the material refractive-index database lookup — those are documented separately.

The method implemented here is a **differential (coupled-wave) method**: the structure is treated as a stack of thin horizontal layers; within each layer, a coupled ordinary differential equation (ODE) system is integrated numerically along the vertical coordinate `y`; the per-layer results are combined into a numerically stable overall response via **S-matrix recursion**. This is closely related to the family of methods used in Rigorous Coupled-Wave Analysis (RCWA) and multilayer optics/ellipsometry.

---

## Table of Contents

0. [Physical Setup](#0-physical-setup)
    1. [Maxwell to Scalar Wavefunction](#01-maxwell-to-scalar-wavefunction)
    2. [Fourier & Floquet Expansion](#02-fourier-and-floquet-expansion)
    3. [Coupled Mode Equation Derivation](#03-coupled-mode-equation-derivation)
    4. [Physical Interpretation](#04-physical-interpretation)
    5. [Geometry Bridge](#05-geometry-bridge)
    6. [Fourier Coefficients of a Piecewise-Constant Layer](#06-fourier-coefficients-of-a-piecewise-constant-layer)
1. [Local Layer Response](#1-local-layer-response)
    1. [From a Second-Order System to a First-Order System](#11-from-a-second-order-system-to-a-first-order-system)
    2. [Dimension of the Solution Space](#12-dimension-of-the-solution-space)
    3. [Homogeneous Media as a Physical Reference Case](#13-homogeneous-media-as-a-physical-reference-case)
    4. [Propagating and Evanescent Orders](#14-propagating-and-evanescent-orders)
    5. [Two Distinct Sets of Asymptotic Parameters](#15-two-distinct-sets-of-asymptotic-parameters)
2. [Boundary Conditions and Modal Amplitudes](#2-boundary-conditions-and-modal-amplitudes)
    1. [The Actual Boundary Conditions](#21-the-actual-boundary-conditions)
    2. [Modal Amplitudes from u and u'](#22-modal-amplitudes-from-u-and-u)
3. [Transfer Matrix](#3-transfer-matrix)
    1. [What a Single Layer Does Physically](#31-what-a-single-layer-does-physically)
    2. [Why a Basis of Trial Solutions Generates the Transfer Matrix](#32-why-a-basis-of-trial-solutions-generates-the-transfer-matrix)
    3. [Block Structure of the Transfer Matrix](#33-block-structure-of-the-transfer-matrix)
    4. [Why Direct Multiplication of Transfer Matrices Is Numerically Problematic](#34-why-direct-multiplication-of-transfer-matrices-is-numerically-problematic)
4. [Scattering Matrix](#4-scattering-matrix)
    1. [Physical Idea of the S-Matrix](#41-physical-idea-of-the-s-matrix)
    2. [Recursion Over Layers](#42-recursion-over-layers)
5. [Layering as a Numerical Discretization](#5-layering-as-a-numerical-discretization)
6. [The Complete Physical Flow Through the Solver](#6-the-complete-physical-flow-through-the-solver)
    1. [Origin of the Single-Order Incident Condition](#61-origin-of-the-single-order-incident-condition)
7. [From Amplitude to Power: Diffraction Efficiency](#7-from-amplitude-to-power-diffraction-efficiency)
    1. [Phase Referencing to a Common Origin](#71-phase-referencing-to-a-common-origin)
    2. [From Field Amplitude to Power](#72-from-field-amplitude-to-power)
    3. [Why Evanescent Orders Automatically Have Zero Efficiency](#73-why-evanescent-orders-automatically-have-zero-efficiency)
8. [Truncation and Its Physical Meaning](#8-truncation-and-its-physical-meaning)
    1. [Finite Mode Number N](#81-finite-mode-number-n)
    2. [Relation Between Evanescent Decay and Discretization Scale](#82-relation-between-evanescent-decay-and-discretization-scale)
9. [The Coating Causal Chain](#9-the-coating-causal-chain)
10. [What This Formalism Does Not Do](#10-what-this-formalism-does-not-do)
11. [Physical Consistency Checks](#11-physical-consistency-checks)
    1. [Homogeneous-Grating Limit](#111-homogeneous-grating-limit)
    2. [Energy Balance](#112-energy-balance)
    3. [Complex Quantities and Absorption](#113-complex-quantities-and-absorption)
12. [Full Equation Chain (Summary)](#12-full-equation-chain-summary)
13. [Conclusion](#13-conclusion)
14. [Mathematical Profs](#14-mathematical-profs)
    1. [From Maxwell's Equations to the Scalar TE Helmholtz Equation](#141-from-maxwells-equations-to-the-scalar-te-helmholtz-equation)
    2. [From Periodicity to the Floquet Expansion](#142-from-periodicity-to-the-floquet-expansion)
    3. [Fourier Expansion of the Material Function](#143-fourier-expansion-of-the-material-function)
    4. [Derivation of the Coupled-Mode Equation](#144-derivation-of-the-coupled-mode-equation)
    5. [Why the Material Fourier Coefficients Couple Different Orders](#145-why-the-material-fourier-coefficients-couple-different-orders)
    6. [Piecewise-Constant Fourier Coefficients](#146-piecewise-constant-fourier-coefficients)
    7. [Homogeneous Medium as a Limiting Case](#147-homogeneous-medium-as-a-limiting-case)
    8. [Derivation of the Modal-Amplitude Transformation](#148-derivation-of-the-modal-amplitude-transformation)
    9. [Propagating and Evanescent Orders](#149-propagating-and-evanescent-orders)
    10. [First-Order Form of the Coupled ODE](#1410-first-order-form-of-the-coupled-ode)
    11. [Fundamental Matrix and Transfer Matrix](#1411-fundamental-matrix-and-transfer-matrix)
    12. [Why the Transfer Matrix Has Dimension $2Q\times2Q$](#1412-why-the-transfer-matrix-has-dimension-2qtimes2q)
    13. [From Transfer Variables to Scattering Variables](#1413-from-transfer-variables-to-scattering-variables)
    14. [Derivation of the Redheffer Star Product](#1414-derivation-of-the-redheffer-star-product)
    15. [Why Scattering-Matrix Composition Is Numerically Better](#1415-why-scattering-matrix-composition-is-numerically-better)
    16. [Incident-Order Selection](#1416-incident-order-selection)
    17. [Derivation of the Diffraction-Efficiency Formula](#1417-derivation-of-the-diffraction-efficiency-formula)
    18. [Why Evanescent Orders Carry Zero Far-Field Power](#1418-why-evanescent-orders-carry-zero-far-field-power)
    19. [Energy Conservation](#1419-energy-conservation)
    20. [Complex Refractive Index and Complex Propagation Constants](#1420-complex-refractive-index-and-complex-propagation-constants)
    21. [Phase Referencing](#1421-phase-referencing)
    22. [Vertical Layering as a Numerical Approximation](#1422-vertical-layering-as-a-numerical-approximation)
    23. [Why the Evanescent Spectrum Controls the Vertical Resolution](#1423-why-the-evanescent-spectrum-controls-the-vertical-resolution)
    24. [Mode Truncation](#1424-mode-truncation)
    25. [The Complete Mathematical Chain](#1425-the-complete-mathematical-chain)
    26. [The Most Important Mathematical Dependencies](#1426-the-most-important-mathematical-dependencies)
    27. [Final Mathematical Interpretation](#1427-final-mathematical-interpretation)

And may see [TMSolver](#TMSolver).
---

## 0. Physical Setup

Overview of the Theoretical Pipeline:

$$\text{Maxwells Equations} \longrightarrow \text{Helmholtz Equation} \longrightarrow \text{Scalar Wave Equation (TE)}$$
$$\downarrow$$
$$\text{Periodic Fourier/Floquet Expansion} \longrightarrow \text{Coupled ODEs in } y \longrightarrow \text{Local Layer Response}$$
$$\downarrow$$
$$\text{Transfer Matrix} \longrightarrow \text{Scattering Matrix} \longrightarrow \text{Reflection Amplitudes} \longrightarrow \text{Diffraction Efficiencies}$$

Parameter Categorization:

- Geometry, Grating, and Material Structure: Defined by the spatial permittivity distribution $\varepsilon_r(x,y)$ and non-magnetic permeability ($\mu = \mu_0$). This encompasses the grating period $d$, layer/coating thicknesses along $y$ and the physical profile/surface topology of the grating.
- Incident Field Parameters: Free-space wavelength $\lambda$ and angle of incidence $\theta$


### 0.1 Maxwell to Scalar Wavefunction

Faradays law of induction and Ampères circuital law:
$$\vec{\nabla} \times \vec{E} = -\frac{\partial \vec{B}}{\partial t}, \quad \vec{\nabla} \times \vec{H} = -\vec{i} + \frac{\partial \vec{D}}{\partial t}$$

with $\vec{D} = \varepsilon_0 \varepsilon_r \vec{E}$ and $\vec{B} = \mu \vec{H}$ (assuming non-magnetic media where $\mu = \mu_0$ and $\vec i = \vec 0$), which can be formulated as:
$$\vec{\nabla} \times \vec{E} = -\mu_0 \frac{\partial \vec{H}}{\partial t}, \quad \vec{\nabla} \times \vec{H} = \varepsilon_0 \varepsilon_r \frac{\partial \vec{E}}{\partial t}$$

Separation of time and space assuming time-harmonic fields $\vec{E} e^{-i\omega t}$ and $\vec{H} e^{-i\omega t}$:
$$\vec{\nabla} \times \vec{E} = i\omega\mu_0\vec{H}, \quad \vec{\nabla} \times \vec{H} = -i \omega \varepsilon_0\varepsilon_r\vec{E}$$

To eliminate $\vec{H}$, take the curl of Faradays law:
$$\vec{\nabla} \times (\vec{\nabla} \times \vec{E}) = i\omega\mu_0 (\vec{\nabla} \times \vec{H})$$

Substituting Ampères law into the right-hand side yields:
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
  $$\frac{\partial^2}{\partial y^2} u(x,y) = \sum_m u_m^{\prime\prime}(y) e^{i\alpha_m x}$$

Combining these gives:

$$\sum_m u_m^{\prime\prime}(y) e^{i\alpha_m x} - \sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x} + k^2(x,y) u(x,y) = 0$$

where $k^2(x,y) \equiv k_0^2 \varepsilon_r(x,y)$ is the periodic wavenumber profile. Due to spatial periodicity, $k^2(x,y)$ can also be expanded into a Fourier series:

$$k^2(x,y) = \sum_p k_p^2(y) e^{i p K x}, \quad K = \frac{2\pi}{d}$$

Inserting this expansion into the material product term gives:

$$\sum_m u_m^{\prime\prime}(y) e^{i\alpha_m x} - \sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x} + \sum_p k_p^2(y) e^{i p K x} \sum_{m} u_m(y) e^{i\alpha_m x} = 0$$

$$\sum_m u_m^{\prime\prime}(y) e^{i\alpha_m x} - \sum_m \alpha_m^2 u_m(y) e^{i\alpha_m x} + \sum_p \sum_m k_p^2(y) u_m(y) e^{i (p K + \alpha_m) x} = 0$$

$$\sum_m \left( u_m^{\prime\prime}(y) - \alpha_m^2 u_m(y) \right) e^{i\alpha_m x} + \sum_p \sum_m k_p^2(y) u_m(y) e^{i (p K + \alpha_m) x} = 0$$

Using the relation for the allowed tangential wavenumbers $\alpha_m = k\sin\theta + m K$ ([Fourier](#02-fourier)), we observe that $p K + \alpha_m = \alpha_{m+p}$. Performing an index substitution by setting $n = m + p$ (or equivalently $p = n - m$), the second summation becomes:

$$\sum_p \sum_m k_p^2(y) u_m(y) e^{i (p K + \alpha_m) x} = \sum_n \left( \sum_m k_{n-m}^2(y) u_m(y) \right) e^{i \alpha_n x}$$

Substituting this back into the full differential equation gives:

$$\sum_m \left( u_m^{\prime\prime}(y) - \alpha_m^2 u_m(y) \right) e^{i\alpha_m x} + \sum_n \left( \sum_m k_{n-m}^2(y) u_m(y) \right) e^{i \alpha_n x} = 0$$

Renaming $m \to n$ in the first sum allows factoring out $e^{i\alpha_n x}$:

$$\sum_n \left( u_n^{\prime\prime}(y) - \alpha_n^2 u_n(y) + \sum_m k_{n-m}^2(y) u_m(y) \right) e^{i \alpha_n x} = 0$$

Since the set of spatial harmonics $\{ e^{i\alpha_n x} \}$ forms an orthogonal basis, this equality must hold independently for each harmonic mode $n$:

$$u_n^{\prime\prime}(y) = \alpha_n^2 u_n(y) - \sum_m k_{n-m}^2(y) u_m(y)$$

By introducing the Kronecker delta:

$$\delta_{nm} = \begin{cases} 1 & \text{if } m = n \\ 0 & \text{if } m \neq n \end{cases}$$

we can rewrite $\alpha_n^2 u_n(y)$ as a sum $\sum_m \alpha_n^2 \delta_{nm} u_m(y)$ and combine both terms:

$$u_n^{\prime\prime}(y) = \sum_m \left( \alpha_n^2 \delta_{nm} u_m(y) \right) - \sum_m \left( k_{n-m}^2(y) u_m(y) \right)$$

Factoring out $u_m(y)$ yields the final system of coupled differential equations:

$$u_n^{\prime\prime}(y) = \sum_m \left[ \alpha_n^2 \delta_{nm} - k^2_{n-m}(y) \right] u_m(y)$$

By defining the matrix operator $M_{nm}(y) \equiv \alpha_n^2 \delta_{nm} - k^2_{n-m}(y)$, this system can be expressed compactly in full vector-matrix notation, this simplifies to the linear system of second-order ordinary differential equations:

$$\mathbf{u}^{\prime\prime}(y) = \mathbf{M}(y) \mathbf{u}(y)$$

### 0.4 Physical Interpretation

The coupling matrix
$$M_{nm}(y) = \alpha_n^2\delta_{nm} - k^2_{n-m}(y)$$
splits into two contributions with distinct physical origins:

- **Diagonal term** $\alpha_n^2\delta_{nm}$ — depends only on the tangential wavenumber of the $n$-th Floquet mode. It is identical in every layer, independent of what material occupies that layer.
- **Off-diagonal term** $-k^2_{n-m}(y)$ — depends only on the spatial material distribution at height $y$, through the Fourier coefficients of $\varepsilon_r(x,y)$.

Hence $M_{nm}(y)$ represents the coupling between Floquet modes induced by the periodic permittivity profile.

**Limiting case (no lateral structure).** If the layer has no lateral structuring, $k_p^2 = 0$ for all $p \neq 0$, and the coupled system decouples into independent equations
$$u_n^{\prime\prime}(y) + \left(k_0^2 - \alpha_n^2\right) u_n(y) = 0,$$
i.e. every Floquet order propagates as an independent plane wave in a homogeneous medium. This is precisely where the physical distinction between a homogeneous slab and a grating layer enters the formalism: **mode coupling exists if and only if the layer has lateral material contrast.**

### 0.5 Geometry Bridge

At a fixed depth $y$, slicing the grating horizontally yields a 1D cross section. Over one period $d$, $k^2(x,y)$ is typically piecewise constant in $x$, taking on the values of whichever materials are present at that height — e.g. vacuum, substrate, one or more coating layers, or several material transitions in sequence.

The geometry therefore supplies, for each $y$, a function
$$x \mapsto k^2(x,y)$$
consisting of constant segments separated by jump discontinuities. Extracting this function is a purely geometric step: for a given $y$, one determines which material occupies each interval and at which $x$-positions the transitions between materials occur. This is the only place where the concrete grating profile enters the formalism — it produces the raw step function that subsequently becomes the input to the Fourier decomposition below. How this step function is obtained algorithmically depends on the profile representation and is deliberately left open here.

### 0.6 Fourier Coefficients of a Piecewise-Constant Layer

Once the step function $k^2(x,y)$ is known for a given $y$, it is not the profile itself that enters the coupled-mode system, but its Fourier decomposition. Let the jump positions be $x_p$ and the jump magnitudes
$$\sigma_p = k^2_{right} - k^2_{left}$$
at each transition. The Fourier coefficients are
$$k_n^2(y) = \frac{1}{d}\int_0^d k^2(x,y)\,e^{-inKx}\,dx, \qquad K = \frac{2\pi}{d}.$$

For $n \neq 0$, since the function is constant between jumps, the integral reduces via integration by parts to a sum over the jump locations alone:
$$k_n^2(y) = -\frac{1}{2\pi i n} \sum_p \sigma_p\, e^{-inKx_p}.$$

For $n = 0$, the coefficient is simply the average value over one period:
$$k_0^2(y) = \frac{1}{d}\int_0^d k^2(x,y)\,dx.$$

The key physical point: a change in coating thickness, profile shape, or material choice only ever enters the ODE system through this route — by changing, at each $y$, the step positions $x_p$ and jump sizes $\sigma_p$, and hence the resulting Fourier coefficients $k_n^2(y)$. No aspect of the physical layer geometry is used by the coupled-mode equations except through these coefficients.

---

## 1. Local Layer Response

### 1.1 From a Second-Order System to a First-Order System

Numerically — and for the general solution theory below — it is convenient to rewrite
$$u^{\prime\prime} = Mu$$
as a system of first order. Defining
$$v = u^{\prime},$$
the second-order system splits into two coupled first-order equations:
$$u^{\prime} = v, \quad v^{\prime} = Mu.$$

Together, these form a single first-order block system:
$$\frac{d}{dy}\begin{pmatrix} u \\ v \end{pmatrix}
= \begin{pmatrix} 0 & I \\ M(y) & 0 \end{pmatrix}
\begin{pmatrix} u \\ v \end{pmatrix}.$$

The state vector at depth $y$ is therefore
$$w = \begin{pmatrix} u \\ u^{\prime} \end{pmatrix},$$
combining both the field amplitudes and their derivatives for every Floquet order. This is the physical reason why the state carries $2(2N+1)$ complex degrees of freedom: $2N+1$ Floquet orders, each contributing one amplitude and one derivative.

### 1.2 Dimension of the Solution Space

There are $2N+1$ Floquet modes retained in the truncated expansion. Each mode obeys a differential equation that is second order in $y$. A second-order linear ODE has a two-dimensional solution space; consequently, the full coupled system possesses
$$2(2N+1)$$
linearly independent solutions. This is not a numerical convenience but a mathematical property of the differential equation itself — it holds regardless of how (or whether) the system is solved numerically.

The full solution space at a given $y$ can therefore be spanned by a set of basis solutions,
$$w^{(1)}(y),\ w^{(2)}(y),\ \dots,\ w^{(2(2N+1))}(y),$$
commonly referred to as **trial solutions**. Any physically admissible field profile within the layer is a linear combination of these trial solutions,
$$w(y) = \sum_{i=1}^{2(2N+1)} c_i\, w^{(i)}(y),$$
for some set of coefficients $c_i$ fixed by the boundary conditions between layers — determining these coefficients is the subject of the next section.

### 1.3 Homogeneous Media as a Physical Reference Case

To understand the boundary conditions used later, consider first a homogeneous medium (no lateral structure). In this case only the zeroth Fourier coefficient of $k^2$ is nonzero, so the coupling matrix becomes diagonal ([Fourier Coefficients of a Piecewise-Constant Layer](#06-fourier-coefficients-of-a-piecewise-constant-layer)):
$$M_{nm} = (\alpha_n^2 - k^2)\,\delta_{nm}.$$
Each Floquet order then obeys an independent equation. Defining
$$\beta_n^2 = k^2 - \alpha_n^2,$$
this becomes
$$u_n''(y) + \beta_n^2\, u_n(y) = 0,$$
with the two independent solutions
$$u_n(y) = A_n e^{i\beta_n y} + B_n e^{-i\beta_n y}, \qquad
u_n'(y) = i\beta_n A_n e^{i\beta_n y} - i\beta_n B_n e^{-i\beta_n y}.$$
This decomposition is the physical basis for the incoming and outgoing amplitudes used at the boundaries of the structure.

### 1.4 Propagating and Evanescent Orders

From $\beta_n^2 = k^2 - \alpha_n^2$, two regimes follow:

**Propagating.** If $|\alpha_n| \le |k|$, then $\beta_n$ is real,
$$\beta_n = \sqrt{k^2 - \alpha_n^2},$$
and $e^{\pm i\beta_n y}$ describes genuine traveling waves.

**Evanescent.** If $|\alpha_n| > |k|$, then $\beta_n = i\gamma_n$ with $\gamma_n > 0$, so
$$e^{i\beta_n y} = e^{-\gamma_n y}, \qquad e^{-i\beta_n y} = e^{+\gamma_n y}.$$
A physically admissible solution in a semi-infinite outer medium must not grow exponentially with distance from the structure. The branch of the square root must therefore be chosen so that the physically allowed direction decays — conventionally by requiring $\operatorname{Im}\beta_n \ge 0$.

### 1.5 Two Distinct Sets of Asymptotic Parameters

The grating is bounded above and below by two, generally different, homogeneous media. This gives two independent sets of propagation constants:
$$\beta_n^{\text{top}} = \sqrt{k_{\text{top}}^2 - \alpha_n^2}, \qquad
\beta_n^{\text{sub}} = \sqrt{k_{\text{sub}}^2 - \alpha_n^2}.$$
These are not two solutions of the same layer — they are the modal parameters of the two asymptotic outer half-spaces, determined solely by the top and substrate refractive indices. No interior layer (e.g. a coating) enters these asymptotic boundary parameters; they depend only on the media the structure is embedded in.

---

## 2. Boundary Conditions and Modal Amplitudes

### 2.1 The Actual Boundary Conditions

The differential equation alone does not determine the solution. In addition, the fields must satisfy the Maxwell boundary conditions at every interface. For TE polarization, this means continuity of the tangential electric field component and of the associated tangential magnetic field component (i.e. $u$ and, proportionally, $u'$).

In the scalar formalism used here, this requirement means that the admissible solutions at the outer boundaries must be expressible as a superposition of incoming and outgoing modes consistent with the asymptotic media on each side — the boundary values cannot be chosen freely; they must match the physically allowed waves in the half-spaces above and below the structure.

### 2.2 Modal Amplitudes from $u$ and $u'$

In a homogeneous medium, evaluating the general solution and its derivative at a boundary gives
$$u_n = A_n + B_n, \qquad u_n' = i\beta_n (A_n - B_n).$$
Inverting this relation:
$$A_n = \frac{1}{2}\left(u_n + \frac{u_n'}{i\beta_n}\right), \qquad
B_n = \frac{1}{2}\left(u_n - \frac{u_n'}{i\beta_n}\right).$$

This is a central step of the formalism:
$$\boxed{(u_n, u_n') \quad \Longleftrightarrow \quad (A_n, B_n)}$$
the two representations — field-and-derivative, and directional amplitude pairs — are fully equivalent. Physically, $A_n$ and $B_n$ correspond to the two opposite propagation directions for propagating modese $n$; which one is "incoming" and which is "outgoing" depends on which outer boundary is considered and on the sign convention chosen for $y$.

---

## 3. Transfer Matrix

### 3.1 What a Single Layer Does Physically

Consider a layer between $y = y_-$ and $y = y_+$. The coupled ODE propagates the full state
$$\mathbf{w}(y_-) = \begin{pmatrix} \mathbf{u}(y_-) \\ \mathbf{u}'(y_-) \end{pmatrix}$$
forward to $\mathbf{w}(y_+)$. Because the governing equation is linear, this propagation defines a linear map
$$\mathbf{w}(y_+) = \mathcal{T}\,\mathbf{w}(y_-),$$
the **transfer matrix $\mathcal{T}$** of the layer. Importantly, $\mathcal{T}$ is not "the physics of the grating" as a whole — it is only the linear map from one layer's state to the next. Think of each layer as a black box: it takes the state entering at one boundary and linearly transforms it into the state leaving at the other, according to its transfer matrix $\mathcal{T_n}$.

### 3.2 Why a Basis of Trial Solutions Generates the Transfer Matrix

Choose a basis $\mathbf{e}_1, \dots, \mathbf{e}_{2Q}$ (with $Q = 2N+1$) of the state space. For each basis vector, the coupled ODE system is solved starting from that vector as initial condition and integrated to $y_+$. The resulting end states form the columns of a fundamental solution matrix $\Phi(y)$. The transfer matrix is then
$$\mathcal{T} = \Phi(y_+)\,\Phi(y_-)^{-1}.$$
This is the mathematical meaning behind solving the system separately for many independent trial solutions, rather than for a single physical solution directly.

### 3.3 Block Structure of the Transfer Matrix

Writing the state in terms of directional amplitudes,

$$\begin{pmatrix} \mathbf{a}_+ \\ \mathbf{b}_+ \end{pmatrix} = \begin{pmatrix} T_{11} & T_{12} \\ T_{21} & T_{22} \end{pmatrix} \begin{pmatrix} \mathbf{a}_- \\ \mathbf{b}_- \end{pmatrix},$$

where $(\mathbf{a}, \mathbf{b})$ denote the appropriate directional amplitude pairs on each side. The precise labeling of the blocks is a matter of convention; what matters physically is that a single layer establishes a linear relation between four sets of amplitudes, two on each side of the layer $\left( R_\pm;T_\pm\right)$.

### 3.4 Why Direct Multiplication of Transfer Matrices Is Numerically Problematic

An evanescent mode can appear across a layer of thickness $h$ with factors as extreme as $e^{+\gamma h}$ and $e^{-\gamma h}$ simultaneously. A single transfer matrix can therefore contain both very large and very small numbers at once. For a stack of many layers,
$$T_{total} = T_M \cdots T_2 T_1,$$
such extreme values compound multiplicatively. Even when the final physical result is moderate, intermediate values in this product can become numerically catastrophic (loss of precision or overflow). This is the physical and numerical motivation for splitting the structure into many thin layers and combining their responses through a scattering-matrix recursion instead of direct transfer-matrix multiplication — a point taken up in the next chapter.

---

## 4. Scattering Matrix

### 4.1 Physical Idea of the S-Matrix

The transfer matrix $\mathcal{T}$ answers: how are the fields/modes on the top and bottom side of a layer related?

The scattering matrix $\mathcal{S}$ instead answers: which outgoing waves result from the incoming waves?
$$\begin{pmatrix} \mathbf{a}_{\text{top,out}} \\ \mathbf{a}_{\text{bot,out}} \end{pmatrix}= S \begin{pmatrix} \mathbf{a}_{\text{top,in}} \\ \mathbf{a}_{\text{bot,in}} \end{pmatrix}.$$

This formulation is physically much better conditioned: growing evanescent solutions never need to be tracked as primary global unknowns, since only physically bounded incoming and outgoing amplitudes appear as variables.

### 4.2 Recursion Over Layers

Suppose a new layer is appended to an already-assembled stack. The new layer contributes a transfer relation; the already-assembled stack below it is represented by its own scattering matrix. The amplitudes at the shared interface between the new layer and the stack are not directly observable — they are internal to the combined structure and can be eliminated algebraically. This elimination is precisely what produces the recursive combination rule for scattering matrices, commonly known as the Redheffer star product [REF!https://en.wikipedia.org/wiki/Redheffer_star_product].

Given two scattering matrices,
$$S_A = \begin{pmatrix} S_{11}^A & S_{12}^A \\ S_{21}^A & S_{22}^A \end{pmatrix}, \qquad
S_B = \begin{pmatrix} S_{11}^B & S_{12}^B \\ S_{21}^B & S_{22}^B \end{pmatrix},$$
representing two adjacent sub-structures (e.g. the stack assembled so far, and the newly added layer), the combined scattering matrix $S = S_A \star S_B$ is

$$S_{11} = S_{11}^A + S_{12}^A\left(I - S_{11}^B S_{22}^A\right)^{-1} S_{11}^B S_{21}^A,$$
$$S_{12} = S_{12}^A\left(I - S_{11}^B S_{22}^A\right)^{-1} S_{12}^B,$$
$$S_{21} = S_{21}^B\left(I - S_{22}^A S_{11}^B\right)^{-1} S_{21}^A,$$
$$S_{22} = S_{22}^B + S_{21}^B\left(I - S_{22}^A S_{11}^B\right)^{-1} S_{22}^A S_{12}^B,$$

where $I$ is an identity matrix, suitable for $S_A$ and $S_B$.

The essential physical point is independent of the exact algebraic form:
$$\boxed{\text{internal waves at the shared interface are eliminated; only the external input/output amplitudes remain.}}$$

The matrix inversions appearing above are, in practice, never carried out as explicit inversions but through numerically stable factorization methods — this is a numerical-stability concern, not a change to the underlying physics.

Think of the transfer matrix $\mathcal{T_n}$ as a black box for each layer $n$, carrying the field state from one boundary to the next, while the scattering matrix $\mathcal{S}$ combines these layers into a global black box that relates only the external incoming and outgoing waves.

---

## 5. Layering as a Numerical Discretization

The real permittivity profile is, in general, continuous in $y$. The solver approximates it by many thin layers, each treated as locally uniform along $y$ so that the coupled ODE can be integrated with locally constant coefficients within that slice.

In the limit $\Delta y \to 0$, the layered decomposition recovers the continuous differential equation exactly. Layering is therefore not an additional physical assumption about the material — it is a numerical discretization in the vertical direction. The physically motivated criterion for how finely to discretize is not primarily the location of material boundaries, but the maximum evanescent growth rate present in the structure and the numerical stability this requires ([Why Direct Multiplication of Transfer Matrices Is Numerically Problematic](#34-why-direct-multiplication-of-transfer-matrices-is-numerically-problematic)).

---

## 6. The Complete Physical Flow Through the Solver

The full solver can now be summarized as a physical chain:

1. **Incident wave.** Given $\lambda$ and $\theta$, compute $k_0 = 2\pi/\lambda$.
2. **Tangential momentum.** The grating periodicity restricts the allowed tangential wavenumbers to $\alpha_n = k_{\text{top}}\sin\theta + nK$.
3. **Local material distribution.** For each $y$, the geometry determines $\varepsilon_r(x,y)$.
4. **Fourier decomposition.** This yields the coefficients $k_n^2(y)$.
5. **Mode coupling.** The Fourier coefficients generate the coupled-mode system $u_n''(y) = \sum_m \left[\alpha_n^2\delta_{nm} - k_{n-m}^2(y)\right] u_m(y).$
6. **Vertical propagation.** This ODE system is integrated through the structure.
7. **Layer response.** Each layer's integration yields a transfer relation between its two boundaries.
8. **Assembly.** The layers are combined via scattering-matrix recursion ([Recursion Over Layers](#42-recursion-over-layers)).
9. **Outer boundary condition.** In the homogeneous medium above the structure, the solution is decomposed into physically outgoing diffraction orders.
10. **Reflection amplitudes.** The resulting coefficients are the reflected diffraction amplitudes, obtained by evaluating the final scattering matrix for the physically prescribed incident condition.

### 6.1 Origin of the Single-Order Incident Condition

There are not arbitrarily many independent inputs to consider. Physically, exactly one incident mode is prescribed: the zeroth diffraction order of the beam incident from above. In the outer region one therefore sets
$$a_0 = 1, \qquad a_n = 0 \ \text{ for all } n \neq 0,$$
i.e. the incident state is the standard basis vector $\mathbf{e}_0$ in mode space.

Multiplying the scattering matrix by this basis vector selects exactly the corresponding column of $S$ — this is why the reflection amplitudes can be read off directly from a single column of the final scattering matrix, rather than requiring a full matrix–vector solve for a general incident condition. With modes indexed symmetrically from $-N$ to $+N$, the order $n=0$ corresponds to the middle position in that indexing, which is why this particular column is singled out.

This is the missing physical link between "the scattering matrix has been assembled" and "the reflection amplitudes are read off from one specific column of it."

---

## 7. From Amplitude to Power: Diffraction Efficiency

### 7.1 Phase Referencing to a Common Origin

A solution of the form $e^{i\beta_n y}$ changes phase under a shift of the coordinate origin. Any amplitude extracted at an internal reference plane $y=a$ (rather than at $y=0$) is therefore not directly the physically meaningful amplitude referenced to the structure's global origin — it must be corrected by a phase factor.

Two effects combine here. First, re-referencing the *outgoing* order $n$ from plane $a$ back to $y=0$ contributes a factor $e^{-i\beta_n a}$. Second, the incident condition itself ($a_0 = 1$) was physically imposed at plane $a$, meaning the incident wave has already accumulated a phase $e^{i\beta_0 a}$ getting there from $y=0$ — this phase must also be removed to reference the whole reflection process consistently to a common origin. Combined, this yields a correction of the form
$$\exp\left[-i(\beta_n+\beta_0)a\right].$$

For the diffraction efficiency itself, this absolute phase is ultimately irrelevant, since only $|B_n|^2$ enters — but it matters if reflection amplitudes (not just efficiencies) are to be compared meaningfully across different reference planes or against other formalisms.

### 7.2 From Field Amplitude to Power

The quantity $|B_n|^2$ is, by itself, only a squared field amplitude — not yet a power. The physically transported power is proportional to the outward normal component of the time-averaged Poynting vector,
$$\mathbf{S} = \frac{1}{2}\operatorname{Re}\left(\mathbf{E}\times\mathbf{H}^*\right).$$
For a plane wave of order $n$ propagating in the outer medium, the power flowing normal to the structure carries a factor proportional to $\operatorname{Re}\beta_n$. Normalizing the outgoing power in order $n$ to the incident power (order $0$) gives the diffraction efficiency
$$\eta_n = |B_n|^2\,\frac{\operatorname{Re}\beta_n}{\operatorname{Re}\beta_0}.$$

### 7.3 Why Evanescent Orders Automatically Have Zero Efficiency

For an evanescent order, $\beta_n = i\gamma_n$ with real $\gamma_n>0$, so $\operatorname{Re}\beta_n = 0$ and therefore
$$\boxed{\eta_n = 0.}$$
This does not mean the evanescent mode is physically absent — it can be significant in the near field close to the structure — but it carries no net normal radiative power into the far field of the homogeneous outer half-space. This is why such modes must still be carried through the full internal calculation (they participate in mode coupling and boundary matching), yet contribute nothing to the far-field efficiency.

---

## 8. Truncation and Its Physical Meaning

### 8.1 Finite Mode Number $N$

The exact theory involves $n \in \mathbb{Z}$, an infinite set of Floquet orders. Any practical solution truncates this to $n=-N,\dots,+N$:
$$\boxed{\text{infinitely many Floquet orders} \ \longrightarrow\ 2N+1 \text{ numerical orders.}}$$
The underlying physics is unchanged by this truncation; only the numerically tracked mode space is made finite. This introduces a convergence question distinct from the physical derivation itself: whether $N$ is large enough that the computed efficiencies are insensitive to further increasing it.

### 8.2 Relation Between Evanescent Decay and Discretization Scale

An evanescent order carries factors $e^{\pm\gamma_n y}$. Across a layer of thickness $\Delta y$, this produces a magnitude factor $e^{|\gamma_n|\Delta y}$. If this factor becomes too large, very large and very small numbers coexist within the same numerical matrices ([Why Direct Multiplication of Transfer Matrices Is Numerically Problematic](#34-why-direct-multiplication-of-transfer-matrices-is-numerically-problematic)), degrading accuracy. Keeping $|\beta_{\max}|\Delta y$ bounded is therefore the physical criterion that ties the required layer resolution to the evanescent decay rates present in the mode spectrum — the near-field character of high orders is directly linked to how finely the vertical direction must be discretized.

---

## 9. The Coating Causal Chain

If a coating changes the local material distribution $\varepsilon_r(x,y)$, this changes $k_p^2(y)$, which changes the coupling matrix $M_{nm}(y)$, which changes the vertical mode profiles $u_n(y)$, which changes the layer transfer relations, which changes the assembled scattering matrix, which changes the reflection amplitudes $B_n$, and therefore the efficiencies $\eta_n$:
$$\boxed{\text{coating geometry} \rightarrow \varepsilon_r(x,y) \rightarrow k_n^2(y) \rightarrow M(y) \rightarrow u_n(y) \rightarrow S \rightarrow B_n \rightarrow \eta_n.}$$
This is the complete causal path by which any change to the layer geometry enters the diffraction efficiencies — there is no other route.

---

## 10. The Central Conceptual Reduction

The physical reduction underlying the entire method is:
$$\boxed{\text{2D electromagnetic problem} \ \Longrightarrow\ \text{1D system of coupled modes.}}$$
The $x$-direction is not eliminated by ignoring it, but by encoding its periodicity exactly into the Fourier/Floquet basis. What remains is a single continuous variable, $y$. This reduction is the conceptual core of the whole formalism.

---

## 11. Physical Consistency Checks

### 11.1 Homogeneous-Grating Limit

A strong test of the full formalism is the limit $\varepsilon_r(x,y) = \text{const}$. Then $k_n^2 = 0$ for all $n\neq0$, all coupling vanishes, and the equations decouple into
$$u_n'' + \beta_n^2 u_n = 0.$$
Only the actually excited incident order remains active, and the system reduces to ordinary plane-wave propagation in a homogeneous medium. A correctly formulated solver must reproduce exactly this trivial physics in this limit.

### 11.2 Energy Balance

For lossless materials, summing over all propagating reflected and transmitted orders should satisfy
$$\boxed{\sum_n R_n + \sum_n T_n = 1}$$
up to numerical error and normalization conventions. For absorbing materials,
$$\sum_n R_n + \sum_n T_n < 1,$$
since part of the energy is dissipated in the material. This energy check is one of the most important independent validation tests for any implementation of the method.

### 11.3 Complex Quantities and Absorption

If the material is absorbing, the refractive index $n = n' + in''$ is complex, and so is $k^2$. Consequently $k_n^2(y)$, $M(y)$, $\beta_n$, $u_n(y)$, and $B_n$ all become complex quantities. Their imaginary parts are not numerical artifacts — they physically encode phase shift and absorption.

---

## 12. Full Equation Chain (Summary)

$$\varepsilon_r(x,y)\ \rightarrow\ k^2(x,y)=k_0^2\varepsilon_r(x,y)\ \rightarrow\ \alpha_n = k_{\text{top}}\sin\theta + n\frac{2\pi}{d}$$
$$u(x,y) = \sum_n u_n(y)\,e^{i\alpha_n x}, \qquad k^2(x,y) = \sum_p k_p^2(y)\,e^{ipKx}$$
$$u_n''(y) = \sum_m\left[\alpha_n^2\delta_{nm} - k_{n-m}^2(y)\right]u_m(y), \qquad \beta_n^2 = k_{\text{medium}}^2 - \alpha_n^2$$
$$\mathbf{w}(y_+) = \mathcal{T}_{\text{layer}}\,\mathbf{w}(y_-), \qquad \mathbf{a}_{\text{out}} = S\,\mathbf{a}_{\text{in}}, \qquad \mathbf{a}_{\text{in}} = \mathbf{e}_0$$
$$B_n = (S\mathbf{e}_0)_n \quad \text{(up to the phase reference of §7.1)}, \qquad \eta_n = |B_n|^2\,\frac{\operatorname{Re}\beta_n}{\operatorname{Re}\beta_0}$$

---

## 13. Conclusion

Physically, this method is not "an algorithm for computing diffraction efficiencies" but the numerical realization of a well-posed electromagnetic boundary-value problem:
$$\boxed{\text{periodic Maxwell problem} \rightarrow \text{coupled Floquet modes} \rightarrow \text{vertical ODE} \rightarrow \text{global scattering matrix} \rightarrow \text{diffracted power.}}$$

Geometry enters the physics exclusively through $\varepsilon_r(x,y) \rightarrow k_n^2(y)$; periodicity generates $\alpha_n = \alpha_0 + nK$; and the coupling of these modes produces diffraction. The two conceptual bridges that tie the whole derivation together are
$$\boxed{\text{Maxwell} \rightarrow \text{scalar TE-Helmholtz} \rightarrow \text{Floquet} \rightarrow \text{coupled ODE}}$$
and
$$\boxed{\text{outer boundary condition} \rightarrow \text{scattering-matrix column} \rightarrow \text{Poynting flux}.}$$

---

## 14. Mathematical Profs

This section fills in those intermediate steps so that the derivation can be followed from one equation to the next without relying on unstated identities or assumptions. The purpose of this section is **not** to introduce additional physics. It is to make the mathematical chain underlying the implementation explicit.

### 14.1 From Maxwell's Equations to the Scalar TE Helmholtz Equation

We start from Maxwell's equations in a linear, isotropic, non-magnetic medium without free currents:

$$ \vec\nabla\times\vec E = -\frac{\partial\vec B}{\partial t}, \qquad \vec\nabla\times\vec H = \frac{\partial\vec D}{\partial t}, $$

with

$$ \vec D=\varepsilon_0\varepsilon_r\vec E, \qquad \vec B=\mu_0\vec H. $$

Assume time-harmonic fields of the form

$$ \vec E(\mathbf r,t) = \vec E(\mathbf r)e^{-i\omega t}, \qquad \vec H(\mathbf r,t) = \vec H(\mathbf r)e^{-i\omega t}. $$

Then

$$ \frac{\partial}{\partial t} \rightarrow -i\omega. $$

Therefore,

$$ \vec\nabla\times\vec E = i\omega\mu_0\vec H $$

and

$$ \vec\nabla\times\vec H = -i\omega\varepsilon_0\varepsilon_r\vec E. $$

Take the curl of Faraday's law:

$$ \vec\nabla\times(\vec\nabla\times\vec E) = i\omega\mu_0 (\vec\nabla\times\vec H). $$

Substituting Ampère's law gives

$$ \vec\nabla\times(\vec\nabla\times\vec E) = i\omega\mu_0 \left( -i\omega\varepsilon_0\varepsilon_r\vec E \right). $$

Since

$$ i(-i)=1, $$

we obtain

$$ \vec\nabla\times(\vec\nabla\times\vec E) = \omega^2\mu_0\varepsilon_0\varepsilon_r\vec E. $$

Define the vacuum wavenumber

$$ k_0 = \omega\sqrt{\mu_0\varepsilon_0} = \frac{\omega}{c} = \frac{2\pi}{\lambda}. $$

Hence

$$ \vec\nabla\times(\vec\nabla\times\vec E) = k_0^2\varepsilon_r\vec E. $$

Now use the vector identity

$$ \boxed{ \vec\nabla\times(\vec\nabla\times\vec E) = \vec\nabla(\vec\nabla\cdot\vec E) - \vec\nabla^2\vec E }. $$

Thus,

$$ \vec\nabla(\vec\nabla\cdot\vec E) - \vec\nabla^2\vec E = k_0^2\varepsilon_r\vec E. $$

#### TE specialization

For the 1D grating considered here, the structure is invariant in the $z$ direction:

$$ \frac{\partial}{\partial z}=0. $$

For TE polarization we choose

$$ \vec E = \begin{pmatrix} 0\\ 0\\ u(x,y) \end{pmatrix}. $$

Consequently,

$$ \vec\nabla\cdot\vec E = \frac{\partial E_z}{\partial z} = 0. $$

The vector equation therefore reduces to

$$ -\vec\nabla^2\vec E = k_0^2\varepsilon_r\vec E. $$

Only the $z$ component is non-zero, so

$$ - \left( \frac{\partial^2u}{\partial x^2} + \frac{\partial^2u}{\partial y^2} \right) = k_0^2\varepsilon_r u. $$

Rearranging,

$$ \boxed{ \frac{\partial^2u}{\partial x^2} + \frac{\partial^2u}{\partial y^2} + k_0^2\varepsilon_r(x,y)u = 0. } $$

This is the scalar Helmholtz equation used by the solver.



### 14.2 From Periodicity to the Floquet Expansion

The grating satisfies

$$ \varepsilon_r(x+d,y) = \varepsilon_r(x,y). $$

Define the grating wavevector

$$ K=\frac{2\pi}{d}. $$

Bloch-Floquet theory states that a solution of a periodic problem can be written as

$$ u(x,y) = e^{i\alpha_0x}p(x,y), $$

where $p$ has the same period as the structure:

$$ p(x+d,y)=p(x,y). $$

Because $p$ is periodic, it has the Fourier expansion

$$ p(x,y) = \sum_{n=-\infty}^{\infty} u_n(y)e^{inKx}. $$

Therefore,

$$ u(x,y) = \sum_n u_n(y)e^{i(\alpha_0+nK)x}. $$

Define

$$ \boxed{ \alpha_n=\alpha_0+nK. } $$

Then

$$ \boxed{ u(x,y) = \sum_nu_n(y)e^{i\alpha_nx}. } $$

For an incident plane wave in the upper medium,

$$ \alpha_0 = k_{\mathrm{top}}\sin\theta. $$

Consequently,

$$ \boxed{ \alpha_n = k_{\mathrm{top}}\sin\theta + n\frac{2\pi}{d}. } $$

The important point is that the infinite set of diffraction orders is a direct consequence of the periodicity. It is not introduced by the numerical method.



### 14.3 Fourier Expansion of the Material Function

Because the permittivity is periodic in $x$, the quantity

$$ k^2(x,y) = k_0^2\varepsilon_r(x,y) $$

is also periodic:

$$ k^2(x+d,y)=k^2(x,y). $$

It therefore has a Fourier expansion

$$ \boxed{ k^2(x,y) = \sum_{p=-\infty}^{\infty} k_p^2(y)e^{ipKx}. } $$

The Fourier coefficients are obtained from the standard orthogonality relation:

$$ \boxed{ k_p^2(y) = \frac{1}{d} \int_0^d k^2(x,y)e^{-ipKx}\,dx. } $$

The normalization follows from

$$ \frac{1}{d} \int_0^d e^{i(n-p)Kx}\,dx = \delta_{np}. $$

Indeed,

$$ \int_0^d e^{i(n-p)Kx}\,dx = \begin{cases} d,&n=p,\\ 0,&n\neq p. \end{cases} $$

Therefore Fourier projection uniquely extracts each coefficient.



### 14.4 Derivation of the Coupled-Mode Equation

Start from

$$ \left( \partial_x^2+\partial_y^2+k^2(x,y) \right)u(x,y)=0. $$

Insert

$$ u(x,y) = \sum_m u_m(y)e^{i\alpha_mx}. $$

#### $x$ derivative

Since

$$ \frac{\partial}{\partial x} e^{i\alpha_mx} = i\alpha_me^{i\alpha_mx}, $$

we have

$$ \frac{\partial^2}{\partial x^2} e^{i\alpha_mx} = -\alpha_m^2e^{i\alpha_mx}. $$

Therefore,

$$ \partial_x^2u = -\sum_m \alpha_m^2u_m(y)e^{i\alpha_mx}. $$

#### $y$ derivative

Because $\alpha_m$ does not depend on $y$,

$$ \partial_y^2u = \sum_m u_m''(y)e^{i\alpha_mx}. $$

#### Material product

Using

$$ k^2(x,y) = \sum_p k_p^2(y)e^{ipKx}, $$

we obtain

$$ k^2u = \left( \sum_p k_p^2e^{ipKx} \right) \left( \sum_m u_me^{i\alpha_mx} \right). $$

Multiplying the two sums,

$$ k^2u = \sum_p\sum_m k_p^2u_m e^{i(pK+\alpha_m)x}. $$

Now use

$$ \alpha_m=\alpha_0+mK. $$

Then

$$ pK+\alpha_m = pK+\alpha_0+mK = \alpha_0+(m+p)K = \alpha_{m+p}. $$

Introduce

$$ n=m+p. $$

Then

$$ p=n-m $$

and therefore

$$ k^2u = \sum_n \left( \sum_m k_{n-m}^2u_m \right) e^{i\alpha_nx}. $$

Putting all terms into the Helmholtz equation gives

$$ \sum_n \left[ u_n'' - \alpha_n^2u_n + \sum_m k_{n-m}^2u_m \right] e^{i\alpha_nx} = 0. $$

Because the functions

$$ e^{i\alpha_nx} $$

are linearly independent, every coefficient must vanish:

$$ u_n'' - \alpha_n^2u_n + \sum_mk_{n-m}^2u_m = 0. $$

Therefore,

$$ \boxed{ u_n'' = \alpha_n^2u_n - \sum_mk_{n-m}^2u_m. } $$

Introduce the Kronecker delta,

$$ \delta_{nm} = \begin{cases} 1,&n=m,\\ 0,&n\neq m. \end{cases} $$

Since

$$ \alpha_n^2u_n = \sum_m \alpha_n^2\delta_{nm}u_m, $$

we obtain

$$ u_n'' = \sum_m \left[ \alpha_n^2\delta_{nm} - k_{n-m}^2 \right]u_m. $$

Define

$$ \boxed{ M_{nm}(y) = \alpha_n^2\delta_{nm} - k_{n-m}^2(y). } $$

Then the complete coupled system becomes

$$ \boxed{ \mathbf u''(y)=M(y)\mathbf u(y). } $$



### 14.5 Why the Material Fourier Coefficients Couple Different Orders

The coupling term is

$$ k_{n-m}^2. $$

Suppose $n=m$. Then

$$ k_{n-m}^2=k_0^2, $$

which corresponds to the spatial average of the material.

If $n\neq m$, then

$$ k_{n-m}^2 $$

is a non-zero Fourier component of the periodic material distribution.

Thus a material Fourier component with index $p$ couples modes satisfying

$$ n-m=p. $$

Equivalently,

$$ \boxed{ n=m+p. } $$

This gives a direct mathematical interpretation of diffraction:

A spatial Fourier component of the grating transfers an integer multiple of the grating momentum $K$ between Floquet orders.

Indeed,

$$ \alpha_n-\alpha_m = (n-m)K = pK. $$

Therefore,

$$ \boxed{ \text{material Fourier order }p \quad\Longleftrightarrow\quad \text{momentum transfer }pK. } $$



### 14.6 Piecewise-Constant Fourier Coefficients

For a fixed $y$, suppose $k^2(x,y)$ is piecewise constant over one period.

Let the discontinuities occur at

$$ x_1,x_2,\ldots,x_P $$

and define the jump at each discontinuity as

$$ \sigma_p = k^2_{\mathrm{right}}-k^2_{\mathrm{left}}. $$

For $n\neq0$,

$$ k_n^2 = \frac{1}{d} \int_0^d k^2(x)e^{-inKx}\,dx. $$

Integration by parts gives

$$ \int_0^d k^2(x)e^{-inKx}\,dx = \left[ \frac{k^2(x)e^{-inKx}}{-inK} \right]_0^d + \frac{1}{inK} \int_0^d \frac{dk^2}{dx} e^{-inKx}\,dx. $$

Because $k^2(x)$ is periodic and

$$ e^{-inKd} = e^{-in2\pi} = 1, $$

the boundary term cancels.

The derivative of a piecewise-constant function consists of delta distributions at the jumps:

$$ \frac{dk^2}{dx} = \sum_p \sigma_p\delta(x-x_p). $$

Therefore,

$$ \int_0^d \frac{dk^2}{dx} e^{-inKx}\,dx = \sum_p \sigma_p e^{-inKx_p}. $$

Hence

$$ k_n^2 = \frac{1}{d} \frac{1}{inK} \sum_p \sigma_p e^{-inKx_p}. $$

Since

$$ K=\frac{2\pi}{d}, $$

we have

$$ dK=2\pi. $$

Thus

$$ \boxed{ k_n^2 = \frac{1}{2\pi in} \sum_p \sigma_p e^{-inKx_p}. } $$

The exact sign depends on the convention used for the jump orientation and Fourier-transform convention. With the convention

$$ \sigma_p=k^2_{\mathrm{right}}-k^2_{\mathrm{left}}, $$

and the Fourier factor $e^{-inKx}$, the expression above follows directly from the distributional derivative. If the implementation defines the jump with the opposite orientation, the equivalent formula appears with a minus sign.

For $n=0$, the formula involving $1/n$ is invalid. Instead,

$$ \boxed{ k_0^2 = \frac{1}{d} \int_0^d k^2(x)\,dx. } $$

Thus the zeroth coefficient is simply the spatial average.



### 14.7 Homogeneous Medium as a Limiting Case

Suppose the medium is homogeneous in $x$:

$$ k^2(x,y)=k^2. $$

Then

$$ k_n^2 = 0 \qquad (n\neq0), $$

while

$$ k_0^2=k^2. $$

Therefore,

$$ k_{n-m}^2 = k^2\delta_{nm}. $$

The coupling matrix becomes

$$ M_{nm} = \alpha_n^2\delta_{nm} - k^2\delta_{nm}. $$

Hence

$$ \boxed{ M_{nm} = (\alpha_n^2-k^2)\delta_{nm}. } $$

The matrix is diagonal, so all modes decouple:

$$ u_n'' = (\alpha_n^2-k^2)u_n. $$

Define

$$ \boxed{ \beta_n^2 = k^2-\alpha_n^2. } $$

Then

$$ \boxed{ u_n''+\beta_n^2u_n=0. } $$

The two independent solutions are

$$ u_n(y) = A_ne^{i\beta_ny} + B_ne^{-i\beta_ny}. $$

Differentiating,

$$ u_n'(y) = i\beta_nA_ne^{i\beta_ny} - i\beta_nB_ne^{-i\beta_ny}. $$

At $y=0$,

$$ u_n=A_n+B_n, $$

and

$$ u_n' = i\beta_n(A_n-B_n). $$

This gives the basic conversion between field variables and directional modal amplitudes.



### 14.8 Derivation of the Modal-Amplitude Transformation

At a reference plane,

$$ u_n=A_n+B_n $$

and

$$ u_n'=i\beta_n(A_n-B_n). $$

Divide the second equation by $i\beta_n$:

$$ \frac{u_n'}{i\beta_n} = A_n-B_n. $$

We therefore have the two equations

$$ u_n=A_n+B_n $$

and

$$ \frac{u_n'}{i\beta_n} = A_n-B_n. $$

Adding them,

$$ u_n+\frac{u_n'}{i\beta_n} = 2A_n. $$

Thus

$$ \boxed{ A_n = \frac12 \left( u_n+\frac{u_n'}{i\beta_n} \right). } $$

Subtracting,

$$ u_n-\frac{u_n'}{i\beta_n} = 2B_n, $$

so

$$ \boxed{ B_n = \frac12 \left( u_n-\frac{u_n'}{i\beta_n} \right). } $$

In matrix form,

$$ \begin{pmatrix} u_n\\ u_n' \end{pmatrix} = \begin{pmatrix} 1&1\\ i\beta_n&-i\beta_n \end{pmatrix} \begin{pmatrix} A_n\\ B_n \end{pmatrix}. $$

The inverse transformation is

$$ \begin{pmatrix} A_n\\ B_n \end{pmatrix} = \frac12 \begin{pmatrix} 1&1/(i\beta_n)\\ 1&-1/(i\beta_n) \end{pmatrix} \begin{pmatrix} u_n\\ u_n' \end{pmatrix}. $$

This is the precise mathematical basis for changing between the state-vector representation and the incoming/outgoing representation.



### 14.9 Propagating and Evanescent Orders

Recall

$$ \beta_n^2 = k^2-\alpha_n^2. $$

#### Propagating case

If

$$ |\alpha_n|<|k|, $$

then $\beta_n^2>0$ for a real positive-index lossless medium, so $\beta_n$ can be chosen real:

$$ \beta_n\in\mathbb R. $$

The solutions are oscillatory:

$$ e^{\pm i\beta_ny}. $$

These represent waves carrying energy in opposite vertical directions.

#### Evanescent case

If

$$ |\alpha_n|>|k|, $$

then

$$ \beta_n^2<0. $$

Write

$$ \beta_n=i\gamma_n, \qquad \gamma_n>0. $$

Then

$$ e^{i\beta_ny} = e^{i(i\gamma_n)y} = e^{-\gamma_ny}, $$

while

$$ e^{-i\beta_ny} = e^{+\gamma_ny}. $$

Thus one solution decays and the other grows.

For a semi-infinite region extending toward positive $y$, physical boundedness requires the growing solution to be absent. For a region extending toward negative $y$, the opposite branch is selected.

This is why the square-root branch cannot be chosen arbitrarily.

A common radiation-condition convention is

$$ \boxed{ \operatorname{Im}\beta_n\ge0. } $$


### 14.10 First-Order Form of the Coupled ODE

Starting from

$$ \mathbf u''(y)=M(y)\mathbf u(y), $$

define

$$ \mathbf v(y)=\mathbf u'(y). $$

Then

$$ \mathbf u'=\mathbf v $$

and

$$ \mathbf v'=M(y)\mathbf u. $$

Define the state vector

$$ \mathbf w = \begin{pmatrix} \mathbf u\\ \mathbf v \end{pmatrix}. $$

Then

$$ \frac{d\mathbf w}{dy} = \begin{pmatrix} 0&I\\ M(y)&0 \end{pmatrix} \mathbf w. $$

Therefore,

$$ \boxed{ \mathbf w'(y)=A(y)\mathbf w(y), } $$

with

$$ \boxed{ A(y) = \begin{pmatrix} 0&I\\ M(y)&0 \end{pmatrix}. } $$

If there are

$$ Q=2N+1 $$

retained diffraction orders, then

$$ \mathbf u\in\mathbb C^Q $$

and

$$ \mathbf w\in\mathbb C^{2Q}. $$

Thus

$$ \boxed{ \dim(\mathbf w)=2(2N+1). } $$

This is the origin of the state-space dimension used by the numerical solver.



### 14.11 Fundamental Matrix and Transfer Matrix

Consider a layer occupying

$$ y_-\le y\le y_+. $$

The first-order system is

$$ \mathbf w'(y)=A(y)\mathbf w(y). $$

A fundamental matrix $\Phi(y)$ is a matrix whose columns are linearly independent solutions:

$$ \Phi(y) = \begin{pmatrix} |&|&&|\\ \mathbf w^{(1)}(y)& \mathbf w^{(2)}(y)&\cdots& \mathbf w^{(2Q)}(y)\\ |&|&&| \end{pmatrix}. $$

Every solution can be written as

$$ \mathbf w(y)=\Phi(y)\mathbf c $$

for some constant coefficient vector $\mathbf c$.

At the lower boundary,

$$ \mathbf w(y_-) = \Phi(y_-)\mathbf c. $$

Therefore,

$$ \mathbf c = \Phi(y_-)^{-1}\mathbf w(y_-). $$

At the upper boundary,

$$ \mathbf w(y_+) = \Phi(y_+)\mathbf c. $$

Substituting the previous expression,

$$ \mathbf w(y_+) = \Phi(y_+) \Phi(y_-)^{-1} \mathbf w(y_-). $$

Hence,

$$ \boxed{ \mathcal T = \Phi(y_+)\Phi(y_-)^{-1} } $$

and

$$ \boxed{ \mathbf w(y_+) = \mathcal T\mathbf w(y_-). } $$

This proves why integrating a complete basis of trial solutions produces the transfer matrix.



### 14.12 Why the Transfer Matrix Has Dimension $2Q\times2Q$

There are $Q=2N+1$ retained Floquet amplitudes.

Each satisfies a second-order differential equation.

Therefore, each order contributes two independent initial conditions:

$$ u_n(y_-), \qquad u_n'(y_-). $$

Hence the total number of independent initial conditions is

$$ 2Q=2(2N+1). $$

The state vector therefore belongs to

$$ \mathbb C^{2Q}. $$

A linear mapping from this space to itself must be represented by a

$$ \boxed{ 2Q\times2Q } $$

matrix.

Thus the dimension of the transfer matrix is not an arbitrary implementation choice.



### 14.13 From Transfer Variables to Scattering Variables

A transfer matrix relates quantities on opposite sides:

$$ \mathbf w_+ = \mathcal T\mathbf w_-. $$

This representation treats both directions symmetrically as state variables.

For scattering problems, however, it is more natural to separate incoming and outgoing waves.

Schematically, write

$$ \begin{pmatrix} \mathbf a_{\mathrm{out}}^{\mathrm{top}}\\ \mathbf a_{\mathrm{out}}^{\mathrm{bot}} \end{pmatrix} = S \begin{pmatrix} \mathbf a_{\mathrm{in}}^{\mathrm{top}}\\ \mathbf a_{\mathrm{in}}^{\mathrm{bot}} \end{pmatrix}. $$

The scattering matrix therefore has the block structure

$$ S= \begin{pmatrix} S_{11}&S_{12}\\ S_{21}&S_{22} \end{pmatrix}. $$

The first index identifies the output side and the second index identifies the input side.

Thus:

- $S_{11}$: top $\rightarrow$ top,
- $S_{12}$: bottom $\rightarrow$ top,
- $S_{21}$: top $\rightarrow$ bottom,
- $S_{22}$: bottom $\rightarrow$ bottom.

The precise naming of $a$ and $b$ variables can differ between implementations, but the physical meaning is always the same: the S-matrix maps incoming channels to outgoing channels.



### 14.14 Derivation of the Redheffer Star Product

Suppose two structures $A$ and $B$ are connected.

Write the scattering equations of $A$ as

$$ \begin{pmatrix} b_1\\ b_2 \end{pmatrix} = \begin{pmatrix} S_{11}^A&S_{12}^A\\ S_{21}^A&S_{22}^A \end{pmatrix} \begin{pmatrix} a_1\\ a_2 \end{pmatrix}. $$

Thus,

$$ b_1 = S_{11}^Aa_1 + S_{12}^Aa_2 $$

and

$$ b_2 = S_{21}^Aa_1 + S_{22}^Aa_2. $$

For structure $B$,

$$ \begin{pmatrix} c_1\\ c_2 \end{pmatrix} = \begin{pmatrix} S_{11}^B&S_{12}^B\\ S_{21}^B&S_{22}^B \end{pmatrix} \begin{pmatrix} d_1\\ d_2 \end{pmatrix}. $$

At the common interface, the internal waves satisfy

$$ a_2=c_1 $$

and

$$ d_1=b_2. $$

Therefore,

$$ a_2 = S_{11}^Bc_? $$

and, after applying the interface identifications consistently,

$$ a_2 = S_{11}^Bb_2 + S_{12}^Bd_2. $$

Using

$$ b_2 = S_{21}^Aa_1 + S_{22}^Aa_2, $$

we obtain

$$ a_2 = S_{11}^B \left( S_{21}^Aa_1 + S_{22}^Aa_2 \right) + S_{12}^Bd_2. $$

Collect the unknown internal amplitude $a_2$:

$$ \left( I-S_{11}^BS_{22}^A \right)a_2 = S_{11}^BS_{21}^Aa_1 + S_{12}^Bd_2. $$

Therefore,

$$ \boxed{ a_2 = \left( I-S_{11}^BS_{22}^A \right)^{-1} \left( S_{11}^BS_{21}^Aa_1 + S_{12}^Bd_2 \right). } $$

The internal amplitude has now been eliminated.

Substituting this result into the expressions for the external outgoing waves yields the Redheffer star-product formulas:

$$ \boxed{ S_{11} = S_{11}^A + S_{12}^A \left( I-S_{11}^BS_{22}^A \right)^{-1} S_{11}^BS_{21}^A } $$

$$ \boxed{ S_{12} = S_{12}^A \left( I-S_{11}^BS_{22}^A \right)^{-1} S_{12}^B } $$

and equivalently,

$$ \boxed{ S_{21} = S_{21}^B \left( I-S_{22}^AS_{11}^B \right)^{-1} S_{21}^A } $$

$$ \boxed{ S_{22} = S_{22}^B + S_{21}^B \left( I-S_{22}^AS_{11}^B \right)^{-1} S_{22}^AS_{12}^B. } $$

The important mathematical operation is the elimination of the internal interface variables.



### 14.15 Why Scattering-Matrix Composition Is Numerically Better

Consider an evanescent mode with

$$ \beta=i\gamma. $$

Propagation through a thickness $h$ produces factors

$$ e^{-\gamma h} $$

and

$$ e^{+\gamma h}. $$

For sufficiently large $\gamma h$,

$$ e^{+\gamma h}\gg1 $$

while

$$ e^{-\gamma h}\ll1. $$

For example, if

$$ \gamma h=100, $$

then

$$ e^{100}\approx2.69\times10^{43}, $$

while

$$ e^{-100}\approx3.72\times10^{-44}. $$

A transfer matrix must represent both scales simultaneously.

After many layers, the numerical dynamic range can become enormous even when the physical reflection and transmission coefficients remain of order unity.

The scattering formulation instead eliminates internal growing amplitudes during the composition process. It therefore avoids carrying the exponentially growing solutions as global external unknowns.

This is the fundamental numerical motivation for the S-matrix formulation.

### 14.16 Incident-Order Selection

The incident field is a single plane wave.

Therefore the incoming amplitude vector has the form

$$ \mathbf a_{\mathrm{in}} = \vec E_0, $$

where $\vec E_0$ is the unit vector corresponding to diffraction order $n=0$.

For orders indexed by

$$ -N,\ldots,-1,0,1,\ldots,N, $$

the position of order zero is

$$ N+1 $$

in one-based indexing.

Thus,

$$ \vec E_0 = \begin{pmatrix} 0\\ \vdots\\ 0\\ 1\\ 0\\ \vdots\\ 0 \end{pmatrix}. $$

If

$$ \mathbf b_{\mathrm{out}} = S\mathbf a_{\mathrm{in}}, $$

then

$$ \mathbf b_{\mathrm{out}} = S\vec E_0. $$

Multiplication by a basis vector selects a matrix column:

$$ \boxed{ S\vec E_0 = S_{(:,\,0)}. } $$

Therefore the reflected diffraction amplitudes are obtained directly from the corresponding column of the final scattering matrix.

No general matrix-vector solve is required for the single incident-order problem.



### 14.17 Derivation of the Diffraction-Efficiency Formula

Consider an outgoing plane-wave order

$$ E_z^{(n)} = B_n e^{i\alpha_nx+i\beta_ny}. $$

For TE polarization, Maxwell's equation gives the corresponding magnetic field.

From

$$ \vec\nabla\times\vec E = i\omega\mu_0\vec H, $$

we obtain

$$ \vec H = \frac{1}{i\omega\mu_0} \vec\nabla\times\vec E. $$

For

$$ \vec E= \hat{\mathbf z} B_ne^{i\alpha_nx+i\beta_ny}, $$

the curl contains the factors $\alpha_n$ and $\beta_n$.

The normal Poynting flux is proportional to

$$ \operatorname{Re} \left( \beta_n \right) |B_n|^2. $$

Consequently, after normalization to the incident order,

$$ \boxed{ \eta_n = |B_n|^2 \frac{\operatorname{Re}\beta_n} {\operatorname{Re}\beta_0}. } $$

For a lossless propagating order, $\beta_n$ is real, so this reduces to

$$ \eta_n = |B_n|^2 \frac{\beta_n}{\beta_0}. $$

The additional factor $\beta_n/\beta_0$ is essential: the squared field amplitude alone is not a power ratio.



### 14.18 Why Evanescent Orders Carry Zero Far-Field Power

For an evanescent mode,

$$ \beta_n=i\gamma_n, \qquad \gamma_n>0. $$

Therefore,

$$ \operatorname{Re}\beta_n=0. $$

The efficiency expression gives

$$ \eta_n = |B_n|^2 \frac{0}{\operatorname{Re}\beta_0} = 0. $$

Thus

$$ \boxed{ \eta_n=0 \qquad \text{for an evanescent far-field order}. } $$

This does **not** mean that $B_n=0$.

An evanescent mode can have a substantial field amplitude near the grating. It simply does not transport net propagating power into the far field.

This distinction is important for understanding why evanescent orders must still be retained in the internal coupled-mode calculation.



### 14.19 Energy Conservation

For a lossless structure, the time-averaged Poynting theorem implies that the total incoming power equals the total outgoing power.

For one incident order,

$$ P_{\mathrm{inc}} = P_{\mathrm{refl}} + P_{\mathrm{trans}}. $$

After normalization by the incident power,

$$ \boxed{ \sum_nR_n+\sum_nT_n=1. } $$

Only propagating orders contribute to these sums.

Numerically, one obtains

$$ \sum_nR_n+\sum_nT_n = 1+\epsilon_{\mathrm{num}}, $$

where $\epsilon_{\mathrm{num}}$ measures numerical error.

Thus energy conservation provides an independent validation of:

- Fourier coefficients,
- modal propagation constants,
- boundary conditions,
- layer integration,
- S-matrix composition,
- amplitude normalization,
- efficiency normalization.

If the materials are absorbing, the outgoing power is smaller:

$$ \boxed{ \sum_nR_n+\sum_nT_n<1. } $$

The missing power corresponds to absorption inside the material.



### 14.20 Complex Refractive Index and Complex Propagation Constants

For an absorbing material,

$$ n=n'+in''. $$

The relative permittivity is related to the refractive index by

$$ \varepsilon_r=n^2. $$

Therefore,

$$ \varepsilon_r = (n'+in'')^2 $$

and hence

$$ \varepsilon_r = (n'^2-n''^2) + 2in'n''. $$

Consequently,

$$ k^2=k_0^2\varepsilon_r $$

is complex.

Its Fourier coefficients are therefore complex:

$$ k_n^2\in\mathbb C. $$

The coupling matrix

$$ M_{nm} = \alpha_n^2\delta_{nm} - k_{n-m}^2 $$

is consequently complex as well.

The vertical propagation constants

$$ \beta_n = \sqrt{k^2-\alpha_n^2} $$

also become complex.

Therefore the field contains both phase accumulation and amplitude attenuation.

The imaginary parts are thus not numerical artifacts. They are part of the physical description of absorption.



### 14.21 Phase Referencing

Suppose an outgoing wave is represented at a reference plane $y=a$:

$$ u_n(a) \propto B_ne^{i\beta_na}. $$

If we want the amplitude referenced to $y=0$, we remove the propagation factor:

$$ B_n^{(0)} = B_n^{(a)}e^{-i\beta_na}. $$

The incident wave has also accumulated a phase between $y=0$ and $y=a$:

$$ A_0(a) = A_0(0)e^{i\beta_0a}. $$

If the incident amplitude is normalized to unity at $y=a$, then the corresponding amplitude at $y=0$ contains the inverse phase factor.

Combining the two reference changes gives a relative phase factor of the form

$$ \boxed{ e^{-i(\beta_n+\beta_0)a}. } $$

The exact sign depends on the chosen propagation-direction convention, but the physical principle is invariant:

> Amplitudes extracted at different reference planes differ by known propagation phases.

Since efficiency depends on

$$ |B_n|^2, $$

a purely real propagation phase does not change the efficiency.



### 14.22 Vertical Layering as a Numerical Approximation

The continuous structure has

$$ M=M(y). $$

The exact problem is

$$ \mathbf w'(y)=A(y)\mathbf w(y). $$

Suppose the vertical coordinate is divided into intervals

$$ [y_j,y_{j+1}] $$

with

$$ \Delta y_j=y_{j+1}-y_j. $$

Within each interval, the material distribution is approximated by a representative value:

$$ M(y)\approx M_j. $$

Then

$$ \mathbf w' = A_j\mathbf w. $$

The exact local solution for a constant $A_j$ is

$$ \mathbf w(y_{j+1}) = e^{A_j\Delta y_j} \mathbf w(y_j). $$

Thus the layer transfer matrix is, formally,

$$ \boxed{ T_j = e^{A_j\Delta y_j}. } $$

For a sufficiently fine discretization,

$$ M(y)\approx M_j $$

within each layer, and the piecewise-constant approximation converges toward the continuous problem as

$$ {\max}_{j} \Delta y_j\rightarrow0. $$

Hence layering is a numerical discretization of the continuous $y$ dependence.



### 14.23 Why the Evanescent Spectrum Controls the Vertical Resolution

For an evanescent mode,

$$ \beta_n=i\gamma_n. $$

Its characteristic variation length is

$$ L_n=\frac{1}{\gamma_n}. $$

A vertical discretization interval $\Delta y$ should therefore resolve the fastest relevant variation.

The propagation factor over one layer is

$$ e^{\gamma_n\Delta y}. $$

The largest retained decay/growth rate,

$$ \gamma_{\max}, $$

sets the most restrictive scale.

A useful dimensionless parameter is

$$ \boxed{ \gamma_{\max}\Delta y. } $$

If this quantity becomes too large, the layer contains an enormous exponential dynamic range.

Therefore, a physically meaningful discretization criterion is to keep

$$ \boxed{ |\beta_{\max}|\Delta y } $$

within a numerically manageable range.

This does not mean that there is one universal value of the allowed product. The appropriate threshold depends on the ODE integrator, floating-point precision, conditioning, mode truncation, and implementation details.



### 14.24 Mode Truncation

The exact Floquet expansion is

$$ u(x,y) = \sum_{n=-\infty}^{\infty} u_n(y)e^{i\alpha_nx}. $$

Numerically, this is replaced by

$$ \boxed{ u(x,y) \approx \sum_{n=-N}^{N} u_n(y)e^{i\alpha_nx}. } $$

The number of retained orders is

$$ Q=2N+1. $$

Consequently,

$$ M\in\mathbb C^{Q\times Q} $$

and the first-order state matrix has dimension

$$ 2Q\times2Q. $$

Increasing $N$ enlarges the represented modal space.

A converged calculation should satisfy

$$ \eta_n(N+\Delta N) \approx \eta_n(N) $$

for all quantities of interest.

Thus two distinct numerical convergence questions must be separated:

1. **Vertical discretization convergence**
   $$ \Delta y\rightarrow0; $$

2. **Fourier/modal truncation convergence**
   $$ N\rightarrow\infty. $$

A result can be converged with respect to one while still being unconverged with respect to the other.



### 14.25 The Complete Mathematical Chain

The entire derivation can now be written without skipping the principal mathematical transitions:

$$ \text{Maxwell} $$

$$ \Downarrow $$

$$ \vec\nabla\times(\vec\nabla\times\vec E) = k_0^2\varepsilon_r\vec E $$

$$ \Downarrow $$

$$ \text{TE specialization} $$

$$ \Downarrow $$

$$ \boxed{ (\partial_x^2+\partial_y^2+k_0^2\varepsilon_r)u=0 } $$

$$ \Downarrow $$

$$ \varepsilon_r(x+d,y)=\varepsilon_r(x,y) $$

$$ \Downarrow $$

$$ \boxed{ u(x,y) = \sum_nu_n(y)e^{i\alpha_nx}, \qquad \alpha_n=\alpha_0+nK } $$

$$ \Downarrow $$

$$ \boxed{ k^2(x,y) = \sum_pk_p^2(y)e^{ipKx} } $$

$$ \Downarrow $$

$$ \boxed{ u_n'' = \sum_m \left[ \alpha_n^2\delta_{nm} - k_{n-m}^2 \right]u_m } $$

$$ \Downarrow $$

$$ \boxed{ \mathbf u''=M(y)\mathbf u } $$

$$ \Downarrow $$

$$ \boxed{ \mathbf w' = \begin{pmatrix} 0&I\\ M&0 \end{pmatrix} \mathbf w } $$

$$ \Downarrow $$

$$ \boxed{ \mathbf w(y_+) = T_{\mathrm{layer}}\mathbf w(y_-) } $$

$$ \Downarrow $$

$$ \boxed{ S_{\mathrm{total}} = S_1\star S_2\star\cdots\star S_M } $$

$$ \Downarrow $$

$$ \boxed{ \mathbf a_{\mathrm{in}}=\vec E_0 } $$

$$ \Downarrow $$

$$ \boxed{ \vec B = S_{\mathrm{total}}\vec E_0 } $$

$$ \Downarrow $$

$$ \boxed{ \eta_n = |B_n|^2 \frac{\operatorname{Re}\beta_n} {\operatorname{Re}\beta_0} } $$

This is the mathematical chain implemented by the solver.



### 14.26 The Most Important Mathematical Dependencies

The derivation can also be viewed as a dependency graph:

$$ \boxed{ \varepsilon_r(x,y) } $$

determines

$$ \boxed{ k^2(x,y)=k_0^2\varepsilon_r(x,y) } $$

which determines its Fourier coefficients

$$ \boxed{ k_n^2(y) } $$

which determine the coupling matrix

$$ \boxed{ M_{nm}(y) = \alpha_n^2\delta_{nm} - k_{n-m}^2(y) } $$

which determines the coupled ODE

$$ \boxed{ \mathbf u''=M\mathbf u } $$

which determines the layer transfer relation

$$ \boxed{ \mathbf w_+ = T_{\mathrm{layer}}\mathbf w_- } $$

which is converted into and combined through scattering matrices

$$ \boxed{ S_{\mathrm{total}} } $$

which, for the incident state $\vec E_0$, gives

$$ \boxed{ \vec B=S_{\mathrm{total}}\vec E_0 } $$

and finally

$$ \boxed{ \eta_n = |B_n|^2 \frac{\operatorname{Re}\beta_n} {\operatorname{Re}\beta_0}. } $$

Thus every numerical quantity in the final diffraction efficiency can be traced back through an explicit mathematical chain to the original spatial permittivity distribution.



### 14.27 Final Mathematical Interpretation

The central mathematical structure of `TESolver` is therefore a sequence of representations of the same electromagnetic boundary-value problem:

$$ \boxed{ \text{Maxwell equations} } $$

are reduced to

$$ \boxed{ \text{scalar PDE} } $$

by the TE symmetry assumption.

The periodicity of the PDE is then represented by

$$ \boxed{ \text{Floquet/Fourier modes} } $$

which transforms the two-dimensional PDE into

$$ \boxed{ \text{a one-dimensional coupled ODE system}. } $$

The ODE system is represented locally by

$$ \boxed{ \text{transfer matrices} } $$

while numerical stability motivates the global representation

$$ \boxed{ \text{scattering matrices}. } $$

Finally, the physical observable is not the complex field amplitude itself, but the normal Poynting flux:

$$ \boxed{ \text{modal amplitude} \rightarrow \text{power} \rightarrow \text{diffraction efficiency}. } $$

The entire method can therefore be summarized mathematically as

$$ \boxed{ \varepsilon_r(x,y) \rightarrow k_n^2(y) \rightarrow M(y) \rightarrow \mathbf u(y) \rightarrow T \rightarrow S \rightarrow B_n \rightarrow \eta_n. } $$

This chain is the mathematical backbone of the implementation. The geometry determines the material Fourier coefficients; the Fourier coefficients determine modal coupling; the coupled modes determine the electromagnetic response; and the resulting outgoing modal amplitudes determine the measurable diffraction efficiencies.

---

# TMSolver

## Purpose and Scope

This section is the TM-polarization counterpart to the TE - it follows the same purpose: explaining the physical and mathematical basis for computing diffraction efficiencies for a single (wavelength, incidence angle) pair — not the concrete implementation, not geometry extraction, not MPI orchestration.

**This document does not repeat derivations that are identical for TE and TM.** The overall solution strategy — Floquet expansion, transfer-matrix propagation through layers, scattering-matrix recursion, mode truncation, layering as a numerical discretization — is polarization-independent machinery. Only the pieces that actually change under $E_z \to H_z$ are re-derived here. Everywhere else, it points to the corresponding section of TE above.

The physical reason a separate derivation is needed at all — rather than simply "the same equation with $H$ instead of $E$" — is that the TE equation has a *constant* coefficient multiplying $u$ ($k_0^2\varepsilon_r$, entering algebraically), whereas the TM equation has $\varepsilon_r$ appearing *inside a derivative* (as $1/\varepsilon_r$, dividing a gradient). This structural difference is the source of every TM-specific result in this document — the coupled-mode equation, the layer boundary conditions, the local first-order system, and the diffraction-efficiency normalization all trace back to it.

---

## Table of Contents

0. [What Carries Over Unchanged from TESolver](#0-what-carries-over-unchanged-from-tesolver)
1. [Physical Setup: TM Wave Equation](#1-physical-setup-tm-wave-equation)
    1. [Maxwell to Scalar TM Wave Equation](#11-maxwell-to-scalar-tm-wave-equation)
    2. [Coupled-Mode Equation for TM](#12-coupled-mode-equation-for-tm)
    3. [The Fourier-Factorization Subtlety](#13-the-fourier-factorization-subtlety)
    4. [Physical Interpretation of the TM Coupling](#14-physical-interpretation-of-the-tm-coupling)
2. [Local Layer Response](#2-local-layer-response)
    1. [First-Order System](#21-first-order-system)
    2. [Homogeneous-Medium Limiting Case](#22-homogeneous-medium-limiting-case)
    3. [Propagating and Evanescent Orders](#23-propagating-and-evanescent-orders)
3. [Boundary Conditions and Modal Amplitudes](#3-boundary-conditions-and-modal-amplitudes)
    1. [The Actual Boundary Conditions for TM](#31-the-actual-boundary-conditions-for-tm)
    2. [Modal Amplitudes from $u$ and $w$](#32-modal-amplitudes-from-u-and-w)
4. [Transfer Matrix and Scattering Matrix](#4-transfer-matrix-and-scattering-matrix)
5. [From Amplitude to Power: TM Diffraction Efficiency](#5-from-amplitude-to-power-tm-diffraction-efficiency)
    1. [The Tangential Electric Field for TM](#51-the-tangential-electric-field-for-tm)
    2. [Efficiency Formula and Its Difference from TE](#52-efficiency-formula-and-its-difference-from-te)
6. [Physical Consistency Checks](#6-physical-consistency-checks)
7. [Mathematical Proofs](#7-mathematical-proofs)
    1. [From Maxwell's Equations to the Scalar TM Wave Equation](#71-from-maxwells-equations-to-the-scalar-tm-wave-equation)
    2. [Derivation of the TM Coupled-Mode Equation](#72-derivation-of-the-tm-coupled-mode-equation)
    3. [Why Naive Termwise Fourier Inversion Is Inconsistent](#73-why-naive-termwise-fourier-inversion-is-inconsistent)
    4. [Reduction to the TE Result in a Homogeneous Layer](#74-reduction-to-the-te-result-in-a-homogeneous-layer)
    5. [Derivation of the TM Boundary-Matching Variable](#75-derivation-of-the-tm-boundary-matching-variable)
    6. [Derivation of the TM Modal-Amplitude Transformation](#76-derivation-of-the-tm-modal-amplitude-transformation)
    7. [Derivation of the TM Diffraction-Efficiency Formula](#77-derivation-of-the-tm-diffraction-efficiency-formula)
8. [Conclusion](#8-conclusion)

---

## 0. What Carries Over Unchanged from TESolver

The following results depend only on periodicity, linearity, and general matrix/ODE structure — not on which field component is scalar. They apply to TM exactly as stated in `TESolver_Algorithm_Reference.md`, with no re-derivation needed:

| Topic | TESolver section |
|---|---|
| Bloch–Floquet origin of $\alpha_n = k_{\text{top}}\sin\theta + nK$ | [Fourier & Floquet Expansion](#02-fourier-and-floquet-expansion), [From Periodicity to the Floquet Expansion](#142-from-periodicity-to-the-floquet-expansion) |
| Geometry bridge (extracting the $x$-cross-section at fixed $y$) | [Geometry Bridge](#05-geometry-bridge) |
| Dimension-counting argument for the state space | [Dimension of the Solution Space](#12-dimension-of-the-solution-space), [Why the Transfer Matrix Has Dimension $2Q\times2Q$](#1412-why-the-transfer-matrix-has-dimension-2qtimes2q) |
| Fundamental-matrix argument: why integrating trial solutions yields $\mathcal T$ | [Why a Basis of Trial Solutions Generates the Transfer Matrix](#32-why-a-basis-of-trial-solutions-generates-the-transfer-matrix), [Fundamental Matrix and Transfer Matrix](#1411-fundamental-matrix-and-transfer-matrix) |
| Why direct transfer-matrix multiplication is numerically unstable | [Why Direct Multiplication of Transfer Matrices Is Numerically Problematic](#34-why-direct-multiplication-of-transfer-matrices-is-numerically-problematic), [Why Scattering-Matrix Composition Is Numerically Better](#1415-why-scattering-matrix-composition-is-numerically-better) |
| Scattering-matrix idea and the Redheffer star product | [Scattering Matrix](#4-scattering-matrix), [From Transfer Variables to Scattering Variables](#1413-from-transfer-variables-to-scattering-variables)–[Derivation of the Redheffer Star Product](#1414-derivation-of-the-redheffer-star-product) |
| Layering as a numerical discretization of continuous $y$-dependence | [Layering as a Numerical Discretization](#5-layering-as-a-numerical-discretization), [Vertical Layering as a Numerical Approximation](#1422-vertical-layering-as-a-numerical-approximation)–[Why the Evanescent Spectrum Controls the Vertical Resolution](#1423-why-the-evanescent-spectrum-controls-the-vertical-resolution) |
| Origin of the single-order incident condition $\mathbf a_{\text{in}} = \mathbf e_0$ | [Origin of the Single-Order Incident Condition](#61-origin-of-the-single-order-incident-condition), [Incident-Order Selection](#1416-incident-order-selection) |
| Phase referencing to a common origin | [Phase Referencing to a Common Origin](#71-phase-referencing-to-a-common-origin), [Phase Referencing](#1421-phase-referencing) |
| Mode truncation ($N\to\infty$ vs. finite $N$) | [Finite Mode Number N](#81-finite-mode-number-n), [Mode Truncation](#1424-mode-truncation) |
| Energy-balance sum rule $\sum_n R_n + \sum_n T_n \le 1$ | [Energy Balance](#112-energy-balance), [Energy Conservation](#1419-energy-conservation) |
| Complex refractive index and absorption entering as complex coefficients | [Complex Quantities and Absorption](#113-complex-quantities-and-absorption), [Complex Refractive Index and Complex Propagation Constants](#1420-complex-refractive-index-and-complex-propagation-constants) |

What changes for TM is confined to: the wave equation itself, the coupling matrix $M(y)$, the first-order state system, the physical quantity that must match at layer boundaries, the modal-amplitude transform, and the diffraction-efficiency normalization. These are derived below.

---

## 1. Physical Setup: TM Wave Equation

### 1.1 Maxwell to Scalar TM Wave Equation

For TM polarization in the same 1D-periodic, $z$-invariant geometry as TE, the roles of $E$ and $H$ swap: the magnetic field is purely out-of-plane, $u(x,y) \equiv H_z(x,y)$, while $E$ lies in the $x$–$y$ plane. Carrying out the same elimination as in the TE case (full derivation in [§7.1](#71-from-maxwells-equations-to-the-scalar-tm-wave-equation)) gives, instead of the constant-coefficient Helmholtz equation, a wave equation with $\varepsilon_r$ *inside* the differential operator:

$$\boxed{\frac{\partial}{\partial x}\!\left(\frac{1}{\varepsilon_r(x,y)}\frac{\partial u}{\partial x}\right) + \frac{\partial}{\partial y}\!\left(\frac{1}{\varepsilon_r(x,y)}\frac{\partial u}{\partial y}\right) + k_0^2\, u(x,y) = 0.}$$

Contrast with the TE equation ([TESolver §0.1](#)):
$$\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + k_0^2\varepsilon_r(x,y)\, u = 0.$$

In TE, $\varepsilon_r$ multiplies $u$ directly — an algebraic, pointwise coefficient. In TM, $1/\varepsilon_r$ sits *inside* a divergence, weighting a gradient before it is differentiated again. This is not a cosmetic difference: it means the Fourier-space representation of the TM equation cannot be obtained by simply replacing $k_0^2\varepsilon_r \to k_0^2/\varepsilon_r$ termwise — the divergence structure has to be expanded explicitly, which is done next.

### 1.2 Coupled-Mode Equation for TM

Introduce the auxiliary field
$$\eta(x,y) \equiv \frac{1}{\varepsilon_r(x,y)},$$
periodic in $x$ with the same period $d$ as $\varepsilon_r$, and expand it in the same Fourier basis used for $k^2(x,y)$ in TE ([TESolver §0.6](#)):
$$\eta(x,y) = \sum_p \eta_p(y)\, e^{ipKx}.$$

Carrying out the Floquet substitution $u(x,y) = \sum_n u_n(y) e^{i\alpha_n x}$ in the divergence-form equation (full derivation in [§7.2](#72-derivation-of-the-tm-coupled-mode-equation)) yields a coupled system that, unlike TE, mixes the mode amplitudes **and** their derivatives through the same convolution structure:

$$\frac{d}{dy}\left[\sum_m H_{nm}(y)\, u_m'(y)\right] - \sum_m \alpha_n\alpha_m H_{nm}(y)\, u_m(y) + k_0^2\, u_n(y) = 0,$$

where $H_{nm}(y) \equiv \eta_{n-m}(y)$ is the Toeplitz (convolution) matrix built from the Fourier coefficients of $1/\varepsilon_r$. In matrix form, defining $D_\alpha = \mathrm{diag}(\alpha_n)$:

$$\boxed{\big(H(y)\,\mathbf u'(y)\big)' = \big(D_\alpha H(y) D_\alpha - k_0^2 I\big)\,\mathbf u(y).}$$

This is the TM analogue of the TE coupled-mode equation $\mathbf u''(y) = M(y)\mathbf u(y)$ ([TESolver §0.3](#)) — but note the derivative now acts on $H(y)\mathbf u'(y)$ as a whole, not on $\mathbf u'(y)$ alone. The physical quantity $H(y)\mathbf u'(y)$ is not a bookkeeping convenience; it is proportional to the tangential electric field ([§5.1](#51-the-tangential-electric-field-for-tm)), which is why it appears as a single differentiated unit.

### 1.3 The Fourier-Factorization Subtlety

A point worth flagging explicitly, because it is a genuine mathematical subtlety specific to TM (and absent from TE): $H(y)$ above is the Toeplitz matrix built from the Fourier coefficients of $\eta = 1/\varepsilon_r$ **computed directly from $\varepsilon_r(x,y)$'s own piecewise-constant profile** — i.e. $\eta_p(y)$ obtained by applying the piecewise-constant Fourier-coefficient formula ([TESolver §0.6](#)) to $1/\varepsilon_r(x,y)$ itself, not by algebraically inverting the Toeplitz matrix of $\varepsilon_r$'s coefficients.

It is a known result in the theory of Fourier-eigenmode methods for discontinuous periodic media that, once the mode sum is truncated to finite $N$, these two routes — (a) Fourier-expand $1/\varepsilon_r$ termwise, versus (b) Fourier-expand $\varepsilon_r$ and invert the resulting truncated Toeplitz matrix — are **not equivalent**, and converge at different rates as $N\to\infty$ at a jump discontinuity of $\varepsilon_r$. This is a consequence of how termwise multiplication of Fourier series behaves at discontinuities versus how matrix inversion of a truncated Toeplitz operator behaves; a short argument for why the two do not commute after truncation is given in [§7.3](#73-why-naive-termwise-fourier-inversion-is-inconsistent). Both routes converge to the same exact physics as $N\to\infty$; they differ only in truncated-$N$ behavior. Which convention is used to build $H(y)$ is therefore a numerical-convergence question, not a physical one — it does not change which equation is being solved, only how quickly the truncated system approximates it. This document states the coupled-mode equation in terms of $H(y)$ without committing to either construction; that choice belongs with the implementation, not the physics.

### 1.4 Physical Interpretation of the TM Coupling

As in TE ([TESolver §0.4](#)), the coupling matrix separates into a part carrying the tangential-momentum structure and a part carrying the material structure — but for TM both roles are played by the *same* matrix $H(y)$, entering twice: once sandwiched between the $\alpha_n,\alpha_m$ factors, and once directly weighting $\mathbf u'$. In a laterally unstructured layer ($\eta_p = 0$ for $p\neq0$), $H$ reduces to $\eta_0(y) I$, both occurrences collapse to a scalar multiple of the identity, and the system decouples mode-by-mode — the same qualitative statement as the TE limiting case ([TESolver §0.4](#)), verified quantitatively in [§2.2](#22-homogeneous-medium-limiting-case).

---

## 2. Local Layer Response

### 2.1 First-Order System

Following the same logic as [TESolver §1.1](#), but now using $w(y) \equiv H(y)\mathbf u'(y)$ as the natural companion field (rather than $\mathbf u'$ itself, since it is $w$, not $\mathbf u'$, that is differentiated as a unit in [§1.2](#12-coupled-mode-equation-for-tm)):

$$\mathbf u'(y) = H(y)^{-1}\, \mathbf w(y), \qquad \mathbf w'(y) = \big(D_\alpha H(y) D_\alpha - k_0^2 I\big)\, \mathbf u(y).$$

The state vector is again $\mathbf w_{\text{full}} = (\mathbf u,\, \mathbf w)$, of the same dimension $2(2N+1)$ as in TE ([TESolver §1.2](#), unchanged — the dimension-counting argument depends only on "second-order ODE per mode," which still holds). The structural difference from TE is that the coupling in the first equation is no longer the identity; $H(y)^{-1}$ appears explicitly. This is the direct algebraic consequence of $1/\varepsilon_r$ sitting inside the divergence in [§1.1](#11-maxwell-to-scalar-tm-wave-equation).

### 2.2 Homogeneous-Medium Limiting Case

In a homogeneous outer medium (no lateral structure), $H_{nm}(y) = \eta_0\,\delta_{nm} = (1/\varepsilon_r)\,\delta_{nm}$, a constant diagonal matrix. Substituting into the first-order system and eliminating $\mathbf w$ (full steps in [§7.4](#74-reduction-to-the-te-result-in-a-homogeneous-layer)) gives

$$u_n''(y) = \big(\alpha_n^2 - k^2\big)\, u_n(y), \qquad k^2 = k_0^2\varepsilon_r,$$

**exactly** the same equation obtained for TE in the homogeneous limit ([TESolver §1.3](#)). This is an essential consistency check: TE and TM must agree in any region without lateral material contrast, since polarization is only a meaningful distinction where there is a preferred in-plane direction — and this is confirmed algebraically here, not merely assumed. Consequently:

$$\boxed{\beta_n^2 = k^2 - \alpha_n^2}$$

is **identical in value and definition** to the TE case. All conclusions in TESolver about propagating/evanescent orders, the outer-medium branch choice, and the two independent asymptotic parameter sets $\beta_n^{\text{top}}, \beta_n^{\text{sub}}$ ([TESolver §1.4–1.5](#)) carry over to TM without modification.

### 2.3 Propagating and Evanescent Orders

Unchanged from TE — see [TESolver §1.4](#). The classification depends only on $\beta_n^2 = k^2-\alpha_n^2$, which by [§2.2](#22-homogeneous-medium-limiting-case) is the same quantity for both polarizations in the (necessarily homogeneous) outer half-spaces.

---

## 3. Boundary Conditions and Modal Amplitudes

### 3.1 The Actual Boundary Conditions for TM

For TE, the Maxwell boundary conditions require continuity of $u = E_z$ and (proportionally) $u'$, because both are tangential-field components at a horizontal interface ([TESolver §2.1](#)).

For TM, the tangential field components at a horizontal interface ($y=\text{const}$) are $H_z = u$ and $E_x$. Continuity of $u$ is unchanged. But $u'$ *by itself* is **not** the continuous quantity anymore — $E_x$ is proportional to $\eta\,\partial u/\partial y$ ([§5.1](#51-the-tangential-electric-field-for-tm)), so in mode space it is $w_n = [H(y)\mathbf u'(y)]_n$, not $u_n'$, that must match across a layer boundary.

$$\boxed{u_n \text{ continuous}, \qquad w_n = [H(y)\,\mathbf u'(y)]_n \text{ continuous}.}$$

This is the precise sense in which the TM boundary conditions differ from TE — not in *which* fields are matched (still one field and one "derivative-like" companion per mode), but in *what that companion physically is*. This is also exactly why $w$, not $\mathbf u'$, was chosen as the first-order companion variable in [§2.1](#21-first-order-system): the state vector $(\mathbf u,\mathbf w)$ is built directly out of the two physically continuous quantities, so the transfer-matrix formalism ([§4](#4-transfer-matrix-and-scattering-matrix)) applies to it without modification.

### 3.2 Modal Amplitudes from $u$ and $w$

In a homogeneous medium, $w_n = \eta_0\, u_n' = (1/\varepsilon_r)\,u_n'$. Repeating the TE derivation ([TESolver §2.2](#)) with this substitution (full steps in [§7.6](#76-derivation-of-the-tm-modal-amplitude-transformation)):

$$\boxed{A_n = \frac12\left(u_n + \frac{\varepsilon_r\, w_n}{i\beta_n}\right), \qquad B_n = \frac12\left(u_n - \frac{\varepsilon_r\, w_n}{i\beta_n}\right).}$$

The extra factor of $\varepsilon_r$ compared to the TE relation ([TESolver §2.2](#)) is the direct trace of $w = (1/\varepsilon_r)u'$ instead of $w=u'$ — inverting that relation to isolate $u'$ reintroduces $\varepsilon_r$.

---

## 4. Transfer Matrix and Scattering Matrix

Nothing about the transfer-matrix or scattering-matrix formalism itself depends on which physical fields make up the state vector — only on the fact that the governing equation is a linear first-order ODE system in a $2Q$-dimensional state, with a well-defined split into "top/bottom" and "in/out" amplitude pairs. Since [§2.1](#21-first-order-system) and [§3.1](#31-the-actual-boundary-conditions-for-tm) establish exactly that structure for TM (with $\mathbf w_{ull}=(\mathbf u,\mathbf w)$ in place of $(\mathbf u,\mathbf u')$), every result in TESolver's transfer- and scattering-matrix chapters transfers verbatim:

- Definition of the transfer matrix via a fundamental-solution basis: [TESolver §3.1–3.2](#), [§14.11](#)
- Block structure $(T_{11},T_{12},T_{21},T_{22})$: [TESolver §3.3](#)
- Numerical breakdown of direct transfer-matrix multiplication for evanescent modes: [TESolver §3.4](#), [§14.15](#)
- Physical idea of the S-matrix and the Redheffer star-product recursion: [TESolver §4](#), [§14.13–14.14](#)

The only substitution needed anywhere in this machinery is: wherever TESolver writes $(\mathbf a,\mathbf b)$ built from $(u_n, u_n')$, read $(\mathbf a,\mathbf b)$ built from $(u_n, w_n)$ per [§3.2](#32-modal-amplitudes-from-u-and-w) above.

---

## 5. From Amplitude to Power: TM Diffraction Efficiency

### 5.1 The Tangential Electric Field for TM

From Ampère's law, $\nabla\times\mathbf H = -i\omega\varepsilon_0\varepsilon_r\mathbf E$, with $\mathbf H = (0,0,u)$:
$$E_x = \frac{i}{\omega\varepsilon_0\varepsilon_r}\,\frac{\partial u}{\partial y}, \qquad E_y = -\frac{i}{\omega\varepsilon_0\varepsilon_r}\,\frac{\partial u}{\partial x}.$$

This confirms the claim used in [§3.1](#31-the-actual-boundary-conditions-for-tm): $E_x \propto \eta\,\partial u/\partial y = w$ (in the homogeneous outer medium where $\eta$ is a scalar), so continuity of the tangential $E$-field is exactly continuity of $w$.

### 5.2 Efficiency Formula and Its Difference from TE

For an outgoing order $u_n = B_n e^{i\alpha_n x + i\beta_n y}$ in a homogeneous outer medium, $\partial u_n/\partial y = i\beta_n u_n$, so from [§5.1](#51-the-tangential-electric-field-for-tm):
$$E_{x,n} \propto \frac{\beta_n}{\varepsilon_r}\,B_n\, e^{i\alpha_n x + i\beta_n y}.$$

The time-averaged normal Poynting flux is proportional to $\operatorname{Re}(E_{x,n} H_{z,n}^*) \propto \operatorname{Re}(\beta_n)/\varepsilon_r \cdot |B_n|^2$ (full derivation in [§7.7](#77-derivation-of-the-tm-diffraction-efficiency-formula), mirroring [TESolver §14.17](#) with $E_x$ and $H_z$ swapped relative to the TE case). Normalizing to the incident order:

$$\boxed{\eta_n^{\text{TM}} = |B_n|^2\,\frac{\operatorname{Re}(\beta_n)/\varepsilon_{r,n}}{\operatorname{Re}(\beta_0)/\varepsilon_{r,0}}.}$$

Compare with the TE result, $\eta_n = |B_n|^2\,\operatorname{Re}(\beta_n)/\operatorname{Re}(\beta_0)$ ([TESolver §7.2](#)). For **reflected** orders, $n$ and the incident order $0$ live in the *same* outer medium, so $\varepsilon_{r,n}=\varepsilon_{r,0}$ and the ratio reduces to the identical TE-looking expression. For **transmitted** orders, $n$ lives in the substrate while the incident order is normalized in the top medium, so $\varepsilon_{r,n}\neq\varepsilon_{r,0}$ in general, and this extra permittivity ratio is a genuine, physically necessary difference from the naive TE-analogous formula — it is not optional or convention-dependent, it follows directly from $E_x \propto \beta_n/\varepsilon_r$ rather than $E_z \propto$ (no $\varepsilon_r$ factor at all, cf. [TESolver §14.17](#)).

Evanescent orders again carry $\operatorname{Re}(\beta_n)=0$ and hence $\eta_n^{TM}=0$ far-field, for the same reason as TE ([TESolver §7.3](#)).

---

## 6. Physical Consistency Checks

All consistency checks from TESolver apply with the substitutions above:

- **Homogeneous-grating limit** ([TESolver §11.1](#)): with $\varepsilon_r(x,y)=\text{const}$, $H(y)$ reduces to a scalar multiple of the identity ([§1.4](#14-physical-interpretation-of-the-tm-coupling)), and the system decouples into plane-wave propagation — same qualitative statement as TE, confirmed quantitatively in [§2.2](#22-homogeneous-medium-limiting-case).
- **Energy balance** ([TESolver §11.2](#)): $\sum_n R_n + \sum_n T_n \le 1$ holds identically for TM, using $\eta_n^{\text{TM}}$ from [§5.2](#52-efficiency-formula-and-its-difference-from-te) in place of the TE efficiency.
- **Complex quantities and absorption** ([TESolver §11.3](#)): unchanged — complex $\varepsilon_r$ makes $H(y)$, $\beta_n$, $u_n(y)$, and $B_n$ complex, with the imaginary parts again encoding absorption rather than being numerical artifacts.
- **TE–TM cross-check specific to this document**: in any layer without lateral structure, TE and TM must produce identical $\beta_n$ and, for reflection, identical efficiencies — this was shown algebraically in [§2.2](#22-homogeneous-medium-limiting-case) and is a strong implementation-independent test unique to having both solvers available.

---

## 7. Mathematical Proofs

This section fills in only the steps that diverge from TESolver's proof chapter (§14). For anything not listed here, see the identically-numbered TESolver derivation (cross-referenced in [§0](#0-what-carries-over-unchanged-from-tesolver)).

### 7.1 From Maxwell's Equations to the Scalar TM Wave Equation

Start again from
$$\nabla\times\mathbf E = i\omega\mu_0\mathbf H, \qquad \nabla\times\mathbf H = -i\omega\varepsilon_0\varepsilon_r\mathbf E,$$
but now with $\mathbf H = (0,0,u(x,y))$ and $\mathbf E$ in the $x$–$y$ plane. From Ampère's law, writing out the curl of a purely-$z$ field:
$$\nabla\times\mathbf H = \left(\frac{\partial u}{\partial y},\,-\frac{\partial u}{\partial x},\,0\right) = -i\omega\varepsilon_0\varepsilon_r\,\mathbf E.$$
Hence
$$E_x = \frac{i}{\omega\varepsilon_0\varepsilon_r}\frac{\partial u}{\partial y}, \qquad E_y = -\frac{i}{\omega\varepsilon_0\varepsilon_r}\frac{\partial u}{\partial x}.$$

Now enforce Faraday's law, $z$-component only (since $\mathbf H$ has only a $z$-component, only the $z$-component of $\nabla\times\mathbf E$ is needed):
$$\frac{\partial E_y}{\partial x} - \frac{\partial E_x}{\partial y} = i\omega\mu_0\, u.$$

Substituting the expressions for $E_x, E_y$:
$$\frac{\partial}{\partial x}\!\left(-\frac{i}{\omega\varepsilon_0\varepsilon_r}\frac{\partial u}{\partial x}\right) - \frac{\partial}{\partial y}\!\left(\frac{i}{\omega\varepsilon_0\varepsilon_r}\frac{\partial u}{\partial y}\right) = i\omega\mu_0\, u.$$

Factor out $-i/(\omega\varepsilon_0)$:
$$-\frac{i}{\omega\varepsilon_0}\left[\frac{\partial}{\partial x}\!\left(\frac{1}{\varepsilon_r}\frac{\partial u}{\partial x}\right) + \frac{\partial}{\partial y}\!\left(\frac{1}{\varepsilon_r}\frac{\partial u}{\partial y}\right)\right] = i\omega\mu_0\, u.$$

Multiply both sides by $i\omega\varepsilon_0$ (using $i\cdot i = -1$ on the left after moving the sign, equivalently multiply by $-i\omega\varepsilon_0/(-i)=i\omega\varepsilon_0$ and simplify):
$$\frac{\partial}{\partial x}\!\left(\frac{1}{\varepsilon_r}\frac{\partial u}{\partial x}\right) + \frac{\partial}{\partial y}\!\left(\frac{1}{\varepsilon_r}\frac{\partial u}{\partial y}\right) = -\omega^2\mu_0\varepsilon_0\, u = -k_0^2\, u,$$
using $k_0^2=\omega^2\mu_0\varepsilon_0$ as in TE ([TESolver §14.1](#)). Rearranged:
$$\boxed{\frac{\partial}{\partial x}\!\left(\frac{1}{\varepsilon_r}\frac{\partial u}{\partial x}\right) + \frac{\partial}{\partial y}\!\left(\frac{1}{\varepsilon_r}\frac{\partial u}{\partial y}\right) + k_0^2 u = 0.}$$

This confirms [§1.1](#11-maxwell-to-scalar-tm-wave-equation). Unlike the TE derivation, no vector identity for $\nabla\times(\nabla\times\mathbf E)$ was needed, because the elimination was carried out on $\mathbf H$ directly via two scalar Maxwell equations rather than by taking a curl of a curl — a direct consequence of $\mathbf H$, not $\mathbf E$, being purely out-of-plane for TM.

### 7.2 Derivation of the TM Coupled-Mode Equation

Insert $u(x,y) = \sum_m u_m(y) e^{i\alpha_m x}$ and $\eta(x,y) = \sum_p \eta_p(y) e^{ipKx}$ into
$$\partial_x(\eta\,\partial_x u) + \partial_y(\eta\,\partial_y u) + k_0^2 u = 0.$$

**First term.** $\partial_x u = \sum_m i\alpha_m u_m(y) e^{i\alpha_m x}$, so
$$\eta\,\partial_x u = \sum_{p,m} \eta_p(y)\, i\alpha_m u_m(y)\, e^{i(pK+\alpha_m)x} = \sum_n\left(\sum_m i\alpha_m\, \eta_{n-m}(y)\, u_m(y)\right) e^{i\alpha_n x}$$
using the same index substitution $n=m+p$ as in TE ([TESolver §14.4](#)). Differentiating once more in $x$ brings down a factor $i\alpha_n$:
$$\partial_x(\eta\,\partial_x u) = \sum_n\left(-\alpha_n \sum_m \alpha_m\, \eta_{n-m}(y)\, u_m(y)\right) e^{i\alpha_n x} = -\sum_n\Big(D_\alpha H(y) D_\alpha\,\mathbf u(y)\Big)_n e^{i\alpha_n x},$$
with $H_{nm}(y)=\eta_{n-m}(y)$ and $D_\alpha=\mathrm{diag}(\alpha_n)$.

**Second term.** By the same Fourier-projection argument,
$$\eta\,\partial_y u = \sum_n\left(\sum_m \eta_{n-m}(y)\, u_m'(y)\right) e^{i\alpha_n x} = \sum_n \big(H(y)\,\mathbf u'(y)\big)_n\, e^{i\alpha_n x}.$$
Differentiating this expression once more in $y$ (the mode index $n$ and the basis $e^{i\alpha_n x}$ do not depend on $y$, so differentiation in $y$ commutes with the mode sum):
$$\partial_y(\eta\,\partial_y u) = \sum_n \frac{d}{dy}\Big[\big(H(y)\mathbf u'(y)\big)_n\Big]\, e^{i\alpha_n x}.$$

**Third term.** Trivially $k_0^2 u = \sum_n k_0^2 u_n(y)\, e^{i\alpha_n x}$.

Collecting all three and using linear independence of $\{e^{i\alpha_n x}\}$ mode-by-mode:
$$\frac{d}{dy}\Big[\big(H(y)\mathbf u'(y)\big)_n\Big] - \big(D_\alpha H(y) D_\alpha\, \mathbf u(y)\big)_n + k_0^2\, u_n(y) = 0,$$
i.e., in matrix form,
$$\boxed{\big(H(y)\,\mathbf u'(y)\big)' = \big(D_\alpha H(y) D_\alpha - k_0^2 I\big)\,\mathbf u(y),}$$
confirming [§1.2](#12-coupled-mode-equation-for-tm).

### 7.3 Why Naive Termwise Fourier Inversion Is Inconsistent

Let $\varepsilon_r(x)$ (fixed $y$) have Fourier coefficients $\varepsilon_p$ and let $\eta(x)=1/\varepsilon_r(x)$ have Fourier coefficients $\eta_p$. Exactly (infinite series),
$$\left(\sum_p \varepsilon_p e^{ipKx}\right)\left(\sum_q \eta_q e^{iqKx}\right) = 1 \quad\Longleftrightarrow\quad \sum_q \varepsilon_{n-q}\,\eta_q = \delta_{n0}\ \ \forall n,$$
i.e. the bi-infinite Toeplitz matrices built from $\{\varepsilon_p\}$ and $\{\eta_p\}$ are exact matrix inverses of one another.

Now truncate both to $|n|,|m|\le N$. Truncating an infinite matrix product and then asking whether the truncated matrices are still inverses of each other are **two operations that do not commute**: $\text{Toeplitz}_N(\{\varepsilon_p\})^{-1} \neq \text{Toeplitz}_N(\{\eta_p\})$ in general, because the truncated product $\sum_{q=-N}^{N}\varepsilon_{n-q}\eta_q$ for $|n|\le N$ still receives contributions in the exact (untruncated) identity from terms with $|q|>N$, which are discarded by truncation. This residual is largest precisely where $\varepsilon_r(x)$ has a jump discontinuity, since Fourier coefficients of a discontinuous function decay only as $1/p$, so the truncation error decays slowly. This is the origin of the factorization ambiguity flagged in [§1.3](#13-the-fourier-factorization-subtlety): "expand $\eta$ termwise" (giving $H=\text{Toeplitz}_N(\{\eta_p\})$) and "expand $\varepsilon_r$ and invert" (giving $H=\text{Toeplitz}_N(\{\varepsilon_p\})^{-1}$) are two different finite-$N$ approximations of the same exact operator, both converging to the same $N\to\infty$ limit but at different rates.

### 7.4 Reduction to the TE Result in a Homogeneous Layer

With $H(y) = \eta_0\, I$ constant (no lateral structure, so $H^{-1}$ is also constant), the first-order system of [§2.1](#21-first-order-system) is
$$\mathbf u' = \eta_0^{-1}\mathbf w = \varepsilon_r\,\mathbf w, \qquad \mathbf w' = \big(\eta_0 D_\alpha^2 - k_0^2 I\big)\mathbf u = \left(\frac{\alpha_n^2}{\varepsilon_r} - k_0^2\right)\delta_{nm}\, u_m.$$

Differentiate the first equation and substitute the second (both hold mode-by-mode, so drop indices):
$$u_n'' = \varepsilon_r\, w_n' = \varepsilon_r\left(\frac{\alpha_n^2}{\varepsilon_r} - k_0^2\right) u_n = \big(\alpha_n^2 - k_0^2\varepsilon_r\big)\, u_n = (\alpha_n^2 - k^2)\, u_n,$$
using $k^2 = k_0^2\varepsilon_r$. This is exactly the TE homogeneous-layer equation ([TESolver §14.7](#)), confirming [§2.2](#22-homogeneous-medium-limiting-case) and, with it, that $\beta_n^2=k^2-\alpha_n^2$ is unchanged between polarizations.

### 7.5 Derivation of the TM Boundary-Matching Variable

The Maxwell tangential-continuity conditions at a horizontal interface ($y$ = const, with unit normal $\hat y$) require the tangential components of $\mathbf E$ and $\mathbf H$ to be continuous: for TM, tangential $\mathbf H$ is just $H_z=u$ (already scalar and automatically tangential), and tangential $\mathbf E$ is $E_x$. By [§5.1](#51-the-tangential-electric-field-for-tm), $E_x \propto \eta\,\partial u/\partial y$, whose mode-space representation, by the same Fourier-projection argument as in [§7.2](#72-derivation-of-the-tm-coupled-mode-equation), is exactly $w_n = [H(y)\mathbf u'(y)]_n$. Hence continuity of tangential $E_x$ is, mode-by-mode, continuity of $w_n$ — establishing [§3.1](#31-the-actual-boundary-conditions-for-tm).

### 7.6 Derivation of the TM Modal-Amplitude Transformation

In a homogeneous medium, from [§7.4](#74-reduction-to-the-te-result-in-a-homogeneous-layer): $u_n=A_n+B_n$, $u_n'=i\beta_n(A_n-B_n)$ (same as TE, since the second-order equation is identical). The TM companion variable is $w_n = \eta_0 u_n' = (1/\varepsilon_r)\, i\beta_n(A_n-B_n)$, i.e.
$$\frac{\varepsilon_r\, w_n}{i\beta_n} = A_n - B_n.$$
Combined with $u_n = A_n+B_n$, adding and subtracting gives
$$A_n = \frac12\left(u_n + \frac{\varepsilon_r w_n}{i\beta_n}\right), \qquad B_n = \frac12\left(u_n - \frac{\varepsilon_r w_n}{i\beta_n}\right),$$
confirming [§3.2](#32-modal-amplitudes-from-u-and-w).

### 7.7 Derivation of the TM Diffraction-Efficiency Formula

For an outgoing order, $H_z = u_n = B_n\, e^{i\alpha_n x + i\beta_n y}$. From [§5.1](#51-the-tangential-electric-field-for-tm),
$$E_{x,n} = \frac{i}{\omega\varepsilon_0\varepsilon_r}\frac{\partial u_n}{\partial y} = \frac{i}{\omega\varepsilon_0\varepsilon_r}\, (i\beta_n)\, B_n\, e^{i\alpha_n x+i\beta_n y} = -\frac{\beta_n}{\omega\varepsilon_0\varepsilon_r}\, B_n\, e^{i\alpha_n x+i\beta_n y}.$$

The time-averaged Poynting vector's $y$-component is proportional to $\operatorname{Re}(E_x H_z^*)$ (up to the sign fixed by $\hat x\times\hat z=-\hat y$, which is common to every order and cancels in the normalization below):
$$\operatorname{Re}(E_{x,n} H_{z,n}^*) = -\frac{1}{\omega\varepsilon_0}\operatorname{Re}\!\left(\frac{\beta_n}{\varepsilon_r}\right)|B_n|^2 \ \overset{\varepsilon_r\ \text{real (lossless)}}{=}\ -\frac{1}{\omega\varepsilon_0}\,\frac{\operatorname{Re}(\beta_n)}{\varepsilon_r}\,|B_n|^2.$$

Normalizing to the incident order (which sees its own outer-medium permittivity $\varepsilon_{r,0}$) cancels the common prefactor $1/(\omega\varepsilon_0)$ and the sign, giving
$$\boxed{\eta_n^{\text{TM}} = |B_n|^2\,\frac{\operatorname{Re}(\beta_n)/\varepsilon_{r,n}}{\operatorname{Re}(\beta_0)/\varepsilon_{r,0}},}$$
confirming [§5.2](#52-efficiency-formula-and-its-difference-from-te). For an absorbing medium, $\varepsilon_r$ is complex and the $\operatorname{Re}(\cdot)$ must in general be taken of the full ratio $\beta_n/\varepsilon_r$ rather than $\beta_n$ and $\varepsilon_r$ separately — the expression above is written for the lossless case for clarity, matching the level of detail at which TESolver's §7.2/§14.17 state the TE formula.

---

## 8. Conclusion

The TM formalism is not an independent theory but the same electromagnetic boundary-value problem, reduced through the same conceptual chain as TE (periodicity → Floquet modes → 1D coupled system → local transfer relation → global scattering matrix → outgoing power), with a single structural change propagating through every step: **the material coefficient sits inside a derivative rather than multiplying the field directly.** This is what turns $M(y)$ from TE's algebraic $\alpha_n^2\delta_{nm}-k_{n-m}^2(y)$ into TM's $H(y)$-sandwiched $D_\alpha H(y) D_\alpha - k_0^2 I$ with a non-trivial first-order coupling $H(y)^{-1}$; what turns the TE boundary pair $(u,u')$ into the TM pair $(u,\,Hu')$; and what turns the TE efficiency's bare $\operatorname{Re}(\beta_n)$ into TM's $\operatorname{Re}(\beta_n)/\varepsilon_{r,n}$. Every other part of the derivation — Floquet expansion, transfer-matrix propagation, S-matrix recursion, layering, mode truncation, energy conservation — is genuinely polarization-independent and is inherited from `TESolver_Algorithm_Reference.md` without change.
