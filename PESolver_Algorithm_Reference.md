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
11. [The Central Conceptual Reduction](#11-the-central-conceptual-reduction)
12. [Physical Consistency Checks](#12-physical-consistency-checks)
    1. [Homogeneous-Grating Limit](#121-homogeneous-grating-limit)
    2. [Energy Balance](#122-energy-balance)
    3. [Complex Quantities and Absorption](#123-complex-quantities-and-absorption)
    4. [y-Dependent Geometry](#124-y-dependent-geometry)
13. [Why the Scattering Matrix Is More Physically Natural Than the Transfer Matrix](#13-why-the-scattering-matrix-is-more-physically-natural-than-the-transfer-matrix)
14. [Full Equation Chain (Summary)](#14-full-equation-chain-summary)
15. [Mental Model](#15-mental-model)
16. [Conclusion](#16-conclusion)

---

## 0. Physical Setup

Overview of the Theoretical Pipeline:

$$\text{Maxwell^{\prime}s Equations} \longrightarrow \text{Helmholtz Equation} \longrightarrow \text{Scalar Wave Equation (TE)}$$
$$\downarrow$$
$$\text{Periodic Fourier/Floquet Expansion} \longrightarrow \text{Coupled ODEs in } y \longrightarrow \text{Local Layer Response}$$
$$\downarrow$$
$$\text{Transfer Matrix} \longrightarrow \text{Scattering Matrix} \longrightarrow \text{Reflection Amplitudes} \longrightarrow \text{Diffraction Efficiencies}$$

Parameter Categorization:

- Geometry, Grating, and Material Structure: Defined by the spatial permittivity distribution $\varepsilon_r(x,y)$ and non-magnetic permeability ($\mu = \mu_0$). This encompasses the grating period $d$, layer/coating thicknesses along $y$ and the physical profile/surface topology of the grating.
- Incident Field Parameters: Free-space wavelength $\lambda$ and angle of incidence $\theta$


### 0.1 Maxwell to Scalar Wavefunction

Faraday^{\prime}s law of induction and Ampère^{\prime}s circuital law:
$$\vec{\nabla} \times \vec{E} = -\frac{\partial \vec{B}}{\partial t}, \quad \vec{\nabla} \times \vec{H} = -\vec{i} + \frac{\partial \vec{D}}{\partial t}$$

with $\vec{D} = \varepsilon_0 \varepsilon_r \vec{E}$ and $\vec{B} = \mu \vec{H}$ (assuming non-magnetic media where $\mu = \mu_0$ and $\vec i = \vec 0$), which can be formulated as:
$$\vec{\nabla} \times \vec{E} = -\mu_0 \frac{\partial \vec{H}}{\partial t}, \quad \vec{\nabla} \times \vec{H} = \varepsilon_0 \varepsilon_r \frac{\partial \vec{E}}{\partial t}$$

Separation of time and space assuming time-harmonic fields $\vec{E} e^{-i\omega t}$ and $\vec{H} e^{-i\omega t}$:
$$\vec{\nabla} \times \vec{E} = i\omega\mu_0\vec{H}, \quad \vec{\nabla} \times \vec{H} = -i \omega \varepsilon_0\varepsilon_r\vec{E}$$

To eliminate $\vec{H}$, take the curl of Faraday^{\prime}s law:
$$\vec{\nabla} \times (\vec{\nabla} \times \vec{E}) = i\omega\mu_0 (\vec{\nabla} \times \vec{H})$$

Substituting Ampère^{\prime}s law into the right-hand side yields:
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

At a fixed depth $y$, slicing the grating horizontally yields a 1D cross section. Over one period $d$, $k^2(x,y)$ is typically piecewise constant in $x$, taking on the values of whichever materials are present at that height - e.g. vacuum, substrate, one or more coating layers, or several material transitions in sequence.

The geometry therefore supplies, for each $y$, a function
$$x \mapsto k^2(x,y)$$
consisting of constant segments separated by jump discontinuities. Extracting this function is a purely geometric step: for a given $y$, one determines which material occupies each interval and at which $x$-positions the transitions between materials occur. This is the only place where the concrete grating profile enters the formalism - it produces the raw step function that subsequently becomes the input to the Fourier decomposition below. How this step function is obtained algorithmically depends on the profile representation and is deliberately left open here.

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

Numerically - and for the general solution theory below - it is convenient to rewrite
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
for some set of coefficients $c_i$ fixed by the boundary conditions between layers - determining these coefficients is the subject of the next section.

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
the two representations — field-and-derivative, and directional amplitude pairs — are fully equivalent. Physically, $A_n$ and $B_n$ correspond to the two opposite propagation directions of mode $n$; which one is "incoming" and which is "outgoing" depends on which outer boundary is considered and on the sign convention chosen for $y$.

## 3. Transfer Matrix

### 3.1 What a Single Layer Does Physically

Consider a layer between $y = y_-$ and $y = y_+$. The coupled ODE propagates the full state
$$\mathbf{w}(y_-) = \begin{pmatrix} \mathbf{u}(y_-) \\ \mathbf{u}'(y_-) \end{pmatrix}$$
forward to $\mathbf{w}(y_+)$. Because the governing equation is linear, this propagation defines a linear map
$$\mathbf{w}(y_+) = \mathcal{T}\,\mathbf{w}(y_-),$$
the **transfer matrix** of the layer. Importantly, $\mathcal{T}$ is not "the physics of the grating" as a whole — it is only the linear map from one layer's state to the next.

### 3.2 Why a Basis of Trial Solutions Generates the Transfer Matrix

Choose a basis $\mathbf{e}_1, \dots, \mathbf{e}_{2Q}$ (with $Q = 2N+1$) of the state space. For each basis vector, the coupled ODE system is solved starting from that vector as initial condition and integrated to $y_+$. The resulting end states form the columns of a fundamental solution matrix $\Phi(y)$. The transfer matrix is then
$$\mathcal{T} = \Phi(y_+)\,\Phi(y_-)^{-1}.$$
This is the mathematical meaning behind solving the system separately for many independent trial solutions, rather than for a single physical solution directly.

### 3.3 Block Structure of the Transfer Matrix

Writing the state in terms of directional amplitudes,
$$\begin{pmatrix} \mathbf{a}_+ \\ \mathbf{b}_+ \end{pmatrix} =
\begin{pmatrix} T_{11} & T_{12} \\ T_{21} & T_{22} \end{pmatrix}
\begin{pmatrix} \mathbf{a}_- \\ \mathbf{b}_- \end{pmatrix},$$
where $(\mathbf{a}, \mathbf{b})$ denote the appropriate directional amplitude pairs on each side. The precise labeling of the blocks is a matter of convention; what matters physically is that a single layer establishes a linear relation between four sets of amplitudes, two on each side of the layer.

### 3.4 Why Direct Multiplication of Transfer Matrices Is Numerically Problematic

An evanescent mode can appear across a layer of thickness $h$ with factors as extreme as $e^{+\gamma h}$ and $e^{-\gamma h}$ simultaneously. A single transfer matrix can therefore contain both very large and very small numbers at once. For a stack of many layers,
$$\mathcal{T}_{\text{total}} = \mathcal{T}_M \cdots \mathcal{T}_2\,\mathcal{T}_1,$$
such extreme values compound multiplicatively. Even when the final physical result is moderate, intermediate values in this product can become numerically catastrophic (loss of precision or overflow). This is the physical and numerical motivation for splitting the structure into many thin layers and combining their responses through a scattering-matrix recursion instead of direct transfer-matrix multiplication — a point taken up in the next chapter.


## 4. Scattering Matrix

### 4.1 Physical Idea of the S-Matrix

The transfer matrix answers: how are the fields/modes on the top and bottom side of a layer related?

The scattering matrix instead answers: which outgoing waves result from the incoming waves?
$$\begin{pmatrix} \mathbf{a}_{\text{top,out}} \\ \mathbf{a}_{\text{bot,out}} \end{pmatrix}
= S \begin{pmatrix} \mathbf{a}_{\text{top,in}} \\ \mathbf{a}_{\text{bot,in}} \end{pmatrix}.$$

This formulation is physically much better conditioned: growing evanescent solutions never need to be tracked as primary global unknowns, since only physically bounded incoming and outgoing amplitudes appear as variables.

### 4.2 Recursion Over Layers

Suppose a new layer is appended to an already-assembled stack. The new layer contributes a transfer relation; the already-assembled stack below it is represented by its own scattering matrix. The amplitudes at the shared interface between the new layer and the stack are not directly observable — they are internal to the combined structure and can be eliminated algebraically. This elimination is precisely what produces the recursive combination rule for scattering matrices, commonly known as the **Redheffer star product**.

Given two scattering matrices,
$$S_A = \begin{pmatrix} S_{11}^A & S_{12}^A \\ S_{21}^A & S_{22}^A \end{pmatrix}, \qquad
S_B = \begin{pmatrix} S_{11}^B & S_{12}^B \\ S_{21}^B & S_{22}^B \end{pmatrix},$$
representing two adjacent sub-structures (e.g. the stack assembled so far, and the newly added layer), the combined scattering matrix $S = S_A \star S_B$ is
$$S_{11} = S_{11}^A + S_{12}^A\left(I - S_{11}^B S_{22}^A\right)^{-1} S_{11}^B S_{21}^A,$$
$$S_{12} = S_{12}^A\left(I - S_{11}^B S_{22}^A\right)^{-1} S_{12}^B,$$
$$S_{21} = S_{21}^B\left(I - S_{22}^A S_{11}^B\right)^{-1} S_{21}^A,$$
$$S_{22} = S_{22}^B + S_{21}^B\left(I - S_{22}^A S_{11}^B\right)^{-1} S_{22}^A S_{12}^B.$$

The essential physical point is independent of the exact algebraic form:
$$\boxed{\text{internal waves at the shared interface are eliminated; only the external input/output amplitudes remain.}}$$

The matrix inversions appearing above are, in practice, never carried out as explicit inversions but through numerically stable factorization methods — this is a numerical-stability concern, not a change to the underlying physics.

## 5. Layering as a Numerical Discretization

The real permittivity profile is, in general, continuous in $y$. The solver approximates it by many thin layers, each treated as locally uniform along $y$ so that the coupled ODE can be integrated with locally constant coefficients within that slice.

In the limit
$$\Delta y \to 0,$$
the layered decomposition recovers the continuous differential equation exactly. Layering is therefore not an additional physical assumption about the material — it is a numerical discretization in the vertical direction. The physically motivated criterion for how finely to discretize is not primarily the location of material boundaries, but the maximum evanescent growth rate present in the structure and the numerical stability this requires (cf. §3.4).

## 6. The Complete Physical Flow Through the Solver

The full solver can now be summarized as a physical chain:

1. **Incident wave.** Given $\lambda$ and $\theta$, compute $k_0 = 2\pi/\lambda$.
2. **Tangential momentum.** The grating periodicity restricts the allowed tangential wavenumbers to $\alpha_n = k_{\text{top}}\sin\theta + nK$.
3. **Local material distribution.** For each $y$, the geometry determines $\varepsilon_r(x,y)$.
4. **Fourier decomposition.** This yields the coefficients $k_n^2(y)$.
5. **Mode coupling.** The Fourier coefficients generate the coupled-mode system
   $$u_n''(y) = \sum_m \left[\alpha_n^2\delta_{nm} - k_{n-m}^2(y)\right] u_m(y).$$
6. **Vertical propagation.** This ODE system is integrated through the structure.
7. **Layer response.** Each layer's integration yields a transfer relation between its two boundaries.
8. **Assembly.** The layers are combined via scattering-matrix recursion (§4.2).
9. **Outer boundary condition.** In the homogeneous medium above the structure, the solution is decomposed into physically outgoing diffraction orders.
10. **Reflection amplitudes.** The resulting coefficients are the reflected diffraction amplitudes, obtained by evaluating the final scattering matrix for the physically prescribed incident condition (§6.1).

### 6.1 Origin of the Single-Order Incident Condition

There are not arbitrarily many independent inputs to consider. Physically, exactly one incident mode is prescribed: the zeroth diffraction order of the beam incident from above. In the outer region one therefore sets
$$a_0 = 1, \qquad a_n = 0 \ \text{ for all } n \neq 0,$$
i.e. the incident state is the standard basis vector $\mathbf{e}_0$ in mode space.

Multiplying the scattering matrix by this basis vector selects exactly the corresponding column of $S$ — this is why the reflection amplitudes can be read off directly from a single column of the final scattering matrix, rather than requiring a full matrix–vector solve for a general incident condition. With modes indexed symmetrically from $-N$ to $+N$, the order $n=0$ corresponds to the middle position in that indexing, which is why this particular column is singled out.

This is the missing physical link between "the scattering matrix has been assembled" and "the reflection amplitudes are read off from one specific column of it."


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

## 8. Truncation and Its Physical Meaning

### 8.1 Finite Mode Number $N$

The exact theory involves $n \in \mathbb{Z}$, an infinite set of Floquet orders. Any practical solution truncates this to $n=-N,\dots,+N$:
$$\boxed{\text{infinitely many Floquet orders} \ \longrightarrow\ 2N+1 \text{ numerical orders.}}$$
The underlying physics is unchanged by this truncation; only the numerically tracked mode space is made finite. This introduces a convergence question distinct from the physical derivation itself: whether $N$ is large enough that the computed efficiencies are insensitive to further increasing it.

### 8.2 Relation Between Evanescent Decay and Discretization Scale

An evanescent order carries factors $e^{\pm\gamma_n y}$. Across a layer of thickness $\Delta y$, this produces a magnitude factor $e^{|\gamma_n|\Delta y}$. If this factor becomes too large, very large and very small numbers coexist within the same numerical matrices (cf. §3.4), degrading accuracy. Keeping $|\beta_{\max}|\Delta y$ bounded is therefore the physical criterion that ties the required layer resolution to the evanescent decay rates present in the mode spectrum — the near-field character of high orders is directly linked to how finely the vertical direction must be discretized.

## 9. The Coating Causal Chain

If a coating changes the local material distribution $\varepsilon_r(x,y)$, this changes $k_p^2(y)$, which changes the coupling matrix $M_{nm}(y)$, which changes the vertical mode profiles $u_n(y)$, which changes the layer transfer relations, which changes the assembled scattering matrix, which changes the reflection amplitudes $B_n$, and therefore the efficiencies $\eta_n$:
$$\boxed{\text{coating geometry} \rightarrow \varepsilon_r(x,y) \rightarrow k_n^2(y) \rightarrow M(y) \rightarrow u_n(y) \rightarrow S \rightarrow B_n \rightarrow \eta_n.}$$
This is the complete causal path by which any change to the layer geometry enters the diffraction efficiencies — there is no other route.

## 10. What This Formalism Does Not Do

For clarity, it is worth stating explicitly what the method does *not* do:

- it does not evaluate a full 2D field on a real-space $(x,y)$ grid as its primary variable;
- it does not trace an independent ray-optics path per diffraction order;
- it does not treat diffraction as a sequence of independent Snell refractions;
- it does not determine the diffraction orders only at the end of the calculation — they are built into the problem from the start via the Floquet basis.

Instead, the fundamental variable throughout is the field superposition
$$u(x,y) = \sum_n u_n(y)\,e^{i\alpha_n x},$$
so diffraction is present in the degrees of freedom from the very first step, not derived afterward.

## 11. The Central Conceptual Reduction

The physical reduction underlying the entire method is:
$$\boxed{\text{2D electromagnetic problem} \ \Longrightarrow\ \text{1D system of coupled modes.}}$$
The $x$-direction is not eliminated by ignoring it, but by encoding its periodicity exactly into the Fourier/Floquet basis. What remains is a single continuous variable, $y$. This reduction is the conceptual core of the whole formalism.

## 12. Physical Consistency Checks

### 12.1 Homogeneous-Grating Limit

A strong test of the full formalism is the limit $\varepsilon_r(x,y) = \text{const}$. Then $k_n^2 = 0$ for all $n\neq0$, all coupling vanishes, and the equations decouple into
$$u_n'' + \beta_n^2 u_n = 0.$$
Only the actually excited incident order remains active, and the system reduces to ordinary plane-wave propagation in a homogeneous medium. A correctly formulated solver must reproduce exactly this trivial physics in this limit.

### 12.2 Energy Balance

For lossless materials, summing over all propagating reflected and transmitted orders should satisfy
$$\boxed{\sum_n R_n + \sum_n T_n = 1}$$
up to numerical error and normalization conventions. For absorbing materials,
$$\sum_n R_n + \sum_n T_n < 1,$$
since part of the energy is dissipated in the material. This energy check is one of the most important independent validation tests for any implementation of the method.

### 12.3 Complex Quantities and Absorption

If the material is absorbing, the refractive index $n = n' + in''$ is complex, and so is $k^2$. Consequently $k_n^2(y)$, $M(y)$, $\beta_n$, $u_n(y)$, and $B_n$ all become complex quantities. Their imaginary parts are not numerical artifacts — they physically encode phase shift and absorption.

### 12.4 $y$-Dependent Geometry

For a strictly rectangular profile, the lateral material distribution within a given layer is constant in $y$, so $\partial M/\partial y = 0$ within that layer. For a sloped or curved profile, the transition positions $x_p(y)$ instead vary continuously with $y$, so $k_n^2(y)$ — and hence $M(y)$ — becomes explicitly $y$-dependent within a layer, not just from layer to layer. Treating $M$ as piecewise constant in $y$ (i.e. neglecting $\partial M/\partial y$ within a layer) is exact only for rectangular profiles; for sloped or curved profiles it is an approximation whose accuracy improves as the layers are made thinner. This does not change the underlying physics — it only affects the accuracy of the derivative information available to the numerical integration.

## 13. Why the Scattering Matrix Is More Physically Natural Than the Transfer Matrix

The transfer matrix tracks a full internal state vector, including combinations of growing and decaying evanescent contributions. The scattering matrix instead directly tracks the relation between the physically and experimentally relevant quantities: incoming amplitudes in, outgoing amplitudes out. For an optical measurement, it is precisely the outgoing diffraction orders that are of interest. The scattering-matrix formulation is therefore not only numerically more stable (§3.4, §4.1) but also conceptually closer to the physically measured quantity.

## 14. Full Equation Chain (Summary)

$$\varepsilon_r(x,y)\ \rightarrow\ k^2(x,y)=k_0^2\varepsilon_r(x,y)\ \rightarrow\ \alpha_n = k_{\text{top}}\sin\theta + n\frac{2\pi}{d}$$
$$u(x,y) = \sum_n u_n(y)\,e^{i\alpha_n x}, \qquad k^2(x,y) = \sum_p k_p^2(y)\,e^{ipKx}$$
$$u_n''(y) = \sum_m\left[\alpha_n^2\delta_{nm} - k_{n-m}^2(y)\right]u_m(y), \qquad \beta_n^2 = k_{\text{medium}}^2 - \alpha_n^2$$
$$\mathbf{w}(y_+) = \mathcal{T}_{\text{layer}}\,\mathbf{w}(y_-), \qquad \mathbf{a}_{\text{out}} = S\,\mathbf{a}_{\text{in}}, \qquad \mathbf{a}_{\text{in}} = \mathbf{e}_0$$
$$B_n = (S\mathbf{e}_0)_n \quad \text{(up to the phase reference of §7.1)}, \qquad \eta_n = |B_n|^2\,\frac{\operatorname{Re}\beta_n}{\operatorname{Re}\beta_0}$$

## 15. Mental Model

The simplest complete picture of this formalism is:

- The grating is a medium that couples different tangential Fourier/Floquet waves to one another.
- Periodicity determines which waves are allowed at all: $\alpha_n = \alpha_0 + nK$.
- The material structure determines how strongly they are coupled: $k_{n-m}^2(y)$.
- Vertical propagation determines how these coupled amplitudes evolve with height: $u'' = M(y)\,u$.
- The outer media determine which of these solutions are propagating or evanescent: $\beta_n^2 = k^2-\alpha_n^2$.
- The scattering matrix determines which outgoing waves result from the one incident beam.
- Their normal power flow yields the measured diffraction efficiency.

## 16. Conclusion

Physically, this method is not "an algorithm for computing diffraction efficiencies" but the numerical realization of a well-posed electromagnetic boundary-value problem:
$$\boxed{\text{periodic Maxwell problem} \rightarrow \text{coupled Floquet modes} \rightarrow \text{vertical ODE} \rightarrow \text{global scattering matrix} \rightarrow \text{diffracted power.}}$$

Geometry enters the physics exclusively through $\varepsilon_r(x,y) \rightarrow k_n^2(y)$; periodicity generates $\alpha_n = \alpha_0 + nK$; and the coupling of these modes produces diffraction. The two conceptual bridges that tie the whole derivation together are
$$\boxed{\text{Maxwell} \rightarrow \text{scalar TE-Helmholtz} \rightarrow \text{Floquet} \rightarrow \text{coupled ODE}}$$
and
$$\boxed{\text{outer boundary condition} \rightarrow \text{scattering-matrix column} \rightarrow \text{Poynting flux}.}$$