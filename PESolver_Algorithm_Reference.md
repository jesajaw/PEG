# PESolver Algorithm Reference

## Purpose and Scope

This document explains how `PESolver` computes the diffraction efficiency of a grating for a single (wavelength, incidence angle) pair. Or to be more acurate, it shows the physical validation how and why you can compute something likes this and also how it is (numerical) implemented. It does not covers the geometry aspect: it consumes this input from `PEG`. It also does not include the MPI orchestration in `mainMPI`, or the material refractive-index database lookup — those are documented separately.

The method implemented here is a **differential (coupled-wave) method**: the structure is treated as a stack of thin horizontal layers; within each layer, a coupled ordinary differential equation (ODE) system is integrated numerically along the vertical coordinate `y`; the per-layer results are combined into a numerically stable overall response via **S-matrix recursion**. This is closely related to the family of methods used in Rigorous Coupled-Wave Analysis (RCWA) and multilayer optics/ellipsometry.

---

## Table of Contents

0. [Physical Setup](#0-physical-setup)
    1. [Maxwell to Scalar Wavefunction](#01-maxwell-to-scalar-wavefunction)
    2. [Fourier & Floquet Expansion](#02-fourier-and-floquet-expansion)
    3. [Coupled Mode Equation Derivation](#03-coupled-mode-equation-derivation)
    4. [Physical Interpretation](#04-physical-interpretation)
    5. [Geometry Bridg](#05-geometry-bridg)
    6. [Fourier Coefficients of a Piecewise-Constant Layer](#06-fourier-coefficients-of-a-piecewise-constant-layer)

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




### 1.1 The grating equation — `computeAlphaAndBeta()`

A grating with period `d` illuminated at incidence angle `theta` can only scatter light into a discrete set of directions, indexed by an integer diffraction order `n` (from `-N` to `N`). This is the grating equation.

Every wave in the problem — incident, or any diffracted order — has a wavevector with a **tangential** component (parallel to the surface, along `x`) and a **normal** component (perpendicular, along `y`). The code:

```cpp
double alpha = k_2 * sin(theta_2) + 2 * M_PI * n / d;
alpha_[i] = alpha;
```

Physically: the incident wave contributes a fixed tangential wavevector `k·sin(theta)`; the grating^{\prime}s periodicity contributes an additional `2π·n/d` per order `n` (the "grating momentum kick"). This is the diffraction equation in wavevector form:

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

The fix: don^{\prime}t integrate the whole height `a` in one shot. Split it into `numLayers_` thinner slices, chosen so that even the *most* evanescent order (largest `|beta_n|`, which occurs at the highest or lowest computed diffraction order — hence checking both `betaM_[0]` and `betaM_[2*N_]`) stays within a safe exponential range over a single layer:

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
