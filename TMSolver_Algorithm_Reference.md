# TMSolver Algorithm Reference

## Purpose and Scope

This document is the TM-polarization counterpart to `TESolver_Algorithm_Reference.md`. It follows the same purpose: explaining the physical and mathematical basis for computing diffraction efficiencies for a single (wavelength, incidence angle) pair — not the concrete implementation, not geometry extraction, not MPI orchestration.

**This document does not repeat derivations that are identical for TE and TM.** The overall solution strategy — Floquet expansion, transfer-matrix propagation through layers, scattering-matrix recursion, mode truncation, layering as a numerical discretization — is polarization-independent machinery. Only the pieces that actually change under $E_z \to H_z$ are re-derived here. Everywhere else, this document points to the corresponding section of `TESolver_Algorithm_Reference.md`.

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
| Bloch–Floquet origin of $\alpha_n = k_{\text{top}}\sin\theta + nK$ | §0.2, §14.2 |
| Geometry bridge (extracting the $x$-cross-section at fixed $y$) | §0.5 |
| Dimension-counting argument for the state space | §1.2, §14.12 |
| Fundamental-matrix argument: why integrating trial solutions yields $\mathcal T$ | §3.2, §14.11 |
| Why direct transfer-matrix multiplication is numerically unstable | §3.4, §14.15 |
| Scattering-matrix idea and the Redheffer star product | §4, §14.13–14.14 |
| Layering as a numerical discretization of continuous $y$-dependence | §5, §14.22–14.23 |
| Origin of the single-order incident condition $\mathbf a_{\text{in}} = \mathbf e_0$ | §6.1, §14.16 |
| Phase referencing to a common origin | §7.1, §14.21 |
| Mode truncation ($N\to\infty$ vs. finite $N$) | §8.1, §14.24 |
| Energy-balance sum rule $\sum_n R_n + \sum_n T_n \le 1$ | §11.2, §14.19 |
| Complex refractive index and absorption entering as complex coefficients | §11.3, §14.20 |

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

Nothing about the transfer-matrix or scattering-matrix formalism itself depends on which physical fields make up the state vector — only on the fact that the governing equation is a linear first-order ODE system in a $2Q$-dimensional state, with a well-defined split into "top/bottom" and "in/out" amplitude pairs. Since [§2.1](#21-first-order-system) and [§3.1](#31-the-actual-boundary-conditions-for-tm) establish exactly that structure for TM (with $\mathbf w_{\text{full}}=(\mathbf u,\mathbf w)$ in place of $(\mathbf u,\mathbf u')$), every result in TESolver's transfer- and scattering-matrix chapters transfers verbatim:

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

Evanescent orders again carry $\operatorname{Re}(\beta_n)=0$ and hence $\eta_n^{\text{TM}}=0$ far-field, for the same reason as TE ([TESolver §7.3](#)).

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
