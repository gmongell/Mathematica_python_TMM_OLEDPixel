# OLED optics, spectral EQE and ensemble-resource derivations

**Research record: 2026-10-16.** This document uses GitHub-supported LaTeX
math, rendered directly in the repository. Open the rendered file on GitHub
or in a Markdown preview with math support, rather than a plain-text viewer.
Equations are restated with project-specific explanations; no third-party
chapter text or figures are reproduced.

## 1. Scope and reading map

This is the consolidated mathematical record for the implemented foundations,
their rationale, numerical transformations and planning estimates.

| Subject | Status | Implementation / companion |
|---|---|---|
| Homogeneous angular emission and one-interface source work | Implemented, benchmarked reference kernels | [dipole_reference.py](../dipole_reference.py), [formulation](DIPOLE_FORMULATION.md) |
| Native optical data, dispersion and layer geometry | Implemented | [spectral_eqe.py](../spectral_eqe.py), [data guide](SPECTRAL_EQE_AND_DATA.md) |
| Dispersive incident-plane-wave R/T/A | Implemented diagnostic, not emitted flux | [spectral_eqe.py](../spectral_eqe.py) |
| Monochromatic/broadband EQE from supplied physical rate budgets | Implemented composition, not automatic stack-derived EQE | [spectral_eqe.py](../spectral_eqe.py) |
| Fixed versus reorienting dipole ensembles | Implemented for monochromatic transition arrays | [spectral_eqe.py](../spectral_eqe.py) |
| Planar workload, finite-pixel FDTD memory/work and optional roofline | Planning models only | [simulation_resources.py](../simulation_resources.py), [resource guide](ENSEMBLE_RESOURCE_ESTIMATES.md) |
| Finite emitting-layer multilayer source/escape calculation | Not implemented | [solver acceptance specification](OLED_SOLVER_SPEC.md) |
| Spatially integrated device EQE, finite-pixel FDTD, GPU backend | Not implemented | Future work; formulas below identify requirements, not delivered solvers |

In particular, source work, external illumination transmission, collected
emitted power and EQE are different observables. A mathematically valid
composition formula does not validate its supplied material data or rates.

## 2. Conventions, units and power normalization

Assume linear, isotropic, nonmagnetic media and harmonic peak amplitudes with
time dependence

$$
\mathbf E(\mathbf r,t)=\operatorname{Re}
 \left[\mathbf E(\mathbf r)e^{-i\omega t}\right],
\qquad
k_0=\frac{2\pi}{\lambda_0}=\frac{\omega}{c}.
$$

The relative index and permittivity are

$$
N_j(\lambda_0)=n_j(\lambda_0)+i\kappa_j(\lambda_0),
\qquad \epsilon_j=N_j^2.
$$

The symbol $\kappa$ denotes extinction coefficient, not a wavevector.
Passive data have $\kappa_j\geq0$ and, in the supported index convention,
$n_j\geq0$. The emitting reference host has real $n_h>0$ and
$k_h=n_hk_0$. Fields propagate as $e^{ik_z z}$ in their outgoing coordinate;
choose longitudinal roots that propagate outward or decay rather than grow.
The API uses vacuum wavelengths and geometric lengths in nm. Ratios and
phases must use consistent length units; SI units are required in power
formulas with $c$, $\epsilon_0$ and a dipole moment in C m.

For a prescribed electric dipole $\mathbf p$, the homogeneous-host power is

$$
P_h=\frac{n_h\omega^4|\mathbf p|^2}{12\pi\epsilon_0c^3}.
\tag{1}
$$

This is a macroscopic lossless-host normalization for a **peak** amplitude,
not an RMS amplitude. No microscopic local-field correction is included.
Using the same dipole and host reference throughout is essential when
combining rates from different orientations or collection channels.

The field radiated back onto the dipole changes source work. With the
outgoing electric dyadic Green tensor defined by
$\mathbf E_{\rm ref}=\mu_0\omega^2\mathbf G_{\rm ref}\mathbf p$,

$$
P_{\rm ref}=\frac{\omega}{2}
 \operatorname{Im}\!\left[\mathbf p^*\cdot\mathbf E_{\rm ref}\right],
\qquad
F_{\rm EM}=1+\frac{6\pi}{k_h}
 \operatorname{Im}\!\left[
 \hat{\mathbf e}^{\,*}\cdot\mathbf G_{\rm ref}\hat{\mathbf e}\right],
\tag{2}
$$

where $\hat{\mathbf e}=\mathbf p/|\mathbf p|$. The coefficient follows by
dividing the source-work expression by (1), using
$\mu_0\epsilon_0c^2=1$. The Green tensor has units of inverse length.
For absorbing surroundings, source work includes dissipated electromagnetic
power as well as outgoing radiation. It excludes intrinsic emitter decay
that is not represented by the electromagnetic environment.

## 3. Homogeneous angular pattern and polarization

For observation direction $\hat{\mathbf r}$, the transverse dipole field is
proportional to
$\hat{\mathbf r}\times(\hat{\mathbf e}\times\hat{\mathbf r})$.
Its squared magnitude is $1-|\hat{\mathbf e}\cdot\hat{\mathbf r}|^2$.
Normalizing its solid-angle integral gives

$$
\frac{1}{P_h}\frac{dP}{d\Omega}
 =\frac{3}{8\pi}
  \left(1-|\hat{\mathbf e}\cdot\hat{\mathbf r}|^2\right),
\qquad
d\Omega=\sin\theta\,d\theta\,d\phi.
\tag{3}
$$

Indeed, $\int d\Omega=4\pi$ and
$\int|\hat{\mathbf e}\cdot\hat{\mathbf r}|^2d\Omega=4\pi/3$,
so the normalized integral is one. Reflection symmetry gives one-half
in either homogeneous hemisphere.

Project the dipole onto the orthogonal far-field polarization axes before
squaring:

$$
\frac{dP_s}{P_h\,d\Omega}=\frac{3}{8\pi}
 |\hat{\mathbf e}\cdot\hat{\boldsymbol\phi}|^2,
\qquad
\frac{dP_p}{P_h\,d\Omega}=\frac{3}{8\pi}
 |\hat{\mathbf e}\cdot\hat{\boldsymbol\theta}|^2.
\tag{4}
$$

| Dipole distribution | $s$ factor, excluding $3/(8\pi)$ | $p$ factor, excluding $3/(8\pi)$ |
|---|---|---|
| Vertical, $z$ | $0$ | $\sin^2\theta$ |
| Horizontal, $x$ | $\sin^2\phi$ | $\cos^2\theta\cos^2\phi$ |
| Horizontal, $y$ | $\cos^2\phi$ | $\cos^2\theta\sin^2\phi$ |
| Uniform random in-plane | $1/2$ | $\cos^2\theta/2$ |
| Uniform isotropic | $1/3$ | $1/3$ |

These averages are averages of **powers**. Averaging unrelated complex
dipole amplitudes first would incorrectly cancel incoherent emission.

## 4. One-interface source-work kernels

Place the interface at $z=0$ and the dipole at $z=d>0$ in medium 1,
with $\epsilon_1=n_h^2$. Medium 2 occupies $z<0$.
For each in-plane wavevector magnitude, define

$$
q=\frac{k_\parallel}{k_h},\qquad
w=\sqrt{1-q^2},\qquad
v=\sqrt{\frac{\epsilon_2}{\epsilon_1}-q^2},
\qquad a=2k_hd=\frac{4\pi n_hd}{\lambda_0}.
\tag{5}
$$

On the real $q$ axis, $w$ is positive real below the host light line
and positive imaginary above it. The substrate root uses the corresponding
outgoing/decaying prescription. Electric-field Fresnel reflection amplitudes
are

$$
r_s=\frac{w-v}{w+v},\qquad
r_p=\frac{\epsilon_2w-\epsilon_1v}
           {\epsilon_2w+\epsilon_1v}.
\tag{6}
$$

Identical media give $r_s=r_p=0$. The explicit perfect-electric-conductor
boundary uses $r_s=-1$, $r_p=1$; it is an ideal boundary, not a finite
tabulated metal index.

The angular-spectrum representation of the reflected Green tensor,
followed by azimuthal integration and division by (1), gives

$$
F_\perp=1+\frac32\operatorname{Re}\int_0^\infty
 \frac{q^3}{w}r_p e^{iaw}\,dq,
\tag{7}
$$

$$
F_\parallel=1+\frac34\operatorname{Re}\int_0^\infty
 \frac{q}{w}\left(r_s-w^2r_p\right)e^{iaw}\,dq.
\tag{8}
$$

Here $\perp$ means dipole normal to the interface, not the $s$ polarization.
The phase is a round trip from source to interface and back.
The **minus sign** in $r_s-w^2r_p$ comes from the reflected $p$-polarization
dyad: incident and reflected polarization vectors have different longitudinal
directions. An intensity average cannot replace this signed field product.
The ideal-mirror tests below independently constrain this sign.

The $q<1$ / $q>1$ split classifies propagation **in the host**.
A host-evanescent component can propagate in a higher-index substrate.
Thus neither the evanescent correction nor $F_{\rm EM}$ alone is an
absorption fraction or an external collection efficiency.

### 4.1 Removing the host light-line singularity

At $q=1$, the $1/w$ factor is integrable but unsuitable for direct point
sampling. In the propagating sector use $u=\sqrt{1-q^2}$, so
$q\,dq=-u\,du$. Reversing the integration limits yields

$$
\Delta F_{\perp,\rm prop}=\frac32\int_0^1
 (1-u^2)\operatorname{Re}[r_pe^{iau}]\,du,
\tag{9}
$$

$$
\Delta F_{\parallel,\rm prop}=\frac34\int_0^1
 \operatorname{Re}[(r_s-u^2r_p)e^{iau}]\,du.
\tag{10}
$$

For the evanescent sector set $t=\sqrt{q^2-1}$, then $w=it$,
$q\,dq=t\,dt$ and $e^{iaw}=e^{-at}$. Finally let $x=at$:

$$
\Delta F_{\perp,\rm evan}=\frac{3}{2a}\int_0^\infty
 \left[1+\left(\frac{x}{a}\right)^2\right]
 e^{-x}\operatorname{Im}r_p\,dx,
\tag{11}
$$

$$
\Delta F_{\parallel,\rm evan}=\frac{3}{4a}\int_0^\infty
 e^{-x}\operatorname{Im}\left[
 r_s+\left(\frac{x}{a}\right)^2r_p\right]dx.
\tag{12}
$$

The coefficients in each transformed integral are evaluated at that integral's
$q$. The Jacobians cancel the host light-line singularity; scaling by $a$
also makes the dominant decay envelope $e^{-x}$ independent of separation.
The full answer is $1+\Delta F_{\rm prop}+\Delta F_{\rm evan}$.

### 4.2 Branches, resonances and error gates

Write $\rho=\epsilon_2/\epsilon_1$. For a lossless dielectric, its light line
appears at

$$
u_b=\sqrt{1-\rho}\quad(0<\rho<1),
\qquad x_b=a\sqrt{\rho-1}\quad(\rho>1).
\tag{13}
$$

Splitting there helps resolve changes in the substrate root. Weak loss moves
these branch locations off the real axis; the integration intervals are split
around their projected centers and widths.

A $p$ surface-mode denominator in (6) gives
$q_{\rm sp}^2=\rho/(\rho+1)$ and hence

$$
t_{\rm sp}=\sqrt{-\frac{1}{\rho+1}},\qquad x_{\rm sp}=at_{\rm sp}.
\tag{14}
$$

For supported lossy-metal cases, intervals around
$\operatorname{Re}x_{\rm sp}$ span geometrically increasing multiples of
$|\operatorname{Im}x_{\rm sp}|$. This resolves narrow resonances and wings
better than a single large adaptive interval.

The implementation integrates to infinity, using an initial tail split
$x=32$ and a refinement split $x=64$. Their quadrature tolerance factors
are $0.05$ and $0.0125$ of the requested tolerances, respectively.
Per-interval absolute budgets share the scaled absolute tolerance.
For each orientation it checks the summed quadrature error estimate and

$$
\delta_{\rm refine}=
 |\Delta F_{\rm prop}^{(2)}-\Delta F_{\rm prop}^{(1)}|
 +|\Delta F_{\rm evan}^{(2)}-\Delta F_{\rm evan}^{(1)}|
 \leq \max(\mathrm{atol},\mathrm{rtol}|F|).
\tag{15}
$$

Failures, nonfinite results and negative total rates raise errors; rates are
not clipped into plausible values. The reported absolute signed tail integral
is a diagnostic, not a rigorous error bound. The exponential integrand is
evaluated as zero beyond $x=700$, outside supported floating-point exponential
resolution at the requested accuracy.

The current kernel rejects $a<10^{-6}$, absorbing hosts and strictly lossless
nonpositive-permittivity substrates. Real-axis poles and ENZ limits need a
separate limiting-absorption prescription. Numerical support does not imply
that local continuum material constants remain valid at molecular distances.

### 4.3 Independent analytic limits

For an ideal mirror, image-dipole interference gives

$$
F_\perp=1+3\frac{\sin a-a\cos a}{a^3},
\qquad
F_\parallel=1-\frac32\left(
 \frac{\sin a}{a}+\frac{\cos a}{a^2}-\frac{\sin a}{a^3}\right).
\tag{16}
$$

Series expansion avoids cancellation near the surface:

$$
F_\perp=2-\frac{a^2}{10}+\frac{a^4}{280}+O(a^6),
\qquad
F_\parallel=\frac{a^2}{5}-\frac{3a^4}{280}+O(a^6).
\tag{17}
$$

The vertical image reinforces the source while the horizontal image cancels
it at vanishing separation. Both enhancements tend to one far from the
interface. These different limits detect sign, phase and prefactor mistakes.

For a lossy half-space at small $a$, large $q$ dominates:
$r_s\to0$, $r_p\to(\epsilon_2-\epsilon_1)/(\epsilon_2+\epsilon_1)$.
Using $\int_0^\infty x^2e^{-x}dx=2$ in (11)-(12) yields

$$
\Delta F_\perp\sim\frac{3}{a^3}
 \operatorname{Im}\frac{\epsilon_2-\epsilon_1}{\epsilon_2+\epsilon_1},
\qquad
\Delta F_\parallel\sim\frac{3}{2a^3}
 \operatorname{Im}\frac{\epsilon_2-\epsilon_1}{\epsilon_2+\epsilon_1}.
\tag{18}
$$

This explains the leading inverse-cube dependence and two-to-one orientation
ratio. It is an asymptotic check, not a substitute for converged integration
or an assertion of material validity at arbitrarily small distances.

## 5. Dispersion, native data and layer thickness

For tabulated data between $\lambda_i$ and $\lambda_{i+1}$, the current
interpolator applies independently to $n$ and $\kappa$:

$$
f(\lambda)=f_i+
 \frac{\lambda-\lambda_i}{\lambda_{i+1}-\lambda_i}(f_{i+1}-f_i),
\qquad f\in\{n,\kappa\}.
\tag{19}
$$

Native database wavelengths in micrometers are converted to nm:
$\lambda_{\rm nm}=1000\lambda_{\rm um}$.
No extrapolation is allowed. Separate $n$ and $\kappa$ tables must both cover
the query, so their usable domain is their intersection. Missing extinction
data require an explicit, recorded lossless assumption; missing data are not
evidence that absorption vanishes.

The supported native Sellmeier records use $\lambda$ in **micrometers**:

$$
\begin{aligned}
\text{formula 1:}\quad
n^2(\lambda)=1+c_0+\sum_m
 \frac{B_m\lambda^2}{\lambda^2-C_m^2},\\
\text{formula 2:}\quad
n^2(\lambda)=1+c_0+\sum_m
 \frac{B_m\lambda^2}{\lambda^2-C_m}.
\end{aligned}
\tag{20}
$$

The two pole conventions must not be interchanged. Evaluation must remain
inside the source range, avoid poles and produce finite positive $n^2$.
Unsupported formulas are rejected, not approximated by these forms.

Each layer has its own $N_j(\lambda)$ and thickness $d_j$. For a plane-wave
component with conserved $k_\parallel$,

$$
k_{z,j}=\sqrt{N_j^2k_0^2-k_\parallel^2},
\qquad \delta_j=k_{z,j}d_j.
\tag{21}
$$

Thickness changes phase and attenuation, while dispersion changes propagation
and Fresnel coefficients. They therefore affect multilayer interference jointly,
not as independent scalar correction factors.

The implemented diagnostic calls the coherent transfer-matrix calculation
for both incident polarizations and averages their **powers**:

$$
R=\frac{R_s+R_p}{2},\qquad
T=\frac{T_s+T_p}{2},\qquad A=1-R-T.
\tag{22}
$$

For the enforced lossless exterior media, power transmission includes
the appropriate flux/admittance ratio, not just $|t|^2$.
This accounts for externally incident illumination. It is **not** the
collected power of a dipole inside the stack, whose source terms and
multiple-reflection feedback differ.

The catalog preserves original source bytes, a pinned source revision,
license/citation/conditions and SHA-256 integrity checks. Checksums identify
changes, not physical accuracy or trusted authorship. The official
refractiveindex.info database is CC0; that does not license unrelated
databases or publication figures/text. See the [data guide](SPECTRAL_EQE_AND_DATA.md)
for the pinned Ag/SiO2 example and provenance details.

## 6. From electromagnetic rates to EQE

### 6.1 Monochromatic transition

Define an intrinsic radiative rate in the homogeneous reference,
$\Gamma_{r,h}$, and a fixed intrinsic nonradiative rate $\Gamma_{nr}$:

$$
q_0=\frac{\Gamma_{r,h}}{\Gamma_{r,h}+\Gamma_{nr}},
\qquad
\Gamma_{nr}=\Gamma_{r,h}\frac{1-q_0}{q_0}
\quad(0<q_0\leq1).
\tag{23}
$$

The optical budgets must share that reference:

$$
F_t=\frac{P_{\rm EM}}{P_h},\qquad
F_e=\frac{P_{\rm collected}}{P_h},\qquad
0\leq F_e\leq F_t.
\tag{24}
$$

Then $\Gamma_{\rm EM}=F_t\Gamma_{r,h}$ and
$\Gamma_{\rm collected}=F_e\Gamma_{r,h}$. Under linear, unsaturated excitation,
the probability of a collected photon per emissive excitation is its collected
rate divided by its entire decay rate. Multiplying by charge-balance
$\eta_c$ and emissive-exciton formation $\eta_x$ gives

$$
\mathrm{EQE}_{\rm mono}(\lambda)=
 \eta_c\eta_x
 \frac{q_0F_e(\lambda)}
 {1-q_0+q_0F_t(\lambda)}.
\tag{25}
$$

The API requires all three efficiencies explicitly. Thickness and optical
index do not determine electrical charge balance, intrinsic molecular yield
or exciton formation; no universal spin factor is inserted.
For $q_0=0$ the implementation returns zero directly. For $q_0=1$,
the optical fraction is $F_e/F_t$ unless $F_t=0$, where steady-state yield
is undefined and raises an error.

Crucially, increasing source work can mean increasing metal dissipation.
A large Purcell/source-work factor does not establish a large EQE.
Collected air escape, substrate emission and detector-NA collection are
different choices of $F_e$.

### 6.2 One broadband emitter population

Let $S_0(\lambda)$ be the intrinsic **photon** spectral density:

$$
S_0(\lambda)\geq0,\qquad
\int S_0(\lambda)\,d\lambda=1.
\tag{26}
$$

For wavelengths in nm its units are $\mathrm{nm}^{-1}$.
The full transition distribution sets one population lifetime:

$$
\Gamma_{\rm EM}=\Gamma_{r,h}\int
 S_0(\lambda)F_t(\lambda)\,d\lambda.
\tag{27}
$$

Consequently,

$$
D=1-q_0+q_0\int S_0(\lambda)F_t(\lambda)\,d\lambda,
\qquad
\frac{d\,\mathrm{EQE}}{d\lambda}
 =\eta_c\eta_x\frac{q_0S_0(\lambda)F_e(\lambda)}{D},
\tag{28}
$$

$$
\mathrm{EQE}_{\rm broadband}
 =\int\frac{d\,\mathrm{EQE}}{d\lambda}\,d\lambda.
\tag{29}
$$

Using a separate denominator from (25) in each bin would assign independent
lifetimes to the bins of this same emitter and changes the physical model.
The implementation uses trapezoidal integration on the supplied grid:

$$
\int f(\lambda)d\lambda\ \simeq\
\sum_{i=0}^{L-2}
\frac{\lambda_{i+1}-\lambda_i}{2}(f_i+f_{i+1}).
\tag{30}
$$

If a measured input is spectral power density $I_\lambda$, first convert
to photon weights:

$$
S_0(\lambda)=
\frac{I_\lambda(\lambda)\lambda/(hc)}
 {\int I_\lambda(\lambda')\lambda'/(hc)\,d\lambda'}.
\tag{31}
$$

Use compatible units when evaluating photon energy. A consistently applied
constant length-conversion factor cancels upon normalization.
The API requires explicit normalization and rejects nonfinite arithmetic;
it does not silently repair spectra. A truncated normalized spectrum models
that truncated distribution, not missing spectral tails.

## 7. Orientation and spatial ensembles

### 7.1 Two-basis power composition

For planar isotropic media and azimuth-integrated collection, symmetry removes
cross terms in the integrated response. A dipole at polar angle $\alpha$
has vertical weight

$$
w_\perp=\frac{|p_z|^2}{|\mathbf p|^2}=\cos^2\alpha,
\qquad
F(w_\perp)=w_\perp F_\perp+(1-w_\perp)F_\parallel.
\tag{32}
$$

This applies separately to total and collected optical rates, provided the
collection geometry has the required azimuthal symmetry.
It does not in general apply to restricted-azimuth collection, anisotropic
media or a coherent phased source array.

Uniform solid-angle sampling is obtained with independent
$u,v\sim U(0,1)$:

$$
\mu=2u-1,\quad\phi=2\pi v,\quad
\hat{\mathbf e}=
 \left(\sqrt{1-\mu^2}\cos\phi,\,
       \sqrt{1-\mu^2}\sin\phi,\,\mu\right).
\tag{33}
$$

Sampling $\alpha$ uniformly instead would oversample the polar directions.
Since $\mu$ is uniform on $[-1,1]$,

$$
\langle w_\perp\rangle=\frac12\int_{-1}^1\mu^2d\mu=\frac13,
\qquad
F_{\rm iso}=\frac13F_\perp+\frac23F_\parallel.
\tag{34}
$$

### 7.2 Fixed dipoles versus rapid reorientation

Write $C=\eta_c\eta_x$. For equally excited fixed dipoles, each orientation
has its own lifetime; therefore average individual yields:

$$
\mathrm{EQE}_{\rm fixed}
 =\frac{C}{N}\sum_{i=1}^N
 \frac{q_0F_{e,i}}{1-q_0+q_0F_{t,i}}.
\tag{35}
$$

If reorientation is rapid relative to decay and samples the stipulated
orientation distribution, average the rates first:

$$
\mathrm{EQE}_{\rm reorient}
 =C\frac{q_0\langle F_e\rangle}
 {1-q_0+q_0\langle F_t\rangle}.
\tag{36}
$$

Averages and nonlinear ratios generally do not commute. For example,
with $q_0=C=1$, equally weighted vertical $(F_t,F_e)=(4,1)$
and horizontal $(1,0.5)$ budgets give

$$
\mathrm{EQE}_{\rm fixed}=\tfrac12(1/4+1/2)=0.375,
\qquad
\mathrm{EQE}_{\rm reorient}=\frac{0.75}{2.5}=0.3.
\tag{37}
$$

These are illustrative budgets, not a measured device.
Both implemented population models exclude correlations, superradiance,
exciton transfer, saturation and orientation-dependent intrinsic emission.
An in-plane ensemble of incoherent dipoles is not a plane wave: a coherent
sum has cross terms in $|\sum_i\mathbf E_i|^2$, whereas the uncorrelated
ensemble average sums powers.

### 7.3 Spatial distributions: formal extension, not implemented device EQE

For a normalized excitation distribution $g(\mathbf r,\hat{\mathbf e})$,
the fixed-emitter monochromatic extension would be

$$
\mathrm{EQE}_{\rm fixed}
 =C\int g(\mathbf r,\hat{\mathbf e})
 \frac{q_0F_e(\lambda;\mathbf r,\hat{\mathbf e})}
 {1-q_0+q_0F_t(\lambda;\mathbf r,\hat{\mathbf e})}
 \,d^3r\,d\Omega_e,
\quad
\int g\,d^3r\,d\Omega_e=1.
\tag{38}
$$

Here the intrinsic/electrical factors are assumed uniform; nonuniform factors
must stay inside the appropriate weighted integrals.
For broadband fixed emitters, construct the integrated denominator (28)
for each population/location before averaging. Equation (38) is a statement
of the needed model, not an implemented spatial-integration API.

In an infinite laterally uniform stack, the basis responses depend on depth
but not lateral source position. For uniform excitation and active density
$\rho_a$ in area $A$ and emitting thickness $h$,

$$
V=Ah,\qquad N_{\rm active}=\rho_aAh.
\tag{39}
$$

Area increases photons and injected carriers together; normalized EQE need
not change. Edges, contacts, lateral scattering, current crowding and heating
break this simplification. Active excitation density is not molecular density.

## 8. Planar CPU work and memory model

Let $L$ be wavelength count, $Z$ sampled depths, $C_d$ independent thickness
configurations, $J$ finite layers and $Q$ quadrature evaluations per orientation,
including refinements. Two reusable bases give

$$
N_{\rm basis}=2LZC_d,\qquad
N_{\rm layer}=2LZC_dQJ,\qquad
W_{\rm planar}=f_\ell N_{\rm layer}.
\tag{40}
$$

$f_\ell$ is an assumed FLOP count per layer evaluation, not measured recursion
throughput. A Cartesian sweep with $m_j$ thickness samples in layer $j$ has
$C_d=\prod_jm_j$, not $\sum_jm_j$.
Area and emitter count do not multiply basis solves under planar symmetry.
The current reference integrator is scalar adaptive SciPy CPU code.

For $b$ bytes per real number, the implemented planning allocation is

$$
\begin{aligned}
B_{\rm workspace}
 &=b\left[16(J+2)+2L(J+2)+8L\right],\\
B_{\rm retained}
 &=4bLZC_d,\\
B_{\rm directions}
 &=3bN\quad\text{if all sampled directions are retained},\\
W_{\rm composition}
 &\simeq20LN\quad\text{operations for one composition workload}.
\end{aligned}
\tag{41}
$$

The factors are explicit budgeting assumptions for streamed layer state,
complex-index tables, reusable spectra and retained rate outputs.
For repeated per-depth/configuration composition, that last workload must
be repeated accordingly. Streaming or analytic orientation distributions
can avoid retaining $N$ directions.
These byte counts exclude interpreter/allocator overhead and full fields;
they are **not measured process peak RAM**.

## 9. Full-pixel FDTD planning model

This section sizes a possible future uniform 3-D Yee-grid solver; no such
solver or GPU backend is delivered. The grid must independently resolve
internal wavelength, film thicknesses, metal skin depth and geometry, and
must pass space/time/PML/source-sampling convergence tests.

For physical dimensions $\ell_x,\ell_y,\ell_z$, cubic cell size $\Delta$
and $p$ PML cells on each face:

$$
n_i=\left\lceil\frac{\ell_i}{\Delta}\right\rceil+2p,\qquad
N_{\rm cell}=n_xn_yn_z.
\tag{42}
$$

With $A_r$ allocated real arrays per cell and $b$ bytes per real:

$$
B_{\rm raw}=N_{\rm cell}A_rb,\qquad
B_{\rm provisioned}=\left\lceil1.3B_{\rm raw}\right\rceil.
\tag{43}
$$

At least six arrays represent the electromagnetic components; the default
$A_r=24$ budgets additional dispersive/PML state rather than deriving it
for a specific algorithm. The 30% headroom is an assumption, not a guarantee.
DFT monitors, snapshots, geometry, MPI halos and device copies are excluded.

For nondispersive media with fastest speed $v_{\max}=c/n_{\min}$,
the cubic-grid Courant bound and chosen timestep are

$$
\Delta t\leq\frac{\Delta}{v_{\max}\sqrt3},
\qquad
\Delta t_{\rm model}=0.99\frac{n_{\min}\Delta}{c\sqrt3},
\qquad
N_t=\left\lceil\frac{T_{\rm sim}}{\Delta t_{\rm model}}\right\rceil.
\tag{44}
$$

The **minimum**, not maximum, nondispersive index sets the fastest-wave
bound. Dispersive auxiliary equations can impose stricter stability limits;
using (44) alone is not a dispersive stability proof.

For $S$ independent source runs, $f_c$ FLOPs and $b_c$ bytes transferred
per cell-step:

$$
U=N_{\rm cell}N_tS,\qquad W_{\rm FDTD}=f_cU,\qquad
B_{\rm traffic}=b_cU.
\tag{45}
$$

Sequential source runs multiply work, not resident grid memory. Concurrent
runs can multiply memory. Twenty sampled positions with three independent
axes give $S=60$, not one solve per physical emitter and not demonstrated
sampling convergence. A suitable broadband pulse can span multiple
wavelengths in one run; do not automatically multiply by planar $L$.

Given **user-supplied sustained** arithmetic throughput $R_F$ and memory
bandwidth $R_B$, the optional ideal roofline is

$$
t_{\rm ideal}\geq
\max\left(\frac{W_{\rm FDTD}}{R_F},
          \frac{B_{\rm traffic}}{R_B}\right).
\tag{46}
$$

The tool reports the maximum as an ideal lower-bound estimate, not predicted
elapsed time. Real performance includes inefficiency, synchronization and
transfers. A bound does not imply feasibility if the grid exceeds RAM/VRAM.
GFLOP/s and GB/s use $10^9$; GiB uses $2^{30}$ bytes.
No sustained GPU rate or GPU runtime was inferred from a device model name.

## 10. Recorded examples and their interpretation

### 10.1 Computed orientation reference

At $\lambda_0=550$ nm, $d=20$ nm, $n_h=1$ and synthetic
$\epsilon_2=-5+0.5i$, the implemented one-interface kernels give:

| Orientation | $P_{\rm EM}/P_h$ |
|---|---:|
| Horizontal plane | 2.07017658 |
| Vertical | 7.75799635 |
| 45-degree tilt | 4.91408646 |
| Analytic isotropic mixture | 3.96611650 |
| 10,000 sampled isotropic directions, seed 42 | 3.96015439 |

The sample has $\langle w_\perp\rangle=0.33228511$, close to $1/3$.
Finite sampling explains the difference; it is not solver error estimation.
These numbers measure source work, **not** EQE or collected photon power.

The recorded five-repeat sequential CPU medians were approximately
0.3857 ms for a matched homogeneous interface and 4.2007 ms for the synthetic
lossy reference, computing both bases per call after warmup.
They are machine/load-specific measurements, not multilayer, FDTD or GPU
speed predictions. No process peak-RAM measurement was performed.

### 10.2 Pixel planning scenarios

Assume $\Delta=10$ nm, physical simulated depth $1$ micrometer,
emitting depth $0.1$ micrometer, $p=10$, $A_r=24$, float32,
$T_{\rm sim}=200$ fs, $n_{\min}=1$, $S=60$,
$f_c=300$ and $b_c=192$ bytes. Double precision needs a separately revised
traffic assumption as well as larger allocated arrays.

| Pixel side (micrometers) | Active volume (cubic micrometers) | Provisioned grid GiB | Total assumed PFLOP |
|---|---:|---:|---:|
| 1 | 0.1 | 0.201 | 0.326 |
| 10 | 10 | 14.511 | 23.576 |
| 50 | 250 | 351.482 | 571.055 |
| 100 | 1,000 | 1,400.333 | 2,275.129 |

PFLOP here is a **work count** of $10^{15}$ operations, not PFLOP/s.
With illustrative active density $1000$ per cubic micrometer, (39) gives
100, 10,000, 250,000 and 1,000,000 active emitters, respectively.
These populations do not establish that 20 sampled positions suffice.

For $L=31$, $J=6$, $Q=2000$, $C_d=1$, the planar sheet has
62 basis integrals; $Z=20$ gives 1,240. Neither count grows with lateral area.
These are workload counts for future finite-stack basis evaluation, not
evidence that its source/escape kernel already exists.

The observed workstation has an Intel Core i9-13900HX, approximately
31.69 GiB system RAM and an RTX 4070 Laptop GPU reporting 8188 MiB VRAM.
Under this grid model the 10-micrometer case exceeds that VRAM but fits
nominal system capacity before other allocations; the two larger cases exceed
both. No FDTD run was performed to validate capacity or elapsed time.

## 11. Reproduction, validation and remaining derivation

Generate a new, non-overwriting report with assumptions and current CPU timing:

```powershell
python -m simulation_resources --output pixel-resource-study.json --benchmark-cpu
```

The reference milestone at code commit `939b689` passed 215 local tests and
all eight Windows/Linux Python 3.11-3.14 CI environments (1,720 executions).
The package build, installed-wheel smoke and runtime dependency audit passed.
See [CI evidence](https://github.com/gmongell/Mathematica_python_TMM_OLEDPixel/actions/runs/38975378348).
Passing tests establish the bounded comparisons below, not independent
experimental validation of full OLED device predictions.

| Mathematical requirement | Named tests |
|---|---|
| Homogeneous normalization and ideal-mirror limits | [reference-kernel tests](../tests/test_dipole_reference.py) |
| Native Sellmeier/index evaluation | `test_rii_native_tables_and_formula` in [spectral tests](../tests/test_spectral_eqe.py) |
| Stack dispersion/thickness and non-EQE diagnostic | `test_stack_dispersion_and_thickness`, `test_plane_wave_thickness_spectrum_is_not_eqe` in [spectral tests](../tests/test_spectral_eqe.py) |
| Equations (25) and (28) | `test_monochromatic_rate_equation_and_yield_boundaries`, `test_broadband_uses_global_rate_denominator` in [spectral tests](../tests/test_spectral_eqe.py) |
| Orientation statistics and noncommuting averages | `test_random_orientations_are_reproducible_and_isotropic`, `test_orientation_averages_powers_before_ratios` in [spectral tests](../tests/test_spectral_eqe.py) |
| Planar workload and population scaling | `test_planar_two_bases_reuse_random_ensemble`, `test_area_volume_population_scaling` in [resource tests](../tests/test_simulation_resources.py) |
| FDTD grid/CFL/memory/roofline | `test_fdtd_grid_cfl_memory_and_roofline_arithmetic`, `test_fdtd_independent_emitters_scale_work_not_serial_memory` in [resource tests](../tests/test_simulation_resources.py) |

The missing central physics is a stable finite-emitting-layer Green/source
calculation between both interfaces, including multiple reflection feedback,
external emitted flux, absorption and guided channels with independent energy
and convergence checks. Only after deriving and validating those budgets
can (25) or (28) produce automatic thickness/index-dependent stack EQE.
Incident plane-wave $T$ and single-interface total source work must not
be substituted for that missing collected-emission calculation.

## 12. Sources and rights

- [B. Auguie, planar decay-rate theory](https://nano-optics.ac.nz/planar/articles/decay_rates_theory.html):
  orientation-resolved planar Green-function formulas and integration changes.
- [Novotny and Hecht, Principles of Nano-Optics, dipole emission near planar interfaces](https://www.cambridge.org/core/books/principles-of-nanooptics/dipole-emission-near-planar-interfaces/059D85E0BA32D70B1B9FCB83A4ABB62A):
  broader theoretical context; citation does not imply all cases are implemented.
- [S. J. Byrnes, Multilayer optical calculations](https://arxiv.org/abs/1603.02720):
  Fresnel/transfer-matrix conventions and power accounting.
- [M. N. Polyanskiy, Scientific Data 11, 94 (2024)](https://doi.org/10.1038/s41597-023-02898-2)
  and the [native database](https://github.com/polyanskiy/refractiveindex.info-database):
  material-data format and provenance. The project example pins database commit
  `c5c2f188e848453def5970e347399d653df2ffc2`.
- [Williams, Waterman and Patterson, Roofline (2009)](https://doi.org/10.1145/1498765.1498785):
  compute/bandwidth performance-bound framework; the workload counts here
  remain project planning assumptions.

Existing project licensing applies; see [NOTICE](../NOTICE) and
[documentation license](../LICENSE-CC-BY-4.0.txt). Preserve attribution
and distinguish measured data from synthetic examples when reusing this record.
