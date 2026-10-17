# Dipole reference kernels: conventions and derivation

For rendered LaTeX equations and the connected EQE/ensemble/resource derivations,
see the [consolidated mathematical record](MATHEMATICAL_DERIVATIONS.md).

Research milestone, 2026-10-16. Implemented in `dipole_reference.py`.
These are homogeneous and **one-interface** reference kernels, not a finite
emitting-layer/multilayer OLED solver, material database, or sale-ready product.

## Geometry and normalization

The interface is z=0. A prescribed harmonic electric dipole is at z=d>0 in
a lossless isotropic nonmagnetic host with real n1>0, epsilon1=n1^2.
The substrate is the z<0 half-space, with relative permittivity epsilon2.
Use exp(-i omega t), SI dipole amplitude p (C m), vacuum wavelength lambda0
(nm at the API), and k1=2 pi n1/lambda0 (inverse nm for phase calculations).
The substrate is passive: Im(epsilon2)>=0. The API also accepts the explicit
ideal-boundary sentinel `"pec"` for a perfect electric conductor.

Homogeneous macroscopic radiated power, for peak (not RMS) dipole amplitude:

$$
P_{\operatorname{host}}=\frac{n_1\omega^4|\mathbf p|^2}{12\pi\epsilon_0c^3}
$$

No microscopic local-field correction is included. The normalized angular
density is dP/(P_host dOmega), with dOmega=sin(theta) dtheta dphi.
For a unit dipole axis e and observation direction r_hat:

$$
\frac{1}{P_{\operatorname{host}}}\frac{dP}{d\Omega}=\frac{3}{8\pi}\left(1-|\hat{\mathbf e}\cdot\hat{\mathbf r}|^2\right)
$$

Project e onto e_phi (s) and e_theta (p) before squaring. In particular:

| Dipole | s density / [3/(8 pi)] | p density / [3/(8 pi)] |
|---|---|---|
| vertical z | 0 | sin(theta)^2 |
| horizontal x | sin(phi)^2 | cos(theta)^2 cos(phi)^2 |
| horizontal y | cos(phi)^2 | cos(theta)^2 sin(phi)^2 |
| random in-plane | 1/2 | cos(theta)^2/2 |
| isotropic ensemble | 1/3 | 1/3 |

Each total density integrates to 1 over the sphere and 1/2 in each hemisphere.
An isotropic ensemble averages **powers**, not amplitudes: vertical weight
1/3 plus random-in-plane weight 2/3.

## Reflected Green tensor and source work

Let q=k_parallel/k1, w=sqrt(1-q^2), and
v=sqrt(epsilon2/epsilon1-q^2). Select the outgoing/decaying branch:
positive real root for propagation, positive imaginary root for evanescence.
The p coefficient here uses the standard electric-field Fresnel convention:

$$
r_s=\frac{w-v}{w+v},\qquad r_p=\frac{\epsilon_2w-\epsilon_1v}{\epsilon_2w+\epsilon_1v}
$$

For `"pec"`, r_s=-1 and r_p=+1. Identical media give both coefficients zero.
The total electromagnetic source-work enhancement relative to P_host is:

$$
\begin{aligned}
F_\perp&=1+\frac32\operatorname{Re}\int_0^\infty\frac{q^3}{w}\,r_p\,e^{2ik_1dw}\,dq\\
F_\parallel&=1+\frac34\operatorname{Re}\int_0^\infty\frac{q}{w}\,(r_s-w^2r_p)\,e^{2ik_1dw}\,dq
\end{aligned}
$$

The **minus** before w^2*r_p in the parallel expression is essential. It
comes from the reflected p dyad, not from a polarization-intensity average.
These expressions follow the reflected dyadic Green-function contribution
to source work (Im[p* . E_reflected]) with homogeneous normalization.
The in-plane orientations are equivalent after azimuthal integration.

These rates exclude intrinsic emitter nonradiative decay. For an absorbing
substrate they include work dissipated in that substrate. They are NOT
far-field escape fractions, photon yield, or external quantum efficiency.
In particular, a host-evanescent component can propagate in a higher-index
lossless substrate; "evanescent correction" does not mean "absorbed power".

## Integration without the host light-line singularity

Write a=2 k1 d. In the q<1 region substitute u=sqrt(1-q^2):

$$
\begin{aligned}
\Delta F_{\perp,\operatorname{prop}}&=\frac32\int_0^1(1-u^2)\operatorname{Re}\left[r_pe^{iau}\right]du\\
\Delta F_{\parallel,\operatorname{prop}}&=\frac34\int_0^1\operatorname{Re}\left[(r_s-u^2r_p)e^{iau}\right]du
\end{aligned}
$$

In q>1 substitute t=sqrt(q^2-1), then x=a*t:

$$
\begin{aligned}
\Delta F_{\perp,\operatorname{evan}}&=\frac{3}{2a}\int_0^\infty\left(1+\frac{x^2}{a^2}\right)e^{-x}\operatorname{Im}r_p\,dx\\
\Delta F_{\parallel,\operatorname{evan}}&=\frac{3}{4a}\int_0^\infty e^{-x}\operatorname{Im}\left[r_s+\frac{x^2}{a^2}r_p\right]dx
\end{aligned}
$$

This removes the integrable q=1 singularity and scales the decaying tail.
For a lossless dielectric, split the substrate light line explicitly:
u=sqrt(1-epsilon2/epsilon1) or x=a*sqrt(epsilon2/epsilon1-1).
For a lossy metal, also split near the complex surface-mode location
t=sqrt(-1/(epsilon2/epsilon1+1)), including its width. Strictly lossless
nonpositive-epsilon substrates are rejected: real-axis poles and ENZ limits
need an explicit limiting-absorption/pole prescription not supplied here.

Weakly absorbing dielectric substrate light lines are also split near their
complex branch location. Surface-mode splits span geometrically increasing
width offsets so narrow resonances and their wings are not left in one large
integration interval. Extremely narrow poles can still fail the accuracy gate.

The full infinite evanescent tail is integrated, not truncated at a chosen
wavevector. Split it initially at x=32, then repeat with x=64 and four-times
tighter quadrature tolerances. Return estimated quadrature error, difference
between refinements, the absolute signed tail contribution, and evaluation count.
The tail diagnostic is not a rigorous error bound. Reject QUADPACK failures,
nonfinite results, negative total rates, or errors/refinement changes above
max(atol, rtol*|F|); never clip a failed rate into a plausible answer.
Very small dimensionless separation a<1e-6 is explicitly unsupported.
Even above that numerical boundary, continuum optical constants may be invalid
at molecular distances. Quadrature estimates are not a proof of convergence
for arbitrarily sharp resonances; independent benchmarks remain necessary.

The raw q-space kernel intentionally rejects q=1 (the singular point).
Use the transformed integrated API, not sampling at the light line.

## Analytic benchmark: ideal mirror

Method of images gives, with a=2 k1 d:

$$
\begin{aligned}
F_\perp&=1+\frac{3(\sin a-a\cos a)}{a^3}\\
F_\parallel&=1-\frac32\left(\frac{\sin a}{a}+\frac{\cos a}{a^2}-\frac{\sin a}{a^3}\right)
\end{aligned}
$$

At small a, the numerically stable limits are
F_vertical=2-a^2/10+a^4/280+... and
F_parallel=a^2/5-3*a^4/280+....
At large separation both tend to 1. These independently derived image-dipole
results detect wrong signs, missing prefactors and phase/unit mistakes.

For a lossy half-space, the quasistatic leading evanescent corrections are:

$$
\Delta F_\perp\sim\frac{3}{a^3}\operatorname{Im}\frac{\epsilon_2-\epsilon_1}{\epsilon_2+\epsilon_1},\qquad \Delta F_\parallel\sim\frac{3}{2a^3}\operatorname{Im}\frac{\epsilon_2-\epsilon_1}{\epsilon_2+\epsilon_1}
$$

This is a mathematical asymptote, not a material/device-validity claim.

## Sources and validation limits

- [Auguie, planar: decay-rate theory](https://nano-optics.ac.nz/planar/articles/decay_rates_theory.html),
  accessed 2026-10-16. Explicit orientation formulas and variable substitutions;
  the implementation here is newly written from the equations, not copied code.
- [Novotny and Hecht, Principles of Nano-Optics, chapter 10](https://www.cambridge.org/core/books/principles-of-nanooptics/dipole-emission-near-planar-interfaces/059D85E0BA32D70B1B9FCB83A4ABB62A).
  Broader planar-interface Green-function framework. This citation is not an
  assertion that all chapter cases are implemented or independently validated.
- [Byrnes, Multilayer optical calculations](https://arxiv.org/abs/1603.02720).
  Background on passive-media conventions and Fresnel power accounting.

Tests compare mirror image-dipole limits, full-sphere homogeneous normalization,
an alternative theta/hyperbolic integration for synthetic dielectric/absorbing
interfaces, dielectric branch cases and a quasistatic asymptote. These are
bounded reference checks, not independent experimental validation. No measured
material constants or redistributed external reference data are supplied.
Finite stacks, emitter-layer feedback, layer-resolved absorption, guided-mode
accounting and air/substrate escape still require subsequent implementation.
