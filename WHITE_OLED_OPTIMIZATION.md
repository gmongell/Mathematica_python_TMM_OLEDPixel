# White-spectrum optical optimization (research prototype)

The module `white_oled_proxy_optimizer.py` performs a reproducible **plane-wave** multilayer thickness optimization with two objectives:

- Maximize a white-spectrum-weighted, projected hemispherical transmittance proxy.
- Minimize the solid-angle-weighted coefficient of variation of spectral transmission versus viewing angle.

The spectrum is an illustrative Gaussian blue/green/red mixture, **not** measured emitter spectral power distribution, D65, or an experimentally validated white point. The Alq3 and silver refractive-index expressions are toy approximations. This model does **not** calculate dipole emission, local density of optical states, internal quantum efficiency, EQE, luminance, CIE 1931 chromaticity, CRI, or viewing-angle color shift. Do not report the objective as OLED efficiency.

## Run

```bash
python -m pip install numpy scipy matplotlib tmm
python -m unittest discover -s tests -v
python white_oled_proxy_optimizer.py --samples 30 --angles 20 --iterations 30 --penalty 0.3 --seed 7
```

The optimizer reports Alq3 and silver thicknesses in **nm**, the objective and transmission/uniformity metrics, and whether differential evolution reported convergence. A finite result is not a physically validated optimum. The output depends on the angle grid, wavelength grid, penalty, toy optical constants, and thickness bounds.

## Roadmap to white OLED efficiency and color uniformity

1. Ingest measured complex refractive-index spectra `n(λ)+ik(λ)` for every layer, plus emitter spectra and photoluminescence quantum yield.
2. Replace incident plane waves with stratified-media electric-dipole Green tensors; resolve horizontal/vertical dipoles, position, polarization, guided/plasmonic losses, and total radiated plus dissipated power.
3. Compute angle-dependent spectral radiant intensity, convert to photometric intensity with the luminous-efficiency function, and calculate CIE XYZ/xy and Δu'v' over viewing angles using tabulated observer data.
4. Optimize external quantum efficiency, luminous efficacy, on-axis white point, and angular color shift with constraints on thickness, optical absorption, and manufacturability; generate Pareto fronts.
5. Validate against independent transfer-matrix/Green-function solvers, limiting cases, energy conservation, and measured angular spectra. Distinguish simulations from experiments.

No proprietary Lumerical/Mathematica solver or experimental dataset was executed to validate this new module.
