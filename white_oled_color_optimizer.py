"""White-light proxy with angular chromaticity optimization.

The model multiplies an illustrative incident spectrum by multilayer plane-wave
transmission. This is NOT an OLED dipole-emission spectrum or external quantum
efficiency. Supply wavelength-matched tabulated CIE CMFs for chromaticity.
"""
import argparse
import json
import numpy as np
from scipy.optimize import differential_evolution
from tmm_OLED9 import calculate_transmission
from white_oled_proxy_optimizer import reference_spectrum, evaluate_proxy
from white_oled_colorimetry import read_cmf_csv, calculate_angular_color


def evaluate_design(alq3_nm, silver_nm, wavelengths_nm, angles_deg, cmf):
    spectrum = reference_spectrum(wavelengths_nm)
    transmission = calculate_transmission(wavelengths_nm, angles_deg, alq3_nm, silver_nm)
    angle_spd = spectrum[:, None] * transmission
    color = calculate_angular_color(wavelengths_nm, angles_deg, angle_spd, cmf)
    proxy = evaluate_proxy(alq3_nm, silver_nm, wavelengths_nm, angles_deg, spectrum)
    return {"proxy": proxy, "color": color}


def optimize_color_uniformity(wavelengths_nm, angles_deg, cmf, angular_penalty=0.2,
                              color_penalty=10.0, maxiter=25, seed=7):
    if not np.isfinite(angular_penalty) or not np.isfinite(color_penalty) or angular_penalty < 0 or color_penalty < 0:
        raise ValueError("Penalties must be nonnegative")
    def objective(x):
        metrics = evaluate_design(x[0], x[1], wavelengths_nm, angles_deg, cmf)
        return (-metrics["proxy"]["hemispherical_transmission_proxy"]
                + angular_penalty * metrics["proxy"]["angular_coefficient_of_variation"]
                + color_penalty * metrics["color"]["max_delta_uv"])
    result = differential_evolution(objective, [(25., 150.), (5., 60.)],
                                    seed=seed, maxiter=maxiter, popsize=6, polish=True)
    metrics = evaluate_design(*result.x, wavelengths_nm, angles_deg, cmf)
    return {"alq3_nm": float(result.x[0]), "silver_nm": float(result.x[1]),
            "objective": float(result.fun), "optimizer_success": bool(result.success),
            "metrics": metrics,
            "warning": "Plane-wave spectral filtering proxy only, not physical OLED efficiency."}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cmf", required=True, help="CSV with wavelength_nm,xbar,ybar,zbar")
    p.add_argument("--iterations", type=int, default=25)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--angular-penalty", type=float, default=0.2)
    p.add_argument("--color-penalty", type=float, default=10.)
    args = p.parse_args()
    if args.iterations < 1:
        p.error("--iterations must be >= 1")
    w, cmf = read_cmf_csv(args.cmf)
    if w[0] < 400 or w[-1] > 700:
        p.error("This illustrative spectral model expects 400–700 nm")
    angles = np.linspace(0, 85, 18)
    result = optimize_color_uniformity(w, angles, cmf, args.angular_penalty,
                                       args.color_penalty, args.iterations, args.seed)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
