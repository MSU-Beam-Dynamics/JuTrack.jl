"""Piwinski Touschek loss rates from JuTrack optics, with explicit SI inputs.

Adapted from Accelerator Toolbox, Copyright 2024 the Accelerator Toolbox
Developers, under Apache-2.0 (licenses/AT-Apache-2.0.txt).
Upstream: atcollab/at, commit d3345747335a2a05bcf0904769ba38fb93396fec,
pyat/at/acceptance/touschek.py. Modified for JuTrack: explicit optical arrays,
complete spatial weights, stable coefficients/Bessel evaluation, and diagnostics.
Physics reference: A. Piwinski, DESY 98-179, arXiv:physics/9903034, Eqs. 31-42.

This is the uncoupled Gaussian-beam model including dispersion in both planes.
It does not compute optics, momentum acceptance, IBS or bunch equilibrium.
"""

import math

import numpy as np
from scipy.constants import c, e, physical_constants
from scipy.integrate import quad
from scipy.special import i0e


_REST_ENERGY_EV = physical_constants["electron mass energy equivalent in MeV"][0] * 1e6
_CLASSICAL_RADIUS_M = physical_constants["classical electron radius"][0]
_TAIL_EXPONENT = 80.0

__all__ = ["touschek_from_optics"]


def _positive(value, name, *, allow_zero=False):
    value = float(value)
    if not math.isfinite(value) or value < 0 or (value == 0 and not allow_zero):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return value


def _array(value, shape, name):
    value = np.asarray(value, dtype=float)
    if value.shape != shape or not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be a finite array of shape {shape}")
    return value


def _coefficients(beta, alpha, dispersion, emits, sigp, beta_rel2, gamma):
    """Stable eigenvalues of Piwinski's positive 2x2 angular-spread matrix."""
    sigma_beta2 = beta * emits
    d = dispersion[[0, 2]]
    dt = alpha * d + beta * dispersion[[1, 3]]
    u = 1.0 / sigp**2 + np.sum(d**2 / sigma_beta2)
    v, w = dt**2 / sigma_beta2
    denominator = u + v + w
    sigma_h2 = 1.0 / denominator
    ax, ay = beta**2 / sigma_beta2
    aa = ax * (u + w) / denominator
    bb = ay * (u + v) / denominator
    cc = sigma_h2 * np.prod(beta * dt / sigma_beta2)
    largest = 0.5 * (aa + bb + math.hypot(aa - bb, 2.0 * cc))
    determinant = ax * ay * u / denominator
    smallest = determinant / largest
    scale = beta_rel2 * gamma**2
    bmin = smallest / scale
    b2 = max(0.0, (largest - smallest) / (2.0 * scale))
    if not all(math.isfinite(x) and x > 0 for x in (sigma_h2, bmin)):
        raise ValueError("Optics/beam inputs do not give finite positive Piwinski coefficients")
    if not math.isfinite(b2):
        raise ValueError("Optics/beam inputs overflow the Piwinski coefficients")
    return sigma_beta2, sigma_h2, bmin, b2


def _integral(tm, bmin, b2, rtol):
    """Eq. 31 integral with tau=tm*exp(u); no Bessel exponential overflow."""
    # Beyond exp(-745) the local rate underflows double precision.
    if bmin * tm > 745.0:
        return 0.0, 0.0
    upper = max(tm * 1.01, _TAIL_EXPONENT / bmin)
    end = math.log(upper / tm)

    def integrand(u):
        tau = tm * math.exp(u)
        log_ratio = u - math.log1p(tau)
        # tau times the original Piwinski bracket; expm1 avoids subtracting
        # nearly equal terms near the scattering threshold.
        bracket = (
            (4.0 * tau + 4.0 + 1.0 / tau) * math.expm1(log_ratio)
            - tau * math.expm1(-0.5 * log_ratio)
            - (2.0 + 0.5 / tau) * log_ratio
        )
        return (
            bracket * math.sqrt(tau / (1.0 + tau))
            * math.exp(-bmin * tau) * float(i0e(b2 * tau))
        )

    transition = math.log(1.0 / (bmin * tm))
    points = [transition] if 0 < transition < end else None
    output = quad(
        integrand, 0.0, end, epsrel=rtol, epsabs=1e-200,
        points=points, limit=250, full_output=1,
    )
    value, error = output[:2]
    if len(output) > 3:
        raise RuntimeError(f"Piwinski quadrature did not converge: {output[3]}")
    if not math.isfinite(value) or not math.isfinite(error) or value < 0:
        raise RuntimeError("Piwinski quadrature returned a nonphysical integral")
    if error > max(1e-200, 5.0 * rtol * abs(value)):
        raise RuntimeError("Piwinski quadrature error exceeds the requested tolerance")
    return value, error


def touschek_from_optics(
    *,
    beta_m,
    alpha,
    dispersion,
    weights_m,
    circumference_m,
    energy_eV,
    emittances_m,
    sigma_z_m,
    sigma_delta,
    momentum_acceptance,
    bunch_charge_C=None,
    bunch_current_A=None,
    repetitions=1,
    rtol=1e-7,
):
    """Calculate lifetime in seconds from prescribed beam parameters and optics.

    Parameters
    ----------
    beta_m, alpha : array, shape (n, 2)
        Horizontal/vertical beta [m] and alpha at the same sample locations.
    dispersion : array, shape (n, 4)
        Dx [m], Dpx, Dy [m], Dpy, differentiated with respect to delta=p/p0-1.
    weights_m : array, shape (n,)
        Positive physical integration lengths represented by each sample. Use
        element lengths for entrance-point quadrature, or nodal quadrature
        weights for a denser profile. These are NOT arbitrary optimizer weights.
    circumference_m : float
        Full ring circumference. Must equal sum(weights_m)*repetitions.
    energy_eV : float
        Total electron energy, not kinetic energy.
    emittances_m : pair
        Prescribed geometric horizontal/vertical emittances [m rad].
    sigma_z_m, sigma_delta : float
        RMS bunch length [m] and relative RMS momentum spread.
    momentum_acceptance : array, shape (n, 2), or pair
        Signed local (negative, positive) acceptance at these same locations,
        e.g. [-0.02, 0.015], NOT percent. A pair explicitly requests a uniform
        acceptance hypothesis. No RF clipping, interpolation or tracking is
        performed. Tracking scan limits are not certified loss boundaries.
    bunch_charge_C, bunch_current_A : float
        Supply exactly one. Current is PER BUNCH, not total ring current.
    repetitions : int
        Explicit identical-cell repetition count. Both optics AND acceptance
        must repeat. Do not use a magnetic period to infer RF-on periodicity.
    rtol : float
        Relative quadrature tolerance (does not control spatial sampling).

    Returns
    -------
    dict
        Python scalars/lists, local fractional loss rates [1/s], separate
        signed-side rates, and spatial contributions to the ring-averaged rate.
        lifetime_s is the inverse initial fractional loss rate, not a constant
        exponential lifetime: with fixed beam shape it is the half-charge time.
        particle_loss_rate_per_s counts lost particles, not collision events.
        Zero charge returns infinite lifetime (not a finite JSON number).

    Notes
    -----
    Assumes an uncoupled Gaussian beam with prescribed beam sizes. Finite
    vertical dispersion is included; general coupled beam covariance is not.
    No equilibrium, IBS, impedance, gas lifetime or acceptance is inferred.
    Uses Piwinski's small-relative-momentum acceptance convention.
    """
    circumference = _positive(circumference_m, "circumference_m")
    energy = _positive(energy_eV, "energy_eV")
    sigz = _positive(sigma_z_m, "sigma_z_m")
    sigp = _positive(sigma_delta, "sigma_delta")
    rtol = _positive(rtol, "rtol")
    if not 50 * np.finfo(float).eps <= rtol < 1:
        raise ValueError("rtol must be at least 50*machine_epsilon and less than 1")
    if isinstance(repetitions, (bool, np.bool_)) or not isinstance(repetitions, (int, np.integer)):
        raise ValueError("repetitions must be a positive integer")
    if repetitions < 1:
        raise ValueError("repetitions must be a positive integer")
    gamma = energy / _REST_ENERGY_EV
    if gamma <= 1:
        raise ValueError("energy_eV must exceed the electron rest energy")
    beta_rel2 = 1.0 - 1.0 / gamma**2

    weights = np.asarray(weights_m, dtype=float)
    if weights.ndim != 1 or not len(weights) or not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("weights_m must be a nonempty vector of finite positive lengths")
    if not math.isclose(float(weights.sum()) * repetitions, circumference, rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError("Spatial coverage: sum(weights_m)*repetitions must equal circumference_m")
    n = len(weights)
    beta = _array(beta_m, (n, 2), "beta_m")
    alphas = _array(alpha, (n, 2), "alpha")
    dispersions = _array(dispersion, (n, 4), "dispersion")
    emits = _array(emittances_m, (2,), "emittances_m")
    if np.any(beta <= 0) or np.any(emits <= 0):
        raise ValueError("beta_m and emittances_m must be positive")
    ma = np.asarray(momentum_acceptance, dtype=float)
    uniform = ma.shape == (2,)
    if uniform:
        ma = np.broadcast_to(ma, (n, 2)).copy()
    ma = _array(ma, (n, 2), "momentum_acceptance")
    if np.any(ma[:, 0] >= 0) or np.any(ma[:, 1] <= 0):
        raise ValueError("momentum_acceptance columns must satisfy negative < 0 < positive")
    if np.any(ma[:, 0] <= -1):
        raise ValueError("Negative momentum acceptance must exceed -1")
    if (bunch_charge_C is None) == (bunch_current_A is None):
        raise ValueError("Supply exactly one of bunch_charge_C or bunch_current_A (per bunch)")
    frev = math.sqrt(beta_rel2) * c / circumference
    if bunch_charge_C is not None:
        charge = _positive(bunch_charge_C, "bunch_charge_C", allow_zero=True)
    else:
        charge = _positive(bunch_current_A, "bunch_current_A", allow_zero=True) / frev
    particles = charge / e

    rates = np.zeros((n, 2))
    errors = np.zeros((n, 2))
    if particles > 0:
        for i in range(n):
            sigma_beta2, sigma_h2, bmin, b2 = _coefficients(
                beta[i], alphas[i], dispersions[i], emits, sigp, beta_rel2, gamma,
            )
            # Eq.31 is algebraically equivalent to AT's prefactor without the
            # ill-conditioned subtraction B1**2-B2**2 in very flat beams.
            prefactor = (
                _CLASSICAL_RADIUS_M**2 * c * particles * np.prod(beta[i])
                * math.sqrt(sigma_h2)
                / (8.0 * math.sqrt(math.pi) * beta_rel2 * gamma**4
                   * np.prod(sigma_beta2) * sigz * sigp)
            )
            for side in range(2):
                if side == 1 and ma[i, 0] == -ma[i, 1]:
                    rates[i, 1], errors[i, 1] = rates[i, 0], errors[i, 0]
                    continue
                try:
                    integral, error = _integral(beta_rel2 * ma[i, side]**2, bmin, b2, rtol)
                except (RuntimeError, OverflowError) as exc:
                    raise RuntimeError(f"Touschek sample {i}, side {side}: {exc}") from exc
                rates[i, side] = prefactor * integral
                errors[i, side] = prefactor * error
    if not np.all(np.isfinite(rates)):
        raise RuntimeError("Nonfinite Touschek loss rate; check beam units and optics")

    local = rates.mean(axis=1)
    normalized_weights = weights * repetitions / circumference
    contributions = local * normalized_weights
    mean_rate = float(contributions.sum())
    quadrature_error = float(np.dot(errors.mean(axis=1), normalized_weights))
    return {
        "lifetime_s": 1.0 / mean_rate if mean_rate > 0 else math.inf,
        "loss_rate_per_s": mean_rate,
        "local_loss_rate_per_s": local.tolist(),
        "side_loss_rates_per_s": rates.tolist(),
        "local_particle_loss_rate_per_s": (local * particles).tolist(),
        "rate_contributions_per_s": contributions.tolist(),
        "particle_loss_rate_per_s": mean_rate * particles,
        "quadrature_error_per_s": quadrature_error,
        "weights_m": weights.tolist(),
        "circumference_m": circumference,
        "repetitions": int(repetitions),
        "bunch_particles": particles,
        "bunch_charge_C": charge,
        "bunch_current_A": charge * frev,
        "energy_eV": energy,
        "emittances_m": emits.tolist(),
        "sigma_z_m": sigz,
        "sigma_delta": sigp,
        "momentum_acceptance": ma.tolist(),
        "acceptance_source": "uniform_assumption" if uniform else "supplied_local_values",
        "method": "Piwinski Eqs.31-42, Gaussian uncoupled beam with both-plane dispersion",
        "numerics": {
            "rtol": rtol,
            "tail_exponent_cutoff": _TAIL_EXPONENT,
            "quadrature_error_excludes_spatial_and_acceptance_errors": True,
        },
    }
