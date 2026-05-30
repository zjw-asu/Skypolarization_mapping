"""
Sky Polarization Simulator — Python conversion of MATLAB codebase.

Computes full-sky Rayleigh scattering polarization maps using the Mueller
matrix (Stokes) formalism for unpolarized sunlight scattered by small
atmospheric aerosol particles.

Functions
---------
compute_sun_position(observer, time_info) -> (zenith, azimuth, elevation)
compute_sigma_s(th_sun, ph_sun, TH_S, PHI_S, mu) -> sigma_S
compute_rayleigh_sky_map(sun_zenith, sun_azimuth, params) -> dict
sky_polarization(observer, time_info, resolution=1.5, wavelength=550e-9) -> dict
"""

import numpy as np


def compute_sun_position(observer, time_info):
    """Compute the sun's position in the sky.

    Parameters
    ----------
    observer : dict with keys 'latitude', 'longitude', 'gmtOffset'
    time_info : dict with keys 'year', 'month', 'day', 'hour', and optional 'dst'

    Returns
    -------
    zenith, azimuth, elevation : float (degrees)
    """
    lat = observer['latitude']
    lon = observer['longitude']
    gmt_offset = observer['gmtOffset']

    year = time_info['year']
    month = time_info['month']
    day = time_info['day']
    hour = time_info['hour']
    dst = time_info.get('dst', 0)

    # Day of year
    days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    is_leap = (year % 4 == 0 and year % 100 != 0) or (year % 400 == 0)
    if is_leap:
        days_per_month[1] = 29

    if month == 1:
        day_of_year = day
    else:
        day_of_year = day + sum(days_per_month[:month - 1])
    n = day_of_year

    # Solar declination (degrees)
    delta = 23.45 * np.sin(np.radians(360 / 365 * (284 + n)))

    # Equation of Time (hours)
    B = 360 / 365 * (n - 81)
    B_rad = np.radians(B)
    EOT = 0.165 * np.sin(2 * B_rad) - 0.126 * np.cos(B_rad) - 0.025 * np.sin(B_rad)

    # Local Solar Time (hours)
    standard_meridian = 15 * gmt_offset
    local_solar_time = hour - (1 / 15) * (standard_meridian - lon) + EOT - dst

    # Solar hour angle (degrees)
    hour_angle = 15 * (local_solar_time - 12)

    # Sun elevation and zenith
    sin_elevation = (np.sin(np.radians(lat)) * np.sin(np.radians(delta)) +
                     np.cos(np.radians(lat)) * np.cos(np.radians(delta)) *
                     np.cos(np.radians(hour_angle)))
    sin_elevation = np.clip(sin_elevation, -1, 1)

    elevation = np.degrees(np.arcsin(sin_elevation))
    zenith = 90 - elevation

    # Sun azimuth (from North, clockwise)
    sun_east = -np.cos(np.radians(delta)) * np.sin(np.radians(hour_angle))
    sun_north = (np.sin(np.radians(delta)) * np.cos(np.radians(lat)) -
                 np.cos(np.radians(delta)) * np.cos(np.radians(hour_angle)) *
                 np.sin(np.radians(lat)))

    azimuth = np.degrees(np.arctan2(sun_east, sun_north)) % 360

    return zenith, azimuth, elevation


def compute_sigma_s(th_sun, ph_sun, TH_S, PHI_S, mu):
    """Rotation angle from scattering plane to local meridian plane.

    Parameters
    ----------
    th_sun : float — Sun zenith angle (degrees)
    ph_sun : float — Sun azimuth angle (degrees)
    TH_S   : ndarray — Sky-point zenith angles (degrees)
    PHI_S  : ndarray — Sky-point azimuth angles (degrees)
    mu     : ndarray — Scattering angles (degrees)

    Returns
    -------
    sigma_S : ndarray — Rotation angles (degrees)
    """
    ph_sun_mod = ph_sun
    if ph_sun > 180:
        ph_sun_mod = ph_sun - 180

    del_phi = PHI_S - ph_sun_mod

    sigma_S = np.zeros_like(TH_S, dtype=float)

    # Identify special cases
    mu_singular = (np.abs(mu) < 1e-10) | (np.abs(mu - 180) < 1e-10)
    thS_zero = TH_S < 1e-10
    thSun_zero = th_sun < 1e-10
    general = ~mu_singular & ~thS_zero & ~thSun_zero

    # Argument for acosd in general formula
    denom = np.sin(np.radians(mu)) * np.sin(np.radians(TH_S))
    denom[denom == 0] = np.nan

    arg = (np.cos(np.radians(th_sun)) -
           np.cos(np.radians(TH_S)) * np.cos(np.radians(mu))) / denom
    arg = np.clip(arg, -1, 1)

    # General case A: del_phi in [0, 180]
    maskA = general & (del_phi >= 0) & (del_phi <= 180)
    if np.any(maskA):
        sigma_S[maskA] = 180 + np.degrees(np.arccos(arg[maskA]))

    # General case B: del_phi outside [0, 180]
    maskB = general & ~maskA
    if np.any(maskB):
        sigma_S[maskB] = 360 - np.degrees(np.arccos(arg[maskB]))

    # Special case: sky point at zenith (TH_S = 0)
    mask_thS0 = ~mu_singular & thS_zero & ~thSun_zero
    if np.any(mask_thS0):
        # Case A
        subA = mask_thS0 & (del_phi >= 0) & (del_phi <= 180)
        if np.any(subA):
            sigma_S[subA] = np.degrees(np.arccos(
                np.cos(np.radians(th_sun)) *
                np.cos(np.radians(PHI_S[subA] - ph_sun_mod))))

        # Case B
        subB = mask_thS0 & ~subA
        if np.any(subB):
            phi_s_vals = PHI_S[subB]
            wrapped = np.zeros_like(phi_s_vals)
            gt = phi_s_vals > ph_sun_mod
            wrapped[gt] = ph_sun_mod + (360 - phi_s_vals[gt])
            wrapped[~gt] = ph_sun_mod - phi_s_vals[~gt]
            sigma_S[subB] = np.degrees(np.arccos(
                np.cos(np.radians(th_sun)) * np.cos(np.radians(wrapped))))

    # Special case: sun at zenith
    sigma_S[~mu_singular & ~thS_zero & thSun_zero] = 0

    # Both at zenith
    sigma_S[thS_zero & thSun_zero] = 0

    # Singular scattering angles
    sigma_S[mu_singular] = 0

    return sigma_S


def compute_rayleigh_sky_map(sun_zenith, sun_azimuth, params):
    """Compute sky polarization via single Rayleigh scattering.

    Parameters
    ----------
    sun_zenith  : float — Sun zenith angle (degrees)
    sun_azimuth : float — Sun azimuth angle (degrees)
    params : dict with keys 'wavelength', 'particleRadius', 'refractiveIndex', 'resolution'

    Returns
    -------
    sky_map : dict with keys 'azimuth', 'zenith', 'DoLP', 'DoCP', 'AoP', 'Q', 'U', 'I', 'mu'
    """
    wavelength = params['wavelength']
    r = params['particleRadius']
    m = params['refractiveIndex']
    resolution = params['resolution']

    th_sun = sun_zenith
    ph_sun = sun_azimuth

    # Build sky grid
    zenith = np.linspace(0, 90, round(90 / resolution) + 1)
    azimuth = np.linspace(0, 360, round(360 / resolution) + 1)

    PHI_S, TH_S = np.meshgrid(azimuth, zenith)

    # Scattering angle mu at each sky point
    mu = np.degrees(np.arccos(np.clip(
        np.sin(np.radians(th_sun)) * np.sin(np.radians(TH_S)) *
        np.cos(np.radians(PHI_S - ph_sun)) +
        np.cos(np.radians(th_sun)) * np.cos(np.radians(TH_S)),
        -1, 1)))

    # Rayleigh Mueller matrix elements
    cosmu = np.cos(np.radians(mu))
    F11 = (1 + cosmu ** 2) / 2
    F12 = -(1 - cosmu ** 2) / 2

    # Rotation angle sigma_S
    sigma_S = compute_sigma_s(th_sun, ph_sun, TH_S, PHI_S, mu)

    # 180-degree symmetry correction
    if ph_sun > 180:
        sigma_S = 180 - sigma_S

    # Scattered Stokes parameters (normalized)
    ratio = F12 / F11
    Q = np.cos(np.radians(2 * sigma_S)) * ratio
    U = -np.sin(np.radians(2 * sigma_S)) * ratio

    # Polarization metrics
    DoLP = np.abs(ratio)
    DoCP = np.zeros_like(mu)
    AoP = np.mod(0.5 * np.degrees(np.arctan2(U, Q)), 180)

    return {
        'azimuth': azimuth,
        'zenith': zenith,
        'DoLP': DoLP,
        'DoCP': DoCP,
        'AoP': AoP,
        'Q': Q,
        'U': U,
        'I': np.ones_like(mu),
        'mu': mu,
    }


def sky_polarization(observer, time_info, resolution=1.5, wavelength=550e-9):
    """Full-sky Rayleigh scattering polarization map.

    Parameters
    ----------
    observer : dict with keys 'latitude', 'longitude', 'altitude', 'gmtOffset'
    time_info : dict with keys 'year', 'month', 'day', 'hour', and optional 'dst'
    resolution : float — Angular resolution in degrees (default 1.5)
    wavelength : float — Wavelength in metres (default 550e-9)

    Returns
    -------
    result : dict — All polarization quantities plus sun position and inputs.
    """
    if 'dst' not in time_info:
        time_info = {**time_info, 'dst': 0}

    # Step 1: Compute sun position
    sun_zenith, sun_azimuth, sun_elevation = compute_sun_position(observer, time_info)

    if sun_elevation <= 0:
        import warnings
        warnings.warn(
            f"Sun is below the horizon (elevation = {sun_elevation:.1f} deg). "
            "Results are physically valid only for sunElevation > 0.")

    # Step 2: Particle and wavelength parameters
    params = {
        'wavelength': wavelength,
        'particleRadius': 100e-9,
        'refractiveIndex': 1.53 + 0.007j,
        'resolution': resolution,
    }

    # Step 3: Compute sky polarization map
    sky_map = compute_rayleigh_sky_map(sun_zenith, sun_azimuth, params)

    # Assemble output
    result = dict(sky_map)
    result['sun'] = {
        'zenith': sun_zenith,
        'azimuth': sun_azimuth,
        'elevation': sun_elevation,
    }
    result['observer'] = observer
    result['time'] = time_info

    return result
