function sigma_S = computeSigmaS(thSun, phSun, TH_S, PHI_S, mu)
%COMPUTESIGMAS  Rotation angle from the scattering plane to the local meridian plane.
%
%   sigma_S = computeSigmaS(thSun, phSun, TH_S, PHI_S, mu)
%
%   sigma_S is the angle by which the Stokes reference frame must be rotated
%   at the scattered-ray direction to transform from the scattering plane
%   coordinate system to the local meridian plane coordinate system.
%
%   The formulas follow the sign convention used in:
%     Bohren & Huffman, "Absorption and Scattering of Light by Small Particles"
%     http://www.oceanopticsbook.info/view/light_and_radiometry/level_2/
%            polarization_scattering_geometry
%
%   INPUTS
%   ------
%   thSun    Sun zenith angle, degrees (scalar)
%   phSun    Sun azimuth angle, degrees (scalar); may be temporarily modified
%            to phSun - 180 when phSun > 180 (handled in computeRayleighSkyMap)
%   TH_S     Sky-point zenith angles,   [Nz x Na] degrees
%   PHI_S    Sky-point azimuth angles,  [Nz x Na] degrees
%   mu       Scattering angles,          [Nz x Na] degrees
%
%   OUTPUT
%   ------
%   sigma_S  [Nz x Na]  Rotation angles, degrees

% When the sun azimuth > 180, the caller passes (phSun - 180) as phSun_mod.
% Here we use phSun as received (the modification is done in computeRayleighSkyMap).
phSun_mod = phSun;
if phSun > 180
    phSun_mod = phSun - 180;
end

del_phi = PHI_S - phSun_mod;   % azimuth difference, used for case selection

sigma_S = zeros(size(TH_S));

% ---- Identify special cases ----------------------------------------------
% Singular scattering: forward (mu=0) or backward (mu=180) scatter.
% Both the incident and scattered rays are collinear, so rotation is undefined.
mu_singular  = (abs(mu) < 1e-10) | (abs(mu - 180) < 1e-10);

% Sky point is at the zenith (TH_S = 0): sin(TH_S) = 0 in the denominator.
thS_zero     = TH_S < 1e-10;

% Sun is at the zenith (thSun = 0): different geometry, sigma_S = 0.
thSun_zero   = thSun < 1e-10;

% General case: no degenerate geometry
general      = ~mu_singular & ~thS_zero & ~thSun_zero;

% ---- Argument for acosd in the general formula ---------------------------
% From spherical trigonometry (law of cosines on the unit sphere):
%
%   cos(sigma_S) = [cos(thSun) - cos(TH_S)*cos(mu)] / [sin(mu)*sin(TH_S)]
%
% The sign of sigma_S (i.e., which quadrant it falls in) is determined by
% whether del_phi falls in [0, 180] or outside that range.
denom = sind(mu) .* sind(TH_S);
denom(denom == 0) = NaN;       % avoid divide-by-zero (handled by masks)

arg = (cosd(thSun) - cosd(TH_S).*cosd(mu)) ./ denom;
arg = max(-1, min(1, arg));    % numerical clamp to [-1, 1]

% ---- General case, Case A: del_phi in [0, 180] ---------------------------
% Sun is in the same or "left" hemisphere relative to the sky point.
maskA = general & (del_phi >= 0) & (del_phi <= 180);
if any(maskA(:))
    sigma_S(maskA) = 180 + acosd(arg(maskA));
end

% ---- General case, Case B: del_phi outside [0, 180] ---------------------
% Sun is in the "right" hemisphere relative to the sky point.
maskB = general & ~maskA;
if any(maskB(:))
    sigma_S(maskB) = 360 - acosd(arg(maskB));
end

% ---- Special case: sky point at zenith (TH_S = 0) -----------------------
% The local meridian plane is ill-defined for a ray aimed straight up.
% Use a limiting form derived from the spherical geometry.
mask_thS0 = ~mu_singular & thS_zero & ~thSun_zero;
if any(mask_thS0(:))

    % Case A: del_phi in [0, 180]
    subA = mask_thS0 & (del_phi >= 0) & (del_phi <= 180);
    if any(subA(:))
        sigma_S(subA) = acosd(cosd(thSun) * cosd(PHI_S(subA) - phSun_mod));
    end

    % Case B: del_phi outside [0, 180]
    subB = mask_thS0 & ~subA;
    if any(subB(:))
        phi_s_vals = PHI_S(subB);
        wrapped    = zeros(size(phi_s_vals));
        gt         = phi_s_vals > phSun_mod;
        wrapped( gt) = phSun_mod + (360 - phi_s_vals( gt));
        wrapped(~gt) = phSun_mod -          phi_s_vals(~gt);
        sigma_S(subB) = acosd(cosd(thSun) .* cosd(wrapped));
    end

end

% ---- Special case: sun at zenith (thSun = 0) ----------------------------
% When the sun is directly overhead, the meridian plane reference is arbitrary.
% sigma_S = 0 by convention (no polarisation rotation needed).
sigma_S(~mu_singular & ~thS_zero & thSun_zero) = 0;

% ---- Special case: both sun and sky point at zenith ----------------------
sigma_S(thS_zero & thSun_zero) = 0;

% ---- Singular scattering angles -----------------------------------------
% Forward and backward scatter produce no net polarisation (DoLP = 0),
% so the rotation angle is irrelevant; set to 0 for clean output.
sigma_S(mu_singular) = 0;

end
