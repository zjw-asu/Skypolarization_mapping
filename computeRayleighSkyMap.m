function skyMap = computeRayleighSkyMap(sunZenith, sunAzimuth, params)
%COMPUTERAYLEIGHSKYMAP  Compute sky polarization via single Rayleigh scattering.
%
%   skyMap = computeRayleighSkyMap(sunZenith, sunAzimuth, params)
%
%   Models the polarization of skylight produced by single Rayleigh scattering
%   of unpolarized sunlight by spherical aerosol particles (particle size << wavelength).
%
%   The calculation uses the Mueller matrix (Stokes) formalism:
%
%     S_out = L(sigma_S) * F_Rayl * L(sigma) * S_in
%
%   where:
%     S_in     = [I, 0, 0, 0]'  (unpolarized sunlight)
%     F_Rayl   = Rayleigh scattering Mueller matrix in the scattering plane
%     L(sigma) = rotation matrix from meridian plane to scattering plane (incident)
%     L(sigma_S) = rotation matrix from scattering plane to meridian plane (scattered)
%
%   Because S_in is unpolarized, L(sigma) * S_in = S_in, so only sigma_S matters.
%
%   INPUTS
%   ------
%   sunZenith   Sun zenith angle, degrees (0 = overhead, 90 = horizon)
%   sunAzimuth  Sun azimuth angle, degrees (from North, clockwise)
%   params      struct with fields:
%     .wavelength       Wavelength in metres
%     .particleRadius   Particle radius in metres
%     .refractiveIndex  Complex refractive index  m = mr + i*mi
%     .resolution       Sky grid angular step size in degrees
%
%   OUTPUT
%   ------
%   skyMap  struct with fields:
%     .azimuth  [1 x Na]   Azimuth angles, degrees
%     .zenith   [1 x Nz]   Zenith angles, degrees
%     .DoLP     [Nz x Na]  Degree of linear polarization [0, 1]
%     .DoCP     [Nz x Na]  Degree of circular polarization (0 for Rayleigh)
%     .AoP      [Nz x Na]  Angle of polarization, degrees [0, 180)
%     .Q        [Nz x Na]  Normalized Stokes Q
%     .U        [Nz x Na]  Normalized Stokes U
%     .I        [Nz x Na]  Normalized intensity (all ones after normalization)
%     .mu       [Nz x Na]  Scattering angle, degrees

% ---- Unpack parameters ---------------------------------------------------
lambda     = params.wavelength;
r          = params.particleRadius;
m          = params.refractiveIndex;
resolution = params.resolution;

thSun = sunZenith;
phSun = sunAzimuth;

% ---- Build sky grid ------------------------------------------------------
% zenith: 0 (overhead) to 90 (horizon); azimuth: 0 to 360 (North, CW)
zenith  = linspace(0, 90,  round(90  / resolution) + 1);
azimuth = linspace(0, 360, round(360 / resolution) + 1);

% Create 2-D grids: each matrix is [Nz x Na]
[PHI_S, TH_S] = meshgrid(azimuth, zenith);

% ---- Rayleigh scattering amplitude scale ---------------------------------
% For a sphere of radius r << lambda and complex refractive index m:
%   size parameter      a   = 2*pi*r / lambda
%   Clausius-Mossotti factor chi = (m^2 - 1) / (m^2 + 2)
%   Scattering intensity scales as  a^6 * |chi|^2
% After normalising scattered Stokes by intensity, this factor cancels.
a            = 2*pi*r / lambda;
chiSquared   = abs((m^2 - 1) / (m^2 + 2))^2;
particleScale = a^6 * chiSquared;  %#ok<NASGU>  (cancels in normalisation)

% ---- Scattering angle mu at each sky point -------------------------------
% mu: angle between the incident ray (from sun) and the scattered ray (to sky point).
% From the dot product of the two unit direction vectors on the unit sphere:
%
%   cos(mu) = sin(thSun)*sin(TH_S)*cos(PHI_S - phSun) + cos(thSun)*cos(TH_S)
%
% mu = 0   : looking directly toward the sun (forward scatter)
% mu = 90  : looking perpendicular to the sun direction (max polarization)
% mu = 180 : looking directly away from the sun (back scatter)
mu = acosd( sind(thSun)*sind(TH_S).*cosd(PHI_S - phSun) + cosd(thSun)*cosd(TH_S) );

% ---- Rayleigh Mueller matrix elements ------------------------------------
% The Rayleigh scattering matrix F is defined in the scattering plane
% (the plane containing both the incident and scattered ray directions).
%
%   F = [ F11  F12   0    0  ]       F11 = (1 + cos^2(mu)) / 2
%       [ F12  F11   0    0  ]       F12 = -(1 - cos^2(mu)) / 2  [always <= 0]
%       [  0    0   F33   0  ]       F33 = cos(mu)
%       [  0    0    0   F33 ]
%
% F12/F11 = -sin^2(mu) / (1 + cos^2(mu))  in [-1, 0]
% DoLP = |F12/F11|                         in [ 0, 1]  (max at mu = 90 deg)
cosmu = cosd(mu);
F11   = (1 + cosmu.^2) / 2;
F12   = -(1 - cosmu.^2) / 2;

% ---- Rotation angle sigma_S ---------------------------------------------
% sigma_S rotates the Stokes reference frame from the scattering plane to
% the local meridian plane at the scattered ray direction (sky point).
% This is required because the scattering plane changes orientation for
% different sky directions, but Stokes parameters must be expressed in a
% common (local meridian) reference frame.
%
% See: Mishchenko et al., "Scattering, Absorption and Emission of Light
%      by Small Particles", Chapter 2.
sigma_S = computeSigmaS(thSun, phSun, TH_S, PHI_S, mu);

% When sun azimuth > 180 deg, apply the 180-degree symmetry correction
% to the rotation matrices (preserves the correct polarization reference frame).
if phSun > 180
    sigma_S = 180 - sigma_S;
end

% ---- Scattered Stokes parameters (normalized) ---------------------------
% For unpolarized input S_in = [I, 0, 0, 0]:
%   L(sigma) * S_in = S_in        (rotation does not change unpolarized light)
%   F_Rayl  * S_in = [F11*I, F12*I, 0, 0]'    (in scattering plane)
%   L(sigma_S) rotates this to the local meridian plane:
%     I_out = F11 * I             (after normalising to 1)
%     Q_out = cos(2*sigma_S) * (F12/F11)
%     U_out = -sin(2*sigma_S) * (F12/F11)
%     V_out = 0                   (Rayleigh scattering preserves zero circular pol.)
ratio = F12 ./ F11;      % always in [-1, 0]; equals -sin^2(mu)/(1+cos^2(mu))

Q = cosd(2*sigma_S) .* ratio;
U = -sind(2*sigma_S) .* ratio;

% ---- Polarization metrics ------------------------------------------------
% Degree of linear polarization: DoLP = sqrt(Q^2 + U^2) = |F12/F11|
DoLP = abs(ratio);    % equivalent to sin^2(mu) / (1 + cos^2(mu))

% Degree of circular polarization: zero for single Rayleigh scattering
DoCP = zeros(size(mu));

% Angle of polarization: orientation of the E-field oscillation [0, 180) deg.
% AoP = 0.5 * atan2(U, Q) maps to (-90, 90]; mod(..., 180) wraps to [0, 180).
AoP = mod(0.5 * atan2d(U, Q), 180);

% ---- Assemble output -----------------------------------------------------
skyMap.azimuth = azimuth;
skyMap.zenith  = zenith;
skyMap.DoLP    = DoLP;
skyMap.DoCP    = DoCP;
skyMap.AoP     = AoP;
skyMap.Q       = Q;
skyMap.U       = U;
skyMap.I       = ones(size(mu));
skyMap.mu      = mu;

end
