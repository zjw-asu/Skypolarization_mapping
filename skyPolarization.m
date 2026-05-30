function result = skyPolarization(observer, time, varargin)
%SKYPOLARIZATION  Full-sky Rayleigh scattering polarization map.
%
%   result = skyPolarization(observer, time)
%   result = skyPolarization(observer, time, 'Name', value, ...)
%
%   Models single Rayleigh scattering of unpolarized sunlight by atmospheric
%   aerosol particles (r << lambda). Uses Mueller matrix formalism with
%   coordinate rotation from the scattering plane to the local meridian plane.
%   The input Stokes vector is S = [I, 0, 0, 0] (fully unpolarized sunlight).
%
%   -------------------------------------------------------------------------
%   INPUTS
%   -------------------------------------------------------------------------
%   observer  struct - Observer's geographic location
%     .latitude   Latitude in degrees  (N positive, S negative)
%     .longitude  Longitude in degrees (E positive, W negative)
%     .altitude   Altitude in metres   (stored in output; not used in physics)
%     .gmtOffset  UTC offset in hours  (e.g. -7 for MST Arizona, +8 for CST China)
%
%   time      struct - Local clock time of observation
%     .year   Four-digit year
%     .month  Month number (1-12)
%     .day    Day of month (1-31)
%     .hour   Local clock time in decimal hours  (e.g. 14.5 = 2:30 PM)
%     .dst    Daylight saving offset in hours    (0 or 1; default 0)
%
%   -------------------------------------------------------------------------
%   OPTIONAL NAME-VALUE PAIRS
%   -------------------------------------------------------------------------
%   'Resolution'   Angular resolution of sky grid in degrees (default: 1.5)
%                  Smaller values give finer maps but increase compute time.
%   'Wavelength'   Wavelength in metres (default: 550e-9, green light)
%
%   -------------------------------------------------------------------------
%   OUTPUT
%   -------------------------------------------------------------------------
%   result  struct
%     .azimuth   [1 x Na]   Azimuth grid, degrees, from North clockwise
%     .zenith    [1 x Nz]   Zenith grid, degrees (0 = overhead, 90 = horizon)
%     .DoLP      [Nz x Na]  Degree of linear polarization [0, 1]
%     .DoCP      [Nz x Na]  Degree of circular polarization (zero for Rayleigh)
%     .AoP       [Nz x Na]  Angle of polarization, degrees [0, 180)
%     .Q         [Nz x Na]  Normalized Stokes Q component
%     .U         [Nz x Na]  Normalized Stokes U component
%     .I         [Nz x Na]  Normalized intensity (1 everywhere for single scatter)
%     .mu        [Nz x Na]  Scattering angle at each sky point, degrees
%     .sun       struct     Sun position
%       .zenith    Sun zenith angle, degrees
%       .azimuth   Sun azimuth angle, degrees (from North CW)
%       .elevation Sun elevation angle, degrees
%     .observer  struct     Copy of input observer
%     .time      struct     Copy of input time
%
%   -------------------------------------------------------------------------
%   EXAMPLE
%   -------------------------------------------------------------------------
%   observer.latitude  = 33.42;    % Tempe, Arizona
%   observer.longitude = -111.93;
%   observer.altitude  = 340;      % metres
%   observer.gmtOffset = -7;
%
%   time.year  = 2024;
%   time.month = 6;
%   time.day   = 21;    % Summer solstice
%   time.hour  = 12.0;  % Local noon
%   time.dst   = 1;     % Daylight saving active
%
%   result = skyPolarization(observer, time);
%   plotSkyPolarization(result);

% ---- Parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p, 'observer', @isstruct);
addRequired(p, 'time',     @isstruct);
addParameter(p, 'Resolution', 1.5,    @(x) isscalar(x) && x > 0 && x <= 15);
addParameter(p, 'Wavelength', 550e-9, @(x) isscalar(x) && x > 0);
parse(p, observer, time, varargin{:});

resolution = p.Results.Resolution;
wavelength = p.Results.Wavelength;

if ~isfield(time, 'dst')
    time.dst = 0;
end

% ---- Step 1: Compute sun position ----------------------------------------
[sunZenith, sunAzimuth, sunElevation] = computeSunPosition(observer, time);

if sunElevation <= 0
    warning('skyPolarization:SunBelowHorizon', ...
        'Sun is below the horizon (elevation = %.1f deg). ', ...
        'Results are physically valid only for sunElevation > 0.', sunElevation);
end

% ---- Step 2: Particle and wavelength parameters --------------------------
% Default: 100 nm radius spherical aerosol, typical for atmospheric haze.
% Refractive index m = 1.53 + 0.007i is standard for continental aerosol.
params.wavelength       = wavelength;
params.particleRadius   = 100e-9;          % metres
params.refractiveIndex  = 1.53 + 0.007i;  % complex index, mr + i*mi
params.resolution       = resolution;

% ---- Step 3: Compute sky polarization map --------------------------------
skyMap = computeRayleighSkyMap(sunZenith, sunAzimuth, params);

% ---- Assemble output struct ----------------------------------------------
result          = skyMap;
result.sun.zenith    = sunZenith;
result.sun.azimuth   = sunAzimuth;
result.sun.elevation = sunElevation;
result.observer      = observer;
result.time          = time;

end
