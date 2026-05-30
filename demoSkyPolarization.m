%DEMOSKYPOLARIZATION  Example script for the sky polarization simulator.
%
%   Demonstrates how to use skyPolarization() to generate and visualise
%   full-sky Rayleigh polarization maps for a given observer location and time.
%
%   Run this script section by section (Ctrl+Enter) or all at once.

clc; clear; close all;

%% =========================================================================
%  SECTION 1: Basic usage - single time snapshot
%  =========================================================================

% --- Observer location ---
% Example: Tempe, Arizona, USA
observer.latitude  =  33.42;    % degrees N
observer.longitude = -111.93;   % degrees E  (negative = West)
observer.altitude  =  340;      % metres (stored in result; not used in physics)
observer.gmtOffset = -7;        % UTC-7  (Mountain Standard Time)

% --- Date and local clock time ---
time.year  = 2024;
time.month = 6;
time.day   = 21;    % Summer solstice
time.hour  = 10.0;  % 10:00 AM local time (decimal hours)
time.dst   = 1;     % Daylight saving is active in summer (+1 hour)

% --- Run simulation ---
fprintf('Computing sky polarization map...\n');
result = skyPolarization(observer, time);

% --- Print summary ---
fprintf('\n--- Sun position ---\n');
fprintf('  Zenith    : %.2f deg\n', result.sun.zenith);
fprintf('  Azimuth   : %.2f deg (from North, clockwise)\n', result.sun.azimuth);
fprintf('  Elevation : %.2f deg\n', result.sun.elevation);

fprintf('\n--- Sky map dimensions ---\n');
fprintf('  Azimuth points : %d  (%.1f deg resolution)\n', length(result.azimuth), ...
    result.azimuth(2) - result.azimuth(1));
fprintf('  Zenith  points : %d\n', length(result.zenith));

fprintf('\n--- Polarization statistics ---\n');
fprintf('  DoLP max  : %.4f  (at scattering angle ~ 90 deg)\n', max(result.DoLP(:)));
fprintf('  DoLP min  : %.4f  (at sun and anti-sun positions)\n', min(result.DoLP(:)));

% --- Visualise ---
plotSkyPolarization(result);

%% =========================================================================
%  SECTION 2: Query specific sky points
%  =========================================================================

% --- Look up polarization at a specific azimuth and zenith ---
query_azimuth = 180;   % due South
query_zenith  = 45;    % 45 deg from overhead

% Find nearest grid indices
[~, ia] = min(abs(result.azimuth - query_azimuth));
[~, iz] = min(abs(result.zenith  - query_zenith));

fprintf('\n--- Polarization at az=%.0f deg, zen=%.0f deg ---\n', ...
    result.azimuth(ia), result.zenith(iz));
fprintf('  Scattering angle (mu) : %.2f deg\n',  result.mu(iz, ia));
fprintf('  DoLP                  : %.4f\n',       result.DoLP(iz, ia));
fprintf('  AoP                   : %.2f deg\n',   result.AoP(iz, ia));
fprintf('  Q (normalised)        : %.4f\n',       result.Q(iz, ia));
fprintf('  U (normalised)        : %.4f\n',       result.U(iz, ia));

%% =========================================================================
%  SECTION 3: Effect of wavelength  (blue vs red sky)
%  =========================================================================

observer_wl = observer;
time_wl     = time;

figure('Color', 'w', 'Name', 'Wavelength comparison');

wavelengths  = [450e-9, 550e-9, 650e-9];   % blue, green, red
wavelabel    = {'450 nm (blue)', '550 nm (green)', '650 nm (red)'};
colors       = {[0.2, 0.2, 1], [0.1, 0.7, 0.1], [1, 0.2, 0.2]};

% DoLP is the same for all wavelengths in Rayleigh theory
% (it depends only on the scattering angle mu, not lambda).
% Intensity, however, scales as lambda^(-4) -- blue is brighter.
fprintf('\n--- DoLP at mu = 90 deg for different wavelengths ---\n');
for k = 1:length(wavelengths)
    r_wl = skyPolarization(observer_wl, time_wl, 'Wavelength', wavelengths(k));
    % DoLP at the anti-solar point (mu ~ 90 deg on the solar-antisolar great circle)
    max_dolp = max(r_wl.DoLP(:));
    intensity_relative = (wavelengths(1) / wavelengths(k))^4;
    fprintf('  %s :  max DoLP = %.4f,  relative intensity = %.2f\n', ...
        wavelabel{k}, max_dolp, intensity_relative);
end
fprintf('(DoLP pattern is wavelength-independent for Rayleigh scattering.)\n');

%% =========================================================================
%  SECTION 4: Diurnal variation of the neutral point pattern
%  =========================================================================

observer_d = observer;
time_d     = time;
time_d.dst = 1;

hours = 6 : 2 : 18;    % 6 AM to 6 PM in 2-hour steps

figure('Color', 'w', 'Name', 'Diurnal variation - DoLP');
n_hours = length(hours);
cols    = ceil(n_hours / 2);

for k = 1:n_hours
    time_d.hour = hours(k);
    r_d = skyPolarization(observer_d, time_d, 'Resolution', 3);

    subplot(2, cols, k);

    [AZ_g, ZE_g] = meshgrid(r_d.azimuth, r_d.zenith);
    x_d =  ZE_g .* sind(AZ_g);
    y_d =  ZE_g .* cosd(AZ_g);

    contourf(x_d, y_d, r_d.DoLP, 50, 'EdgeColor', 'none');
    colormap(gca, parula(256));
    caxis([0, 1]); axis equal; axis off;

    % Horizon and sun
    th_c = linspace(0, 360, 360);
    hold on;
    plot(90*sind(th_c), 90*cosd(th_c), 'k-', 'LineWidth', 1);
    sunX_d =  r_d.sun.zenith * sind(r_d.sun.azimuth);
    sunY_d =  r_d.sun.zenith * cosd(r_d.sun.azimuth);
    plot(sunX_d, sunY_d, 'w*', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold off;

    title(sprintf('%02d:00  (elv=%.0f°)', hours(k), r_d.sun.elevation), 'FontSize', 9);
end
sgtitle(sprintf('DoLP throughout the day  |  %04d-%02d-%02d  |  %.1f°N %.1f°E', ...
    time.year, time.month, time.day, observer.latitude, observer.longitude), ...
    'FontSize', 10, 'Interpreter', 'none');

%% =========================================================================
%  SECTION 5: Access raw data for further processing
%  =========================================================================

% The result struct gives you direct access to all computed quantities.
%
% result.azimuth   [1 x Na]   azimuth grid, degrees
% result.zenith    [1 x Nz]   zenith grid, degrees
% result.DoLP      [Nz x Na]  degree of linear polarization
% result.AoP       [Nz x Na]  angle of polarization, degrees [0, 180)
% result.Q         [Nz x Na]  normalised Stokes Q
% result.U         [Nz x Na]  normalised Stokes U
% result.mu        [Nz x Na]  scattering angle at each sky point, degrees
% result.sun       struct     .zenith, .azimuth, .elevation

% Example: find the great circle of maximum polarization (mu = 90 deg)
result2 = skyPolarization(observer, time, 'Resolution', 1.0);
mask_90 = abs(result2.mu - 90) < 2;   % pixels within 2 deg of mu=90

fprintf('\n--- Pixels near mu = 90 deg (maximum DoLP arc) ---\n');
fprintf('  Mean DoLP  : %.4f\n', mean(result2.DoLP(mask_90)));
fprintf('  Mean AoP   : %.2f deg\n', mean(result2.AoP(mask_90)));
fprintf('  Count      : %d pixels\n', sum(mask_90(:)));
