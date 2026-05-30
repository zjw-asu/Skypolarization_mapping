function plotSkyPolarization(result)
%PLOTSKYPOLARIZATION  Visualise sky polarization in 2D fisheye and 3D sphere views.
%
%   plotSkyPolarization(result)
%
%   Produces TWO figures from the struct returned by skyPolarization():
%
%   Figure 1 – 2D Fisheye Polar Projection
%     Three side-by-side panels: DoLP | AoP | DoP
%     Centre = zenith overhead.  Edge = horizon (90° from zenith).
%     AoP panel includes polarization tick marks whose length encodes DoLP
%     and whose angle encodes the E-field oscillation direction.
%     Dashed rings mark 30° and 60° zenith distances.
%
%   Figure 2 – 3D Hemisphere Projection
%     Same three quantities mapped as colour onto a 3D dome.
%     Sun position is marked with a yellow sphere.
%     All panels can be rotated interactively (click-and-drag).
%
%   INPUT
%   -----
%   result   Struct returned by skyPolarization().
%
%   EXAMPLE
%   -------
%   result = skyPolarization(observer, time);
%   plotSkyPolarization(result);

% =========================================================================
%  Pre-compute shared quantities
% =========================================================================

% Degree of polarization  DoP = sqrt(Q^2 + U^2 + V^2) / I
% For single Rayleigh scatter V=0, so DoP equals DoLP numerically.
% We compute it from the Stokes components to make the formula explicit.
DoP  = sqrt(result.Q.^2 + result.U.^2);   % [Nz x Na]
DoLP = result.DoLP;                         % [Nz x Na]
AoP  = result.AoP;                          % [Nz x Na], degrees [0, 180)

% Grid: azimuth along columns (Na), zenith along rows (Nz)
[AZ_grid, ZE_grid] = meshgrid(result.azimuth, result.zenith);

% ---- 2D fisheye Cartesian (x = East, y = North) -------------------------
x2d =  ZE_grid .* sind(AZ_grid);
y2d =  ZE_grid .* cosd(AZ_grid);

% ---- 3D hemisphere unit-sphere coordinates -------------------------------
% x = East,  y = North,  z = Up
X3d = sind(ZE_grid) .* sind(AZ_grid);
Y3d = sind(ZE_grid) .* cosd(AZ_grid);
Z3d = cosd(ZE_grid);

% ---- Sun position --------------------------------------------------------
sunAboveHorizon = result.sun.elevation > 0;

% 2D fisheye
sunX2d =  result.sun.zenith * sind(result.sun.azimuth);
sunY2d =  result.sun.zenith * cosd(result.sun.azimuth);

% 3D hemisphere
sunX3d = sind(result.sun.zenith) * sind(result.sun.azimuth);
sunY3d = sind(result.sun.zenith) * cosd(result.sun.azimuth);
sunZ3d = cosd(result.sun.zenith);

% ---- Shared title string -------------------------------------------------
obs = result.observer;
tim = result.time;
mainTitle = sprintf( ...
    '%04d-%02d-%02d   %05.2f h local  |  %.2f°%s  %.2f°%s  |  Sun: zen=%.1f°  az=%.1f°  elv=%.1f°', ...
    tim.year, tim.month, tim.day, tim.hour, ...
    abs(obs.latitude),  char('N' + ('S'-'N')*(obs.latitude  < 0)), ...
    abs(obs.longitude), char('E' + ('W'-'E')*(obs.longitude < 0)), ...
    result.sun.zenith, result.sun.azimuth, result.sun.elevation);

% =========================================================================
%  FIGURE 1 – 2D Fisheye Polar
% =========================================================================
figure('Color', 'w', 'Name', '2D Fisheye Polar – Sky Polarization', ...
       'Position', [40, 120, 1380, 520]);

panel2d = { DoLP, [0 1],    parula(512), 'DoLP',       'Degree of Linear Polarization',  false; ...
            AoP,  [0 180],  hsv(512),    'AoP  (deg)', 'Angle of Polarization',           true;  ...
            DoP,  [0 1],    parula(512), 'DoP',        'Degree of Polarization',          false };
%          data   clim      colormap     cbar-label    panel-title                        ticks?

th_ring = linspace(0, 360, 361);   % used to draw circles

for k = 1:3
    subplot(1, 3, k);

    data      = panel2d{k, 1};
    climits   = panel2d{k, 2};
    cmap      = panel2d{k, 3};
    cbLabel   = panel2d{k, 4};
    panTitle  = panel2d{k, 5};
    showTicks = panel2d{k, 6};

    % -- filled colour map in fisheye coordinates --------------------------
    contourf(x2d, y2d, data, 256, 'EdgeColor', 'none');
    colormap(gca, cmap);
    caxis(climits);
    cb = colorbar('Location', 'southoutside');
    cb.Label.String  = cbLabel;
    cb.Label.FontSize = 9;

    hold on;

    % -- reference rings: 30° and 60° zenith distance ----------------------
    for zring = [30, 60]
        plot(zring*sind(th_ring), zring*cosd(th_ring), '--', ...
             'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
        text(0, -zring - 1.5, sprintf('%d°', zring), ...
             'FontSize', 7, 'Color', [0.45 0.45 0.45], ...
             'HorizontalAlignment', 'center');
    end

    % -- horizon ring (90°) ------------------------------------------------
    plot(90*sind(th_ring), 90*cosd(th_ring), 'k-', 'LineWidth', 1.8);

    % -- cardinal direction labels -----------------------------------------
    text( 0,   98, 'N', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    text( 0,  -98, 'S', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    text( 98,   0, 'E', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);
    text(-98,   0, 'W', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 11);

    % -- AoP tick marks (on AoP panel only) --------------------------------
    if showTicks
        drawPolarizationTicks(result.azimuth, result.zenith, AoP, DoLP);
    end

    % -- sun marker --------------------------------------------------------
    if sunAboveHorizon
        plot(sunX2d, sunY2d, 'o', ...
             'MarkerSize', 16, 'MarkerFaceColor', [1.0, 0.92, 0.0], ...
             'MarkerEdgeColor', [0.55, 0.40, 0.0], 'LineWidth', 2);
        text(sunX2d + 5, sunY2d + 5, 'Sun', 'FontSize', 8, ...
             'Color', [0.3 0.3 0.3]);
    else
        % Show where sun would be below horizon
        plot(sunX2d, sunY2d, 'v', ...
             'MarkerSize', 10, 'MarkerFaceColor', [0.5, 0.5, 0.5], ...
             'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
    end

    axis equal; axis off;
    xlim([-107, 107]); ylim([-107, 107]);
    title(panTitle, 'FontSize', 12, 'FontWeight', 'bold');
    hold off;
end

sgtitle(mainTitle, 'FontSize', 9.5, 'Interpreter', 'none', 'FontWeight', 'normal');
drawnow;

% =========================================================================
%  FIGURE 2 – 3D Hemisphere
% =========================================================================
figure('Color', 'w', 'Name', '3D Hemisphere – Sky Polarization', ...
       'Position', [40, 80, 1380, 600]);

panel3d = { DoLP, [0 1],   parula(512), 'DoLP',       'Degree of Linear Polarization'; ...
            AoP,  [0 180], hsv(512),    'AoP  (deg)', 'Angle of Polarization';          ...
            DoP,  [0 1],   parula(512), 'DoP',        'Degree of Polarization'         };

% Azimuth lines for compass grid (4 lines at 0/90/180/270)
azGrid = linspace(0, 360, 361);
th_zen = linspace(0, 90, 91);          % zenith 0 → 90 (horizon)

for k = 1:3
    ax3 = subplot(1, 3, k);

    data     = panel3d{k, 1};
    climits  = panel3d{k, 2};
    cmap     = panel3d{k, 3};
    cbLabel  = panel3d{k, 4};
    panTitle = panel3d{k, 5};

    % -- hemisphere surface ------------------------------------------------
    hs = surf(X3d, Y3d, Z3d, data, ...
              'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceAlpha', 1.0);
    colormap(ax3, cmap);
    caxis(climits);
    cb3 = colorbar('Location', 'southoutside');
    cb3.Label.String  = cbLabel;
    cb3.Label.FontSize = 9;

    hold on;

    % -- ground disc -------------------------------------------------------
    th_g = linspace(0, 2*pi, 361);
    fill3(cos(th_g), sin(th_g), zeros(1, 361), [0.88, 0.93, 0.88], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.5);

    % -- horizon ring at z = 0 ---------------------------------------------
    plot3(cos(th_g), sin(th_g), zeros(1, 361), 'k-', 'LineWidth', 1.5);

    % -- elevation reference rings at 30° and 60° elevation ---------------
    for elv_ring = [30, 60]
        zen_ring = 90 - elv_ring;
        xr = sind(zen_ring) * sind(azGrid);
        yr = sind(zen_ring) * cosd(azGrid);
        zr = cosd(zen_ring) * ones(size(azGrid));
        plot3(xr, yr, zr, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end

    % -- cardinal meridian lines (N-S, E-W) --------------------------------
    for az_line = [0, 90, 180, 270]
        xl = sind(th_zen) .* sind(az_line);
        yl = sind(th_zen) .* cosd(az_line);
        zl = cosd(th_zen);
        plot3(xl, yl, zl, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    end

    % -- cardinal direction labels at horizon level ------------------------
    text( 0,    1.13, 0.0, 'N', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text( 0,   -1.13, 0.0, 'S', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text( 1.13,  0,   0.0, 'E', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');
    text(-1.13,  0,   0.0, 'W', 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');

    % -- elevation ring labels ---------------------------------------------
    text(sind(60)*sind(45), sind(60)*cosd(45), cosd(60) + 0.03, ...
         '30° elv', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
    text(sind(30)*sind(45), sind(30)*cosd(45), cosd(30) + 0.03, ...
         '60° elv', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);

    % -- sun marker --------------------------------------------------------
    if sunAboveHorizon
        % Raised slightly above the surface to be visible
        offset = 0.04;
        plot3(sunX3d*(1+offset), sunY3d*(1+offset), sunZ3d*(1+offset), 'o', ...
              'MarkerSize', 18, 'MarkerFaceColor', [1.0, 0.92, 0.0], ...
              'MarkerEdgeColor', [0.55, 0.40, 0.0], 'LineWidth', 2.5);
        % Vertical line from ground to sun position for depth reference
        plot3([sunX3d, sunX3d], [sunY3d, sunY3d], [0, sunZ3d], ...
              ':', 'Color', [0.7 0.6 0.1], 'LineWidth', 1.2);
    end

    % -- axis, view, lighting ----------------------------------------------
    axis equal; axis off;
    xlim([-1.3, 1.3]); ylim([-1.3, 1.3]); zlim([-0.05, 1.15]);

    % Perspective view: slightly south-west and 30° above horizon
    view(-30, 30);

    % Subtle diffuse lighting for depth cues without distorting colour
    light('Position', [0, 0, 3], 'Style', 'infinite');
    lighting flat;
    material([0.85, 0.15, 0, 1, 0]);   % mostly ambient, slight diffuse, no specular

    title(panTitle, 'FontSize', 12, 'FontWeight', 'bold');

    % Enable interactive rotation
    rotate3d(ax3, 'on');

    hold off;
end

sgtitle(mainTitle, 'FontSize', 9.5, 'Interpreter', 'none', 'FontWeight', 'normal');
drawnow;

end   % end of plotSkyPolarization


% =========================================================================
%  LOCAL HELPER – polarization tick marks on the 2D fisheye panel
% =========================================================================
function drawPolarizationTicks(azimuth, zenith, AoP, DoLP)
%DRAWPOLARIZATIONTICKS  Overlay E-field direction segments on the AoP panel.
%
%   Segments are drawn at a subsampled grid of sky points.
%   Length  ∝ DoLP  (longer = more polarized)
%   Angle   = AoP   (orientation of the E-field, measured from North CW)
%
%   Convention (x=East, y=North):
%     AoP =   0°  →  N-S segment  → (dx, dy) = (0,  1)
%     AoP =  90°  →  E-W segment  → (dx, dy) = (1,  0)
%     AoP =  45°  →  NE-SW        → (dx, dy) = (sin45, cos45)

nAz = length(azimuth);
nZe = length(zenith);

% Aim for about 18 evenly spaced ticks along each axis
stepAz = max(1, round(nAz / 18));
stepZe = max(1, round(nZe / 18));

TICK_SCALE = 4.2;   % max half-length of a tick in the polar units (degrees)

az_idx = 1 : stepAz : nAz;
ze_idx = 2 : stepZe : nZe;   % skip index 1 (zenith = 0, degenerate position)

for ia = az_idx
    for iz = ze_idx
        zen = zenith(iz);
        az  = azimuth(ia);
        aop = AoP(iz, ia);
        len = DoLP(iz, ia) * TICK_SCALE;   % half-length

        % Centre of tick in fisheye (x = East, y = North)
        xc =  zen * sind(az);
        yc =  zen * cosd(az);

        % Direction of E-field oscillation:
        %   (dx, dy) = (sin(AoP), cos(AoP))  in (East, North)
        dx = sind(aop) * len;
        dy = cosd(aop) * len;

        plot([xc - dx, xc + dx], [yc - dy, yc + dy], '-', ...
             'Color', [0.08, 0.08, 0.08], 'LineWidth', 0.9);
    end
end

end   % end of drawPolarizationTicks
