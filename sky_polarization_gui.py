"""
Sky Polarization GUI — Interactive 2D fisheye and 3D hemisphere visualization.

Provides a tkinter-based GUI to configure observer location and time, compute
Rayleigh scattering polarization maps.

Usage
-----
    python sky_polarization_gui.py
"""

import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from sky_polarization import sky_polarization


# ---------------------------------------------------------------------------
#  MATLAB parula-like colormap (256-sample approximation)
# ---------------------------------------------------------------------------
_PARULA_DATA = [
    (0.2422, 0.1504, 0.6603), (0.2504, 0.1650, 0.7076),
    (0.2578, 0.1818, 0.7511), (0.2647, 0.1978, 0.7952),
    (0.2706, 0.2147, 0.8364), (0.2751, 0.2342, 0.8710),
    (0.2783, 0.2559, 0.8991), (0.2803, 0.2782, 0.9221),
    (0.2813, 0.3006, 0.9414), (0.2810, 0.3228, 0.9579),
    (0.2795, 0.3447, 0.9717), (0.2760, 0.3667, 0.9829),
    (0.2699, 0.3892, 0.9906), (0.2602, 0.4123, 0.9952),
    (0.2440, 0.4358, 0.9988), (0.2206, 0.4603, 0.9973),
    (0.1963, 0.4847, 0.9892), (0.1834, 0.5074, 0.9798),
    (0.1786, 0.5289, 0.9698), (0.1764, 0.5499, 0.9594),
    (0.1687, 0.5703, 0.9492), (0.1540, 0.5902, 0.9385),
    (0.1460, 0.6091, 0.9259), (0.1380, 0.6276, 0.9125),
    (0.1248, 0.6459, 0.8979), (0.1113, 0.6635, 0.8828),
    (0.0952, 0.6798, 0.8668), (0.0689, 0.6948, 0.8492),
    (0.0297, 0.7082, 0.8305), (0.0036, 0.7197, 0.8106),
    (0.0067, 0.7293, 0.7886), (0.0433, 0.7374, 0.7646),
    (0.0964, 0.7441, 0.7394), (0.1408, 0.7493, 0.7139),
    (0.1717, 0.7532, 0.6891), (0.1938, 0.7560, 0.6651),
    (0.2161, 0.7575, 0.6395), (0.2470, 0.7574, 0.6105),
    (0.2906, 0.7556, 0.5780), (0.3406, 0.7522, 0.5432),
    (0.3909, 0.7473, 0.5078), (0.4378, 0.7410, 0.4733),
    (0.4813, 0.7334, 0.4405), (0.5221, 0.7248, 0.4087),
    (0.5612, 0.7152, 0.3779), (0.5992, 0.7047, 0.3480),
    (0.6363, 0.6932, 0.3190), (0.6727, 0.6811, 0.2905),
    (0.7082, 0.6682, 0.2629), (0.7426, 0.6546, 0.2364),
    (0.7752, 0.6405, 0.2117), (0.8065, 0.6256, 0.1891),
    (0.8361, 0.6100, 0.1695), (0.8637, 0.5934, 0.1543),
    (0.8884, 0.5760, 0.1450), (0.9096, 0.5577, 0.1440),
    (0.9261, 0.5385, 0.1527), (0.9380, 0.5185, 0.1680),
    (0.9449, 0.4975, 0.1877), (0.9471, 0.4755, 0.2094),
    (0.9450, 0.4528, 0.2310), (0.9394, 0.4297, 0.2509),
    (0.9310, 0.4063, 0.2678), (0.9765, 0.9831, 0.0538),
]

def _make_parula():
    """Build a parula-like colormap by interpolation."""
    # Use a well-known 64-stop parula approximation
    colors_64 = [
        [0.2081, 0.1663, 0.5292], [0.2116, 0.1898, 0.5777],
        [0.2123, 0.2138, 0.6270], [0.2081, 0.2386, 0.6771],
        [0.1959, 0.2645, 0.7279], [0.1707, 0.2919, 0.7792],
        [0.1253, 0.3242, 0.8303], [0.0591, 0.3598, 0.8683],
        [0.0117, 0.3875, 0.8820], [0.0060, 0.4086, 0.8828],
        [0.0165, 0.4266, 0.8786], [0.0329, 0.4430, 0.8720],
        [0.0498, 0.4586, 0.8641], [0.0629, 0.4737, 0.8554],
        [0.0723, 0.4887, 0.8467], [0.0779, 0.5040, 0.8384],
        [0.0793, 0.5200, 0.8312], [0.0749, 0.5375, 0.8263],
        [0.0641, 0.5570, 0.8240], [0.0488, 0.5772, 0.8228],
        [0.0343, 0.5966, 0.8199], [0.0265, 0.6137, 0.8135],
        [0.0239, 0.6287, 0.8038], [0.0231, 0.6418, 0.7913],
        [0.0228, 0.6535, 0.7768], [0.0267, 0.6642, 0.7607],
        [0.0384, 0.6743, 0.7436], [0.0590, 0.6838, 0.7254],
        [0.0843, 0.6928, 0.7062], [0.1133, 0.7015, 0.6859],
        [0.1453, 0.7098, 0.6646], [0.1801, 0.7177, 0.6424],
        [0.2178, 0.7250, 0.6193], [0.2586, 0.7317, 0.5954],
        [0.3022, 0.7376, 0.5712], [0.3482, 0.7424, 0.5473],
        [0.3953, 0.7459, 0.5244], [0.4420, 0.7481, 0.5033],
        [0.4871, 0.7491, 0.4840], [0.5300, 0.7491, 0.4661],
        [0.5709, 0.7485, 0.4494], [0.6099, 0.7473, 0.4337],
        [0.6473, 0.7456, 0.4188], [0.6834, 0.7435, 0.4044],
        [0.7184, 0.7411, 0.3905], [0.7525, 0.7384, 0.3768],
        [0.7858, 0.7356, 0.3633], [0.8185, 0.7327, 0.3498],
        [0.8507, 0.7299, 0.3360], [0.8824, 0.7274, 0.3217],
        [0.9139, 0.7258, 0.3063], [0.9450, 0.7261, 0.2886],
        [0.9739, 0.7314, 0.2634], [0.9938, 0.7455, 0.2403],
        [0.9990, 0.7653, 0.2164], [0.9955, 0.7861, 0.1967],
        [0.9880, 0.8066, 0.1794], [0.9789, 0.8271, 0.1633],
        [0.9697, 0.8481, 0.1475], [0.9626, 0.8705, 0.1309],
        [0.9589, 0.8949, 0.1132], [0.9598, 0.9218, 0.0948],
        [0.9661, 0.9514, 0.0755], [0.9763, 0.9831, 0.0538],
    ]
    return LinearSegmentedColormap.from_list('parula', colors_64, N=512)

PARULA = _make_parula()


# ---------------------------------------------------------------------------
#  Plotting helpers
# ---------------------------------------------------------------------------

def _draw_polarization_ticks(ax, azimuth, zenith, AoP, DoLP):
    """Overlay E-field direction tick marks on the AoP panel."""
    nAz = len(azimuth)
    nZe = len(zenith)
    stepAz = max(1, round(nAz / 18))
    stepZe = max(1, round(nZe / 18))
    TICK_SCALE = 4.2

    az_idx = range(0, nAz, stepAz)
    ze_idx = range(1, nZe, stepZe)  # skip zenith=0

    for ia in az_idx:
        for iz in ze_idx:
            zen = zenith[iz]
            az = azimuth[ia]
            aop = AoP[iz, ia]
            length = DoLP[iz, ia] * TICK_SCALE

            xc = zen * np.sin(np.radians(az))
            yc = zen * np.cos(np.radians(az))

            dx = np.sin(np.radians(aop)) * length
            dy = np.cos(np.radians(aop)) * length

            ax.plot([xc - dx, xc + dx], [yc - dy, yc + dy], '-',
                    color=(0.08, 0.08, 0.08), linewidth=0.9)


def plot_2d_fisheye(fig, result):
    """Draw the 2D fisheye polar projection (3 panels) on *fig*."""
    fig.clear()

    DoP = np.sqrt(result['Q'] ** 2 + result['U'] ** 2)
    DoLP = result['DoLP']
    AoP = result['AoP']

    AZ_grid, ZE_grid = np.meshgrid(result['azimuth'], result['zenith'])
    x2d = ZE_grid * np.sin(np.radians(AZ_grid))
    y2d = ZE_grid * np.cos(np.radians(AZ_grid))

    sun = result['sun']
    sun_above = sun['elevation'] > 0
    sunX2d = sun['zenith'] * np.sin(np.radians(sun['azimuth']))
    sunY2d = sun['zenith'] * np.cos(np.radians(sun['azimuth']))

    obs = result['observer']
    tim = result['time']
    lat_ch = 'N' if obs['latitude'] >= 0 else 'S'
    lon_ch = 'E' if obs['longitude'] >= 0 else 'W'
    main_title = (
        f"{tim['year']:04d}-{tim['month']:02d}-{tim['day']:02d}   "
        f"{tim['hour']:05.2f} h local  |  "
        f"{abs(obs['latitude']):.2f}\u00b0{lat_ch}  "
        f"{abs(obs['longitude']):.2f}\u00b0{lon_ch}  |  "
        f"Sun: zen={sun['zenith']:.1f}\u00b0  "
        f"az={sun['azimuth']:.1f}\u00b0  "
        f"elv={sun['elevation']:.1f}\u00b0"
    )

    panels = [
        (DoLP, (0, 1), PARULA, 'DoLP', 'Degree of Linear Polarization', False),
        (AoP, (0, 180), 'hsv', 'AoP  (deg)', 'Angle of Polarization', True),
        (DoP, (0, 1), PARULA, 'DoP', 'Degree of Polarization', False),
    ]

    th_ring = np.linspace(0, 360, 361)

    for k, (data, clim, cmap, cb_label, pan_title, show_ticks) in enumerate(panels):
        ax = fig.add_subplot(1, 3, k + 1)

        cf = ax.contourf(x2d, y2d, data, levels=256, vmin=clim[0], vmax=clim[1])
        cf.set_cmap(cmap if isinstance(cmap, str) else cmap)
        ax.set_facecolor('white')

        # Colorbar below each panel
        cb = fig.colorbar(cf, ax=ax, orientation='horizontal', fraction=0.06, pad=0.08)
        cb.set_label(cb_label, fontsize=9)

        # Reference rings at 30 and 60 degrees
        for zring in [30, 60]:
            ax.plot(zring * np.sin(np.radians(th_ring)),
                    zring * np.cos(np.radians(th_ring)),
                    '--', color=(0.55, 0.55, 0.55), linewidth=0.7)
            ax.text(0, -zring - 1.5, f'{zring}\u00b0',
                    fontsize=7, color=(0.45, 0.45, 0.45), ha='center')

        # Horizon ring
        ax.plot(90 * np.sin(np.radians(th_ring)),
                90 * np.cos(np.radians(th_ring)),
                'k-', linewidth=1.8)

        # Cardinal direction labels
        ax.text(0, 98, 'N', ha='center', fontweight='bold', fontsize=11)
        ax.text(0, -98, 'S', ha='center', fontweight='bold', fontsize=11)
        ax.text(98, 0, 'E', ha='center', fontweight='bold', fontsize=11)
        ax.text(-98, 0, 'W', ha='center', fontweight='bold', fontsize=11)

        # AoP tick marks
        if show_ticks:
            _draw_polarization_ticks(ax, result['azimuth'], result['zenith'], AoP, DoLP)

        # Sun marker
        if sun_above:
            ax.plot(sunX2d, sunY2d, 'o', markersize=14,
                    markerfacecolor=(1.0, 0.92, 0.0),
                    markeredgecolor=(0.55, 0.40, 0.0), markeredgewidth=2)
            ax.text(sunX2d + 5, sunY2d + 5, 'Sun', fontsize=8, color=(0.3, 0.3, 0.3))
        else:
            ax.plot(sunX2d, sunY2d, 'v', markersize=10,
                    markerfacecolor=(0.5, 0.5, 0.5),
                    markeredgecolor='k', markeredgewidth=1.2)

        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim(-107, 107)
        ax.set_ylim(-107, 107)
        ax.set_title(pan_title, fontsize=12, fontweight='bold')

    fig.suptitle(main_title, fontsize=9.5, fontweight='normal')
    fig.tight_layout(rect=[0, 0, 1, 0.95])


def plot_3d_hemisphere(fig, result):
    """Draw the 3D hemisphere projection (3 panels) on *fig*."""
    fig.clear()

    DoP = np.sqrt(result['Q'] ** 2 + result['U'] ** 2)
    DoLP = result['DoLP']
    AoP = result['AoP']

    AZ_grid, ZE_grid = np.meshgrid(result['azimuth'], result['zenith'])

    X3d = np.sin(np.radians(ZE_grid)) * np.sin(np.radians(AZ_grid))
    Y3d = np.sin(np.radians(ZE_grid)) * np.cos(np.radians(AZ_grid))
    Z3d = np.cos(np.radians(ZE_grid))

    sun = result['sun']
    sun_above = sun['elevation'] > 0
    sunX3d = np.sin(np.radians(sun['zenith'])) * np.sin(np.radians(sun['azimuth']))
    sunY3d = np.sin(np.radians(sun['zenith'])) * np.cos(np.radians(sun['azimuth']))
    sunZ3d = np.cos(np.radians(sun['zenith']))

    obs = result['observer']
    tim = result['time']
    lat_ch = 'N' if obs['latitude'] >= 0 else 'S'
    lon_ch = 'E' if obs['longitude'] >= 0 else 'W'
    main_title = (
        f"{tim['year']:04d}-{tim['month']:02d}-{tim['day']:02d}   "
        f"{tim['hour']:05.2f} h local  |  "
        f"{abs(obs['latitude']):.2f}\u00b0{lat_ch}  "
        f"{abs(obs['longitude']):.2f}\u00b0{lon_ch}  |  "
        f"Sun: zen={sun['zenith']:.1f}\u00b0  "
        f"az={sun['azimuth']:.1f}\u00b0  "
        f"elv={sun['elevation']:.1f}\u00b0"
    )

    panels = [
        (DoLP, (0, 1), PARULA, 'DoLP', 'Degree of Linear Polarization'),
        (AoP, (0, 180), 'hsv', 'AoP  (deg)', 'Angle of Polarization'),
        (DoP, (0, 1), PARULA, 'DoP', 'Degree of Polarization'),
    ]

    az_grid_line = np.linspace(0, 360, 361)
    th_zen = np.linspace(0, 90, 91)

    for k, (data, clim, cmap_name, cb_label, pan_title) in enumerate(panels):
        ax = fig.add_subplot(1, 3, k + 1, projection='3d')

        # Normalize data for color mapping
        if isinstance(cmap_name, str):
            cmap_obj = cm.get_cmap(cmap_name, 512)
        else:
            cmap_obj = cmap_name

        norm_data = (data - clim[0]) / (clim[1] - clim[0])
        norm_data = np.clip(norm_data, 0, 1)
        facecolors = cmap_obj(norm_data)

        ax.plot_surface(X3d, Y3d, Z3d, facecolors=facecolors,
                        rstride=1, cstride=1, shade=False, antialiased=False)

        # Ground disc (filled polygon)
        th_g = np.linspace(0, 2 * np.pi, 361)
        verts = [list(zip(np.cos(th_g), np.sin(th_g), np.zeros(361)))]
        ground = Poly3DCollection(verts, facecolors=(0.88, 0.93, 0.88),
                                  alpha=0.5, edgecolors='none')
        ax.add_collection3d(ground)

        # Horizon ring
        ax.plot(np.cos(th_g), np.sin(th_g), np.zeros(361), 'k-', linewidth=1.5)

        # Elevation reference rings at 30 and 60 degrees
        for elv_ring in [30, 60]:
            zen_ring = 90 - elv_ring
            xr = np.sin(np.radians(zen_ring)) * np.sin(np.radians(az_grid_line))
            yr = np.sin(np.radians(zen_ring)) * np.cos(np.radians(az_grid_line))
            zr = np.cos(np.radians(zen_ring)) * np.ones_like(az_grid_line)
            ax.plot(xr, yr, zr, '--', color=(0.6, 0.6, 0.6), linewidth=0.8)

        # Cardinal meridian lines
        for az_line in [0, 90, 180, 270]:
            xl = np.sin(np.radians(th_zen)) * np.sin(np.radians(az_line))
            yl = np.sin(np.radians(th_zen)) * np.cos(np.radians(az_line))
            zl = np.cos(np.radians(th_zen))
            ax.plot(xl, yl, zl, '-', color=(0.6, 0.6, 0.6), linewidth=0.6)

        # Cardinal direction labels
        ax.text(0, 1.13, 0, 'N', fontweight='bold', fontsize=10, ha='center')
        ax.text(0, -1.13, 0, 'S', fontweight='bold', fontsize=10, ha='center')
        ax.text(1.13, 0, 0, 'E', fontweight='bold', fontsize=10, ha='center')
        ax.text(-1.13, 0, 0, 'W', fontweight='bold', fontsize=10, ha='center')

        # Elevation ring labels
        ax.text(np.sin(np.radians(60)) * np.sin(np.radians(45)),
                np.sin(np.radians(60)) * np.cos(np.radians(45)),
                np.cos(np.radians(60)) + 0.03,
                '30\u00b0 elv', fontsize=7, color=(0.4, 0.4, 0.4))
        ax.text(np.sin(np.radians(30)) * np.sin(np.radians(45)),
                np.sin(np.radians(30)) * np.cos(np.radians(45)),
                np.cos(np.radians(30)) + 0.03,
                '60\u00b0 elv', fontsize=7, color=(0.4, 0.4, 0.4))

        # Sun marker
        if sun_above:
            offset = 0.04
            ax.plot([sunX3d * (1 + offset)], [sunY3d * (1 + offset)],
                    [sunZ3d * (1 + offset)], 'o', markersize=14,
                    markerfacecolor=(1.0, 0.92, 0.0),
                    markeredgecolor=(0.55, 0.40, 0.0), markeredgewidth=2.5)
            ax.plot([sunX3d, sunX3d], [sunY3d, sunY3d], [0, sunZ3d],
                    ':', color=(0.7, 0.6, 0.1), linewidth=1.2)

        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_zlim(-0.05, 1.15)
        ax.view_init(elev=30, azim=-30)
        ax.axis('off')
        ax.set_title(pan_title, fontsize=12, fontweight='bold')

        # Add a ScalarMappable for colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap_obj,
                                   norm=plt.Normalize(vmin=clim[0], vmax=clim[1]))
        sm.set_array([])
        cb = fig.colorbar(sm, ax=ax, orientation='horizontal',
                          fraction=0.06, pad=0.05, shrink=0.8)
        cb.set_label(cb_label, fontsize=9)

    fig.suptitle(main_title, fontsize=9.5, fontweight='normal')
    fig.tight_layout(rect=[0, 0, 1, 0.95])


# ---------------------------------------------------------------------------
#  GUI Application
# ---------------------------------------------------------------------------

class SkyPolarizationApp:
    def __init__(self, root):
        self.root = root
        self.root.title('Sky Polarization Simulator')
        self.root.state('zoomed')  # maximize on Windows

        # ---- Control panel (left) ------------------------------------------
        ctrl_frame = ttk.LabelFrame(root, text='Parameters', padding=10)
        ctrl_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)

        # Observer
        ttk.Label(ctrl_frame, text='Observer Location',
                  font=('Segoe UI', 10, 'bold')).grid(
            row=0, column=0, columnspan=2, sticky='w', pady=(0, 5))

        fields_obs = [
            ('Latitude (\u00b0N):', '33.42'),
            ('Longitude (\u00b0E):', '-111.93'),
            ('Altitude (m):', '340'),
            ('GMT Offset (h):', '-7'),
        ]
        self.obs_vars = {}
        for i, (label, default) in enumerate(fields_obs):
            ttk.Label(ctrl_frame, text=label).grid(row=i + 1, column=0, sticky='e', padx=(0, 4))
            var = tk.StringVar(value=default)
            ttk.Entry(ctrl_frame, textvariable=var, width=12).grid(row=i + 1, column=1, sticky='w')
            self.obs_vars[label] = var

        # Time
        ttk.Label(ctrl_frame, text='Observation Time',
                  font=('Segoe UI', 10, 'bold')).grid(
            row=6, column=0, columnspan=2, sticky='w', pady=(12, 5))

        fields_time = [
            ('Year:', '2024'),
            ('Month:', '6'),
            ('Day:', '21'),
            ('Hour (decimal):', '10.0'),
            ('DST (0 or 1):', '1'),
        ]
        self.time_vars = {}
        for i, (label, default) in enumerate(fields_time):
            ttk.Label(ctrl_frame, text=label).grid(row=i + 7, column=0, sticky='e', padx=(0, 4))
            var = tk.StringVar(value=default)
            ttk.Entry(ctrl_frame, textvariable=var, width=12).grid(row=i + 7, column=1, sticky='w')
            self.time_vars[label] = var

        # Options
        ttk.Label(ctrl_frame, text='Options',
                  font=('Segoe UI', 10, 'bold')).grid(
            row=13, column=0, columnspan=2, sticky='w', pady=(12, 5))

        ttk.Label(ctrl_frame, text='Resolution (\u00b0):').grid(row=14, column=0, sticky='e', padx=(0, 4))
        self.res_var = tk.StringVar(value='1.5')
        ttk.Entry(ctrl_frame, textvariable=self.res_var, width=12).grid(row=14, column=1, sticky='w')

        ttk.Label(ctrl_frame, text='Wavelength (nm):').grid(row=15, column=0, sticky='e', padx=(0, 4))
        self.wl_var = tk.StringVar(value='550')
        ttk.Entry(ctrl_frame, textvariable=self.wl_var, width=12).grid(row=15, column=1, sticky='w')

        # Buttons
        btn_frame = ttk.Frame(ctrl_frame)
        btn_frame.grid(row=16, column=0, columnspan=2, pady=(15, 5))

        self.compute_btn = ttk.Button(btn_frame, text='Compute & Plot',
                                      command=self.compute_and_plot)
        self.compute_btn.pack(fill=tk.X, pady=2)

        # Status
        self.status_var = tk.StringVar(value='Ready.')
        ttk.Label(ctrl_frame, textvariable=self.status_var,
                  foreground='gray').grid(
            row=17, column=0, columnspan=2, sticky='w', pady=(10, 0))

        # Sun info display
        self.sun_info_var = tk.StringVar(value='')
        ttk.Label(ctrl_frame, textvariable=self.sun_info_var,
                  foreground='#336', wraplength=200, justify='left').grid(
            row=18, column=0, columnspan=2, sticky='w', pady=(5, 0))

        # ---- Plot area (right) — notebook with tabs -----------------------
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Tab 1: 2D Fisheye
        tab_2d = ttk.Frame(self.notebook)
        self.notebook.add(tab_2d, text='  2D Fisheye Polar  ')

        self.fig_2d = plt.Figure(figsize=(14, 5.2), dpi=100, facecolor='white')
        self.canvas_2d = FigureCanvasTkAgg(self.fig_2d, master=tab_2d)
        toolbar_2d = NavigationToolbar2Tk(self.canvas_2d, tab_2d)
        toolbar_2d.update()
        self.canvas_2d.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab 2: 3D Hemisphere
        tab_3d = ttk.Frame(self.notebook)
        self.notebook.add(tab_3d, text='  3D Hemisphere  ')

        self.fig_3d = plt.Figure(figsize=(14, 6), dpi=100, facecolor='white')
        self.canvas_3d = FigureCanvasTkAgg(self.fig_3d, master=tab_3d)
        toolbar_3d = NavigationToolbar2Tk(self.canvas_3d, tab_3d)
        toolbar_3d.update()
        self.canvas_3d.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Auto-compute with default parameters on startup
        self.root.after(100, self.compute_and_plot)

    def _read_params(self):
        """Read GUI fields and return observer/time dicts and options."""
        observer = {
            'latitude': float(self.obs_vars['Latitude (\u00b0N):'].get()),
            'longitude': float(self.obs_vars['Longitude (\u00b0E):'].get()),
            'altitude': float(self.obs_vars['Altitude (m):'].get()),
            'gmtOffset': float(self.obs_vars['GMT Offset (h):'].get()),
        }
        time_info = {
            'year': int(self.time_vars['Year:'].get()),
            'month': int(self.time_vars['Month:'].get()),
            'day': int(self.time_vars['Day:'].get()),
            'hour': float(self.time_vars['Hour (decimal):'].get()),
            'dst': int(self.time_vars['DST (0 or 1):'].get()),
        }
        resolution = float(self.res_var.get())
        wavelength = float(self.wl_var.get()) * 1e-9  # nm -> m
        return observer, time_info, resolution, wavelength

    def compute_and_plot(self):
        """Run computation and update both plot tabs."""
        self.status_var.set('Computing...')
        self.compute_btn.config(state='disabled')
        self.root.update_idletasks()

        try:
            observer, time_info, resolution, wavelength = self._read_params()
            result = sky_polarization(observer, time_info,
                                      resolution=resolution,
                                      wavelength=wavelength)

            sun = result['sun']
            self.sun_info_var.set(
                f"Sun position:\n"
                f"  Zenith:    {sun['zenith']:.2f}\u00b0\n"
                f"  Azimuth:   {sun['azimuth']:.2f}\u00b0\n"
                f"  Elevation: {sun['elevation']:.2f}\u00b0\n\n"
                f"DoLP max: {result['DoLP'].max():.4f}\n"
                f"DoLP min: {result['DoLP'].min():.4f}")

            # Plot 2D
            plot_2d_fisheye(self.fig_2d, result)
            self.canvas_2d.draw()

            # Plot 3D
            plot_3d_hemisphere(self.fig_3d, result)
            self.canvas_3d.draw()

            self.status_var.set('Done.')

        except Exception as e:
            self.status_var.set(f'Error: {e}')
            self.sun_info_var.set('')
            raise

        finally:
            self.compute_btn.config(state='normal')


def main():
    root = tk.Tk()
    app = SkyPolarizationApp(root)
    root.mainloop()


if __name__ == '__main__':
    main()
