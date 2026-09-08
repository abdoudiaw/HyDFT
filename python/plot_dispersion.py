"""Plot dispersion relations produced by hydft_linear and, optionally, overlay
the frequency/damping fitted from a single-mode simulation probe file.

    python plot_dispersion.py out_dispersion.dat [out_probe.dat q]
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from hydft_io import read_dispersion, read_probe, fit_damped_mode

disp = read_dispersion(sys.argv[1])
fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
ax[0].plot(disp['q'], disp['re_w'], label='Re ω (with viscosity, relaxation)')
ax[0].plot(disp['q'], disp['w0'], '--', label='q sqrt(n0 K/m) (undamped)')
ax[0].set_xlabel('q = k a'); ax[0].set_ylabel('ω/ω_p'); ax[0].legend()
ax[1].plot(disp['q'], -disp['im_w'])
ax[1].set_xlabel('q = k a'); ax[1].set_ylabel('-Im ω/ω_p')
if len(sys.argv) >= 4:
    t, dn, _ = read_probe(sys.argv[2])
    q = float(sys.argv[3])
    w, g = fit_damped_mode(t, dn.real, t_start=0.1*t[-1])
    ax[0].plot([q], [w], 'ko', label='simulation'); ax[1].plot([q], [-g], 'ko')
    print(f'simulation: q = {q:.4f}  Re w = {w:.6f}  Im w = {g:.6f}')
plt.tight_layout()
plt.savefig(sys.argv[1].replace('.dat', '.png'), dpi=150)
print('wrote', sys.argv[1].replace('.dat', '.png'))
