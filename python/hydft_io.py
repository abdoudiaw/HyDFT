"""Readers for HyDFT output files.

    diag   = read_diag('out_diag.dat')        # dict of columns
    probe  = read_probe('out_probe.dat')      # t, complex dn_k/n0, complex mom_k
    fld    = read_fields('out_fld_000100.bin')  # dict with n, mom, pi, grid info
    disp   = read_dispersion('out_dispersion.dat')
    w, S   = read_dsf('out_dsf.dat')          # S has one column per q
"""
import io
import re
import numpy as np


def _load(path):
    # Fortran prints three-digit exponents without the 'E' (1.0-100); restore it.
    with open(path) as f:
        text = re.sub(r'(\d)([+-]\d{3})\b', r'\1E\2', f.read())
    return np.loadtxt(io.StringIO(text), comments='#')


def read_diag(path):
    d = np.atleast_2d(_load(path))
    names = ['step', 't', 'dt', 'mass', 'momx', 'momy', 'momz', 'ekin', 'free_energy', 'etot', 'nmin', 'nmax']
    return {k: d[:, i] for i, k in enumerate(names)}


def read_probe(path):
    d = np.atleast_2d(_load(path))
    return d[:, 0], d[:, 1] + 1j*d[:, 2], d[:, 4] + 1j*d[:, 5]


def read_dispersion(path):
    d = np.atleast_2d(_load(path))
    return {'q': d[:, 0], 're_w': d[:, 1], 'im_w': d[:, 2], 'w0': d[:, 3], 'S': d[:, 4]}


def read_dsf(path):
    with open(path) as f:
        header = f.readline()
    qs = [float(s.split()[0]) for s in header.split('q=')[1:]]
    d = np.atleast_2d(_load(path))
    return d[:, 0], d[:, 1:], qs


def read_fields(path):
    with open(path, 'rb') as f:
        nx, ny, nz, ndim, npi = np.fromfile(f, dtype=np.int32, count=5)
        t, lx, ly, lz = np.fromfile(f, dtype=np.float64, count=4)
        n = np.fromfile(f, dtype=np.float64, count=nx*ny*nz).reshape((nz, ny, nx)).transpose()
        mom = np.fromfile(f, dtype=np.float64, count=3*nx*ny*nz).reshape((3, nz, ny, nx)).transpose()
        pi = None
        if npi > 0:
            pi = np.fromfile(f, dtype=np.float64, count=npi*nx*ny*nz).reshape((npi, nz, ny, nx)).transpose()
    x = np.arange(nx)*lx/nx
    return {'t': t, 'nx': nx, 'ny': ny, 'nz': nz, 'ndim': ndim, 'lx': lx, 'ly': ly, 'lz': lz,
            'x': x, 'n': n, 'mom': mom, 'pi': pi}


def fit_damped_mode(t, y, t_start=0.0):
    """Frequency and damping of y(t) ~ A exp(g t) cos(w t + phi), t >= t_start
    (zero crossings + peak envelope; same method as the Fortran tests)."""
    m = t >= t_start
    t, y = t[m], y[m]
    s = np.sign(y)
    idx = np.where(s[:-1]*s[1:] < 0)[0]
    tz = t[idx] + y[idx]/(y[idx] - y[idx+1])*(t[idx+1] - t[idx])
    if len(tz) < 3:
        return np.nan, np.nan
    w = np.pi*(len(tz) - 1)/(tz[-1] - tz[0])
    tp, lp = [], []
    for a, b in zip(tz[:-1], tz[1:]):
        mm = (t > a) & (t < b)
        if mm.any():
            j = np.argmax(np.abs(y[mm]))
            tp.append(t[mm][j]); lp.append(np.log(np.abs(y[mm][j])))
    g = np.polyfit(tp, lp, 1)[0]
    return w, g
