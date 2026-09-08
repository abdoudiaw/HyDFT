"""Plot S(q, w)/n0 from hydft_linear:  python plot_dsf.py out_dsf.dat"""
import sys
import matplotlib.pyplot as plt
from hydft_io import read_dsf

w, S, qs = read_dsf(sys.argv[1])
for j, q in enumerate(qs):
    plt.plot(w, S[:, j]/S[:, j].max(), label=f'q = {q:.2f}')
plt.xlabel('ω/ω_p'); plt.ylabel('S(q,ω) / max'); plt.legend()
plt.tight_layout()
plt.savefig(sys.argv[1].replace('.dat', '.png'), dpi=150)
print('wrote', sys.argv[1].replace('.dat', '.png'))
