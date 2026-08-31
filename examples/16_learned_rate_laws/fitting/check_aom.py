#!/usr/bin/env python3
"""
check_aom.py -- hold the shipped aom.gnn up against the law it was fitted from.

aom_samples.csv came from a dual-Monod rate law, so unlike a network fitted to real data this one
has a right answer to be compared with.  This draws 300 fresh samples the network never saw and
reports, per species, how well it tracks that law.

    python examples/check_aom.py
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# walk up until we find tools/, so this works wherever the example folder is placed
ROOT = HERE
while not os.path.isdir(os.path.join(ROOT, "tools")) and os.path.dirname(ROOT) != ROOT:
    ROOT = os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(ROOT, "tools"))
import train_graphnet as tg

VMAX, K_CH4, K_SO4 = 4.0e-3, 1.0e-3, 5.0e-4      # mol/L/h, mol/L, mol/L

net, M = tg.read_gnn(os.path.join(HERE, "aom.gnn"))
rng = np.random.default_rng(3)
C = np.column_stack([rng.uniform(1e-5, 5e-3, 300), rng.uniform(1e-5, 8e-3, 300),
                     rng.uniform(0.0, 2e-3, 300),  rng.uniform(1e-4, 4e-3, 300)])
truth = VMAX * C[:, 0]/(K_CH4 + C[:, 0]) * C[:, 1]/(K_SO4 + C[:, 1]) / 3600.0   # per second
pred = tg.predict(net, M, C)

print("aom.gnn against the dual-Monod law it was fitted from, on 300 unseen samples")
worst = 1.0
for k, (nm, sgn) in enumerate([("CH4", -1), ("SO4", -1), ("HS", 1), ("HCO3", 1)]):
    r = np.corrcoef(pred[:, k], sgn * truth)[0, 1]
    worst = min(worst, r)
    print("  %-5s R = %.4f" % (nm, r))
r = np.corrcoef(pred[:, 4], truth)[0, 1]
worst = min(worst, r)
print("  growth R = %.4f" % r)
print("  peak rate %.3e mol/L/s, the law gives %.3e" % (np.abs(pred[:, 0]).max(), truth.max()))
print("  CH4 consumed in %d/300, HS produced in %d/300"
      % (int((pred[:, 0] < 0).sum()), int((pred[:, 2] > 0).sum())))
print("\n%s" % ("the example network reproduces its source law"
                if worst > 0.9 else "** the example network does not track its source law **"))
sys.exit(0 if worst > 0.9 else 1)
