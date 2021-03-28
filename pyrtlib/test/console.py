import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
np.seterr('raise')
import sys
sys.path.append("/Users/slarosa/dev/pyrtlib")
from pyrtlib.atmp import AtmosphericProfiles as atmp
from pyrtlib.main import tb_cloud_rte
from pyrtlib.utils import ppmv2gkg, mr2rh

fig, ax = plt.subplots(1, 1)
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('BT (K)')
atm = ['Tropical',
       'Midlatitude Summer',
       'Midlatitude Winter',
       'Subarctic Summer',
       'Subarctic Winter',
       'U.S. Standard']
# img = plt.imread("/Users/slarosa/dev/pyrtlib/resources/logo/logo_large.png")

for i in range(1, 7):
    z, p, d, t, md, gasids = atmp.gl_atm(atmp.TROPICAL)

    gkg = ppmv2gkg(md[:, 0], gasids[0])
    rh = mr2rh(p, t, gkg)[0] / 100

    ang = np.array([90.])
    frq = np.arange(20, 201, 1)
    nf = len(frq)

    denliq = np.zeros(z.shape)
    denice = np.zeros(z.shape)
    cldh = np.zeros((2, 0))

    df = tb_cloud_rte(z, p, t, rh, denliq, denice, cldh, frq, ang, absmdl='rose19sd', ray_tracing=True, from_sat=True)
    df = df.set_index(frq)
    # np.savetxt('/Users/slarosa/Downloads/adry_{}.txt'.format(atm[i-1].replace(' ', '_')), a.reshape(181, 50))
    # np.savetxt('/Users/slarosa/Downloads/awet_{}.txt'.format(atm[i - 1].replace(' ', '_')), b.reshape(181, 50))

    df.tbtotal.plot(ax=ax, label=atm[i-1])
    df.to_csv("/Users/slarosa/Downloads/out_{}.csv".format(atm[i-1].replace(' ', '_')))


# x0, x1 = plt.xlim()
# y0, y1 = plt.ylim()
# ax.imshow(img, zorder=0, extent=[x0, x1, y0, y1], aspect='equal', alpha=0.1)
ax.grid(which='both')
ax.legend()

plt.title("Model")

fig.savefig("/Users/slarosa/Downloads/spectrum.png")
