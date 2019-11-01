import numpy as np

def find_alpha_and_phi(dpx, dpy):

    phi = np.sqrt(dpx ** 2 + dpy ** 2) / 2.0

    if dpy>=0.:
        if dpx>=0:
            # First quadrant
            if np.abs(dpx) >= np.abs(dpy):
                # First octant
                alpha = 1.
            else:
                # Second octant
                alpha = 2.
        else: #dpx<0
            # Second quadrant
            if np.abs(dpx) <  np.abs(dpy):
                # Third octant
                alpha = 3.
            else:
                # Forth  octant
                alpha = 4.
    else: #dpy<0
        if dpx<=0:
            # Third quadrant
            if np.abs(dpx) >= np.abs(dpy):
                # Fifth octant
                alpha = 5.
            else:
                # Sixth octant
                alpha = 6.
        else: #dpx>0
            # Forth quadrant
            if np.abs(dpx) <= np.abs(dpy):
                # Seventh octant
                alpha = 7.
            else:
                # Eighth octant
                alpha = 8.

#     if phi < 1e-20:
#         alpha = 0.0
# 
# 
#     elif np.abs(dpx) >= np.abs(dpy):
#         alpha = np.arctan(dpy / dpx)
#         if dpx < 0:
#             phi = -phi
#     else:
#         alpha = np.sign(dpy) * (np.pi / 2 - np.abs(np.arctan(dpx / dpy)))
#         if dpy < 0:
#             phi = -phi

    return alpha, phi





XX, YY = meshgrid(linspace(-1., 1, 100), linspace(-1, 1, 101))
AA, PP = np.vectorize(find_alpha_and_phi)(XX, YY)

import matplotlib.pyplot as plt
plt.close('all')

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
ax1.pcolormesh(XX, YY, AA)
plt.colorbar()

plt.show()
