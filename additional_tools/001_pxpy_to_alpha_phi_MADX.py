import numpy as np

from cpymad.madx import Madx

madx = Madx()

def find_alpha_and_phi_MADX(deltax,deltay):

    madx.input(f'''
    option,-info;
    dlt_px={deltax};
    dlt_py={deltay};
    ''')

    madx.input('''

    bbbmacros_absphi = sqrt(dlt_px^2 + dlt_py^2) / 2.0;

    if (bbbmacros_absphi < 1e-20)
        {xang_labelNIR_nn = bbbmacros_absphi;xplane_labelNIR_nn = 0.0;}
    else
        {if (dlt_py>=0)
            {if (dlt_px>=0) ! First QUADRANT
                {if (abs(dlt_px) >= abs(dlt_py)) ! First OCTANT
                    {xang_labelNIR_nn=bbbmacros_absphi;
                     xplane_labelNIR_nn = atan(dlt_py/dlt_px);}
                else                             ! Second OCTANT
                    {
                     xang_labelNIR_nn=bbbmacros_absphi;
                     xplane_labelNIR_nn = 0.5*pi - atan(dlt_px/dlt_py);}
                }
            else !dlt_px<0  ! Second QUADRANT
                {if (abs(dlt_px) < abs(dlt_py))  ! Third OCTANT
                    {
                     xang_labelNIR_nn=bbbmacros_absphi;
                     xplane_labelNIR_nn = 0.5*pi - atan(dlt_px/dlt_py);}
                else                             ! Fourth OCTANT
                    {
                     xang_labelNIR_nn=-bbbmacros_absphi;
                     xplane_labelNIR_nn = atan(dlt_py/dlt_px);}
                }
            }
        else !dlt_py<0
            {if (dlt_px<=0) ! Third QUADRANT
                {if (abs(dlt_px) >= abs(dlt_py)) ! Fifth OCTANT
                    {
                     xang_labelNIR_nn=-bbbmacros_absphi;
                     xplane_labelNIR_nn = atan(dlt_py/dlt_px);}
                else                             ! Sixth OCTANT
                    {
                     xang_labelNIR_nn=-bbbmacros_absphi;
                     xplane_labelNIR_nn = 0.5*pi - atan(dlt_px/dlt_py);}
                }
            else !dlt_px>0  ! Fourth QUADRANT
                {if (abs(dlt_px) < abs(dlt_py))  ! Seventh OCTANT
                    {
                     xang_labelNIR_nn=-bbbmacros_absphi;
                     xplane_labelNIR_nn = 0.5*pi - atan(dlt_px/dlt_py);}
                else                             ! Eighth  OCTANT
                    {
                     xang_labelNIR_nn=bbbmacros_absphi;
                     xplane_labelNIR_nn = atan(dlt_py/dlt_px);}
                }
            }
        }
    ''')

    return madx.eval('xplane_labelNIR_nn'),madx.eval('xang_labelNIR_nn')



XX, YY = np.meshgrid(np.linspace(-1., 1, 100), np.linspace(-1, 1, 101))
AA, PP = np.vectorize(find_alpha_and_phi_MADX)(XX, YY)

import matplotlib.pyplot as plt
plt.close('all')

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
mpbl = ax1.pcolormesh(XX, YY, AA)
plt.colorbar(mpbl)
fig1.suptitle('alpha')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
mpbl = ax2.pcolormesh(XX, YY, np.sign(PP))
plt.colorbar(mpbl)
fig2.suptitle('sign(phi)')

plt.show()
