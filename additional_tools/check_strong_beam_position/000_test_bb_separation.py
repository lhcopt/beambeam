import numpy as np
import pysixtrack
import helpers as hp
import matplotlib.pyplot as plt
import pdb;
import textwrap
#pdb.set_trace()

plt.close('all')

def write_fort3(filename='./fort.3', partnum=1.e11, emitnx_um=5., emitny_um=5., sigz_m= 0.00001,\
        sige=0., ibeco=1, ibsix=2, sigma_xx_mm2=1.,sigma_yy_mm2=1.,xang_rad=0.,xplane_rad=0.,h_sep_mm=0.,v_sep_mm=0.):
        '''
        Write of the fort.3 file. Look at https://sixtrack.web.cern.ch/SixTrack/docs/user_manual.pdf for the variable meanings (Section 6.6, November 2019).
        '''
        myString=f'''\
                GEOME-STRENG TITLE:simple_cell
                PRINTOUT OF INPUT PARAMETERS--------------------------------------------
                NEXT
                TRACKING PARAMETERS-----------------------------------------------------
                1 0 1 0.00 0.00 0 1
                1 1 0 1 2
                0 0 1 1 1000 50000 2
                NEXT
                INITIAL COORDINATES-----------------------------------------------------
                        2 0. 0. 0.999999 0
                        0.
                        0.
                        0. 
                        0.
                        0.
                        0.
                        0.0
                        0.0
                        0.0 
                        0.0
                        0.0
                        0.0002
                        7000000.
                        7000000.
                        7000000.
                NEXT
                SYNC
                        1  0.030031     5.000 0.     20.000000    938.272081 1
                1.        1.
                NEXT
                DUMP
                ALL 1 664 101 dump3.dat
                NEXT
                LINEAR OPTICS
                ELEMENT  0 1 1 2.5 2.5
                NEXT
                BEAM
                EXPERT
                {partnum}  {emitnx_um}  {emitny_um}  {sigz_m}  {sige}  {ibeco} 0 0 0 
                bb  {ibsix}  {xang_rad} {xplane_rad} {h_sep_mm}  {v_sep_mm}  
                {sigma_xx_mm2}  0.  1.  {sigma_yy_mm2}  0.
                1.  0. 0. 0. 0.  1.
                bb4d  0  {sigma_xx_mm2}  {sigma_yy_mm2}  {h_sep_mm}  {v_sep_mm}  1.  0.
                NEXT
                ITERATION-ACCURACY------------------------------------------------------
                50 1D-12 1D-15
                10 1D-10 1D-10
                10 1D-6  1D-6
                        1D-6  1D-9  1D-9
                NEXT
                ENDE====================================================================
                '''
        with open(filename,'w+') as fid: 
                fid.write(textwrap.dedent(myString))
        return

def BBKick(x,Np,gamma_r, sigma):
    r_p=1.534698e-18
    if x==0:
        return(0)
    else:
        return 2*r_p*Np/gamma_r/x*(1-np.exp(-x**2/(2*sigma**2)))
BBKickVector=np.vectorize(BBKick)

# Check with no separation
if 0:
        my_partnum=1.e11
        my_xang_rad=0.
        my_xplane_rad=0.
        my_h_sep_mm=1.e-14 # for the 4D bb needs to be non zero (otherwise nan)
        my_v_sep_mm=0.
        my_sigma_xx_mm2=1.
        write_fort3(filename='./fort.3',partnum=my_partnum, xang_rad=my_xang_rad, xplane_rad=my_xplane_rad, h_sep_mm=my_h_sep_mm, v_sep_mm=my_v_sep_mm, sigma_xx_mm2=my_sigma_xx_mm2)

        part_on_CO = pysixtrack.Particles(p0c=7e12, x=0.)
        
        # along the x-axis        

        x_test= np.linspace(-2e-2, 2e-2, 100)
        y_test = 0*np.linspace(-2e-2, 2e-2, 100)

        x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
                Dx_wrt_CO_m=x_test,
                Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
                Dsigma_wrt_CO_m=0, Ddelta_wrt_CO=0, n_turns=2,
                mode = 'ebe')

        plt.figure()

        plt.subplot(221)
        plt.plot(x_test, px_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.title('The 6D kick')
        plt.plot(x_test,BBKickVector(x_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(222)
        plt.plot(x_test, px_tbt[3, :]-px_tbt[2, :],label='SixTrack results')
        plt.plot(x_test,BBKickVector(x_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.title('The 4D kick')
        plt.legend(loc='best')

        plt.subplot(223)
        plt.plot(x_test, py_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('py')
        plt.plot(x_test, 0*x_test,'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(224)
        plt.plot(x_test, py_tbt[3, :]-py_tbt[2, :],label='SixTrack results')
        plt.plot(x_test, 0*x_test,'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.legend(loc='best')

        # along the y-axis        
        plt.figure()

        x_test= 0*np.linspace(-2e-2, 2e-2, 100)
        y_test = np.linspace(-2e-2, 2e-2, 100)

        x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
                Dx_wrt_CO_m=x_test,
                Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
                Dsigma_wrt_CO_m=0, Ddelta_wrt_CO=0, n_turns=2,
                mode = 'ebe')

        plt.subplot(221)
        plt.plot(y_test, px_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('px')
        plt.title('The 6D kick')
        plt.plot(y_test,0*y_test,'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(222)
        plt.plot(y_test, px_tbt[3, :]-px_tbt[2, :],label='SixTrack results')
        plt.plot(y_test,0*y_test,'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('px')
        plt.title('The 4D kick')
        plt.legend(loc='best')

        plt.subplot(223)
        plt.plot(y_test, py_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('py')
        plt.plot(y_test, BBKickVector(y_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(224)
        plt.plot(y_test, py_tbt[3, :]-py_tbt[2, :],label='SixTrack results')
        plt.plot(y_test, BBKickVector(y_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('px')
        plt.legend(loc='best')

# Check with h-separation
if 1:
        my_partnum=1.e11
        my_xang_rad=0.
        my_xplane_rad=0.
        my_h_sep_mm=1. # for the 4D bb needs to be non zero (otherwise nan)
        my_v_sep_mm=0.
        my_sigma_xx_mm2=1.
        write_fort3(filename='./fort.3',partnum=my_partnum, xang_rad=my_xang_rad, xplane_rad=my_xplane_rad, h_sep_mm=my_h_sep_mm, v_sep_mm=my_v_sep_mm, sigma_xx_mm2=my_sigma_xx_mm2)

        part_on_CO = pysixtrack.Particles(p0c=7e12, x=0.)
        
        # along the x-axis        

        x_test= np.linspace(-2e-2, 2e-2, 100)
        y_test = 0*np.linspace(-2e-2, 2e-2, 100)

        x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
                Dx_wrt_CO_m=x_test,
                Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
                Dsigma_wrt_CO_m=0, Ddelta_wrt_CO=0, n_turns=2,
                mode = 'ebe')

        plt.figure()

        plt.subplot(221)
        plt.plot(x_test, px_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.title('The 6D kick')
        plt.plot(x_test,BBKickVector(x_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(222)
        plt.plot(x_test, px_tbt[3, :]-px_tbt[2, :],label='SixTrack results')
        plt.plot(x_test,BBKickVector(x_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.title('The 4D kick')
        plt.legend(loc='best')

        plt.subplot(223)
        plt.plot(x_test, py_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('py')
        plt.plot(x_test, 0*x_test,'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(224)
        plt.plot(x_test, py_tbt[3, :]-py_tbt[2, :],label='SixTrack results')
        plt.plot(x_test, 0*x_test,'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.legend(loc='best')

        # along the y-axis        
        plt.figure()

        x_test= 0*np.linspace(-2e-2, 2e-2, 100)
        y_test = np.linspace(-2e-2, 2e-2, 100)

        x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
                Dx_wrt_CO_m=x_test,
                Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
                Dsigma_wrt_CO_m=0, Ddelta_wrt_CO=0, n_turns=2,
                mode = 'ebe')

        plt.subplot(221)
        plt.plot(y_test, px_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('px')
        plt.title('The 6D kick')
        plt.plot(y_test,0*y_test,'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(222)
        plt.plot(y_test, px_tbt[3, :]-px_tbt[2, :],label='SixTrack results')
        plt.plot(y_test,0*y_test,'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('px')
        plt.title('The 4D kick')
        plt.legend(loc='best')

        plt.subplot(223)
        plt.plot(y_test, py_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('py')
        plt.plot(y_test, BBKickVector(y_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.subplot(224)
        plt.plot(y_test, py_tbt[3, :]-py_tbt[2, :],label='SixTrack results')
        plt.plot(y_test, BBKickVector(y_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('y [m]')
        plt.ylabel('px')
        plt.legend(loc='best')
        plt.show()
quit();
# Check horizontal separation
if 0:
        my_h_sep_mm=5.
        my_v_sep_mm=0.
        x_test= np.linspace(-10e-2, 10e-2, 1000)
        y_test = 0.

        x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
                Dx_wrt_CO_m=x_test,
                Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
                Dsigma_wrt_CO_m=0, Ddelta_wrt_CO=0, n_turns=10,
                mode = 'ebe')

        plt.figure()
        plt.plot(x_test, px_tbt[1, :],label='SixTrack results')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.title('The 6D kick')
        plt.plot(x_test,BBKickVector(x_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.legend(loc='best')

        plt.figure()
        plt.plot(x_test, px_tbt[3, :]-px_tbt[2, :],label='SixTrack results')
        plt.plot(x_test,BBKickVector(x_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
                -BBKickVector(0+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
        plt.grid()
        plt.xlabel('x [m]')
        plt.ylabel('px')
        plt.title('The 4D kick')
        plt.legend(loc='best')


# Check vertical separation
x_test= 0
y_test = np.linspace(-10e-2, 10e-2, 1000)
my_h_sep_mm=0.
my_v_sep_mm=5.

x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
        Dx_wrt_CO_m=x_test,
        Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
        Dsigma_wrt_CO_m=0, Ddelta_wrt_CO=0, n_turns=10,
        mode = 'ebe')

plt.figure()
plt.plot(y_test, py_tbt[1, :],label='SixTrack results')
plt.grid()
plt.xlabel('y [m]')
plt.ylabel('py')
plt.title('The 6D kick')
plt.plot(y_test,BBKickVector(y_test+my_v_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
        -BBKickVector(0+my_v_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
plt.legend(loc='best')

plt.figure()
plt.plot(y_test, py_tbt[3, :]-py_tbt[2, :],label='SixTrack results')
plt.plot(y_test,BBKickVector(y_test+my_h_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3)\
        -BBKickVector(0+my_v_sep_mm, my_partnum, part_on_CO.gamma0, np.sqrt(my_sigma_xx_mm2)*1e-3),'.r',label='analytical reference')
plt.grid()
plt.xlabel('y [m]')
plt.ylabel('py')
plt.title('The 4D kick')
plt.legend(loc='best')

# Plot separations
z_test = np.linspace(-10e-2, 10e-2, 1000)
y_test = 0.
x_test = 0

x_tbt, px_tbt, y_tbt, py_tbt, sigma_tbt, delta_tbt = hp.track_particle_sixtrack(partCO=part_on_CO.to_dict(),
        Dx_wrt_CO_m=x_test,
        Dpx_wrt_CO_rad=0, Dy_wrt_CO_m=y_test, Dpy_wrt_CO_rad=0,
        Dsigma_wrt_CO_m=z_test, Ddelta_wrt_CO=0, n_turns=10,
        mode = 'ebe')

plt.figure()
plt.plot(z_test,  px_tbt[1, :])
plt.show()






