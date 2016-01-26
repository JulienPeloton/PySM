import numpy as np
import healpy as hp
from pysm import *
import ConfigParser

def main():

#Read config file.
	Config = ConfigParser.ConfigParser()
	Config.read('main_config.ini')
	out = output(Config._sections['GlobalParameters'])
	cmb_order = Config.getint('CMB', 'order')
	cmb_specs  = Config.get('CMB', 'input_spectra')
	cmb_output = Config.get('CMB','output_dir')
	compute_lensing = Config.getboolean('CMB','compute_lensed_cmb')

####This is edited from taylens code.
	
	if compute_lensing == True:
		
		print('Using taylens to compute temperature map.')
		print '----------------------------------------------------- \n'
		synlmax = 8*out.nside #this used to be user-defined.
		data = np.transpose(np.loadtxt(cmb_specs))
		lmax_cl = len(data[0])+1
		l = np.arange(int(lmax_cl+1))
		synlmax = min(synlmax, l[-1])

#Reading input spectra in CAMB format.  CAMB outputs l(l+1)/2pi hence the corrections.

		cl_tebp_arr=np.zeros([10,lmax_cl+1])
		cl_tebp_arr[0,2:]=2*np.pi*data[1]/(l[2:]*(l[2:]+1))      #TT
		cl_tebp_arr[1,2:]=2*np.pi*data[2]/(l[2:]*(l[2:]+1))      #EE
		cl_tebp_arr[2,2:]=2*np.pi*data[3]/(l[2:]*(l[2:]+1))      #BB
		cl_tebp_arr[3,2:]=2*np.pi*data[5]/(l[2:]*(l[2:]+1))**2   #PP
		cl_tebp_arr[4,2:]=2*np.pi*data[4]/(l[2:]*(l[2:]+1))      #TE
		cl_tebp_arr[5,:] =np.zeros(lmax_cl+1)                    #EB
		cl_tebp_arr[6,:] =np.zeros(lmax_cl+1)                    #BP
		cl_tebp_arr[7,:] =np.zeros(lmax_cl+1)                    #TB
		cl_tebp_arr[8,2:]=2*np.pi*data[7]/(l[2:]*(l[2:]+1))**1.5 #EP
		cl_tebp_arr[9,2:]=2*np.pi*data[6]/(l[2:]*(l[2:]+1))**1.5 #TP

# Coordinates of healpix pixel centers
		ipos = np.array(hp.pix2ang(out.nside, np.arange(12*(out.nside**2))))

# Simulate a CMB and lensing field
		cmb, aphi = simulate_tebp_correlated(cl_tebp_arr,out.nside,synlmax)
		
		if cmb.ndim == 1: cmb = np.reshape(cmb, [1,cmb.size])

# Compute the offset positions
		phi, phi_dtheta, phi_dphi = hp.alm2map_der1(aphi,out.nside,lmax=synlmax)

		del aphi
		
		opos, rot = offset_pos(ipos, phi_dtheta, phi_dphi, pol=True, geodesic=False) #geodesic used to be used defined.
		del phi, phi_dtheta, phi_dphi

# Interpolate maps one at a time
		maps  = []
		for comp in cmb:
			for m in taylor_interpol_iter(comp, opos, cmb_order, verbose=False, lmax=None): #lmax here needs to be fixed
				pass
			maps.append(m)
		del opos, cmb
#save the map computed for future referemce.
		rm = apply_rotation(maps, rot)
		np.save(cmb_output+'cmb_spec',cmb_specs)
		hp.write_map(cmb_output+'lensed_cmb.fits',rm)
	else:
#option to use an already-computed lensed cmb map.	
		print('Scaling a lensed cmb temperature map.')
		print '----------------------------------------------------- \n'
		rm = hp.read_map(Config.get('CMB','lensed_cmb'),field=(0,1,2),verbose=False)

	for i in out.output_frequency:
		map_cmb= tuple([convert_units(['u','K_CMB'],['u','K_RJ'],i)*m for m in rm])
		hp.write_map(out.output_dir+'pysm_run/'+'lensed_cmb_%d.fits'%(i),map_cmb,coord='G',column_units=out.output_units)

