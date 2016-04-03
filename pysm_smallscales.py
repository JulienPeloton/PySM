import healpy as hp
import numpy as np
from scipy.optimize import minimize
import scipy.constants as const
from pysm import read_map_wrapped


def safe_log10(x,minval=0.0000000001):
    return np.log10(x.clip(min=minval))

def b_l(ell,sigma):
    return np.exp(-0.5*ell*(ell+1)*sigma**2)

def sigma(fwhm):
    return fwhm*(2.*np.pi/(60.*360.))/(np.sqrt(8.*np.log(2.)))

def fit_cl_conv(map,l_fit,theta):
    cl = hp.anafast(map)
    ell = np.arange(cl.size)
    coeffs = np.polyfit(safe_log10(ell[l_fit[0]:l_fit[1]]),safe_log10(cl[l_fit[0]:l_fit[1]]),1)
    cl_pl = np.power(10,coeffs[1])*np.power(ell,coeffs[0])
    print('Coefficients of power law fit:', coeffs)
    return cl_pl*(1-b_l(ell,sigma(theta*60.))**2)

def ss_model(gamma,cl_pl,map_in,nside_o):
    g = hp.synfast(cl_pl[1:],nside=nside_o,verbose=False,new=True)
    ss = (g-np.mean(g))/np.sqrt(np.var(g))*(np.abs(map_in)**gamma)
    return ss

def chi_sq(gamma,map_in,cl_pl,pl_fit,nside_o):
    model = ss_model(gamma,cl_pl,map_in,nside_o)
    cl_model = hp.anafast(model)
    ell = np.arange(cl_model.size)
    chi_sq_arg = (cl_model[1:]/cl_pl[1:]-1)**2*(2*ell[1:]+1)/2
    return np.sum(chi_sq_arg[pl_fit[0]:pl_fit[1]])

def generate_ss_map(map_in,nside_o,l_fit,pl_fit,theta):
    print('Fitting power law.')
    cl_pl = fit_cl_conv(map_in,l_fit,theta)
    print('Minimizing chi squared at large l.')
    res = minimize(chi_sq,1.3,args=(map_in,cl_pl,pl_fit,nside_o),method='Powell')
    return ss_model(res.x,cl_pl,map_in,nside_o)

""" 
OPTIONS:

Details of the map to be extended to small scales.  The map will be
 read in and degraded to the specified nside_in before being used.  Therefore 
if the template is nside 1024 this must be specified as:
 
>nside_in = 1024

The output is generated at nside_out.

The other specified information is:

l_fit:  the multipole range over which the power law is fitted.
theta_res:  the angular scale (in deg) of the limiting resolution of the input.
pl_fit: range of multipoles over which the chi_sq fits.

"""

nside_in = 256
nside_out = 256

l_fit = [50.,100.]
pl_fit = [200.,600.]

theta_res = 1.   

file_in = './Ancillaries/ThermalDust/mbb/smoothed_dust_polq_em.fits'
file_out = './Ancillaries/ThermalDust/mbb/smoothed_dust_polq_em_ss.fits'


"""
MAIN CODE
"""

map_in = read_map_wrapped(file_in, nside_in, 0)
map_ss = generate_ss_map(map_in,nside_out,l_fit,pl_fit,theta_res)

print('Writing new small scale map to '+file_out)

hp.write_map(file_out,map_in+map_ss)





