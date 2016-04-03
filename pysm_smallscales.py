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

def local_mean(m,pix,radius): #radius is in radians
    nside= hp.get_nside(m)
    pos = hp.pix2vec(nside,pix)
    inds = hp.query_disc(nside,pos,radius)
    disc_mean = np.mean(m[inds])
    return disc_mean

def mod_gaussian_model(gamma,cl_pl,map_in,nside_o):
    g = hp.synfast(cl_pl[1:],nside=nside_o,verbose=False,new=True)
    ss = (g-np.mean(g))/np.sqrt(np.var(g))*(np.abs(map_in)**gamma)
    return ss

def lognormal_model(gamma,map_in,m_0,R_mean,d_g,d_g_var):
    norm = np.min(map_in)/m_0
    x = R_mean**gamma*d_g-0.5*d_g_var
    return norm*(np.exp(x)-1.)

def chi_sq(gamma,map_in,cl_pl,pl_fit,nside_o,model,m_0=None,R_mean=None,d_g=None,d_g_var=None):
    if model == 'mod_gaussian': 
        model = mod_gaussian_model(gamma,cl_pl,map_in,nside_o)
    if model == 'lognormal':
        model = lognormal_model(gamma,map_in,m_0,R_mean,d_g,d_g_var)
    cl_model = hp.anafast(model)
    ell = np.arange(cl_model.size)
    chi_sq_arg = (cl_model[1:]/cl_pl[1:]-1)**2*(2*ell[1:]+1)/2
    print gamma
    return np.sum(chi_sq_arg[pl_fit[0]:pl_fit[1]])

def mod_gaussian_ss_map(map_in,nside_o,l_fit,pl_fit,theta):
    print('Fitting power law.')
    cl_pl = fit_cl_conv(map_in,l_fit,theta)
    print('Minimizing chi squared at large l.')
    res = minimize(chi_sq,1.3,args=(map_in,cl_pl,pl_fit,nside_o,'mod_gaussian'),method='Powell')
    return mod_gaussian_model(res.x,cl_pl,map_in,nside_o)

def lognormal_ss_map(map_in,nside_o,l_fit,pl_fit,theta):
    m_0 = np.mean(map_in);npix=map_in.size
    d_ls = (map_in-m_0)/m_0
    print('Fitting power law.')
    cl_pl = fit_cl_conv(d_ls,l_fit,theta)
    print('Computing local scaling.')
    R_mean = np.zeros(npix)
    for i in np.arange(npix):
        R_mean[i] = local_mean(map_in,i,(np.pi/180.)*4.)/m_0
    print('Minimizing chi squared.')
    d_g = hp.synfast(cl_pl[1:],nside=nside_o,new=True)
    d_g_var = np.var(d_g)
    res = minimize(chi_sq,.6,args=(map_in,cl_pl,pl_fit,nside_o,'lognormal',m_0,R_mean,d_g,d_g_var),method='Powell')
    d_ss = lognormal_model(res.x,map_in,m_0,R_mean,d_g,d_g_var)
    return m_0*d_ss

def generate_ss_map(map_in,nside_o,l_fit,pl_fit,theta,method):
    if method=='mod_gaussian': return mod_gaussian_ss_map(map_in,nside_o,l_fit,pl_fit,theta)
    if method=='lognormal': return lognormal_ss_map(map_in,nside_o,l_fit,pl_fit,theta)



""" 
OPTIONS:

method    : options are 'lognormal' or 'mod_gaussian'.  'lognormal' can be used to generate
            intensity map small scales as it ensures that the final map is positive.  
            'mod_gaussian' can be used to add small scales to Q and U maps.  

nside_in  : nside of input map. 
nside_out : nside of output map.
theta_res : the angular scale (in deg) of the limiting resolution of the input.

l_fit     : the multipole range over which the power law is fitted.
pl_fit    : range of multipoles over which the chi_sq fits.

file_in   : address of input map. Should not be a multiple-map fits file.
file_out  : where to put the final map with added small scales.
"""

method = 'lognormal'  #options are 'lognormal' or 'mod_gaussian'

nside_in = 256
nside_out = 256
theta_res = 1.

l_fit = [20.,50.]
pl_fit = [200.,600.]

file_in = './Ancillaries/Synchrotron/mamd2008/haslam408_dsds_Remazeilles2014_8p33_mono_sub.fits'
file_out = './Ancillaries/Synchrotron/mamd2008/haslam408_dsds_Remazeilles2014_8p33_mono_sub_ss.fits'


"""
MAIN CODE
"""

map_in = read_map_wrapped(file_in, nside_in, 0)
map_ss = generate_ss_map(map_in,nside_out,l_fit,pl_fit,theta_res,method)

print('Writing new small scale map to '+file_out)

hp.write_map(file_out,map_in+map_ss)





