import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

nside = 256
npix = hp.nside2npix(nside)
fsky_pix = 1./npix

ndeg = (180./np.pi)**2*(4.*np.pi)

beta_mean = 1.59
beta_sigma_deg = 0.3

beta_sigma_pix = np.sqrt(beta_sigma_deg**2*(npix/ndeg))

beta_d = np.random.normal(beta_mean,beta_sigma_pix,npix)

hp.write_map('model2_dust_beta.fits',beta_d)
