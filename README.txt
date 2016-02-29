PySM v-0.1 
Beta version.
Most recent version available at: https://github.com/bthorne93/PySM

Authors: Ben Thorne, David Alonso, Sigurd Naess, Jo Dunkley
Contact: ben.thorne@physics.ox.ac.uk
------------------------------------------------------------------------

This code generates full-sky simulations of Galactic foregrounds in intensity and 
polarization relevant for CMB experiments. The components simulated are: thermal dust, 
synchrotron, AME, and CMB. Free-free will be added soon.

Current version v-0.1 generates degree-scale smoothed maps, at Healpix Nside=256. There
exist options to integrate over a top hat bandpass, and to add instrument noise in 
intensity and polarization.
Smaller-scales and different Nside available soon.

There will be scope for a few options for the model for each component, attempting to 
be consistent with current data. The current v-0.1 version has a single nominal option 
for each component. 

This code is based on the large-scale Galactic part of Planck Sky Model code and uses 
some of its inputs (http://www.apc.univ-paris7.fr/~delabrou/PSM/psm.html, 
astro-ph/1207.3675).

If you use the code for any publications, please acknowledge it.

-----------------------------------------------------------------------

Dependencies:
This code uses python, and needs the healpy, numpy, scipy and astropy libraries. 
Versions of those that it is known to work with are:

    - python 2.7.6
    - healpy 1.9.1     
    - numpy 1.8.1
    - scipy 0.14.0
    - astropy 1.1.1

It requires at least: 

    - healpy 1.9.1
   	    
Note that the healpy.write_map function will not work properly with outdated versions 
of healpy.  healpy.write_map will also throw a warning when run with the most recent 
versions of healpy and astropy because healpy uses a deprecated astropy function. 
This does not affect the outcome of the code.

--------------------------------------------------------------------------
 
To run the code, in the directory containing main.py run:

    > python main.py main_config.ini

The default outputs are Healpix maps, at the specified frequencies, of the 
summed emission of all the chosen components. The default output directory is './Output/'.

To change the parameters of the simulation edit the 'main_config.ini' file (or 
create a separate configuration file). The different parameters are described 
in the comments of this ini file as well as the individual model config files 
in './ConfigFiles/<model>_config.ini'. 

--------------------------------------------------------------------------

The nominal models used for the components are: 

'dust1' = Thermal dust: Thermal dust is modelled as a single-component modified 
 black body (mbb).  We use dust templates for emission at 545 GHz in intensity and 
 353 GHz in polarisation from the Planck-2015 analysis, and scale these to different 
 frequencies with a mbb spectrum using the spatially varying temperature and spectral 
 index obtained from the Planck data using the Commander code (Planck Collaboration 
 2015, arXiv:1502.01588). Note that it therefore assumes the same spectral index for
 polarization as for intensity.  All input templates provided with the code have 
 already been degraded to Nside=256 and smoothed to degree scale. 

'synch1' = Synchrotron:  A power law scaling is used for the synchrotron emission, with 
a spatially varying spectral index.  The emission templates are smoothed to degree scale
and are the Haslam 408 MHz data reprocessed by Remazeilles et al 2015 MNRAS 451, 4311, 
and the WMAP 7-year 23 GHz Q/U maps (Jarosik et al 2011 ApJS, 192, 14J), smoothed to 3 
degree FWHM and with smaller scales added using the PSM code (Delabrouille et al. A&A 
553, A96, 2013). The spectral index map was derived using a combination of the Haslam 
408 MHz data and WMAP 23 GHz 7-year data (Miville-Deschenes, M.-A. et al., 2008, A&A, 490, 1093). 
The same scaling is used for intensity and polarization.  This is the same prescription 
as used in the Planck Sky Model's v1.7.8 'power law' option (Delabrouille et al. A&A 553, 
A96, 2013), but with the Haslam map updated to the Remazeilles version. A 'curved power 
law' model is also supported with a single isotropic curvature index.

'spdust1' = Spinning Dust: We model the AME as a sum of two spinning dust populations 
based on the Commander code (Planck Collaboration 2015, arXiv:1502.01588). A component 
is defined by an emission template at a reference frequency and a peak frequency of the 
emission law. Both populations have a spatially varying emission template, one 
population has a spatially varying peak frequency, and the other population has a 
spatially constant peak frequency.  The emission law is generated using the SpDust2 code 
(Ali-Haimoud 2008, http://arxiv.org/abs/0812.2904)  

'cmb1' = CMB: A lensed CMB realisation is computed using Taylens, a code to compute 
a lensed CMB realisation using nearest-neighbour Taylor interpolation 
(https://github.com/amaurea/taylens; Naess, S. K. and Louis, T. JCAP 09 001, 2013, 
astro-ph/1307.0719). This code takes, as an input, a set of unlensed Cl's generated 
using CAMB (http://www.camb.info/). The params.ini is in the Ancillary directory.

----------------------------------------------------------------------------


		
