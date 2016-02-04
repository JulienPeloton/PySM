PySM v-0.1 
Most recent version available at: https://github.com/bthorne93/PySM

-----------------------------------------------------------------------------------------

This code generates full-sky simulations of galactic foregrounds in intensity and polarization up to 857GHz. The components simulated are: thermal dust, synchrotron, and CMB, with spinning dust and free-free to be added soon.  

-----------------------------------------------------------------------------------------

This code is known to work with: 

    - python 2.7.6
    - healpy 1.9.1     
    - numpy 1.8.1
    - scipy 0.14.0
    - astropy 1.1.1

It requires at least: 

    - healpy 1.9.1
   	    
Note that the healpy.write_map function will not work properly with outdated versions of healpy.  healpy.write_map will also throw a warning when run with the most recent versions of healpy and astropy because healpy uses a deprecated astropy function. This does not affect the outcome of the code.

-----------------------------------------------------------------------------------------
 
To run the code, in the directory containing main.py type:

    > python main.py main_config.ini

The default outputs are maps, at the specified frequencies, of the summed emission of all the specified components. The default output directory is './Output/'.

To change the parameters of the simulation edit the 'main_config.ini' file (or create a separate configuration file). The different parameters are fully described in the comments of this file as well as the individual model config files in './ConfigFiles/<model>_config.ini'. 

-----------------------------------------------------------------------------------------

The nominal models used for the components are: 

    Thermal dust: Thermal dust is modelled as a single-component modified black body (mbb).  We use templates for emission at 545GHz in intensity and 353GHz in polarisation, and scale these to different frequencies with a mbb spectrum using the spatially varying temperature and spectral index obtained from the Planck data using the Commander code (Planck Collaboration, Adam, R., Ade, P. A. R., et al. 2015, arXiv:1502.01588). Note that it therefore assumes the same spectral index for polarization as for intensity.  All input templates provided with the code have already been degraded to N_side 256 and smoothed to degree scale. 

    Synchrotron:  A power law scaling is used for the synchrotron emission, with a spatially varying spectral index.  The emission template and spectral index map were derived using a combination of the Haslam 408 MHz data and WMAP 23 GHz data (Miville-Deschenes, M.-A. et al., 2008, A&A, 490, 1093). The same scaling is used for intensity and polarization.  This is the same prescription as used in the Planck Sky Model's 'power law' option (Delabrouille et al. A&A 553, A96, 2013). A 'curved power law' model is also supported with a single isotropic curvature index.

    CMB: A lensed CMB realisation is computed using taylens, a code to compute a lensed CMB realisaion using nearest-neighbour Taylor interpolation (Naess, S. K. and Louis, T. Journal of Cosmology and Astroparticle Physics September 2013). This code takes, as an input, a set of unlensed Cl's generated using CAMB (http://www.camb.info/).

-----------------------------------------------------------------------------------------


		


