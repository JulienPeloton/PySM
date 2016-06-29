import ConfigParser, os
import pysm_synchrotron,pysm_thermaldust,pysm_cmb,pysm_spinningdust, pysm_noise, pysm_freefree
from pysm import output,config2list
import healpy as hp
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Code to simulate galactic foregrounds.')
parser.add_argument('config_file', help='Main configuration file.')


##Get the output directory in order to save the configuration file.
Config = ConfigParser.ConfigParser()
Config.read(parser.parse_args().config_file)
out = output(Config._sections['GlobalParameters'])

if out.debug == True:

##Print information about the run:
    print '----------------------------------------------------- \n'
    print ''.join("%s: %s \n" % item   for item in vars(out).items())
    print '-----------------------------------------------------'

#Save the configuration file.
if not os.path.exists(out.output_dir): os.makedirs(out.output_dir)
with open(out.output_dir+out.output_prefix+'main_config.ini','w') as configfile: Config.write(configfile)

sky = np.zeros(hp.nside2npix(out.nside))
print '----------------------------------------------------- \n'
#Create synchrotron, dust, AME, freefree,  and cmb maps at output frequencies then add noise.
if 'synchrotron' in out.components:
    sky = pysm_synchrotron.main(parser.parse_args().config_file)

if 'thermaldust' in out.components:
    sky = sky + pysm_thermaldust.main(parser.parse_args().config_file)

if 'spinningdust' in out.components:
    sky = sky + pysm_spinningdust.main(parser.parse_args().config_file)

if 'freefree' in out.components:
    sky = sky + pysm_freefree.main(parser.parse_args().config_file)

if 'cmb' in out.components:
    sky = sky + pysm_cmb.main(parser.parse_args().config_file)

if out.instrument_noise == True:
    sky = sky + pysm_noise.instrument_noise(parser.parse_args().config_file)

#Smooth sky maps by defined FWHM.
import matplotlib.pyplot as plt

#Change to orering: (frequency, stokes param, pixels)
sky = np.swapaxes(sky,0,1)

#Smooth maps if asked for.
if out.smoothing:
    print 'Smoothing output maps.'
    print '----------------------------------------------------- \n'
    for i,fwhm in enumerate(out.fwhm):
        sky[i,...] = hp.smoothing(sky[i,...],fwhm=(np.pi/180.)*fwhm,verbose=False)

#Write maps to ouput directory.
comps = str()
for i in sorted(out.components): comps = comps + i[0:5] + '_'
if out.instrument_noise: comps = comps + 'noisy_'
fname = list()

for i,freq in enumerate(out.output_frequency): 
    
    fname = ''.join([out.output_prefix,comps, str(freq).replace('.', 'p'),'_', str(out.nside), '.fits'])
    path = os.path.join(out.output_dir, fname)

    hp.write_map(path, hp.ud_grade(sky[i,...], nside_out=out.nside), coord='G', column_units = ''.join(out.output_units), column_names = None, extra_header = config2list(Config))

    if out.debug: print 'Written to %s'%path

print '-----------------------------------------------------\n'
print 'PySM completed successfully. \n' 
print '-----------------------------------------------------'

