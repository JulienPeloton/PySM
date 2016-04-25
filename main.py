import ConfigParser, os
import pysm_synchrotron,pysm_thermaldust,pysm_cmb,pysm_spinningdust, pysm_noise
from pysm import output
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
#Create synchrotron, dust, AME,  and cmb maps at output frequencies then add noise.
if 'synchrotron' in out.components:
    sky = pysm_synchrotron.main(parser.parse_args().config_file)

if 'thermaldust' in out.components:
    sky = sky + pysm_thermaldust.main(parser.parse_args().config_file)

if 'spinningdust' in out.components:
    sky = sky + pysm_spinningdust.main(parser.parse_args().config_file)

if 'cmb' in out.components:
    sky = sky + pysm_cmb.main(parser.parse_args().config_file)

if out.instrument_noise == True:
    sky = sky + pysm_noise.instrument_noise(parser.parse_args().config_file)



comps =str()
for i in sorted(out.components): comps = comps+i[0:4]+'_'
if out.instrument_noise == True: comps = comps + 'noisy_'
fname = list()
for i in range(len(out.output_frequency)): 
    
    fname.append(comps+str(out.output_frequency[i]).replace('.','p')+'_'+str(out.nside)+'.fits')
    hp.write_map(out.output_dir+out.output_prefix+fname[i],hp.ud_grade(sky[:,i,:],nside_out=out.nside),coord='G', column_units=out.output_units[0]+out.output_units[1])

    if out.debug == True:
        print 'Written to '+out.output_dir+out.output_prefix+fname[i]

print '-----------------------------------------------------\n'
print 'PySM completed successfully. \n' 
print '-----------------------------------------------------'

