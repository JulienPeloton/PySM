import ConfigParser, os
import pysm_synchrotron,pysm_thermaldust, pysm_cmb
from pysm import output
import healpy as hp
import numpy as np

##Get the output directory in order to save the configuration file.
Config = ConfigParser.ConfigParser()
Config.read('main_config.ini')
out = output(Config._sections['GlobalParameters'])

##Print information about the run:
print '----------------------------------------------------- \n'
print 'PYSM \n'
print '----------------------------------------------------- \n'
print ''.join("%s: %s \n" % item   for item in vars(out).items())
print '-----------------------------------------------------'

#Save the configuration file.
if not os.path.exists(out.output_dir): os.makedirs(out.output_dir)
with open(out.output_dir+'main_config.ini','a') as configfile: Config.write(configfile)

sky = np.zeros(hp.nside2npix(out.nside))

#Create synchrotron, dust, and cmb maps at output frequencies.
if 'synchrotron' in out.components:
    sky = pysm_synchrotron.main()
if 'thermaldust' in out.components:
    sky = sky + pysm_thermaldust.main()
if 'cmb' in out.components:
    sky = sky + pysm_cmb.main()


comps =str()
for i in sorted(out.components): comps = comps+i[0]+'_'
fname = list()
for i in range(len(out.output_frequency)): 
    fname.append(comps+str(out.output_frequency[i])+'.fits')
    hp.write_map(out.output_dir+fname[i],sky[:,i,:],coord='G',column_units=out.output_units)
    print 'Written to '+out.output_dir+fname[i]

print '-----------------------------------------------------\n'
print 'PySM completed successfully. \n' 
print '-----------------------------------------------------'

exit()
