import numpy as np
import healpy as hp
from pysm import scale_freqs, convert_units, component, output
import ConfigParser

def main():

#Read in configuration file to classes.
	Config = ConfigParser.ConfigParser()
	Config.read('main_config.ini')
	out = output(Config._sections['GlobalParameters'])

	Config.read('./ConfigFiles/'+Config.get('Synchrotron','model')+'_config.ini')
	synch = component(Config._sections['Synchrotron'])

	print('Computing synchrotron map with parameters:')
	print '----------------------------------------------------- \n'
	print ''.join("%s: %s \n" % item   for item in vars(synch).items())
        print '----------------------------------------------------- \n'

#The unit conversion takes care of the scaling being done in MJysr. After scaling we convert to whatever the output units are.
	conv1 = convert_units(synch.template_units, ['M','Jysr'], synch.freq_ref)
	conv2 = convert_units(['M','Jysr'],out.output_units,out.output_frequency)
	unit_conversion = conv1*conv2.reshape((len(out.output_frequency),1))
#Do the scaling.

	scaled_map_synch = scale_freqs(synch, out, pol=False)*synch.em_template*unit_conversion
	scaled_map_synch_pol = scale_freqs(synch, out, pol=True)[np.newaxis,...]*np.array([synch.polq_em_template,synch.polu_em_template])[:,np.newaxis,:]*unit_conversion

	if out.debug == True:
		syn = np.concatenate([scaled_map_synch[np.newaxis,...],scaled_map_synch_pol])
		for i in range(0,len(out.output_frequency)):
			hp.write_map(out.output_dir+'synch_%d.fits'%(out.output_frequency[i]),syn[:,i,:],coord='G',column_units=out.output_units)

	return np.concatenate([scaled_map_synch[np.newaxis,...],scaled_map_synch_pol])



