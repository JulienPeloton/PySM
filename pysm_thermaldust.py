import numpy as np
import healpy as hp
import ConfigParser
from pysm import scale_freqs, convert_units, bandpass_integrated, output, component

def main():
	
	print('Computing dust map using modified black body for a single component.')
	print '----------------------------------------------------- \n'
#Read configuration into classes
	Config = ConfigParser.ConfigParser()
	Config.read('main_config.ini')
	out = output(Config._sections['GlobalParameters'])

	Config.read('./ConfigFiles/'+Config.get('ThermalDust','model')+'_config.ini')
	dust = component(Config._sections['ThermalDust'])

#In this case the scaling is done in uK_RJ, so the unit conversionn is different to synchrotron.
	conv1 = convert_units(dust.template_units, ['u','K_RJ'], dust.freq_ref)
	conv2 = convert_units(['u','K_RJ'],out.output_units,out.output_frequency)
	unit_conversion = conv1*conv2.reshape((len(out.output_frequency),1))

#Do the scaling.
	scaled_map_dust = scale_freqs(dust,out,pol=False)*dust.em_template*unit_conversion
	scaled_map_dust_pol = scale_freqs(dust,out,pol=True)[np.newaxis,...]*np.array([dust.polq_em_template,dust.polu_em_template])[:,np.newaxis,:]*unit_conversion
	
	if out.debug == True:
		dus = np.concatenate([scaled_map_dust[np.newaxis,...],scaled_map_dust_pol])
		for i in range(0,len(out.output_frequency)):
			hp.write_map(out.output_dir+'dust_%d.fits'%(out.output_frequency[i]),dus[:,i,:],coord='G',column_units=out.output_units)

       	return np.concatenate([scaled_map_dust[np.newaxis,...],scaled_map_dust_pol])



