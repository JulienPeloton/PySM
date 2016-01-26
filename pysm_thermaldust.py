import numpy as np
import healpy as hp
import ConfigParser
from pysm import scale_freqs, convert_units, bandpass_integrated, output, component

def main():
	
	print('Computing dust map using modified black body for a single component.')
	print '----------------------------------------------------- \n'
#Read configuration into classes
	Config = ConfigParser.ConfigParser()
	Config.read('pysm_config.ini')
	out = output(Config._sections['GlobalParameters'])
	dust = component(Config._sections['ThermalDust'])

#In this case the scaling is done in uK_RJ, so the unit conversionn is different to synchrotron.
	conv1 = convert_units(dust.template_units, ['u','K_RJ'], dust.freq_ref)
	conv2 = convert_units(['u','K_RJ'],out.output_units,out.output_frequency)
	unit_conversion = conv1*conv2.reshape((len(out.output_frequency),1))

#Do the scaling.
	scaled_map_dust = scale_freqs(dust.model,out.output_frequency,dust.beta_template,dust.freq_ref,None,None,dust.temp_template)*dust.em_template*unit_conversion
	scaled_map_dust_polq = scale_freqs(dust.model,out.output_frequency,dust.beta_template,dust.pol_freq_ref,None,None,dust.temp_template)*dust.polq_em_template*unit_conversion
	scaled_map_dust_polu = scale_freqs(dust.model,out.output_frequency,dust.beta_template,dust.pol_freq_ref,None,None,dust.temp_template)*dust.polu_em_template*unit_conversion

	for i in range(0,len(out.output_frequency)):
		d = [scaled_map_dust[i],scaled_map_dust_polq[i],scaled_map_dust_polu[i]]
		hp.write_map(dust.output_dir+'pysm_run/'+'dust_%d.fits'%(out.output_frequency[i]),d,coord='G',column_units=out.output_units)

	del dust, d
