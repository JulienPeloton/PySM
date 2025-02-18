import numpy as np
import healpy as hp
from pysm import scale_freqs, convert_units, component, output
import ConfigParser

def scale_synch_pop(i_pop,npop,out,Config) :

	if npop==1 : secname='Synchrotron'
	else : secname='Synchrotron_%d'%i_pop
	synch = component(Config._sections[secname],out.nside)
	print('Computing synchrotron maps (%d-th component).'%i_pop)
	print '----------------------------------------------------- \n'
	if out.debug == True:
		print ''.join("%s: %s \n" % item   for item in vars(synch).items())
		print '----------------------------------------------------- \n'

#The unit conversion takes care of the scaling being done in uK_RJ. After scaling we convert to whatever the output units are.
	conv_I = convert_units(synch.template_units, ['u','K_RJ'], synch.freq_ref)
	conv_pol = convert_units(synch.template_units, ['u','K_RJ'], synch.pol_freq_ref) 
	conv2 = convert_units(['u','K_RJ'],out.output_units,out.output_frequency)
	unit_conversion_I = conv_I*conv2.reshape((len(out.output_frequency),1))
	unit_conversion_pol = conv_pol*conv2.reshape((len(out.output_frequency),1))
#Do the scaling.

	scaled_map_synch = scale_freqs(synch, out, pol=False)*synch.em_template*unit_conversion_I
	scaled_map_synch_pol = scale_freqs(synch, out, pol=True)[np.newaxis,...]*np.array([synch.polq_em_template,synch.polu_em_template])[:,np.newaxis,:]*unit_conversion_pol

#This section forces P/I<0.75. This is done using the same procedure as the PSM 1.7.8 psm_synchrotron.pro.

	P = np.sqrt(scaled_map_synch_pol[0,:,:]**2+scaled_map_synch_pol[1,:,:]**2)/scaled_map_synch
	F = 0.75*np.tanh(P/0.75)/P
	scaled_map_synch_pol[0,:,:]=F*scaled_map_synch_pol[0,:,:]
	scaled_map_synch_pol[1,:,:]=F*scaled_map_synch_pol[1,:,:]

#-------

	if out.debug == True:
		syn = np.concatenate([scaled_map_synch[np.newaxis,...],scaled_map_synch_pol])
		for i in range(0,len(out.output_frequency)):
			hp.write_map(out.output_dir+out.output_prefix+'synch%d'%i_pop+'_%d'%(out.output_frequency[i])+'_'+str(out.nside)+'.fits',syn[:,i,:],coord='G',column_units=out.output_units)

	return np.concatenate([scaled_map_synch[np.newaxis,...],scaled_map_synch_pol])

def main(fname_config):

#Read in configuration file to classes.
	Config = ConfigParser.ConfigParser()
	Config_model = ConfigParser.ConfigParser()

	Config.read(fname_config)
	out = output(Config._sections['GlobalParameters'])

	a=Config_model.read('./ConfigFiles/'+Config.get('Synchrotron','model')+'_config.ini')
	if a==[] :
		print 'Couldn\'t find file'+' ./ConfigFiles/'+Config.get('Synchrotron','model')+'_config.ini'
		exit(1)
	
	npops = len(Config_model.sections())
		
	with open(out.output_dir+out.output_prefix+'synchrotron_config.ini','w') as configfile: Config_model.write(configfile)
	
	synch_out = 0.
	
	for i_pop in np.arange(npops)+1 :
		synch_out += scale_synch_pop(i_pop,npops,out,Config_model)

	return synch_out
