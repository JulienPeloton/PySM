import numpy as np
import healpy as hp
import ConfigParser
from pysm import scale_freqs, convert_units, output, component

def scale_dust_pop(i_pop,npop,out,Config):

	if npop==1 : secname='ThermalDust'
	else : secname='ThermalDust_%d'%i_pop
	dust = component(Config._sections[secname],out.nside)
	print('Computing dust maps (%d-th component)'%i_pop)
	print '----------------------------------------------------- \n'
	if out.debug == True: 
		print ''.join("%s: %s \n" % item   for item in vars(dust).items())
		print '----------------------------------------------------- \n'

	conv_I = convert_units(dust.template_units, ['u','K_RJ'], dust.freq_ref)
        conv_pol =  convert_units(dust.template_units, ['u','K_RJ'], dust.pol_freq_ref)
        conv2 = convert_units(['u','K_RJ'],out.output_units,out.output_frequency)
        unit_conversion_I = conv_I*conv2.reshape((len(out.output_frequency),1))
        unit_conversion_pol = conv_pol*conv2.reshape((len(out.output_frequency),1))

	scaled_map_dust = scale_freqs(dust,out,pol=False)*dust.em_template*unit_conversion_I
	scaled_map_dust_pol = scale_freqs(dust,out,pol=True)[np.newaxis,...]*np.array([dust.polq_em_template,dust.polu_em_template])[:,np.newaxis,:]*unit_conversion_pol


	if out.debug == True:
                dus = np.concatenate([scaled_map_dust[np.newaxis,...],scaled_map_dust_pol])
                for i in range(0,len(out.output_frequency)):
                        hp.write_map(out.output_dir+out.output_prefix+'dust%d'%i_pop+'_%d'%(out.output_frequency[i])+'_'+str(out.nside)+'.fits',dus[:,i,:],coord='G',column_units=out.output_units)

	return np.concatenate([scaled_map_dust[np.newaxis,...],scaled_map_dust_pol])

def main(fname_config):
	
#Read configuration into classes
	Config = ConfigParser.ConfigParser()
	Config_model = ConfigParser.ConfigParser()

	Config.read(fname_config)
	out = output(Config._sections['GlobalParameters'])

	a=Config_model.read('./ConfigFiles/'+Config.get('ThermalDust','model')+'_config.ini')
	if a==[] :
		print 'Couldn\'t find file '+'./ConfigFiles/'+Config.get('ThermalDust','model')+'_config.ini'
		exit(1)
	npops = len(Config_model.sections())

	with open(out.output_dir+out.output_prefix+'thermaldust_config.ini','w') as configfile: Config_model.write(configfile)

	dust_out = 0.

	for i_pop in np.arange(npops)+1 : 
		dust_out += scale_dust_pop(i_pop,npops,out,Config_model)

	return dust_out


