import ConfigParser, os
import pysm_synchrotron, pysm_thermaldust, pysm_cmb,pysm_spinningdust, pysm_noise, pysm_freefree
from pysm import output, config2list#, run_pysm_comp, file_path, smooth_write
import healpy as hp
import numpy as np
import argparse
from multiprocessing import Manager, Process

pysm_submods = {
    'thermaldust':pysm_thermaldust,
    'synchrotron':pysm_synchrotron,
    'cmb':pysm_cmb,
    'noise':pysm_noise,
    'freefree':pysm_freefree,
    'spinningdust':pysm_spinningdust
}

def run_pysm_comp(pysm_comp,config_file,result):
    result[pysm_comp] = pysm_submods[pysm_comp].main(config_file)

def file_path(o,freq):
    comps = str()
    for i in sorted(o.components): comps='_'.join([comps,i[0:5]])
    fname = ''.join([o.output_prefix,comps, str(freq).replace('.', 'p'),'_', str(o.nside), '.fits'])
    path = os.path.join(o.output_dir, fname)
    return path

def smooth_write(sky_freq,o,freq,fwhm):
    if o.smoothing: sky_freq = hp.smoothing(sky_freq,fwhm=(np.pi/180.)*fwhm,verbose=False)

    path = file_path(o,freq)
    hp.write_map(path, hp.ud_grade(sky_freq, nside_out=o.nside), coord='G', column_units = ''.join(o.output_units), column_names = None, extra_header = config2list(Config))
    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Code to simulate galactic foregrounds.')
    parser.add_argument('config_file', help='Main configuration file.')
    
##Get the output directory in order to save the configuration file.
    Config = ConfigParser.ConfigParser()
    Config.read(parser.parse_args().config_file)
    out = output(Config._sections['GlobalParameters'])

    if out.debug:

##Print information about the run:
        print '----------------------------------------------------- \n'
        print ''.join("%s: %s \n" % item   for item in vars(out).items())
        print '-----------------------------------------------------'

#Save the configuration file.
    if not os.path.exists(out.output_dir): os.makedirs(out.output_dir)
    with open(out.output_dir+out.output_prefix+'main_config.ini','w') as configfile: Config.write(configfile)


    print '----------------------------------------------------- \n'
#Create synchrotron, dust, AME, freefree,  and cmb maps at output frequencies then add noise.
    if out.instrument_noise: out.components.append('noise')

    result = Manager().dict()
    
    processes = [Process(target=run_pysm_comp,args=(comp,parser.parse_args().config_file,result)) for comp in out.components]
    for p in processes: p.start()
    for p in processes: p.join()
    sky = sum(result.values())

#Smooth sky maps by defined FWHM.
#Change to orering: (frequency, stokes param, pixels)
    sky = np.swapaxes(sky,0,1)


    if out.smoothing:
        print 'Smoothing output maps.'
        print '----------------------------------------------------- \n'

    processes = [Process(target=smooth_write,args=(sky[i,...],out,freq,fwhm)) for i,(freq,fwhm) in enumerate(zip(out.output_frequency,out.fwhm))]
    for p in processes: p.start()
    for p in processes: p.join()

    print '-----------------------------------------------------\n'
    print 'PySM completed successfully. \n'
    print '-----------------------------------------------------'




