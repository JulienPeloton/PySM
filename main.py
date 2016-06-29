import ConfigParser, os
import pysm_synchrotron,pysm_thermaldust,pysm_cmb,pysm_spinningdust, pysm_noise, pysm_freefree
from pysm import output,config2list
import healpy as hp
import numpy as np
import argparse
from multiprocessing import Manager, Process, Queue

def run_pysm_comp(pysm_comp,config_file,result):
    if 'synchrotron' == pysm_comp:
        result[pysm_comp] = pysm_synchrotron.main(config_file)

    if 'thermaldust' == pysm_comp:
        result[pysm_comp] =  pysm_thermaldust.main(config_file)

    if 'spinningdust' == pysm_comp:
        result[pysm_comp] = pysm_spinningdust.main(config_file)

    if 'freefree' == pysm_comp:
        result[pysm_comp] = pysm_freefree.main(config_file)

    if 'cmb' == pysm_comp:
        result[pysm_comp] = pysm_cmb.main(config_file)

    if 'noise' == pysm_comp:
        result[pysm_comp] == pysm_noise.main(config_file)

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


#    result = Queue()
    result = Manager().dict()
    
    processes = [Process(target=run_pysm_comp,args=(comp,parser.parse_args().config_file,result)) for comp in out.components]
    for p in processes: p.start()
    for p in processes: p.join()
    sky = sum(result.values())

#Smooth sky maps by defined FWHM.
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
    for i in sorted(out.components): comps='_'.join([comps,i[0:5]])

    for i,freq in enumerate(out.output_frequency): 
    
        fname = ''.join([out.output_prefix,comps, str(freq).replace('.', 'p'),'_', str(out.nside), '.fits'])
        path = os.path.join(out.output_dir, fname)
        
        hp.write_map(path, hp.ud_grade(sky[i,...], nside_out=out.nside), coord='G', column_units = ''.join(out.output_units), column_names = None, extra_header = config2list(Config))

        if out.debug: print 'Written to %s'%path

    print '-----------------------------------------------------\n'
    print 'PySM completed successfully. \n' 
    print '-----------------------------------------------------'

