import numpy as np
import scipy as sp
import healpy as hp, sys, time
from scipy.misc import factorial, comb

constants = {
    'T_CMB':2.7255,
    'h':6.62607004e-34,
    'k_B':1.38064852e-23,
    'c':2.99792e8
}

units = {
    'n':1.e-9,
    'u':1.e-6,
    'm':1.e-3,
    '1': 1.,
    'k':1.e3,
    'M':1.e6,
    'G':1.e9,
    'K_RJ': lambda x: 3.072387e4*x**2,
    'K_CMB': lambda x:  (0.017608676*x/(np.exp(0.017608676*x)-1))**2*np.exp(0.017608676*x)*3.072387e4*x**2,
    'Jysr': lambda x: 1.
}

class component(object):

    def __init__(self,cdict):
        keys = cdict.keys()
        if 'pol' in keys:
            self.pol = cdict['pol']
        if 'model' in keys:
            self.spectral_model = cdict['spectral_model']
        if 'em_template' in keys:
            self.em_template = hp.read_map(cdict['em_template'],verbose=False)
        if 'beta_template' in keys:
            self.beta_template = hp.read_map(cdict['beta_template'],verbose=False)
        if 'temp_template' in keys:
            self.temp_template = hp.read_map(cdict['temp_template'],verbose=False)
        if 'curvefreq' in keys:
            self.curvefreq = float(cdict['curvefreq'])
        if 'beta_curve' in keys:
            self.beta_curve = float(cdict['beta_curve'])
        if 'polq_em_template' in keys: 
            self.polq_em_template = hp.read_map(cdict['polq_em_template'],verbose=False)
        if 'polu_em_template' in keys:
            self.polu_em_template = hp.read_map(cdict['polu_em_template'],verbose=False)
        if 'freq_ref' in keys:
            self.freq_ref = float(cdict['freq_ref'])
        if 'pol_freq_ref' in keys:
            self.pol_freq_ref = float(cdict['pol_freq_ref'])
        if 'template_units' in keys:
            self.template_units = [cdict['template_units'][0],cdict['template_units'][1:]]
        if 'output_dir' in keys:
            self.output_dir = cdict['output_dir']

class output(object):
    def __init__(self, config_dict):
        self.debug = 'True' in config_dict['debug']
        self.components = [i for i in config_dict['components'].split()]
        self.output_frequency = [float(i) for i in config_dict['output_frequency'].split()]
        self.output_units = [config_dict['output_units'][0],config_dict['output_units'][1:]]
        self.nside = int(config_dict['nside'])
        self.output_dir = config_dict['output_dir']
        self.cmb_seed = config_dict['cmb_seed']

def convert_units(u_from, u_to, freq): #freq in GHz
    return units[u_from[0]]*units[u_from[1]](np.asarray(freq))/(units[u_to[0]]*units[u_to[1]](np.asarray(freq)))

def scale_freqs(component, output, pol, samples=100., width=10.):
    
     freq = np.asarray(output.output_frequency)
     if pol == False: freq_ref = component.freq_ref
     if pol ==True: freq_ref = component.pol_freq_ref

     if component.spectral_model=="curvepowerlaw": return (freq[...,np.newaxis]/freq_ref)**(component.beta_template+component.beta_curve*np.log10(freq[...,np.newaxis]/component.freq_curve))
     elif component.spectral_model=="powerlaw": return (freq[...,np.newaxis]/freq_ref)**component.beta_template
     elif component.spectral_model=="thermaldust":
        exponent=(constants['h']/constants['k_B'])*(freq[...,np.newaxis]*1.e9/component.temp_template)
        exponent_ref=(constants['h']/constants['k_B'])*(freq_ref*1.e9/component.temp_template)
        return (freq[...,np.newaxis]/freq_ref)**(component.beta_template+1)*((np.exp(exponent_ref)-1.)/(np.exp(exponent)-1.))
     else:
        print('No law selected')

        exit()

def bandpass_integrated(component,out,width,samples=10.):

    band_freqs = np.asarray(output_frequency)[...,np.newaxis]+np.linspace(-width/2.,width/2.,num=samples)

    return np.sum(scale_freqs(model,band_freqs,beta,freq_ref,freq_curve,beta_curve,temp),axis=np.ndim(band_freqs)-1)/samples

#This code is edited from the taylens code: Naess, S. K. and Louis, T. 2013 'Lensing simulations by Taylor expansion -  not so inefficient after all'  Journal of Cosmology and Astroparticle Physics September 2013
#Available at: https://github.com/amaurea/taylens

# This generates correlated T,E,B and Phi maps                                             
def simulate_tebp_correlated(cl_tebp_arr,nside,lmax) :
        alms=hp.synalm(cl_tebp_arr,lmax=lmax,new=True)
        aphi=alms[-1]
        acmb=alms[0:-1]
#Set to zero above map resolution to avoid aliasing                                        
        beam_cut=np.ones(3*nside)
        for ac in acmb :
                hp.almxfl(ac,beam_cut,inplace=True)
        cmb=np.array(hp.alm2map(acmb,nside,pol=True,verbose=False))

        return cmb,aphi

# This function is the core of Taylens.                                                              
def taylor_interpol_iter(m, pos, order=3, verbose=False, lmax=None):
        """Given a healpix map m[npix], and a set of positions                                       
        pos[{theta,phi},...], evaluate the values at those positions                                 
        using harmonic Taylor interpolation to the given order (3 by                                 
        default). Successively yields values for each cumulative order                               
        up to the specified one. If verbose is specified, it will print                              
        progress information to stderr."""
        nside = hp.npix2nside(m.size)
        if lmax is None: lmax = 3*nside
        # Find the healpix pixel centers closest to pos,                                             
        # and our deviation from these pixel centers.                                                
        ipos = hp.ang2pix(nside, pos[0], pos[1])
        pos0 = np.array(hp.pix2ang(nside, ipos))
        dpos = pos[:2]-pos0
        # Take wrapping into account                                                                 
        bad = dpos[1]>np.pi
        dpos[1,bad] = dpos[1,bad]-2*np.pi
        bad = dpos[1]<-np.pi
        dpos[1,bad] = dpos[1,bad]+2*np.pi

        # Since healpix' dphi actually returns dphi/sintheta, we choose                              
        # to expand in terms of dphi*sintheta instead.                                               
        dpos[1] *= np.sin(pos0[0])
        del pos0

        # We will now Taylor expand our healpix field to                                             
        # get approximations for the values at our chosen                                            
        # locations. The structure of this section is                                                
        # somewhat complicated by the fact that alm2map_der1 returns                                 
        # two different derivatives at the same time.                                                
        derivs = [[m]]
        res = m[ipos]
        yield res
        for o in range(1,order+1):
                # Compute our derivatives                                                            
                derivs2 = [None for i in range(o+1)]
                used    = [False for i in range(o+1)]
                # Loop through previous level in steps of two (except last)                          
                if verbose: tprint("order %d" % o)
                for i in range(o):
                        # Each alm2map_der1 provides two derivatives, so avoid                       
                        # doing double work.                                                         
                        if i < o-1 and i % 2 == 1:
                                continue
                        a = hp.map2alm(derivs[i], use_weights=True, lmax=lmax, iter=0)
                        derivs[i] = None
                        dtheta, dphi = hp.alm2map_der1(a, nside, lmax=lmax)[-2:]
                        derivs2[i:i+2] = [dtheta,dphi]
                        del a, dtheta, dphi
                        # Use these to compute the next level                                        
                        for j in range(i,min(i+2,o+1)):
                                if used[j]: continue
                                N = comb(o,j)/factorial(o)
                                res += N * derivs2[j][ipos] * dpos[0]**(o-j) * dpos[1]**j
                                used[j] = True
                                # If we are at the last order, we don't need to waste memory         
                                # storing the derivatives any more                                   
                                if o == order: derivs2[j] = None
                derivs = derivs2
                yield res

# The following functions are support routines for reading                                           
# input data and preparing it for being lensed. Most of them                                 
# are only needed to take care of tiny, curvature-related                                            
# effects that can be safely ignored.                                                                
def readspec(fname):
        """Read a power spectrum with columns [l,comp1,comp2,....]                                   
        into a 2d array indexed by l. Entries with missing data are                                  
        filled with 0."""
        tmp = np.loadtxt(fname).T
        l, tmp = tmp[0], tmp[1:]
        res = np.zeros((len(tmp),np.max(l)+1))
        res[:,np.array(l,dtype=int)] = tmp
        return res

def offset_pos(ipos, dtheta, dphi, pol=False, geodesic=False):
        """Offsets positions ipos on the sphere by a unit length step                                
        along the gradient dtheta, dphi/sintheta, taking the curvature                               
        of the sphere into account. If pol is passed, also computes                                  
        the cos and sin of the angle by which (Q,U) must be rotated to                               
        take into account the change in local coordinate system.                                     
                                                                                                     
        If geodesic is passed, a quick and dirty, but quite accurate, approximation                  
        is used.                                                                                     
                                                                                                     
        Uses the memory of 2 maps (4 if pol) (plus that of the input maps)."""
        opos = np.zeros(ipos.shape)
        if pol and not geodesic: orot = np.zeros(ipos.shape)
        else: orot = None
        if not geodesic:
                # Loop over chunks in order to conserve memory                                       
                step = 0x10000
                for i in range(0, ipos.shape[1], step):
                        small_opos, small_orot = offset_pos_helper(ipos[:,i:i+step], dtheta[i:i+step\
], dphi[i:i+step], pol)
                        opos[:,i:i+step] = small_opos
                        if pol: orot[:,i:i+step] = small_orot
        else:
                opos[0] = ipos[0] + dtheta
                opos[1] = ipos[1] + dphi/np.sin(ipos[0])
                opos = fixang(opos)
        return opos, orot

def offset_pos_helper(ipos, dtheta, dphi, pol):
        grad = np.array((dtheta,dphi))
        dtheta, dphi = None, None
        d = np.sum(grad**2,0)**0.5
        grad  /= d
        cosd, sind = np.cos(d), np.sin(d)
        cost, sint = np.cos(ipos[0]), np.sin(ipos[0])
        ocost  = cosd*cost-sind*sint*grad[0]
        osint  = (1-ocost**2)**0.5
        ophi   = ipos[1] + np.arcsin(sind*grad[1]/osint)
        if not pol:
                return np.array([np.arccos(ocost), ophi]), None
        A      = grad[1]/(sind*cost/sint+grad[0]*cosd)
        nom1   = grad[0]+grad[1]*A
        denom  = 1+A**2
        cosgam = 2*nom1**2/denom-1
        singam = 2*nom1*(grad[1]-grad[0]*A)/denom
        return np.array([np.arccos(ocost), ophi]), np.array([cosgam,singam])

def fixang(pos):
        """Handle pole wraparound."""
        a = np.array(pos)
        bad = np.where(a[0] < 0)
        a[0,bad] = -a[0,bad]
        a[1,bad] = a[1,bad]+np.pi
        bad = np.where(a[0] > np.pi)
        a[0,bad] = 2*np.pi-a[0,bad]
        a[1,bad] = a[1,bad]+np.pi
        return a

def apply_rotation(m, rot):
        """Update Q,U components in polarized map by applying                                        
        the rotation rot, representat as [cos2psi,sin2psi] per                                       
        pixel. Rot is one of the outputs from offset_pos."""
        if len(m) < 3: return m
        if rot is None: return m
        m = np.asarray(m)
        res = m.copy()
        res[1] = rot[0]*m[1]-rot[1]*m[2]
        res[2] = rot[1]*m[1]+rot[0]*m[2]
        return m

# Set up progress prints                                                                             
t0 = None
def silent(msg): pass
def tprint(msg):
        global t0
        if t0 is None: t0 = time.time()
        print >> sys.stderr, "%8.2f %s" % (time.time()-t0,msg)

