# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 09:52:05 2017
 
@author: useradmin
 
This 
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import yaml
from skimage.transform import warp, AffineTransform
import NanoImagingPack as nip
# own functions


is_debug = False 

class MuScatModel(object):
    def __init__(self, size=(256,256,256), lambda0=405e-9, nEmbb=1.33, dxy=1, dz=1, NAo=.5, NAc=.3, dn=0.1):
 
        ''' Create Multiple Scattering Class;
        INPUTS:
        my_mat_paras - parameter file from MATLAB containing the experimental details - has to be load previously!
        '''
        self.lambda0 = lambda0
        self.nEmbb = nEmbb
        self.lambdaM = self.lambda0/self.nEmbb
        self.dn = dn
        self.size = size
        self.dxy = dxy
        self.dz = dz
        self.NAo = NAo
        self.NAc = NAc
        
        
    def repmat4d(self,inputarr, n4dim):
        return np.tile(np.reshape(inputarr, [inputarr.shape[0], inputarr.shape[1], 1, 1]), [1, 1, 1, n4dim])
    
    
    #@define_scope
    def computesys(self, obj=None, is_tomo = False, is_compute_psf='BORN', is_mictype='BF'):
        # is_compute_psf: 'BORN', 'sep'
        """ This computes the FWD-graph of the Q-PHASE microscope;
        1.) Compute the physical dimensions
        2.) Compute the sampling for the waves
        3.) Create the illumination waves depending on System's properties
 
        # 'is_mictype' gives the possibility to switch between different microscope types like Brightfield, Phasecontrast etc. 
        BF - Brightfield 
        PC - Phase-Contrast (Zernike)
        DIC - Differential Interference Contrast
        DPC - Differential Phase Contrast (Oblique Illumination)
        AIP - Annular intensity projection      
        
        
        
        ##### IMPORTANT! ##### 
        The ordering of the channels is as follows:
            Nillu, Nz, Nx, Ny
        """
        # define whether we want to pad the experiment 
        self.is_mictype = is_mictype
        self.is_tomo = is_tomo
        self.is_compute_psf = is_compute_psf
            
        # Allocate memory for the object 
        if obj is None:
            self.obj = np.zeros(self.size)
        else:
            self.obj = obj
            
                
        ## Establish normalized coordinates.
        #-----------------------------------------
        vxx = nip.xx(self.size[1:],freq='ftfreq') * self.lambdaM * self.nEmbb / (self.dxy * self.NAo)
        vyy = nip.yy(self.size[1:],freq='ftfreq') * self.lambdaM * self.nEmbb / (self.dxy * self.NAo)
         
        # AbbeLimit=lambda0/NAo;  # Rainer's Method
        # RelFreq = rr(size,'freq')*AbbeLimit/dx;  # Is not generally right (dx and dy)
        self.RelFreq = np.sqrt(nip.abssqr(vxx) + nip.abssqr(vyy));    # spanns the frequency grid of normalized pupil coordinates
        self.Po = (self.RelFreq < 1.0);   # Create the pupil of the objective lens        
        self.Po = np.expand_dims(self.Po,0)
        '''
        if(1):
            print('lens shadig (e.g. vignatting is applied!!)')
            self.vignetting = np.exp(-nip.rr(size=self.Po.shape,freq='ftfreq')**2/0.125)
            #self.vignetting -= np.min(self.vignetting)
            self.Po *= self.vignetting
            self.Po /= np.max(self.Po)
        '''

        # Prepare the normalized spatial-frequency grid
        self.S = self.NAc/self.NAo   # Coherence factor
 
        # eventually introduce a phase factor to approximate the experimental data better
        print('----------> Be aware: We are taking aberrations into account!')
        # Assuming: System has coma along X-direction
        #self.myaberration = np.sum(self.zernikefactors * self.myzernikes, axis=2)
        #self.Po *= np.exp(1j*self.myaberration)
        
        # do individual pupil functoins according to is_mictype
        if(self.is_mictype=='BF' or self.is_mictype=='DF' or self.is_mictype=='AIP'):
            # Brightfield/Generic case
            self.Po = self.Po
            
        elif(self.is_mictype=='PC'):
            # Anullar Phase-plate with generic absorption
            print('We are taking a phase-ring in the Pupil plane into account!')
            p_o = (self.NAc*.9)/self.NAo;   # Coherence factor
            p_i = (self.NAci*1.1)/self.NAo;   # Coherence factor
            self.myphaseplate = (1.*(self.RelFreq < p_o) * 1.*(self.RelFreq > p_i))>0 # Create the pupil of the condenser plane
            self.Po *= np.exp(1j*np.pi/2*self.myphaseplate)
        elif(self.is_mictype=='DIC'):
            # DIC Phase mask from: https://github.com/mattersoflight/microlith/blob/master/%40microlith/computesys.m
            shearangle = 45
            shear = .48/2
            bias = 25
            freqGridShear=vxx*np.cos(shearangle*np.pi/180)+vyy*np.sin(shearangle*np.pi/180);
            halfshear=shear/2;
            halfbias=0.5*bias*np.pi/180;
            self.Po *= 1j*np.sin(2*np.pi*freqGridShear*halfshear-halfbias)
                    
        #self.Po = np.fft.fftshift(self.Po)# Need to shift it before using as a low-pass filter    Po=np.ones((np.shape(Po)))
        # do individual illumination sources according to is_mictype
        if(self.is_mictype=='BF' or self.is_mictype == 'DIC'):
            # Brightfield/Generic case
            self.Ic = self.RelFreq <= self.S
        elif(self.is_mictype=='PC' or self.is_mictype=='DF'):
        #if hasattr(self, 'NAci'):
            # Anullar illumination e.g. for PC or DF 
            if self.NAci == None or self.NAci < 0:
                self.NAci = self.NAc - .1
                print('I detected a darkfield illumination aperture, but value was not set! ')
            self.S_o = self.NAc/self.NAo;   # Coherence factor
            self.S_i = self.NAci/self.NAo;   # Coherence factor
            self.Ic = (1.*(self.RelFreq < self.S_o) * 1.*(self.RelFreq > self.S_i))>0 # Create the pupil of the condenser plane
        elif(self.is_mictype=='AIP'):
            self.Ic = (self.RelFreq <= self.S)
            self.Ic = np.roll(self.Ic,self.size[1]//4-4) # oblique illumination
             
        ## Forward propagator  (Ewald sphere based) DO NOT USE NORMALIZED COORDINATES HERE
        self.kxysqr= (nip.abssqr(nip.xx((self.size[1], self.size[2]), freq='ftfreq') / self.dxy) + 
                      nip.abssqr(nip.yy((self.size[1], self.size[2]), freq='ftfreq') / self.dxy)) + 0j
        self.k0=1/self.lambdaM
        self.kzsqr= nip.abssqr(self.k0) - self.kxysqr
        self.kz=np.sqrt(self.kzsqr)
        self.kz[self.kzsqr < 0]=0
        self.dphi = 2*np.pi*self.kz*self.dz # exp(1i*kz*dz) would be the propagator for one slice

        self.Nc=np.sum(self.Ic>0)
        print('Number of Illumination Angles / Plane waves: '+str(self.Nc))
        

        ## Get a list of vector coordinates corresponding to the pixels in the mask
        xfreq= nip.xx((self.size[1], self.size[2]),freq='ftfreq');
        yfreq= nip.yy((self.size[1], self.size[2]),freq='ftfreq');

        # Calculate the computatonal grid/sampling
        self.kxcoord = np.reshape(xfreq[self.Ic>0],[1, 1, 1, self.Nc]);    # NA-positions in condenser aperture plane in x-direction
        self.kycoord = np.reshape(yfreq[self.Ic>0],[1, 1, 1, self.Nc]);    # NA-positions in condenser aperture plane in y-direction
        self.RefrCos = np.reshape(self.k0/self.kz[self.Ic>0],[1, 1, 1, self.Nc]);   # 1/cosine used for the application of the refractive index steps to acount for longer OPD in medium under an oblique illumination angle
            
        ## Generate the illumination amplitudes
        self.intensityweights = self.Ic[self.Ic>0]
        self.A_input = self.intensityweights *np.exp((2*np.pi*1j) *
            (self.kxcoord * self.repmat4d(nip.xx((self.size[1], self.size[2])), self.Nc) 
            + self.kycoord * self.repmat4d(nip.yy((self.size[1], self.size[2])), self.Nc))) # Corresponds to a plane wave under many oblique illumination angles - bfxfun
            
        ## propagate field to z-stack and sum over all illumination angles
        self.Alldphi = -(np.reshape(np.arange(0, self.size[0], 1), [1, 1, self.size[0]]))*np.repeat(self.dphi[:, :, np.newaxis], self.size[0], axis=2)
            
        # Ordinary backpropagation. This is NOT what we are interested in:
        self.myAllSlicePropagator = np.transpose(np.exp(1j*self.Alldphi) * (np.repeat(self.dphi[:, :, np.newaxis], self.size[0], axis=2) >0), [2, 0, 1]);  # Propagates a single end result backwards to all slices

        ## propagate the field through the entire object for all angles simultaneously
        self.A_prop = np.transpose(self.A_input,[3, 0, 1, 2])  # ??????? what the hack is happening with transpose?!

        # Precalculate the oblique effect on OPD to speed it up
        print('--------> ATTENTION: I added no pi factor - is this correct?!')
        self.RefrEffect = np.reshape(np.squeeze(1j * 2* np.pi/self.lambda0* self.dz * self.RefrCos), [self.Nc, 1, 1]) # pi-factor
        
        myprop = np.exp(1j * self.dphi) * (self.dphi > 0);  # excludes the near field components in each step
        myprop = self.repmat4d(myprop, self.Nc)
        self.myprop = np.transpose(myprop, [3, 0, 1, 2])  # ??????? what the hack is happening with transpose?!
        
    
 
            
    #@define_scope
    def computemodel(self, is_allslices=False):
        ''' Perform Multiple Scattering here
        1.) Multiple Scattering is perfomed by slice-wise propagation the E-Field throught the sample
        2.) Each Field has to be backprojected to the BFP
        3.) Last step is to ceate a focus-stack and sum over all angles
 
        This is done for all illumination angles (coming from illumination NA
        simultaneasly)'''
        
        print("Buildup Q-PHASE Model ")
        ###### make sure, that the first dimension is "batch"-size; in this case it is the illumination number
        # @beniroquai It's common to have to batch dimensions first and not last.!!!!!
        # the following loop propagates the field sequentially to all different Z-planes
 
         
        if(self.is_tomo):
            print('Experimentally using the tomographic scheme!')
            self.A_prop = np.conj(self.A_prop)
             
        
        if self.is_compute_psf is 'BPM':
            A_input=self.A_prop
                
        # Initiliaze memory
        self.allSumAmp = 0
        
        # compute multiple scattering only in higher Born order        
        if self.is_compute_psf is 'BPM':
            A_prop = self.propslices(is_allslices)
        
        if self.is_compute_psf=='BPM':
            # The input field of the PSF calculation is the BFP of the objective lens
            A_prop = nip.ift(nip.ft(A_prop, axes=[1,2])*np.expand_dims(self.Po,-1), axes=[1,2])
            
        # Experimenting with pseudo tomographic data? No backpropgation necessary!
        if self.is_tomo:
            return self.computetomo()         
                
        if self.is_compute_psf is 'BPM':
            # now backpropagate the multiple scattered slice to a focus stack - i.e. convolve with PTF
            self.allSumAmp = self.propangles(A_prop)
            self.allSumAmp = self.allSumAmp*np.exp(-1j*np.pi/2) # align it with the born model
        else:
            if self.is_compute_psf=='BORN': 
                self.TF_ASF, self.TF_ATF = self.compute_psf()
            if self.is_compute_psf=='3QDPC': 
                self.TF_ASF, self.TF_ATF = self.compute_psf_3dqdpc()
            elif self.is_compute_psf=='sep': #
                self.TF_ASF, self.TF_ATF = self.compute_psf_sepang()

            # ASsign dummy variable- not used
            self.TF_allSumAmp = None 
            
            

        return self.allSumAmp

    def propslices(self,is_allslices=False):
        ''' Here we do the multiple scattering e.g. the forward propagation of 
        all light modes through the 3D sample. It's basically the split-step 
        fourier method (SSFM) beam propagation method (BPM). The last slice 
        holds all scattering of the sample which can be used to generate a 
        focus stack'''
        
        # simulate multiple scattering through object
    
        #print('---------ATTENTION: We are inverting the RI!')
        self.allslices=[]
        for pz in range(0, self.size[0]):
            f = np.exp(-self.RefrEffect*self.obj[-pz,:,:])
            self.A_prop = self.A_prop * np.expand_dims(f,-1)  # refraction step
            self.A_prop = nip.ift(nip.ft(self.A_prop,axes=[1,2]) * self.myprop,axes=[1,2]) # diffraction step 
            if is_allslices:
                self.allslices.append(self.A_prop)
        return self.A_prop
        
    def propangles(self, A_prop):
        ''' Here we generate a focus stack for all illumination modes independtly. 
        it follows the routine "propslices". '''
        
        # create Z-Stack by backpropagating Information in BFP to Z-Position
        self.kzcoord = np.reshape(self.kz[self.Ic>0], [1, 1, 1, self.Nc]);
 
        # self.mid3D = ([np.int(np.ceil(self.A_input.shape[0] / 2) + 1), np.int(np.ceil(self.A_input.shape[1] / 2) + 1), np.int(np.ceil(self.size[0] / 2) + 1)])
        self.mid3D = ([np.int32(self.size[0]//2), np.int32(self.A_input.shape[0] // 2), np.int32(self.A_input.shape[1]//2)])
        allSumAmp=[]
        for pillu in range(0, self.Nc):
            OneAmp = np.expand_dims(A_prop[pillu, :, :], 0)
            # Fancy backpropagation assuming what would be measured if the sample was moved under oblique illumination:
            # The trick is: First use conceptually the normal way
            # and then apply the XYZ shift using the Fourier shift theorem (corresponds to physically shifting the object volume, scattered field stays the same):
            self.AdjustKXY = np.squeeze(np.conj(self.A_input[:,:,:,pillu])) # Maybe a bit of a dirty hack, but we first need to shift the zero coordinate to the center
            self.AdjustKZ = np.exp(np.transpose(
                2 * np.pi * 1j * self.dz * np.reshape(np.arange(0, self.size[0], 1), # We want to start from first Z-slice then go to last which faces the objective lens
                        [1, 1, self.size[0]]) * self.kzcoord[:, :, :,pillu], [2, 1, 0]))
            self.allAmp = nip.ift(nip.ft2d(np.squeeze(OneAmp)) * self.myAllSlicePropagator,axes=[1,2]) * self.AdjustKZ * self.AdjustKXY  # * (TF_AdjustKZ);  # 2x bfxfun.  Propagates a single amplitude pattern back to the whole stack
            #tf_global_phase = tf.cast(tf.angle(self.TF_allAmp[self.mid3D[0],self.mid3D[1],self.mid3D[2]]), tf.complex64)
            #tf_global_phase = tf.cast(np.random.randn(1)*np.pi,tf.complex64)
            #self.TF_allAmp = self.TF_allAmp * tf.exp(1j*tf_global_phase) # Global Phases need to be adjusted at this step!  Use the zero frequency
            
            
    
            '''
            self.allAmp_3dft = np.fft.fftn(np.expand_dims(self.allAmp, axis=0))
            #tf_global_phase = tf.angle(self.TF_allAmp_3dft[0,0,0,0])#tf.angle(self.TF_allAmp_3dft[0, self.mid3D[2], self.mid3D[1], self.mid3D[0]])
            global_phase = np.angle(self.allAmp_3dft[0,0,0,0])#tf.angle(self.TF_allAmp_3dft[0, self.mid3D[2], self.mid3D[1], self.mid3D[0]])
            global_phase = f_global_phase, tf.complex64)

            self.TF_allAmp = self.TF_allAmp * tf.exp(-1j * tf_global_phase);  # Global Phases need to be adjusted at this step!  Use the zero frequency
            #print('Global phase: '+str(tf.exp(1j*tf.cast(tf.angle(self.TF_allAmp[self.mid3D[0],self.mid3D[1],self.mid3D[2]]), tf.complex64).eval()))
            '''
        
            allSumAmp.append(self.allAmp) #/ self.intensityweights[pillu];  # Superpose the Amplitudes
            
            # Normalize the image such that the values do not depend on the fineness of
            # the source grid.
            #TF_allSumAmp = TF_allSumAmp/tf.cast(np.sum(self.Ic), tf.complex64) # tf.cast(tf.reduce_max(tf.abs(self.TF_allSumAmp)), tf.complex64) # self.Nc #/

        # Normalize along Z to account for energy conservation
        #TF_mynorm = tf.cast(tf.sqrt(tf_helper.tf_abssqr(tf.reduce_sum(TF_allSumAmp , axis=(1,2))))/np.prod(self.size[1:3]),tf.complex64)
        
        allSumAmp = np.sum(np.array(allSumAmp), 0)
        
        mynorm = np.sqrt(nip.abssqr(np.mean(allSumAmp , axis=(1,2))))
        mynorm = np.expand_dims(np.expand_dims(mynorm,1),1)
        allSumAmp = allSumAmp/mynorm;
        print('BPM Normalization accounts for ENERGY conservation!!')

        # Following is the normalization according to Martin's book. It ensures
        # that a transparent specimen is imaged with unit intensity.
        # normfactor=abs(Po).^2.*abs(Ic); We do not use it, because it leads to
        # divide by zero for dark-field system. Instead, through normalizations
        # perfomed above, we ensure that image of a point under matched
        # illumination is unity. The brightness of all the other configurations is
        # relative to this benchmark.
        return allSumAmp
     
    def computetomo(self):
        print('Only Experimental! Tomographic data?!')
        # Bring the slice back to focus - does this make any sense?! 
        print('----------> Bringing field back to focus')
        TF_centerprop = tf.exp(-1j*tf.cast(self.Nz/2*tf.angle(self.TF_myprop), tf.complex64))
        self.TF_A_prop = tf.ifft2d(tf.fft2d(self.TF_A_prop) * TF_centerprop) # diffraction step
        return self.TF_A_prop


    def addRegularizer(self, is_tv, is_gr, is_pos):
        print('Do stuff')
 
    def defineWriter(self, logs_path = '/tmp/tensorflow_logs/example/'):
        # create writer object
        self.logs_path = logs_path 
        self.writer = tf.summary.FileWriter(logs_path, graph=tf.get_default_graph())
         
     
    def writeParameterFile(self, mylr, mytv, myeps, filepath = '/myparameters.yml'):
        ''' Write out all parameters to a yaml file in case we need it later'''
        mydata = dict(
                shiftIcX = float(self.shiftIcX),
                shiftIcY = float(self.shiftIcY),                
                NAc = float(self.NAc),
                NAo = float(self.NAo), 
                Nc = float(self.Nc), 
                Nx = float(self.Nx), 
                Ny = float(self.Ny),
                Nz = float(self.Nz),
                dx = float(self.dxy),
                dy = float(self.dy),
                dz = float(self.dz),
                dn = float(self.dn),
                lambda0 = float(self.lambda0),
                lambdaM = float(self.lambdaM),
                learningrate = mylr, 
                lambda_reg = mytv, 
                eps_tv = myeps) 
                #zernikfactors = float(self.zernikefactors))
 
        with open(filepath, 'w') as outfile:
                yaml.dump(mydata, outfile, default_flow_style=False)
                 
    
            
            
         
    def saveFigures_list(self, savepath, myfwdlist, mylosslist, myfidelitylist, myneglosslist, mytvlosslist, result_phaselist, result_absorptionlist, 
                              globalphaselist, globalabslist, mymeas, figsuffix):
        ''' We want to save some figures for debugging purposes'''

        # get latest resutl from lists
        myfwd = myfwdlist[-1]
        my_res = result_phaselist[-1]
        my_res_absorption = result_absorptionlist[-1]
        
             
        plt.figure()
        plt.subplot(231), plt.title('REAL XZ'),plt.imshow(np.real(myfwd)[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('REAL YZ'),plt.imshow(np.real(myfwd)[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('REAL XY'),plt.imshow(np.real(myfwd)[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(234), plt.title('Imag XZ'),plt.imshow(np.imag(myfwd)[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Imag XZ'),plt.imshow(np.imag(myfwd)[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Imag XY'),plt.imshow(np.imag(myfwd)[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.savefig(savepath+'/res_myfwd_realimag'+figsuffix+'.png'), plt.show()

        if(0):
            plt.figure()    
            plt.subplot(231), plt.title('ABS XZ'),plt.imshow(np.abs(myfwd)[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
            plt.subplot(232), plt.title('ABS YZ'),plt.imshow(np.abs(myfwd)[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
            plt.subplot(233), plt.title('ABS XY'),plt.imshow(np.abs(myfwd)[myfwd.shape[0]//2 ,:,:]), plt.colorbar()# plt.show()
            #myfwd=myfwd*np.exp(1j*2)
            plt.subplot(234), plt.title('Angle XZ'),plt.imshow(np.angle(myfwd)[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
            plt.subplot(235), plt.title('Angle YZ'),plt.imshow(np.angle(myfwd)[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
            plt.subplot(236), plt.title('Angle XY'),plt.imshow(np.angle(myfwd)[myfwd.shape[0]//2 ,:,:]), plt.colorbar(), plt.show()
            plt.savefig(savepath+'/res_myfwd_ampphase'+figsuffix+'.png'), plt.show()
     
        # This is the measurment
        plt.figure()
        plt.subplot(231), plt.title('REAL XZ'),plt.imshow(np.real(mymeas)[:,mymeas.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('REAL YZ'),plt.imshow(np.real(mymeas)[:,:,mymeas.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('REAL XY'),plt.imshow(np.real(mymeas)[mymeas.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(234), plt.title('Imag XZ'),plt.imshow(np.imag(mymeas)[:,mymeas.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Imag XZ'),plt.imshow(np.imag(mymeas)[:,:,mymeas.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Imag XY'),plt.imshow(np.imag(mymeas)[mymeas.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.savefig(savepath+'/res_mymeas'+figsuffix+'.png'), plt.show()
     
        # This is the residual
        myresi = tf_helper.abssqr(myfwd-mymeas)
        plt.figure()
        plt.subplot(331), plt.title('Residual REAL XZ'),plt.imshow((np.real(myfwd)-np.real(mymeas))[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(332), plt.title('Residual REAL YZ'),plt.imshow((np.real(myfwd)-np.real(mymeas))[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(333), plt.title('Residual REAL XY'),plt.imshow((np.real(myfwd)-np.real(mymeas))[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(334), plt.title('Residual Imag XZ'),plt.imshow((np.imag(myfwd)-np.imag(mymeas))[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(335), plt.title('Residual Imag XZ'),plt.imshow((np.imag(myfwd)-np.imag(mymeas))[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(336), plt.title('Residual Imag XY'),plt.imshow((np.imag(myfwd)-np.imag(mymeas))[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(337), plt.title('Residual abssqr XZ'), plt.imshow((((myresi)**.1))[:,myresi.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(338), plt.title('Residual abssqr XY'), plt.imshow((((myresi)**.1))[:,:,myresi.shape[2]//2]), plt.colorbar()#, plt.show()    
        plt.subplot(339), plt.title('Residual abssqr Yz'), plt.imshow((((myresi)**.1))[myresi.shape[0]//2,:,:]), plt.colorbar()#, plt.show()    
        

        plt.savefig(savepath+'/res_myresidual'+figsuffix+'.png'), plt.show()
         
        # diplay the error over time
        plt.figure()
        plt.subplot(231), plt.title('Error/Cost-function'), plt.semilogy((np.array(mylosslist)))#, plt.show()
        plt.subplot(232), plt.title('Fidelity-function'), plt.semilogy((np.array(myfidelitylist)))#, plt.show()
        plt.subplot(233), plt.title('Neg-loss-function'), plt.plot(np.array(myneglosslist))#, plt.show()
        plt.subplot(234), plt.title('TV-loss-function'), plt.semilogy(np.array(mytvlosslist))#, plt.show()
        plt.subplot(235), plt.title('Global Phase'), plt.plot(np.array(globalphaselist))#, plt.show()
        plt.subplot(236), plt.title('Global ABS'), plt.plot(np.array(globalabslist))#, plt.show()
        plt.savefig(savepath+'/myplots'+figsuffix+'.png'), plt.show()
         
        # Display RI result
        plt.figure()
        plt.subplot(231), plt.title('Result Phase: XZ'),plt.imshow(my_res[:,my_res.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('Result Phase: XZ'),plt.imshow(my_res[:,:,my_res.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('Result Phase: XY'),plt.imshow(my_res[my_res.shape[0]//2,:,:]), plt.colorbar()
        plt.subplot(234), plt.title('Result Abs: XZ'),plt.imshow(my_res_absorption[:,my_res.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Result abs: XZ'),plt.imshow(my_res_absorption[:,:,my_res.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Result abs: XY'),plt.imshow(my_res_absorption[my_res.shape[0]//2,:,:]), plt.colorbar()
        plt.savefig(savepath+'/RI_abs_result'+figsuffix+'.png'), plt.show()
         
                 
    def saveFigures(self, sess, savepath, tf_fwd, np_meas, mylosslist, myfidelitylist, myneglosslist, mytvlosslist, globalphaselist, globalabslist, 
                    result_phaselist=None, result_absorptionlist=None, init_guess=None, figsuffix=''):
        ''' We want to save some figures for debugging purposes'''
        # This is the reconstruction
        if(init_guess is not None):
            myfwd, mymeas, my_res, my_res_absorption, myzernikes = sess.run([tf_fwd, self.tf_meas, self.TF_obj, self.TF_obj_absorption, self.TF_zernikefactors], 
                    feed_dict={self.tf_meas:np_meas, self.TF_obj:np.real(init_guess), self.TF_obj_absorption:np.imag(init_guess)})
        else:
            myfwd, mymeas, my_res, my_res_absorption, myzernikes = sess.run([tf_fwd, self.tf_meas, self.TF_obj, self.TF_obj_absorption, self.TF_zernikefactors], feed_dict={self.tf_meas:np_meas})
             
        plt.figure()
        plt.subplot(231), plt.title('ABS XZ'),plt.imshow(np.abs(myfwd)[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('ABS YZ'),plt.imshow(np.abs(myfwd)[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('ABS XY'),plt.imshow(np.abs(myfwd)[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(234), plt.title('Angle XZ'),plt.imshow(np.angle(myfwd)[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Angle XZ'),plt.imshow(np.angle(myfwd)[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Angle XY'),plt.imshow(np.angle(myfwd)[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.savefig(savepath+'/res_myfwd'+figsuffix+'.png'), plt.show()
     
        # This is the measurment
        plt.figure()
        plt.subplot(231), plt.title('ABS XZ'),plt.imshow(np.abs(mymeas)[:,mymeas.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('ABS YZ'),plt.imshow(np.abs(mymeas)[:,:,mymeasomcputeshape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('ABS XY'),plt.imshow(np.abs(mymeas)[mymeas.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(234), plt.title('Angle XZ'),plt.imshow(np.angle(mymeas)[:,mymeas.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Angle XZ'),plt.imshow(np.angle(mymeas)[:,:,mymeas.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Angle XY'),plt.imshow(np.angle(mymeas)[mymeas.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.savefig(savepath+'/res_mymeas'+figsuffix+'.png'), plt.show()
     
        # This is the residual
        plt.figure()
        plt.subplot(231), plt.title('Residual ABS XZ'),plt.imshow((np.abs(myfwd)-np.abs(mymeas))[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('Residual ABS YZ'),plt.imshow((np.abs(myfwd)-np.abs(mymeas))[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('Residual ABS XY'),plt.imshow((np.abs(myfwd)-np.abs(mymeas))[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.subplot(234), plt.title('Residual Angle XZ'),plt.imshow((np.angle(myfwd)-np.angle(mymeas))[:,myfwd.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Residual Angle XZ'),plt.imshow((np.angle(myfwd)-np.angle(mymeas))[:,:,myfwd.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Residual Angle XY'),plt.imshow((np.angle(myfwd)-np.angle(mymeas))[myfwd.shape[0]//2,:,:]), plt.colorbar()#, plt.show()
        plt.savefig(savepath+'/res_myresidual'+figsuffix+'.png'), plt.show()
         
        # diplay the error over time
        plt.figure()
        plt.subplot(231), plt.title('Error/Cost-function'), plt.semilogy((np.array(mylosslist)))#, plt.show()
        plt.subplot(232), plt.title('Fidelity-function'), plt.semilogy((np.array(myfidelitylist)))#, plt.show()
        plt.subplot(233), plt.title('Neg-loss-function'), plt.plot(np.array(myneglosslist))#, plt.show()
        plt.subplot(234), plt.title('TV-loss-function'), plt.semilogy(np.array(mytvlosslist))#, plt.show()
        plt.subplot(235), plt.title('Global Phase'), plt.plot(np.array(globalphaselist))#, plt.show()
        plt.subplot(236), plt.title('Global ABS'), plt.plot(np.array(globalabslist))#, plt.show()
        plt.savefig(savepath+'/myplots'+figsuffix+'.png'), plt.show()
         
        # Display RI result
        plt.figure()
        plt.subplot(231), plt.title('Result Phase: XZ'),plt.imshow(my_res[:,my_res.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(232), plt.title('Result Phase: XZ'),plt.imshow(my_res[:,:,my_res.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(233), plt.title('Result Phase: XY'),plt.imshow(my_res[my_res.shape[0]//2,:,:]), plt.colorbar()
        plt.subplot(234), plt.title('Result Abs: XZ'),plt.imshow(my_res_absorption[:,my_res.shape[1]//2,:]), plt.colorbar()#, plt.show()
        plt.subplot(235), plt.title('Result abs: XZ'),plt.imshow(my_res_absorption[:,:,my_res.shape[2]//2]), plt.colorbar()#, plt.show()
        plt.subplot(236), plt.title('Result abs: XY'),plt.imshow(my_res_absorption[my_res.shape[0]//2,:,:]), plt.colorbar()
        plt.savefig(savepath+'/RI_abs_result'+figsuffix+'.png'), plt.show()
         
        # Display recovered Pupil
        plt.figure()
        myshiftX = sess.run(self.TF_shiftIcX)
        myshiftY = sess.run(self.TF_shiftIcY)
        plt.subplot(131), plt.title('Po Phase'), plt.imshow(np.angle(sess.run(self.TF_Po_aberr))), plt.colorbar()
        plt.subplot(132), plt.title('Ic, shiftX: '+str(myshiftX)+' myShiftY: '+str(myshiftY)), plt.imshow(np.abs(sess.run(self.TF_Po_aberr))), plt.colorbar()
        plt.subplot(133), plt.bar(np.linspace(1, np.squeeze(myzernikes.shape), np.squeeze(myzernikes.shape)), myzernikes, align='center', alpha=0.5)
        plt.ylabel('Zernike Values')
        plt.title('Zernike Coefficients (Noll)')
        plt.savefig(savepath+'/recovered_pupil'+figsuffix+'.png'), plt.show()

        # Eventually write H5 stacks to disc
        if(result_phaselist is not None): data.export_realdatastack_h5(savepath+'/myrefractiveindex'+figsuffix+'.h5', 'temp', np.array(result_phaselist))
        if(result_absorptionlist is not None): data.export_realdatastack_h5(savepath+'/myrefractiveindex_absorption'+figsuffix+'.h5', 'temp', np.array(result_absorptionlist))
        #myobj = my_res+1j*my_res_absorption
        #np.save('my_res_cmplx', myobj)         
                 
    def initObject():
        print('To Be done')
                    
