from unidec_modules import mzmlparse_auto as automzml
import unidec
from unidec_modules import unidectools as ud
from unidec_modules.ChromEng import *
from UniChrom2 import *
from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
from metaunidec.meta_import_wizard.meta_import_wizard import *  
from metaunidec.mudeng import *

import plate_map as pm

class MegaUniDec():
    def __init__(self):
        # initiate UniChrome and MetaUniDec
        self.chrom = ChromEngine()
        self.meta = MetaUniDec()
        self.peaks = []

        # Set default params
        # TIC processing
        self.ticlb = 40
        self.ticub = 0 
        self.fwhmfrac = 0.7 # amount of each peak to select as a fraction of its FWHM
        self.fwhmoffset = -0.1 # % offset of peak window selection (earlier tends to be better)
        self.peakwindow = 10 # detection window(local max plus or minus window)
        self.peakthresh = 0.4 # above threshold of peakthresh * max data intensity
        # self.ticpeaks = None

        # manual peak selection params
        self.totalpeaks = None
        self.firstpeak = None
        self.peakspacing = None
        self.peaktotal = None

        # plate map 
        self.header_names = {'Well ID': {'dtype':str, 'long':True, 'short_row': False, 'short_col':False},
                        'Type': {'dtype':str, 'long':True, 'short_row': True, 'short_col':True},
                        'Protein Name': {'dtype':str, 'long':True, 'short_row': True, 'short_col':True},
                        'Protein Concentration': {'dtype':float, 'long':True, 'short_row': True, 'short_col':True},
                        'Reagent Name': {'dtype':str, 'long':True, 'short_row': True, 'short_col':True},
                        'Reagent Concentration': {'dtype':float, 'long':True, 'short_row': True, 'short_col':True},
                        'Concentration Units':{'dtype':str, 'long':True, 'short_row': True, 'short_col':True},
                        'Time':{'dtype':float, 'long':True, 'short_row': True, 'short_col':True},
                            }
        self.plate_size = 96
        self.plate_map_bool = False

    def import_mzml(self, path, show_tic = True):
        
        self.path = path 
        self.dirname, self.filename = os.path.split(path)
        self.chrom.load_mzml(path, load_hdf5=False) # change this eventually

        self.tic = self.chrom.tic

        # make new folder to house megaunidec files 
        self.folder = self.dirname+"\\megaunidecfiles_" + self.filename
        Path(self.folder).mkdir(parents = True, exist_ok = True)

        if show_tic == True:
            self.plot_tic()     

    def plot_tic(self, show_peaks = False, show_windows = False, *args, **kwargs):
        plt.plot(self.tic[:, 0], self.tic[:, 1])

        if show_peaks == True:
            plt.scatter(self.ticpeaks[:, 0], self.ticpeaks[:, 1], marker = 'x', *args, **kwargs)
        if show_windows == True:
            for w in self.chrompeaks_tranges:
                plt.axvspan(w[0], w[1], alpha = 0.3, color = 'orange') 


        plt.ylabel("Intensity / a.u.")
        plt.xlabel("retention time / mins")
        plt.show() 

    def update_peaks(self):
        self.peaktimes = [self.firstpeak + i*self.peakspacing for i in range(self.peaktotal)]
        print(self.peaktimes)

    def define_peaks(self, firstpeak, spacing, total, ll = 5, ul = 5):
        """Manually defines peak windows in TIC"""
        
        self.firstpeak, self.peakspacing, self.peaktotal = firstpeak, spacing, total
        self.peakwindowll = ll
        self.peakwindowul = ul

        self.update_peaks()

        self.chrompeaks_tranges = [[(i-self.peakwindowll)/60, (i+self.peakwindowul)/60] for i in self.peaktimes]
        self.chrom.chrompeaks_tranges = self.chrompeaks_tranges # just to be safe and if chom eng is needed again

    def pick_peaks(self, show_tic = False):

        """Picks peaks in TIC and calculates their FWHM"""
        self.ticpeaks =  ud.peakdetect(self.tic, window = self.peakwindow, threshold = self.peakthresh)

        self.ticFWHM = []
        for p in self.ticpeaks:
            self.ticFWHM.append(ud.calc_FWHM(p[0], self.tic))


        if show_tic == True:
            self.plot_tic(show_peaks=True)

    def get_peak_windows(self, show_tic = False):

        self.pick_peaks()

        self.chrompeaks_tranges = []
        for t in self.ticFWHM:
            frac = (t[0] - self.fwhmfrac*t[0])/2 # get fraction of FWHM to add/subtract to x1 and x2
            lb = t[1][0] + frac  # get new x1 and x2 values for peak window
            ub = t[1][1] - frac

            self.chrompeaks_tranges.append([lb, ub])
        self.chrom.chrompeaks_tranges = self.chrompeaks_tranges

        if show_tic == True:
            self.plot_tic(show_peaks = True, show_windows = True)

    def extract_peaks(self, show_spectra = False):
        """Extract defined windows and export to text files"""
        
        for t in self.chrompeaks_tranges:
            data = self.chrom.get_data_from_times(t[0], t[1])
            filename = "{:.3f}-{:.3f}-mean.txt".format(t[0], t[1])
            np.savetxt(self.folder+"\\"+filename, data, delimiter = "\t")

            if show_spectra == True:
                pass
    
    def import_plate_map(self, platemap, size = 96, map_type = 'long'):
        
        if map_type == 'short':
            platemap = pm.short_map(platemap, size = size)

        if map_type == 'long':
            platemap = pm.plate_map(platemap, size = size)

        # remove empty wells
        filt = platemap['Type'] != 'empty'
        used = platemap[filt] # get used wells
        self.wells = list(used.index) # get well ID's
        self.time = np.array(used['Time'], dtype = float)

        self.plate_map_bool = True


    def to_hdf5(self, folder = None, hdf5_name = None, var1 = 0, var2 = 0,
                filenames = None, overwrite = True):
        """Exports text files to hdf5"""      

        if hdf5_name == None:
            if self.filename[-5:].lower() == '.mzml':
                self.hdf5_name = self.dirname + "\\" + self.filename[:-5] + ".hdf5"
            else:
                self.hdf5_name = self.dirname + "\\" + self.filename + ".hdf5"
            
        else:
            self.hdf5_name = hdf5_name

        if overwrite == True: # stops some errors
            if os.path.exists(self.hdf5_name) == True:
                os.remove(self.hdf5_name)
        
        if folder == None:
            folder = self.folder

        self.meta.data.new_file(self.hdf5_name)

        for i, filename in enumerate(os.listdir(folder)): # add option to specify files to include?
            if filename[-8:] == 'mean.txt':
                peakll, peakull, temp = filename.split('-')

                path = folder+"/"+filename

                self.meta.data.add_file(path = path)
                typelist = [float, int] # if vars not float/int, assume to be iterable TODO: add try except clause 
                if type(var1) not in typelist:
                    self.meta.data.spectra[-1].var1 = var1[i]
                else:
                    self.meta.data.spectra[-1].var1 = var1
                if type(var2) not in typelist:
                    self.meta.data.spectra[-1].var2 = var2[i]
                else:
                    self.meta.data.spectra[-1].var2 = var2
                if filenames != None:
                    self.meta.data.spectra[-1].name = filenames[i]
                else:
                    self.meta.data.spectra[-1].name = filename 
        
        self.meta.data.export_hdf5()
        print("{} created".format(self.hdf5_name))

    def to_meta(self, path = None):
        """uploads HDF5 file to metaunidec for processing and deconvolution"""

        # interpretation of UniDecApp.on_open()
        self.meta.config.dirname = self.dirname
        if path == None:
            self.meta.config.hdf_file = self.hdf5_name
        else:
            self.meta.config.hdf_file = path


        self.meta.open(self.meta.config.hdf_file)

        # set new default config 
        # -- Data processing 
        self.meta.config.subtype = 2 # background subtraction - subtract curved
        self.meta.config.subbuff = 100 # background subtraction amount(subtract curved) 0 = 0ff, 100 = good amount when on
        self.meta.config.datanorm = 0 # turn off data normalisation

        # -- Deconvolution
        self.meta.config.numit = 100 # number of iterations

        # mass range (default = 5000.0 to 500000.0 Da)
        self.meta.config.massub = 15000 # upper 
        self.meta.config.masslb = 11000 # lower

        self.meta.config.massbins = 0.1 # sample mass every 0.1 Da

        # FWHM 
        # eng.get_auto_peak_width()
        self.meta.config.mzsig = 0 

        # charge range
        self.meta.config.startz = 1
        self.meta.config.endz = 30

        # smoothing 
        self.meta.config.zzsig = 1 # charge smooth width (smooth charge state distributions)
        self.meta.config.psig = 1 # smooth nearby points (point smooth width, some = 1)
        self.meta.config.beta = 0 # suppress artifacts (beta, some = 50)

        self.meta.config.psfun = 0 # Peak shape function (gaussian, lorentzian, split G/L)

        # -- Peak Selection and plotting
        self.meta.config.peaknorm = 0 # Normalise peaks (0 = off)
        self.meta.config.datanorm = 0
        self.meta.config.peakwindow = 100 # peak window / Da
        self.meta.config.exnorm = 0 # extract normalisation
        self.meta.config.peakthresh = 0.05
        self.meta.config.nativeub = 10
        self.meta.config.nativelb = -10
        self.meta.config.write_hdf5() # update config to HDF5

        return self.meta


    def to_unidec(self):
        pass
        


        








        

        





