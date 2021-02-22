from unidec_modules import mzmlparse_auto as automzml
import unidec
from unidec_modules.ChromEng import *
from UniChrom2 import *
from pathlib import Path
import os
import numpy as np
import matplotlib.pyplot as plt
from metaunidec.meta_import_wizard.meta_import_wizard import *  
from metaunidec.mudeng import *

class MegaUniDec():
    def __init__(self):
        # initiate UniChrome and MetaUniDec
        self.chrom = ChromEngine()
        self.meta = MetaUniDec()

    def import_mzml(self, path, show_tic = True):
        
        self.path = path 
        self.dirname, self.filename = os.path.split(path)
        self.chrom.load_mzml(path, load_hdf5=False) # change this eventually

        self.tic = self.chrom.tic

        if show_tic == True:

            plt.plot(self.tic[:, 0], self.tic[:, 1])
            plt.ylabel("Intensity / a.u.")
            plt.xlabel("retention time / s")
            plt.show()


    def update_peaks(self):
        self.peaktimes = [self.firstpeak + i*self.peakspacing for i in range(self.totalpeaks)]
        print(self.peaktimes)

    def define_peaks(self, firstpeak, spacing, total, ll = 5, ul = 5):
        """Manually select peaks in TIC"""
        self.firstpeak, self.peakspacing, self.totalpeaks = firstpeak, spacing, total
        self.peakwindowll = ll
        self.peakwindowul = ul

        self.update_peaks()

        self.chrompeaks_tranges = [[(i-self.peakwindowll)/60, (i+self.peakwindowul)/60]for i in self.peaktimes]
        self.chrom.chrompeaks_tranges = self.chrompeaks_tranges # just to be safe and if chom eng is needed again


    def extract_peaks(self, show_spectra = False):
        """Extract defined windows and export to text files"""
        # make new folder 
        self.folder = self.dirname+"\\megaunidecfiles_" + self.filename
        Path(self.folder).mkdir(parents = True, exist_ok = True)
        for i, t in enumerate(self.chrompeaks_tranges):
            data = self.chrom.get_data_from_times(t[0], t[1])
            filename = "{:.3f}-{:.3f}-mean.txt".format(t[0], t[1])
            np.savetxt(self.folder+"\\"+filename, data, delimiter = "\t")

            if show_spectra == True:
                pass
    
    def to_hdf5(self, folder = None, hdf5_name = None, plate_map = None, var1 = 0, var2 = 0,
                overwrite = True):
        """Exports text files to hdf5"""

        if plate_map == None: # TO DO: Add plate map stuff
            from_pm = 0 

        if hdf5_name == None:
            if self.filename[-5:].lower() == '.mzml':
                self.hdf5_name = self.dirname + "\\" + self.filename[:-5] + ".hdf5"
            else:
                self.hdf5_name = self.dirname + "\\" + self.filename + ".hdf5"
            
        else:
            self.hdf5_name = hdf5_name

        if overwrite == True:
            if os.path.exists(self.hdf5_name) == True:
                os.remove(self.hdf5_name)
        
        if folder == None:
            folder = self.folder

        self.meta.data.new_file(self.hdf5_name)
        for filename in os.listdir(folder): # add option to specify files to include?
            if filename[-8:] == 'mean.txt':
                peakll, peakull, temp = filename.split('-')

                v1, v2, path = var1, from_pm, folder+"/"+filename # change to add option to specify unique variables in each spectrum uploaded to HDF5

                self.meta.data.add_file(path = path)
                self.meta.data.spectra[-1].var1 = v1
                self.meta.data.spectra[-1].var2 = v2
                self.meta.data.spectra[-1].name = filename 
        
        self.meta.data.export_hdf5()
        print("{} created".format(self.hdf5_name))

    def to_meta(self, path = None):
        """uploads HDF5 file to metaunidec for processing and deconvolution"""

        # (interpretation of UniDecApp.on_open())
        
        
        self.meta.config.dirname = self.dirname
        if path == None:
            self.meta.config.hdf_file = self.hdf5_name
        else:
            self.meta.config.hdf_file = path


        self.meta.open(self.meta.config.hdf_file)

        return self.meta


    def to_unidec(self):
        pass
        


        








        

        





