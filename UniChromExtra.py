import plate_map as pm 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from unidec_modules import unidectools as ud
from copy import deepcopy

rheaders = {"Reaction":str, "Species":str, "Concentration":float,
           "Units":str, "Mass":float, "Reagent Type":str, "Sequence":str}

pheaders = {"Well ID":{'dtype':str, 'long':True, 'short_row': False, 'short_col':False}, 
            "Type":{'dtype':str, 'long':True, 'short_row': True, 'short_col':True}, 
            "Reaction":{'dtype':str, 'long':True, 'short_row': True, 'short_col':True}, 
            "Time":{'dtype':str, 'long':True, 'short_row': True, 'short_col':True}}



def read_in_long(map_path, sheet_names = ['plate map', 'species map']):
    """

    Note - plate map must come first in sheet_names. 

    """
    rmap = pd.read_excel(map_path, 
                        dtype = rheaders,
                        sheet_name = sheet_names[1])

    pdata_types = {i[0]: i[1]['dtype'] for i in pheaders.items()} 
    pmap = pd.read_excel(map_path, dtype = pdata_types, index_col = "Well ID", skiprows = 1)

    # add valid column
    pheaders2 = [x for x in pheaders.keys() if pheaders[x]['long']]
    pmapdf = pm.empty_map(size = 6, header_names = pheaders)
    pmapdf.update(pmap)

    return rmap, pmapdf

# TODO: short map 

def update_config(eng):
    eng.config.subtype = 2 # background subtraction - subtract curved
    eng.config.subbuff = 100 # background subtraction amount(subtract curved) 0 = 0ff, 100 = good amount when on
    eng.config.datanorm = 0 # turn off data normalisation

    # -- Deconvolution
    eng.config.numit = 100 # number of iterations

    # mass range (default = 5000.0 to 500000.0 Da)
    eng.config.massub = 15000 # upper 
    eng.config.masslb = 11000 # lower

    eng.config.massbins = 0.1 # sample mass every 0.1 Da

    # FWHM 
    # eng.get_auto_peak_width()
    eng.config.mzsig = 0 

    # charge range
    eng.config.startz = 1
    eng.config.endz = 30

    # smoothing 
    eng.config.zzsig = 1 # charge smooth width (smooth charge state distributions)
    eng.config.psig = 1 # smooth nearby points (point smooth width, some = 1)
    eng.config.beta = 0 # suppress artifacts (beta, some = 50)

    eng.config.psfun = 0 # Peak shape function (gaussian, lorentzian, split G/L)

    # -- Peak Selection and plotting
    eng.config.peaknorm = 0 # Normalise peaks (0 = off)
    eng.config.datanorm = 0
    eng.config.peakwindow = 10 # peak window / Da
    eng.config.exnorm = 0 # extract normalisation
    eng.config.peakthresh = 0.05
    # eng.config.nativeub = 10
    # eng.config.nativelb = -10
    eng.data.export_hdf5()

    return eng

def update_vars(eng, pmap, skip_empty = False, groupby = 'Time'):
    """Updates vars for each spectra with well ID (starting from top) and stpulated 'groupby' parameter (usually 'Time')"""
    spectra = eng.data.spectra

    # filter pmap 
    if skip_empty == True:
        pmap = pmap[pmap['Type'] != 'empty']
    if len(spectra) == len(pmap):
        
        for i, s in enumerate(spectra):
            
            well_id = pmap.index[i]
            groupbyvar = pmap[groupby].iloc[i]
            s.attrs['Variable 1'] = well_id
            s.var1 = well_id
            s.attrs['Variable 2'] = groupbyvar
            s.var2 = groupbyvar

    else:
        raise Exception("Check plate map and TIC - are lengths the same?")
    return eng

def plot_tic(eng, peak_windows = False, *args, **kwargs):

    plt.figure(*args, **kwargs)
    plt.plot(eng.tic[:, 0], eng.tic[:, 1])

    if peak_windows == True:
        for w in eng.chrompeaks_tranges:
            plt.axvspan(w[0], w[1], alpha = 0.3, color = 'orange') 

    plt.ylabel("Intensity / a.u.")
    plt.xlabel("retention time / mins")

    plt.show() 

class Species:
    def __init__(self, dictionary, name, rmap = None, pmap = None, peak = None, integral = []):
        self.__dict__.update(dictionary)
        # self.__name__ = name
        self.name = name
        self.peak = peak
        self.rmap = rmap
        self.pmap = pmap
        self.integral = []

    def __repr__(self):
        keys = [key for key, val in self.__dict__.items()]
        vals = [val for key, val in self.__dict__.items()]
        return "<"+self.name+ "("+", ".join("{} = {}".format(*t) for t in zip(keys, vals))+")>"        

class Time:
    """Class that contains the data sorted by a stipulated groupby variable (normally time). """
    def __init__(self, time, species, coord = None, name = 'Time', rmap = None, pmap = None):

        self.species = species # list of species in each time point
        self.time = time
        self.coord = coord
        self.__name__= name
        self.spectra = None
        self.thresh = 1
        self.rmap = rmap
        self.pmap = pmap
    
    def __repr__(self):
        species_str = "species = (" + ", ".join(s.__name__ for s in self.species)+")"
        if self.spectra == None:
            return self.__name__+ "(time = {}, coord = {}, ".format(self.time, self.coord) + species_str + ")" 
        else:
            return self.__name__+ "(time = {}, coord = {}, spectra = {}, ".format(self.time, self.coord, self.spectra)+ species_str + ")" 
        
    def extract_masses(self):
        self.theory_masses = np.array([sp.Mass for sp in self.species])
        self.species_name = np.array([sp.name for sp in self.species], dtype = str)
        self.data_masses= np.array([p.mass for p in self.spectra.pks.peaks])
        self.pks = np.array([p for p in self.spectra.pks.peaks])

def update_spectra (reactions_groupby, eng):
    spectra = {s.var1:s for s in eng.data.spectra}    
    for rkey, rval in reactions_groupby.items():
            for t in rval:
                t.spectra = spectra[t.coord]
                t.extract_masses()
    return reactions_groupby

def process(rmap, pmap, eng = None, groupby = 'Time', skip_empty = True):
    reactions_species = {}
    reactions_groupby = {}

    eng = update_vars(eng, pmap, skip_empty)
    

    for rkey, rval in rmap.groupby(['Reaction']):
        reactions_species[rkey] = [Species(spval.to_dict('records')[0], name = spkey) for spkey, spval in rval.groupby('Species')]

    #------------------------------------------------------------------------
    for rkey, rval in pmap.groupby('Reaction'):
        # update each species with time points and coordinates
        for species in reactions_species[rkey]:
            setattr(species, groupby, rval[groupby].values)
            setattr(species, "wells", np.array(rval.index))
        # ----------------------------------------------------
        time = rval[groupby].values # change??
        coords = np.array(rval[groupby].index)
        reactions_groupby[rkey] = []
        for t, c in zip(time, coords):
            # make new species instance for each time point - avoids some memory issues
            species = [deepcopy(s) for s in reactions_species[rkey]]
                
            reactions_groupby[rkey].append(Time(species = species, time = t, coord = c))
    # -----------------------------------------------------------------------
    if eng != None:

        reactions_groupby = update_spectra(reactions_groupby, eng)

    return reactions_groupby

def integrate_all(eng, int_range = None):
    """Creates new spectra attribute and stores areas for each peak there. (eventually put this in UniChrom)"""
    spectra = eng.data.spectra
    
    if int_range == None:
        lb, ub = -eng.config.peakwindow, eng.config.peakwindow
    elif type(int_range) == float:
        lb, ub = int_range[0], int_range[1]
    else:
        lb, ub = -int_range, int_range
        
    for s in spectra:
        peak_ints = []

        for p in s.pks.peaks:
                p.integralrange = [p.mass+lb, p.mass+ub]
                ints = (ud.integrate(s.massdat, p.integralrange[0], p.integralrange[1]))
                p.integral = ints
                peak_ints.append((ints))
        
        s.integrals = peak_ints

    return eng
    
def set_spectra_colors(eng, cmap = 'rainbow'):
    
    cmap = plt.get_cmap(cmap)
    colors = cmap(np.linspace(0, 1, len(eng.data.spectra)))
    for i, s in enumerate(eng.data.spectra):
        s.color = colors[i]
    return eng


def plot_all(eng, show_ints = True, xlim = [], combine = False, cmap = 'Set1'):
    """Plots each spectra stored in Unichrom class"""
    eng = set_spectra_colors(eng, cmap)
    spectra = eng.data.spectra
    xcounter = 0
    ycounter = 0

    if combine == True:
        fig, ax = plt.subplots(dpi = 100)

    for s in spectra:
        if combine == False:
            fig, ax = plt.subplots()
            ax.plot(s.massdat[:, 0], s.massdat[:, 1], color = s.color, linewidth = 0.5)

        if combine == True:        
            ax.plot(s.massdat[:, 0]+xcounter, s.massdat[:, 1]+ycounter, color = s.color, linewidth = 0.5, label = s.var1)
            
        ax.set_xlabel('Mass / Da')
        ax.set_ylabel('Intensity')
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        for i, p in enumerate(s.pks.peaks):
            if combine == False:
                ax.scatter(p.mass, p.height, marker = '^', color = p.color, s=10)
            if xlim != []:
                ax.set_xlim(xlim[0], xlim[1])
                if combine == True:
                    ax.set_xlim(xlim[0], xlim[1]+xcounter)
            if show_ints == True:
                ints = s.integrals[i][1]
                if combine == False:
                    ax.fill_between(ints[:, 0], ints[:, 1], color = p.color, alpha = 0.3)
                else:
                    ax.fill_between(ints[:, 0]+xcounter, ints[:, 1]+ycounter, ycounter, color = p.color, alpha = 0.25)
                    ax.legend()
        xcounter+=s.massdat[:, 0].max()*0.05
        ycounter+=s.massdat[:, 1].max()*0.05

    plt.show()

def match_peaks(eng, reactions):

    for rkey, rval in reactions.items():
        for t in rval:
            tm, dm = np.meshgrid(t.theory_masses, t.data_masses)
            diff = abs(tm - dm)
            diff[diff>10] = np.nan
            for i, d in enumerate(diff):
                if np.isnan(d).all()==False:
                    minimum = np.nanargmin(d)
                    data_peak = t.data_masses[i]
                    t.species[minimum].peak = t.pks[i]
                    t.species[minimum].integral = t.pks[i].integral[0]
    return reactions

def get_data_from_dct(dct, groupby = 'time', data = 'integral'):
    
    time = []
    speciesdct = {}
    speciestimedct = {}
    for t in dct:
        tattr = getattr(t, groupby)
        time.append(tattr)
        for s in t.species:
            
            if s.name in speciesdct:
                speciesdct[s.name].append(getattr(s, data))
                speciestimedct[s.name].append(tattr)
            else:
                speciesdct[s.name] = [getattr(s, data)]
                speciestimedct[s.name] = [t.time]

    return speciesdct, speciestimedct

def plot_data(speciesdct, speciestimedct, combine = True, *args, **kwargs):

    fig, ax = plt.subplots()
    for key, val in speciesdct.items():
        if combine == False:
            fig, ax = plt.subplots()
        ax.plot(speciestimedct[key], val, marker = 'x', label = key, *args, **kwargs)
        ax.set_xlabel('time / s')
        ax.set_ylabel('AUC')
        y_labels = ax.get_yticks()
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
        ax.legend(bbox_to_anchor=(1, 1), loc = 'upper left', frameon = False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    plt.show()
    
# TODO: Set unique colour for each species, choice of cmap 
# TODO: Compare relative AUC's
# TODO: Add functionality to UniChrom API
