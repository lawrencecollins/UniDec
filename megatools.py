from unidec_modules.ChromEng import *
from unidec_modules import unidectools as ud



def get_fwhm(spectra):

    spectradict = {}

    for s in spectra.data.spectra:
        height = []
        mass = []
        spectradict[s.name] = []
        for p in s.pks.peaks:
            height.append(p.height)
            mass.append(p.mass)

            FWHM = ud.calc_FWHM(p.mass, s.massdat) 

            spectradict[s.name].append([FWHM, p.mass, p.height])
            
    return spectradict
    