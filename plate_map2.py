import plate_map as pm 
import pandas as pd


rheaders = {"Reaction Name":str, "Species":str, "Concentration":float,
           "Units":str, "Mass":float, "Reagent Type":str, "Sequence":str}

pheaders = {"Well ID":{'dtype':str, 'long':True, 'short_row': False, 'short_col':False}, 
            "Type":{'dtype':str, 'long':True, 'short_row': True, 'short_col':True}, 
            "Reaction Name":{'dtype':str, 'long':True, 'short_row': True, 'short_col':True}, 
            "Time":{'dtype':str, 'long':True, 'short_row': True, 'short_col':True}}

map_name = "plate map uni chrom update.xlsx"


def read_in(map_path, map_type = 'long'):

    if map_type == 'long':
        rmap = pd.read_excel(map_path, 
                            dtype = rheaders, index_col = "Species",
                            sheet_name = "reaction map")

        pdata_types = {i[0]: i[1]['dtype'] for i in pheaders.items()} 
        pmap = pd.read_excel(map_path, dtype = pdata_types, index_col = "Well ID", skiprows = 1)

        # add valid column
        pheaders2 = [x for x in pheaders.keys() if pheaders[x]['long']]
        pmapdf = pm.empty_map(size = 6, header_names = pheaders)
        pmapdf.update(pmap)

    return rmap, pmapdf

def process_data(rmap, pmapdf, spectra = None, groupby = 'Time', maptogroupby = 'plate_map'):
    """Collects metadata and data for each reaction on plate map into single dictionary of dictionaries"""
    # filter by reactions 
    reactions = {key:val for key, val in rmap.groupby('Reaction Name')}
    wells = {key: val for key, val in pmapdf.groupby('Reaction Name')} 
    

    # TODO: add check that reaction names (.unique()) appear in both reagent/reaction map and plate map

    if maptogroupby == "plate_map": # NOTE: Duplicate rows (i.e. multiple wells for same time point) will be merged into single container
        
        platemapdict = {}
        for key, val in wells.items():
            grouped = val.groupby(groupby) 
            platemapdict[key] = {}
            # unpack each time point
            platemapdict[key][groupby] = {float(k):[v] for k, v in grouped} # where v is the dataframe row containing the metadata from the platemap

            # if spectra != None:
            platemapdict[key]['spectra'] = "insert spectra here"

    for key, val in reactions.items(): 
        platemapdict[key]['metadata'] = val
    return platemapdict


