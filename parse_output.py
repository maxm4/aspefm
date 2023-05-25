# -*- coding: utf-8 -*-

"""
parse_output.py

Parses clingo output into Pickle and/or JSON

"""

import numpy as np
import pandas, re
from argparse import ArgumentParser
import pickle, sys, warnings

def get_flux_vectors(fname):
    """
    Reads the Clingo[LP] output and gets the flux vector for each solution
    TODO: Not use this function and get answers with clingo API instead
    
    Params: 
        fname: input file name, output from Clingo[LP]

    Returns:
        fluxes: list of flux vectors
        names: list of reaction or enzyme subset names
    """      
    f = open(fname)
    lines = f.readlines()
    fluxes = []
    names = []
    # uses the fact that the min function is 0 and always returns 0
    # after beginning brace index
    bb_idx = len('(0.0, {')
    # this filtering is a bit too specific to be reused...
    lines_flux = filter(lambda w: w.startswith("(0.0,"), lines)
    spl = None
    for l in lines_flux:
        # end brace index
        eb_idx = l.rindex('})')
        # remove beginning and ending braces
        l = l[bb_idx:eb_idx]
        spl = l.split(",")
        ls = list(map(lambda w: float(w.split(":")[1]), spl))
        fluxes.append(np.array(ls))
    if spl is None:
        raise ValueError('unsat' if any('UNSATISFIABLE' in x for x in lines) else 'error') # Clingo LP output format error
    names = spl
    
    def get_word_in_quotes(e):
        res = re.compile(r'flux\(("?.*)"?\)').search(e)
        res = "None" if res is None else res.group(1)
        if res.startswith('"') and res.endswith('"'): res = res[1:-1]
        return res
    
    names = list(map(get_word_in_quotes, names))
    f.close()
    return fluxes, names



def write_flux_vectors(inname, pklfile, jsonfile, csvfile, split_revs, normalize_max, normalize_biom):
    """
    Writes the flux vectors retrieved from a clingo[LP] output
    
    Params: 
        inname: input file name, output from Clingo[LP]
        pklfile: file to store flux modes in pickle format
        jsonfile: file to store flux modes in json format
        csvfile: file to store flux modes in csv format
        split_revs: do not merge reversible reaction fluxes
        normalize_max: if toggled, normalize output flux vectors by their maximum flux value
        normalize_biom: biomass reaction name, if exists, normalize output flux vectors by their biomass flux value
    """      
    fluxes, names = get_flux_vectors(inname)
    fvs = np.vstack(fluxes)
    pfvs = pandas.DataFrame(data=fvs, columns=names)
    pfvs = pfvs.reindex(sorted(pfvs.columns), axis=1)

    if normalize_biom is not None:
        pfvs = pfvs.apply(lambda x: x/x[normalize_biom], axis=1) # note: NaN producing if flux(biom) == 0
    elif normalize_max:
        pfvs = pfvs.apply(lambda x: x/x.max(), axis=1)

    if not split_revs:
        names = pfvs.columns[~pfvs.columns.str.endswith('rev')] 
        revs = pfvs.columns[pfvs.columns.str.endswith('rev')]
        drev = pfvs.copy()[revs] # type: pandas.DataFrame
        drev = -drev # type: ignore
        drev.columns = drev.columns.map(lambda x: x[:-4])  
        pfvs = pfvs[names]
        pfvs[drev.columns] += drev 
            
    if pklfile is not None:
        with open(pklfile, 'wb') as pkl:
            pickle.dump(pfvs, pkl)

    if jsonfile is not None:
        with open(jsonfile, 'w') as js:
            pfvs.to_json(path_or_buf=js, orient='records')

    if csvfile is not None:
        with open(csvfile, 'w') as csv:
            pfvs.to_csv(path_or_buf=csv)
            
    if pklfile is None and jsonfile is None and csvfile is None:
        warnings.warn('Please specify a format to convert to (pickle, json, csv) along with a filename')

    
  
if __name__== "__main__":
    parser = ArgumentParser()
    parser.add_argument('infile', metavar='clingoLP.file', help='Input file name')
    parser.add_argument('--pickle', default=None, metavar='Pickle.file', help='Store flux modes with pickle')
    parser.add_argument('--json', default=None, metavar='JSON.file', help='Store flux modes in JSON format')
    parser.add_argument('--csv', default=None, metavar='CSV.file', help='Store flux modes in CSV format')
    parser.add_argument('--normalize-max', action='store_true', help='Normalize reaction fluxes by their maximum flux (toggle)')
    parser.add_argument('--normalize-biom', default=None, help='Normalize reaction fluxes by their biomass flux (specify name)')  
    parser.add_argument('--split-revs', action='store_true', help='Merges reversible reaction fluxes')
    opts = parser.parse_args()
    write_flux_vectors(opts.infile, opts.pickle, opts.json, opts.csv, opts.split_revs, opts.normalize_max, opts.normalize_biom)
