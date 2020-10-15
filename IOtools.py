#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 16:40:27 2020

@author: nico
"""

import json 
import numpy as np
import os
import re

def load_js(fname):
    """
    Parameters
    ----------
    fname: str
    
    Returns
    -------
    obj
        content of the json file, generally dict
    """
    with open(fname,"r") as f:
        jsdict = json.load(f)
    return jsdict

def load_pointers(d,path=None):
    """
    Parameters
    ----------
    d: dict
        the data dictionary to process. Supposes d["prop"] = [[val1/pointer1, line], [val2, line], ..]
    path: str/None
        the path where the file in pointer is. Default is "."
        
    Does
    ----
    Replaces any pointer with the corresponding values
    """
    npz = {}
    if path is None:
        path = os.getcwd()
    for k,v in d.items():
        if type(v[0][0]) == str and re.match(".+npz", v[0][0]):
            for n,i in enumerate(v):
                if v[n][0] not in npz.keys():
                    npz[v[n][0]] = np.load(os.path.join(path,v[n][0]),allow_pickle=True)
                d[k][n][0] = npz[v[n][0]][k][n]
                
def load_data(fname):
    """
    Parameters
    ----------
    fname: str

    Returns
    -------
    dict
        the dictionary with values instead of pointers
    """
    splt = os.path.split(fname)
    path = splt[0] if splt[0] else os.getcwd()
    jsd = load_js(fname)
    load_pointers(jsd,path=path)
    return jsd