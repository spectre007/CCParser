#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 23:06:40 2017

@author: alex
"""
import re
import numpy as np
import pandas as pd
import xarray as xr
from .ParserData import MolecularOrbitals

#class QCMethod(object):
#    pass

class SCF(object):
    def __init__(self):
        self.hooks = {"orb_energies": "Orbital Energies (a.u.)"}
    
    def orb_energies(self, l_index, data):
        """ Parse Hartree-Fock molecular orbital energies """
        if self.hooks["orb_energies"] in data[l_index]:
            n = 0
            while True:
                if "------------" in data[l_index+2+n]:
                    n += 1
                    break
                else:
                    n += 1
            s = "".join(data[l_index+2:l_index+2+n])
#            occpatt = r"-- Occupied --\s+(-?\d+\.\d+\s*){1,}\s+"
            occp = r"-- Occupied --([^A-Z]*)-- Virtual --"
            virt = r"-- Virtual --([^A-Z]*)[-]{4,}"
#            rem = re.search(occpatt,s,re.MULTILINE)
            rem2 = re.search(occp,s,re.MULTILINE)
            rem3 = re.search(virt,s,re.MULTILINE)
            if rem2:  
                a_occ = re.findall("-?\d+\.\d+",rem2.group(0),re.M)
            if rem3:
                a_virt = re.findall("-?\d+\.\d+",rem3.group(0),re.M)
            alpha = MolecularOrbitals(a_occ,a_virt)
            alpha.RVS(50)
            return alpha

class ADC(object):
    def __init__(self):
        #self.type = ["Excited States","Perturbation Theory"]
        self.hooks = {"scf_energy": "HF Summary",
                      "mp_energy": r"(RI-)?MP\([2-3]\)",
                      "exc_energies": "Excitation energy:",
                      "osc_strength": "Osc. strength:"}
        # have rawData as attribute??
    def exc_energies(self, l_index, data):
        """ Parse excitation energies in eV """
        if self.hooks["exc_energies"] in data[l_index]:
            return float(data[l_index].split()[-2])
        
    def osc_strength(self, l_index, data):
        """ Parse oscillator strengths """
        if self.hooks["osc_strength"] in data[l_index]:
            return float(data[l_index].split()[-1])
        
    def scf_energy(self, l_index, data):
        """ Parse SCF energy from adcman """
        if self.hooks["scf_energy"] in data[l_index]:
            return float(data[l_index+2].split()[-2])
        
    def mp_energy(self, l_index, data):
        """ Parse MP() reference state energy from adcman """
        match = re.search(self.hooks["mp_energy"], data[l_index])
        if match:
            return match.group(0)
        