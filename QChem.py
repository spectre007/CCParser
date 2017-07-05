#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 23:06:40 2017

@author: alex
"""
import re
#import numpy as np
#import pandas as pd
#import xarray as xr
from .ParserData import MolecularOrbitals
from .QCBase import QCMethod, VarNames as V


class SCF(QCMethod):
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"scf_energy" : "SCF   energy in the final basis set",
                      "orb_energies": "Orbital Energies (a.u.)"}
    
    def scf_energy(self, l_index, data):
        """ Get Hartree-Fock energy [a.u.] from scfman """
        self.add_variable(self.func_name(), V.scf_energy)
        if self.hooks["scf_energy"] in data[l_index]:
            return float(data[l_index].split()[-1])
    
    def orb_energies(self, l_index, data):
        """ Parse Hartree-Fock molecular orbital energies """
        #self.map["mos"] = "orb_energies"
        self.add_variable(self.func_name(), V.mo_energies)
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
            #alpha.RVS(50)
            return alpha

class ADC(QCMethod):
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        #self.type = ["Excited States","Perturbation Theory"]
        self.hooks = {"scf_energy": "HF Summary",
                      "mp_energy": r"(RI-)?MP\([2-3]\) Summary",
                      "exc_energies": "Excitation energy:",
                      "osc_strength": "Osc. strength:",
                      "amplitudes": "Important amplitudes:"}

    def exc_energies(self, l_index, data):
        """ Parse excitation energies [eV] """
        self.add_variable(self.func_name(), V.exc_energy_rel)
        if self.hooks["exc_energies"] in data[l_index]:
            return float(data[l_index].split()[-2])
        
    def osc_strength(self, l_index, data):
        """ Parse oscillator strengths """
        self.add_variable(self.func_name(), V.osc_str)
        if self.hooks["osc_strength"] in data[l_index]:
            return float(data[l_index].split()[-1])
        
    def scf_energy(self, l_index, data):
        """ Parse SCF energy [a.u.] from adcman """
        self.add_variable(self.func_name(), V.scf_energy)
        if self.hooks["scf_energy"] in data[l_index]:
            return float(data[l_index+2].split()[-2])
        
    def mp_energy(self, l_index, data):
        """ Parse MP() reference state energy [a.u.] from adcman """
        self.add_variable(self.func_name(), V.mp_energy)
        match = re.search(self.hooks["mp_energy"], data[l_index])
        if match:
            return float(data[l_index+3].split()[2])
    
    def amplitudes(self, l_index, data):
        """ Parse occ -> virt amplitudes """
        self.add_variable(self.func_name(), V.amplitudes)
        if self.hooks["amplitudes"] in data[l_index]:
            try:
                # regex has the awesome captures feature
                import regex
                expr = r"(?:(?P<orb>\d+) [^0-9]+)+ (?P<ampl>-?\d+\.\d+)"
                srch = regex.search
                have_regex = True
            except ImportError:
                 # explicit expression for orbital transitions
                 expr = r"(\d+) [^0-9]+ ((\d+) [^0-9]+)? (\d+) [^0-9]+ ((\d+) [^0-9]+)? (-?\d+\.\d+)"
                 srch = re.search
                 have_regex = False
        
            idx = l_index+3
            amplist = []
            while "--------" not in data[idx]:
                match = srch(expr, data[idx])
                if match:
                    if have_regex:
                        orb = list(map(int, match.captures("orb")))
                        amp = list(map(float, match.captures("ampl")))
                        amplist.append(orb + amp)
                    else:
                        amplist.append([int(match.group(x)) for x in [1,3,4,6] if (match.group(x) != None)]+[float(match.group(7))])
                idx += 1
        return amplist
                
        