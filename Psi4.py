#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 1 15:35:15 2018

@author: alex
"""
import re
import numpy as np
import logging
#from .ParserData import MolecularOrbitals, Amplitudes
from .QCBase import QCMethod, VarNames as V

# create module logger
mLogger = logging.getLogger("CCParser.Psi4")


class SAPT(QCMethod):
    """Parse quantities related to Symmetry-Adapted Perturbation Theory."""
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"sapt_components" : "SAPT Results",
                      "sapt0_total" : "Total SAPT0"}

    def sapt_components(self, n, data):
        """ Parse SAPT components in [mEh] from SAPT results section.
        
        Returns
        -------
        l : list
            List of SAPT components in form [elst, exch, indc, disp]
        """
        self.add_variable(self.func_name(), V.sapt_components)
        keys = ["Electrostatics", "Exchange", "Induction", "Dispersion"]
        l = [0,0,0,0]
        i = 0
        if self.hooks["sapt_components"] in data[n]:
            while True:#dirrrty
                curr_line = data[n+i+1].split()
                if len(curr_line) > 0:
                    if "Total" in data[n+i+1].split()[0]:
                        break
                    descriptor = data[n+i+1].split()[0]
                    if descriptor in keys:
                        l[keys.index(descriptor)] = data[n+i+1].split()[1]
                i += 1
        mLogger.info("SAPT components in [mEh]",
                     extra={"Parsed":V.sapt_components})
        return l
        
    def sapt0_total(self, n, data):
        """ Parse total interaction energy at SAPT0 level.
        
        Returns
        -------
        sapt0 : float
            Total interaction energy at SAPT0 level in 
        """
        self.add_variable(self.func_name(), V.sapt0_total)
        if self.hooks["sapt0_total"] in data[n]:
            sapt0 = data[n].split()[2]
            mLogger.info("SAPT0 interaction energy in [mEh]",
                         extra={"Parsed":V.sapt0_total})
            return sapt0