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
from .QCBase import var_tag

# create module logger
mLogger = logging.getLogger("CCParser.Psi4")

class General(QCMethod):
    """Parse quantities related to Symmetry-Adapted Perturbation Theory."""
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        # [!!!!] just temporary fix since I can't make out a has_finished hook!
        self.hooks = {"has_finished" : "Warning! sapt0 does not have an associated derived wavefunction."}

    @var_tag(V.has_finished)
    def has_finished(self, i, data):
        """ Parse final statement that indicates if Gaussian finished
        without errors. """
        mLogger.info("whether Psi4 has finished successfully",
                     extra={"Parsed": V.has_finished})
        return True

class SAPT(QCMethod):
    """Parse quantities related to Symmetry-Adapted Perturbation Theory."""
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"sapt_components" : "SAPT Results",
                      "sapt0_total" : "Total SAPT0"}

    @var_tag(V.sapt_components)
    def sapt_components(self, i, data):
        """ Parse SAPT components in [mEh] from SAPT results section.
        
        Returns
        -------
        l : list
            List of SAPT components in form [elst, exch, indc, disp]
        """
        keys = ["Electrostatics", "Exchange", "Induction", "Dispersion"]
        l = [0,0,0,0]
        n = 0
        if self.hooks["sapt_components"] in data[n]:
            while True:#dirrrty
                curr_line = data[i+n+1].split()
                if len(curr_line) > 0:
                    if "Total" in data[i+n+1].split()[0]:
                        break
                    descriptor = data[i+n+1].split()[0]
                    if descriptor in keys:
                        l[keys.index(descriptor)] = float(data[i+n+1].split()[1])
                n += 1
        mLogger.info("SAPT components in [mEh]",
                     extra={"Parsed":V.sapt_components})
        return l

    @var_tag(V.sapt0_total)
    def sapt0_total(self, n, data):
        """ Parse total interaction energy at SAPT0 level.
        
        Returns
        -------
        sapt0 : float
            Total interaction energy at SAPT0 level in 
        """
        sapt0 = float(data[n].split()[2])
        mLogger.info("SAPT0 interaction energy in [mEh]",
                     extra={"Parsed":V.sapt0_total})
        return sapt0
