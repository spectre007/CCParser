#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 13:46:56 2018

@author: alex
"""

import re
import numpy as np
import logging
from .ParserData import MolecularOrbitals, Amplitudes
from .QCBase import QCMethod, VarNames as V
from .QCBase import var_tag
from .ParserTools import is_float

# create module logger
mLogger = logging.getLogger("CCParser.Gaussian")

def parse_tddft_list(i, data, columns=[0], asmatrix=False):
    """General function to parse tables from the TDDFT output.
    
    Parameters
    ----------
    i : int
        Line number of hook.
    data : list
        Readlines container.
    columns : list
        Indices of columns to parse (python counting).
    asmatrix : bool
        Whether or not to return a numpy.matrix object from parsed data.
    
    Returns
    -------
    the_list : array_like
        List or np.matrix containing the requested values.
    """
    index_line = i+2
    the_list = []
    while True:
        line = data[index_line].split()
        ncol = len(line)
        if ncol == 0:
            break
        elif not is_float(line[0]):
            break
        the_list.append(list(map(float, [line[i] for i in columns])))
        index_line += 1
        
    if asmatrix:
        return np.asmatrix(the_list)
    else:
        return the_list

class General(QCMethod):
    """Quantities that are not related to any method."""
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"has_finished" : "Normal termination of Gaussian"}

    @var_tag(V.has_finished)
    def has_finished(self, i, data):
        """ Parse final statement that indicates if Gaussian finished
        without errors. """
        mLogger.info("whether Gaussian has finished successfully",
                     extra={"Parsed": V.has_finished})
        return True

class TDDFT(QCMethod):
    """Parse TDDFT output."""
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"transition_dipole":("Ground to excited state transition"
                          " electric dipole moments"),
                      "exc_energy_rel": (r"Excited State\s+\d+:\s+(?P<Label>"
                          r"[A-Za-z0-9_?-]+)\s+(?P<ExcEV>\d+\.\d+) eV\s+"
                          r"(?P<ExcNM>\d+\.\d+) nm\s+f=(?P<Osc>\d+\.\d+)\s+"
                          r"<S\*\*2>=(?P<S2>\d+\.\d+)\s*"),
                      "state_label": (r"Excited State\s+\d+:\s+(?P<Label>"
                          r"[A-Za-z0-9_?-]+)\s+(?P<ExcEV>\d+\.\d+) eV\s+"
                          r"(?P<ExcNM>\d+\.\d+) nm\s+f=(?P<Osc>\d+\.\d+)\s+"
                          r"<S\*\*2>=(?P<S2>\d+\.\d+)\s*"),
                      "oscillator_strength": (r"Excited State\s+\d+:\s+(?P<Label>"
                          r"[A-Za-z0-9_?-]+)\s+(?P<ExcEV>\d+\.\d+) eV\s+"
                          r"(?P<ExcNM>\d+\.\d+) nm\s+f=(?P<Osc>\d+\.\d+)\s+"
                          r"<S\*\*2>=(?P<S2>\d+\.\d+)\s*"),
                      "amplitudes": r"Excited State\s+\d+:\s+"}

    @var_tag(V.transition_dipole)
    def transition_dipole(self, i, data):
        """Parse transition dipole moment in [a.u.]."""
        mLogger.info("transition dipole moment",
                     extra={"Parsed": V.transition_dipole})
        return parse_tddft_list(i, data, columns=[1,2,3])

    @var_tag(V.exc_energy_rel)
    def exc_energy_rel(self, i , data):
        """Parse excitation energy in [eV]."""
        match = re.search(self.hooks[self.func_name()], data[i])
        mLogger.info("relative TDDFT excitation energy/-ies [eV]",
                     extra={"Parsed": V.exc_energy_rel})
        return float(match.group("ExcEV"))

    @var_tag(V.state_label)
    def state_label(self, i , data):
        """Parse state label (multiplicity and symmetry group)."""
        match = re.search(self.hooks[self.func_name()], data[i])
        mLogger.info("state multiplicity and symmetry",
                     extra={"Parsed": V.state_label})
        return match.group("Label")

    @var_tag(V.osc_str)
    def oscillator_strength(self, i , data):
        """Parse oscillator strength."""
        match = re.search(self.hooks[self.func_name()], data[i])
        mLogger.info("oscillator strength",
                     extra={"Parsed": V.osc_str})
        return float(match.group("Osc"))

    @var_tag(V.amplitudes)
    def amplitudes(self, i, data):
        """ Parse occ -> virt amplitudes """
        pattern = r"\s+(?P<occ>\d+) -> (?P<virt>\d+)\s+(?P<v>-?\d+\.\d+)\s*"
        j = 0 # line counter
        amplist = []
        while True:
            m = re.search(pattern, data[i+1+j])
            if m:
                amplist.append([int(m.group("occ")), int(m.group("virt")),
                                float(m.group("v"))])
                j += 1
            else:
                break
        mLogger.info("TDDFT amplitudes", extra={"Parsed":V.amplitudes})
        return Amplitudes.from_list(amplist, factor=2.0)

class Freq(QCMethod):
    """Parse frequency output."""
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"vibrational_freq" : "and normal coordinates:",
                      "infrared_intensity" : "and normal coordinates:"}

    @var_tag(V.vib_freq)
    def vibrational_freq(self, i , data):
        """Parse vibrational frequencies in [cm-1]."""
        n = 0
        freqs = []
        while "------" not in data[i+n]:
            if "Frequencies --" in data[i+n]:
                freqs += data[i+n].split()[2:]
            n += 1
        freqs = list(map(float, freqs))
        mLogger.info("vibrational frequencies in [cm-1]",
                     extra={"Parsed":V.vib_freq})
        return freqs

    @var_tag(V.vib_intensity)
    def infrared_intensity(self, i, data):
        """Parse IR intensity in [km/mol]."""
        n = 0
        intensity = []
        while "------" not in data[i+n]:
            if "IR Inten    --" in data[i+n]:
                intensity += data[i+n].split()[3:]
            n += 1
        intensity = list(map(float, intensity))
        mLogger.info("IR intensities in [km/mol]",
                     extra={"Parsed":V.vib_intensity})
        return intensity
