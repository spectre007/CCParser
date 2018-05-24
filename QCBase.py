#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import logging
#import numpy as np
"""
Created on Thu Jun 15 04:48:19 2017

@author: alex
"""

class VarNames(object):
    """ General mapping of member variable names to parsing function names.
    After adding a new parsing function, add its result name here. It will
    then be callable as for instance ccp.Parser.results.my_new_var[0] or
    if it was parsed several times as ccp.Parser.results.my_new_var.get_last()"""
    # -- Energies --
    scf_energy = "scf_energy"
    mp_energy = "mp_energy"
    exc_energy_abs = "exc_energy_abs"
    exc_energy_rel = "exc_energy_rel"
    mo_energies = "mo_energies"
    nuc_repulsion = "nuc_rep"
    sapt0_total = "sapt0_total"
    sapt_components = "sapt_comp"
    
    # -- Excited States --
    osc_str = "osc_str"
    amplitudes = "ampl"
    diff_dens_anl = "diff_dens_anl"
    trans_dens_anl = "trans_dens_anl"
    total_dipole = "tot_dip"
    dipole_moment = "dip_moment"
    
    # -- General --
    xyz_coordinates = "xyz"
    pt_order = "pt_order"
    overlap_matrix = "S"
    orthonorm_matrix = "X"
    mo_coefficients = "C"
    alpha_dens_mat = "P_alpha"
    n_bas = "nbas"
    mulliken = "mulliken"
    chelpg = "chelpg"
    
    # -- Issues & Warnings --
    has_converged = "state_cnvgd"
    
    # -- Frozen-Density Embedding
    fde_omega_ref = "fde_omega_ref"
    fde_omega_I = "fde_omega_I"
    fde_trust = "fde_trust"
    fde_delta_lin = "fde_delta_lin"
    fde_electrostatic = "fde_elstat"
    fde_non_electrostatic = "fde_non_elstat"
    fde_timing = "fde_times"
    fde_scf_vemb = "fde_scf_vemb"
    fde_scf_vemb_components = "fde_scf_comp"
    fde_expansion = "fde_expansion"
    fde_method_rhoB = "fde_method_B"
    fde_isA_imported = "fde_import_A"
    fde_isB_imported = "fde_import_B"
    fde_Tfunc = "fde_Tfunc"
    fde_XCfunc = "fde_XCfunc"
    fde_non_elstat_ref = "fde_non_elstat_ref"
    
    
    
    


class QCMethod(object):
    """ Base class for Quantum Chemistry methods """
    def __init__(self):
        self.map = {}
        
    def func_name(self):
        """
        :return: name of caller
        """
        return sys._getframe(1).f_code.co_name
    
    def add_variable(self, func_name, var_name):
        self.map[func_name] = var_name
        
    def print_parsed(self, var_name, description):
        """Print short parsing information
        
        Parameters
        ----------
        var_name : string
            Variable name from the VarNames class
        description : string
            Short (!) description of the object that was parsed.
        """
        print("[results.{0:}] Parsed {1:}.".format(var_name, description))

class GenFormatter(logging.Formatter):
    """ Generic Formatter for logging.Handler
    
    (see https://stackoverflow.com/questions/1343227/
    can-pythons-logging-format-be-modified-depending-on-the-message-log-level)
    """
    default_formatter = logging.Formatter('%(levelname)s in %(name)s: %(message)s')

    def __init__(self, formats):
        """ formats is a dict { loglevel : logformat } """
        self.formatters = {}
        for loglevel in formats:
            self.formatters[loglevel] = logging.Formatter(formats[loglevel])

    def format(self, record):
        formatter = self.formatters.get(record.levelno, self.default_formatter)
        return formatter.format(record)
# TODO: finish this
#class Printer(object):
#    """ Regulates output of the CompChemParser module """
#    out_basename = "CCParser"
#    out_extension = "log"
#    
#    def __init__(self):
#        pass
#    
#    def write(self, string):
#        pass
    
class AtomicBasis(object):
    def __init__(self, xyz, expo, coef, name=""):
        #self.format = None
        self.name = name
        self.center = xyz
        self.exponents = expo
        self.coefficients = coef

#    def __init__(self, name=None, xyz=None, expo=None, coeff=None):
#        if name == None:
#            pass#do something
#        else:
#            pass
        
    @classmethod
    def from_file(cls, atom):
        pass

    def parse_basis(self):
        pass


class BasisSet(object):
    """ General class to deal with basis set issues unrelated to the
    program-specific options """
    def __init__(self):
        self.format = None
        
