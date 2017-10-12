#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
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
    
    # -- Excited States --
    osc_str = "osc_str"
    amplitudes = "ampl"
    diff_dens_anl = "diff_dens_anl"
    trans_dens_anl = "trans_dens_anl"
    
    # -- General --
    xyz_coordinates = "xyz"
    pt_order = "pt_order"
    overlap_matrix = "S"
    orthonorm_matrix = "X"
    mo_coefficients = "C"
    n_bas = "nbas"
    
    # -- Issues & Warnings --
    has_converged = "state_cnvgd"
    
    # -- Frozen-Density Embedding
    fde_trust_first = "fde_trust_init"
    fde_trust = "fde_trust"
    fde_delta_lin = "fde_delta_lin"
    fde_electrostatic = "fde_elstat_int"
    fde_timing = "fde_times"
    fde_scf_vemb = "fde_scf_vemb"
    
    
    
    


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
        
#    def generate_map(self, func_name, var_name):
#        self.var[var_name] = func_name
        
#    def gen_map(self):
#        func_list = [func for func in dir(Foo) if callable(getattr(Foo, func)) and not func.startswith("__")]
#        for func in func_list:
#            self.map[func] = 0
