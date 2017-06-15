#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
"""
Created on Thu Jun 15 04:48:19 2017

@author: alex
"""

class VarNames(object):
    # -- Energies --
    scf_energy = "scf_energy"
    mp_energy = "mp_energy"
    exc_energy_abs = "exc_energy_abs"
    exc_energy_rel = "exc_energy_rel"
    mo_energies = "mo_energies"
    
    # -- Transition properties --
    osc_str = "osc_str"
    
    # -- Misc --
    pt_order = "pt_order"
    
    
    


class QCMethod(object):
    """ Base class for Quantum Chemistry methods """
    def __init__(self):
        self.map = {}
        self.var = {
                "scf_energy" : 0,
                "exc_energy" : 0,
                "osc_str" : 0,
                "mos" : 0
                }
        
    def func_name(self):
        """
        :return: name of caller
        """
        return sys._getframe(1).f_code.co_name
    
    def add_variable(self, func_name, var_name):
        self.map[func_name] = var_name
        
#    def generate_map(self, func_name, var_name):
#        self.var[var_name] = func_name
        
#    def gen_map(self):
#        func_list = [func for func in dir(Foo) if callable(getattr(Foo, func)) and not func.startswith("__")]
#        for func in func_list:
#            self.map[func] = 0