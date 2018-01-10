#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 23:06:40 2017

@author: alex
"""
import re
import numpy as np
#import pandas as pd
#import xarray as xr
from .ParserData import MolecularOrbitals, Amplitudes
from .QCBase import QCMethod, VarNames as V

def parse_symmetric_matrix(readlin, n, asmatrix=True):
    """Parse a symmetric matrix printed columnwise 
    
    Parameters
    ----------
    readlin : list
        Readlines list object
    n : int
        Line number of identifier
    asmatrix : bool
        Whether to return a numpy.matrix object or not
    
    Returns
    -------
    numpy.matrix
        Parsed AO matrix as numpy.matrix object if asmatrix=True
    list
        Parsed AO matrix as list of lists if asmatrix=False
    """
    matrix = []
    cols = 0
    index_line = n+1
    first_batch = True
    stop_signals = ["Gap", "=", "eV", "Convergence", "criterion"]
    while not any(stop in readlin[index_line].split() for stop in stop_signals):#loop over blocks
        ncol = len(readlin[index_line].split())
        # adding rows scheme -> take line split as is
        j = 0
        if cols > 0:
            first_batch = False
        while True: #loop over lines in one block
            if len(readlin[index_line+j+1].split()) != ncol+1 \
            or not readlin[index_line+j+1].split()[0].isdigit() \
            or any(stop in readlin[index_line+j+1].split() for stop in stop_signals):
                break
            if first_batch:
                matrix.append([])
            matrix[j] += list(map(float, readlin[index_line+j+1].split()[1:]))
            j += 1
        index_line += j+1
        cols += ncol
    if asmatrix: # return np.matrix object
        return np.asmatrix(matrix)
    else: # return list of lists
        return matrix

def parse_AO_matrix(readlin, n, asmatrix=True):
    """Parse a matrix that uses AO descriptors, e.g. MO coefficients 
    
    The matrix is printed columnwise (here hardcoded to max. 6 columns) and has
    AO descriptors on the left side, for instance "25 C 1   pz". The head of
    each column contains the column number (Fortran counting) and the
    eigenvalue, both of which are not parsed.
    
    Parameters
    ----------
    readlin : list
        Readlines list object
    n : int
        Line number of identifier
    asmatrix : bool
        Whether to return a numpy.matrix object or not
    
    Returns
    -------
    numpy.matrix
        Parsed AO matrix as numpy.matrix object if asmatrix=True
    list
        Parsed AO matrix as list of lists if asmatrix=False
    """
    index_line = n+1
    cols = 0
    pattern = r"(\d{1,4}?)\s+(?P<atom>\w{1,3}?\s?(\d+)?)\s+(?!wall)(\w+\s?\d?)\s+(?P<coeffs>(-?\d+\.\d+\s+){1,6})"#hardcoded up to 6 columns
    matrix = []
    while True: # loop over block
        ncol = len(readlin[index_line].split())
        try:
            readlin[index_line].split()[0].isdigit()
        except IndexError:
            break
        if ncol == 0:
            break
        matrix += [[] for _dummy in range(ncol)]
        nbas = 0
        while True: # loop over lines
            match = re.search(pattern, readlin[index_line+nbas+2])
            if match:
                nbas += 1
                tmp = match.group("coeffs").split()
                for k in range(ncol):
                    matrix[cols+k].append(float(tmp[k]))
            else:
                break
        index_line += nbas+2
        cols += ncol
    if asmatrix:
        return np.asmatrix(matrix)
    else:
        return matrix


class General(QCMethod):
    """Parse general information like basis set, number of atoms, etc. """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"xyz_coordinates" : "Standard Nuclear Orientation (Angstroms)",
                      "electrons" : r"There are \w+(?P<alpha>\d+) alpha and \w+(?P<beta>\d+) beta electrons",
                      "nuc_rep" : "Nuclear Repulsion Energy =",
                      "n_basis" : r"There are (?P<shells>\d+) shells and (?P<bsf>\d+) basis functions",
                      "mulliken": "Ground-State Mulliken Net Atomic Charges",
                      "chelpg": "Ground-State ChElPG Net Atomic Charges"}
    
    def xyz_coordinates(self, l_index, data):
        """ Get XYZ coordinates """
        self.add_variable(self.func_name(), V.xyz_coordinates)
        if self.hooks[self.func_name()] in data[l_index]:
            xyz_dat, n = [], 0
            while True:
                if "------------" in data[l_index+3+n]:
                    break
                else:
                    xyz_dat.append(data[l_index+3+n].split()[1:])
                    n += 1
            xyz_dat = [[x[0]]+list(map(float,x[1:])) for x in xyz_dat]
            self.print_parsed(V.xyz_coordinates, "XYZ coordinates")
            return xyz_dat

    def nuc_rep(self, l_index, data):
        """ Get nuclear repulsion energy [a.u.] """
        self.add_variable(self.func_name(), V.nuc_repulsion)
        if self.hooks[self.func_name()] in data[l_index]:
            self.print_parsed(V.nuc_repulsion, "nuclear repulsion energy")
            return float(data[l_index].split()[-2])
    
    def n_basis(self, l_index, data):
        """ Parse number of basis functions """
        self.add_variable(self.func_name(), V.n_bas)
        match = re.search(self.hooks["n_basis"], data[l_index])
        if match:
            self.print_parsed(V.n_bas, "number of basis functions")
            return int(match.group("bsf"))
    
    def mulliken(self, l_index, data):
        """ Parse ground-state mulliken charges """
        self.add_variable(self.func_name(), V.mulliken)
        if self.hooks[self.func_name()] in data[l_index]:
            chg, n = [], 0
            while True:
                if "-------" in data[l_index+4+n]:
                    break
                else:
                    chg.append(data[l_index+4+n].split()[1:])
                    n += 1
            chg = [[x[0]]+[float(x[1])] for x in chg]
            self.print_parsed(V.mulliken, "Ground-State Mulliken charges")
            return chg
        
    def chelpg(self, l_index, data):
        """ Parse ground-state ChElPG charges """
        self.add_variable(self.func_name(), V.chelpg)
        if self.hooks[self.func_name()] in data[l_index]:
            chg, n = [], 0
            while True:
                if "-------" in data[l_index+4+n]:
                    break
                else:
                    chg.append(data[l_index+4+n].split()[1:])
                    n += 1
            chg = [[x[0]]+[float(x[1])] for x in chg]
            self.print_parsed(V.chelpg, "Ground-State ChElPG charges")
            return chg
            

class SCF(QCMethod):
    """ Parse scfman related quantities """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"scf_energy" : "SCF   energy in the final basis set",
                      "mo_energies": "Orbital Energies (a.u.)",
                      "overlap_matrix" : " Overlap Matrix",
                      "orthonorm_matrix" : " Orthonormalization Matrix",
                      "alpha_density_matrix" : " Alpha Density Matrix",
                      "mo_coefficients_r" : "RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS"}
    
    def scf_energy(self, l_index, data):
        """ Get Hartree-Fock energy [a.u.] from scfman """
        self.add_variable(self.func_name(), V.scf_energy)
        if self.hooks["scf_energy"] in data[l_index]:
            self.print_parsed(V.scf_energy, "SCF energy")
            return float(data[l_index].split()[-1])
    
    def mo_energies(self, l_index, data):
        """ Parse Hartree-Fock molecular orbital energies """
        #self.map["mos"] = "orb_energies"
        self.add_variable(self.func_name(), V.mo_energies)
        if self.hooks["mo_energies"] in data[l_index]:
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
            self.print_parsed(V.mo_energies, "molecular orbital energies")
            return alpha
        
    def overlap_matrix(self, l_index, data):
        """ Parse overlap matrix S """
        self.add_variable(self.func_name(), V.overlap_matrix)
        self.print_parsed(V.overlap_matrix, "overlap matrix")
        return parse_symmetric_matrix(data, l_index)
    
    def orthonorm_matrix(self, l_index, data):
        """ Parse orthonormalization matrix X """
        self.add_variable(self.func_name(), V.orthonorm_matrix)
        self.print_parsed(V.orthonorm_matrix, "orthonormalization matrix")
        return parse_symmetric_matrix(data, l_index)
    
    def alpha_density_matrix(self, l_index, data):
        """ Parse alpha density matrix P_alpha """
        self.add_variable(self.func_name(), V.alpha_dens_mat)
        self.print_parsed(V.alpha_dens_mat, "SCF alpha density matrix")
        return parse_symmetric_matrix(data, l_index)
    
    def mo_coefficients_r(self, l_index, data):
        """ Parse MO coefficients C for restricted SCF"""
        self.add_variable(self.func_name(), V.mo_coefficients)
        C = parse_AO_matrix(data, l_index)
        self.print_parsed(V.mo_coefficients, "{0:} molecular orbital coefficients".format(C.shape[0]))
        return C
        

class ADC(QCMethod):
    """ Parse adcman related quantities """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        #self.type = ["Excited States","Perturbation Theory"]
        self.hooks = {"scf_energy": "HF Summary",
                      "mp_energy": r"(RI-)?MP\([2-3]\) Summary",
                      "exc_energies": "Excitation energy:",
                      "osc_strength": "Osc. strength:",
                      "has_converged" : r"Excited state \d+ \(.*?\)\s+\[(.*?)\]",
                      "amplitudes": "Important amplitudes:",
                      "total_dipole" : "Total dipole [Debye]",
                      "diff_dens_anl": "Exciton analysis of the difference density matrix",
                      "mulliken_adc": "Mulliken Population Analysis"}

    def exc_energies(self, l_index, data):
        """ Parse excitation energies [eV] """
        self.add_variable(self.func_name(), V.exc_energy_rel)
        if self.hooks["exc_energies"] in data[l_index]:
            self.print_parsed(V.exc_energy_rel, "relative ADC(x) excitation energy/-ies in [eV]")
            return float(data[l_index].split()[-2])
        
    def osc_strength(self, l_index, data):
        """ Parse oscillator strengths """
        self.add_variable(self.func_name(), V.osc_str)
        if self.hooks["osc_strength"] in data[l_index]:
            self.print_parsed(V.osc_str, "ADC oscillator strength/s")
            return float(data[l_index].split()[-1])
        
    def scf_energy(self, l_index, data):
        """ Parse SCF energy [a.u.] from adcman """
        self.add_variable(self.func_name(), V.scf_energy)
        if self.hooks["scf_energy"] in data[l_index]:
            self.print_parsed(V.scf_energy, "SCF energy in [a.u.] (adcman)")
            return float(data[l_index+2].split()[-2])
        
    def mp_energy(self, l_index, data):
        """ Parse MP() reference state energy [a.u.] from adcman """
        self.add_variable(self.func_name(), V.mp_energy)
        match = re.search(self.hooks["mp_energy"], data[l_index])
        if match:
            self.print_parsed(V.mp_energy, "MP(x) energy in [a.u.] (adcman)")
            return float(data[l_index+3].split()[2])
    def has_converged(self, l_index, data):
        """ Parse if state has converged. """
        self.add_variable(self.func_name(), V.has_converged)
        match = re.search(self.hooks["has_converged"], data[l_index])
        if match:
            self.print_parsed(V.has_converged, "if states converged (adcman)")
            return True if match.group(1) == "converged" else False

    def amplitudes(self, l_index, data):
        """ Parse occ -> virt amplitudes """
        self.add_variable(self.func_name(), V.amplitudes)
        if self.hooks["amplitudes"] in data[l_index]:
            try:
                # regex has the awesome captures feature
                import regex
                expr = r"(?:(?P<orb>\d+) [^0-9]+)+ (?P<ampl>-?\d+\.\d+)" # FIXME: [^0-9] might be a problem if symmetry is on
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
        self.print_parsed(V.amplitudes, "ADC amplitudes")
        return Amplitudes.from_list(amplist, factor=2.0)

    def total_dipole(self, l_index, data):
        """ Parse total dipole moment [Debye] for HF, MP2 and ADC states
        
        Returns
        -------
        float
            The total dipole moment [Debye]"""
        self.add_variable(self.func_name(), V.total_dipole)
        if self.hooks["total_dipole"] in data[l_index]: 
            self.print_parsed(V.total_dipole, "Total dipole")
            return float(data[l_index].split()[-1])
    
    def diff_dens_anl(self, l_index, data):
        """ Parse difference density matrix analysis block """
        self.add_variable(self.func_name(), V.diff_dens_anl)
        return 0
#        if self.hooks[S"diff_dens_anl"] in data[l_index]:
#            d = {}
#            p = [r"\[(-?\d+\.\d+), (-?\d+\.\d+), (-?\d+\.\d+)\]",
#                 r"(-?\d+\.\d+)"]
#            pattern = "|".join(p)
#            i = 1
#            s = ""
#            while True:
#                if "Transition density matrix analysis:" in data[l_index+i]:
#                    print("diff_dens_anl: combined ",i," lines.")
#                    break
#                else:
#                  s += data[l_index+i]
#                  i += 1
#            
#            # check if f has list elements
##            d["hole_pos"] = f[0]
##            d["elec_pos"] = f[1]
##            d["eh_dist"] = f[2]
##            d["hole_size"] = f[3]
##            d["hole_comp"] = f[4]
##            d["elec_size"] = f[5]
##            d["elec_comp"] = f[6]
#            return f
    def mulliken_adc(self, l_index, data):
        """ Parse MP(x) and ADC(x) mulliken charges """
        self.add_variable(self.func_name(), V.mulliken)
        if self.hooks[self.func_name()] in data[l_index]:
            chg, n = [], 0
            while True:
                if "-------" in data[l_index+3+n]:
                    break
                else:
                    chg.append(data[l_index+3+n].split()[1:])
                    n += 1
            chg = [[x[0]]+[float(x[1])] for x in chg]
            self.print_parsed(V.mulliken, "MP(x)/ADC(x) Mulliken charges")
            return chg
                
    
class FDE_ADC(QCMethod):
    """ Parsing related to FDE-ADC implementation in Q-Chem """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"fde_trust_first": "FDE control parameter",
                      "fde_electrostatic": "rho_A <-> rho_B",
                      "fde_trust": "lambda(FDE)",
                      "fde_delta_lin": "Delta_Lin:",
                      "fde_timing": "FDE timings",
                      "fde_scf_vemb": "Integrated total embedding potential"}

        
    def fde_trust_first(self, l_index, data):
        """ Parse FDE trust parameter (before construction of embedding potential) [ppm] """
        self.add_variable(self.func_name(), V.fde_trust_first)
        if self.hooks["fde_trust_first"] in data[l_index]:
            self.print_parsed(V.fde_trust_first, "a priori FDE overlap parameter Lambda")
            return float(data[l_index].split()[5])
    
    def fde_trust(self, l_index, data):
        """ Parse FDE trust parameter [ppm] from final FDE output """
        self.add_variable(self.func_name(), V.fde_trust)
        if self.hooks["fde_trust"] in data[l_index]:
            self.print_parsed(V.fde_trust, "state specific FDE overlap parameter Lambda")
            return float(data[l_index].split()[1])
        
    def fde_delta_lin(self, l_index, data):
        """ Parse 1st order term [eV] of LinFDET approximation """
        self.add_variable(self.func_name(), V.fde_delta_lin)
        if self.hooks["fde_delta_lin"] in data[l_index]:
            self.print_parsed(V.fde_delta_lin, "1st order term of linearized FDET")
            return float(data[l_index].split()[1])
    
    
    def fde_electrostatic(self, l_index, data):
        """ Parse state-specific electrostatic interactions [a.u.] from FDE summary
        all at once. The order is rho_A<->rho_B, rho_A<->Nuc_B, rho_B<->Nuc_A, Nuc_A<->Nuc_B
        
        Returns
        -------
        list
            A list containing electrostatic in the same order as in the output (see above) """
        self.add_variable(self.func_name(), V.fde_electrostatic)
        if self.hooks["fde_electrostatic"] in data[l_index]:
            elstat = []
            for i in range(4):
                elstat.append(float(data[l_index+i].split()[3]))
            self.print_parsed(V.fde_electrostatic, "FDE electrostatic contributions [rho_A<->rho_B, rho_A<->Nuc_B, rho_B<->Nuc_A, Nuc_A<->Nuc_B]")
            return elstat

    def fde_timing(self, l_index, data):
        """ Parses FDE timings from the FDE summary. The order is FDE-ADC, RhoA_ref, RhoB, v_emb.

        Returns
        -------
        list
            A list of tuples (CPU, wall) containing the times in seconds
            according to the above mentioned order. """
        self.add_variable(self.func_name(), V.fde_timing)
        if self.hooks["fde_timing"] in data[l_index]:
            times_list = [0 for i in range(4)]
            order_dict = {"FDE-ADC": 0, "RhoA_ref generation": 1, "RhoB generation": 2, "v_emb": 3}
            pattern = r"(?P<label>\b.+\b)\s+(?P<cpu>\d+\.\d+)\s+\(.+\)\s+(?P<wall>\d+\.\d+)"
            for i in range(4):
                match = re.search(pattern, data[l_index+4+i])
                if match:
                    lbl = match.group("label")
                    times_list[order_dict[lbl]] = (float(match.group("cpu")),
                        float(match.group("wall")))
            self.print_parsed(V.fde_timing, "final FDE timings")
            return times_list
    
    def fde_scf_vemb(self, l_index, data):
        """ Fetches HF expectation value of the embedding potential in atomic units.
        
        Returns
        -------
        float
            Energy in a.u.
        """
        self.add_variable(self.func_name(), V.fde_scf_vemb)
        if self.hooks["fde_scf_vemb"] in data[l_index]:
            self.print_parsed(V.fde_scf_vemb, "HF expectation value of the embedding potential")
            return float(data[l_index].split()[4])
        
        
