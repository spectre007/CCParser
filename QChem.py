#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 23:06:40 2017

@author: alex
"""
import re
import numpy as np
import logging
from .ParserData import MolecularOrbitals, Amplitudes
from .QCBase import QCMethod, VarNames as V
from .QCBase import var_tag
from .ParserTools import is_square
from .constants import eV2Hartree

# create module logger
mLogger = logging.getLogger("CCParser.QChem")

def extract_floats(p_string):
    """ Extracts floats from convoluted string.

    Example: '1.23345-412.2451515'. The function will
    return a list [1.23345, -412.2451515].
    """
    pattern = r"[-+]?\d+\.\d+"
    floats = re.findall(pattern, p_string)
    # if len(floats) == 0:
    #     floats = [p_string]
    # else:
    floats = list(map(float, floats))
    return floats

def clean_line_split(line_split):
    """ Clean line split if floats are glued together """
    # take care of index glued to float --> 1st element
    # e.g. 1-627.1223
    if len(line_split) == 0:
        return []
    first = extract_floats(line_split[0])
    match = re.search(r"[+-]?\d+\.\d+", first[0])
    if len(first) == 1 and match:
        line_split[0] = line_split[0].replace(first[0], '')
        line_split.insert(1, first[0])

    # if two matrix elements are stuck together
    multi_list = [extract_floats(elem) for elem in line_split]
    flattened  = [item for sublist in multi_list for item in sublist]
    return flattened

def parse_symmetric_matrix(n, readlin, asmatrix=True):
    """Parse a symmetric matrix printed columnwise

    Parameters
    ----------
    n : int
        Line number of identifier
    readlin : list
        Readlines list object
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
    ncol_ref = len(readlin[index_line].split())
    while True:#loop over blocks
        index_ls = readlin[index_line].split()
        ncol = len(index_ls)
        # abort if more or zero columns or due to stop signal
        if (ncol > ncol_ref or ncol == 0) \
        or any(stop in index_ls for stop in stop_signals):
            break
        # abort if first element is not a number
        if not index_ls[0].isdigit():
            break
        # adding rows scheme -> take line split as is
        j = 0
        if cols > 0:
            first_batch = False
        while True: #loop over lines in one block
            line_split = readlin[index_line+j+1].split()
            # abort if no elements in line split
            if len(line_split) == 0:
                break
            # abort if first element is not a number
            if not line_split[0].isdigit():
                break
            # if needed clean the split from glued floats
            if len(line_split) == ncol and \
                    len(line_split[0].split('-')) > 1:
                line_split = clean_line_split(line_split)
            # abort if incorrect number of elements or signal
            # if len(line_split) != ncol+1 \
            if any(stop in line_split for stop in stop_signals):
                break
            if first_batch:
                matrix.append([])
            # everything's fine? let's make a matrix
            matrix[j] += list(map(float, line_split[1:]))
            j += 1
        index_line += j+1#update index line
        cols += ncol#update total number of columns processed
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

def parse_STS_table(i, data, cols=[0,1], fmt=[int, float]):
    """ Parse table in STSman format.

       Electron Dipole Moments of Singlet Excited State
     -----------------------------------------------------
        State     X           Y           Z(a.u.)
     -----------------------------------------------------
           1     6.500046    1.125952    1.637921
    """
    # clean up input
    if type(cols) == int:
        cols = [cols]
    if type(fmt) != list:
        fmt = [fmt]

    # find out units from header
    i_header = i+2 #line containing table header
    # header = data[i_header].split()
    # units = []
    # mark = (-1, "n/a")
    # for p, head in enumerate(reversed(header)):
    #     if head == "State":
    #         units.append("n/a")
    #     elif head == "States":
    #         units += ["n/a", "n/a"]
    #     elif "(eV)" in head:
    #         mark = (p, "eV")
    #         units.append("eV")
    #     elif "(a.u.)":
    #         mark = (p, "n/a")
    #         units.append("n/a")
    #     elif mark[0] >= 0:
    #         units.append(mark[1])
    # units = [u for u in reversed(units)]

    i_a      = i_header+2 #first line of table
    j = 0
    values = []
    while True:
        if "-----" in data[i_a+j]:
            break
        line_split = data[i_a+j].split()
        # select columns of interest
        selection  = [line_split[m] for m in cols]
        # sel_units  = [units[m] for m in cols]
        # convert string to specified format
        values.append([conv(selection[m]) for m, conv in enumerate(fmt)])
        # if necessary convert unit to atomic units
        # for m, u in enumerate(sel_units):
        #     if u == "n/a":
        #         continue
        #     elif u == "eV":
        #         values[j][m] = values[j][m]*eV2Hartree
        j += 1
    return values

# def parse_inline_vec(line, asarray=True):
#     """ Extracts a vector of the format
#     '[ 1.000, 2.000, 3.000]' from the current line"""
#     pattern = r"[+-]?\d+\.\d*"
#     match = re.findall(pattern, line)
#     if len(match) > 0:
#         match = list(map(float, match))
#         if asarray:
#             return np.asarray(match)
#         else:
#             return match

def parse_libwfa_vec(hook_string, data, i):
    """ Extract vector from libwfa block
    which ends with an empty line."""
    j = 1
    vec = []
    while True:
        ls = data[i+j].split()
        if len(ls) == 0:
            break
        elif hook_string in data[i+j]:
            vec = extract_floats(data[i+j])
            break
        j += 1
    return vec

def parse_libwfa_float(hook_string, data, i):
    """ Extracts a float from libwfa block
    which ends with an empty line."""
    j = 1
    dist = 0.0
    while True:
        ls = data[i+j].split()
        if len(ls) == 0:
            break
        elif hook_string in data[i+j]:
            dist = float(ls[-1])
            break
        j += 1
    return dist

def parse_converged_genscfman(i, data, offset=0):
    """ Extract the total energy from converged SCF cycles."""
    j = offset
    patt = r"\s+\d+\s+(-?\d+\.\d+)\s+(\d+\.\d+[e]-\d+)\s+\d{5}\s(Convergence criterion met)"
    end_string = "Timing for Total SCF:"
    SCF_conv = 0
    while True:
        if end_string in data[i+j]:
            break
        match = re.search(patt, data[i+j])
        if match:
            SCF_conv = match.group(1)
            break
        j += 1
    return SCF_conv

###############################################################################
class Input(QCMethod):
    """Parse input section. """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {'molecule' : '$molecule',
        }

    @var_tag(V.molecule)
    def molecule(self, i, data):
        """ Get molecular/fragment configuration """
        mLogger.info("input molecule(s)", extra={"Parsed" : V.molecule})
        mol, n = {}, 0
        ifrag = 0
        curr_mol = 'total'
        mol = {'atoms' : [], 'xyz' : [], 'charge' : [], 'multiplicity' : []}
        while True:
            cwline = data[i+1+n]
            if len(cwline.split()) == 2:#may be not strict enough
                elconf = list(map(int, cwline.split()))
                mol['charge'].append(elconf[0])
                mol['multiplicity'].append(elconf[1])
            elif "$end" in cwline:
                break
            elif '--' in cwline:
                ifrag += 1
                mol['atoms'].append([])
                mol['xyz'].append([])
            elif "read" in cwline:
                mol = 'read'
            else:
                if ifrag == 0:
                    mol['atoms'].append(cwline.split()[0])
                    mol['xyz'].append(list(map(float, cwline.split()[1:])))
                else:
                    mol['atoms'][ifrag-1].append(cwline.split()[0])
                    mol['xyz'][ifrag-1].append(list(map(float, cwline.split()[1:])))
            n += 1
        return mol

class General(QCMethod):
    """Parse general information like basis set, number of atoms, etc. """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {'version' : r'Q-Chem\s*(\d+\.\d+), Q-Chem, Inc.,',
                "xyz_coordinates" : "Standard Nuclear Orientation (Angstroms)",
                "electrons" : r"There are \s+(?P<alpha>\d+) alpha and \s+(?P<beta>\d+) beta electrons",
                "nuc_rep" : "Nuclear Repulsion Energy =",
                "n_basis" : r"There are (?P<shells>\d+) shells and (?P<bsf>\d+) basis functions",
                "mulliken": "Ground-State Mulliken Net Atomic Charges",
                "chelpg": "Ground-State ChElPG Net Atomic Charges",
                "has_finished": "Thank you very much for using Q-Chem.",
                "basis_name": "Requested basis set is"}

    @var_tag(V.version)
    def version(self, i, data):
        """ Extract version number of Q-Chem """
        mLogger.info("Q-Chem version number", extra={"Parsed": V.version})
        match = re.search(self.hooks['version'], data[i])
        if match:
            return match.group(1)

    @var_tag(V.xyz_coordinates)
    def xyz_coordinates(self, i, data):
        """ Get XYZ coordinates """
        mLogger.info("XYZ coordinates", extra={"Parsed":V.xyz_coordinates})
        xyz_dat, n = [], 0
        while True:
            if "------------" in data[i+3+n]:
                break
            else:
                xyz_dat.append(data[i+3+n].split()[1:])
                n += 1
        xyz_dat = [[x[0]]+list(map(float,x[1:])) for x in xyz_dat]
        return xyz_dat

    @var_tag(V.n_occ)
    def electrons(self, i, data):
        """ Parse number of occupied alpha and beta electrons"""
        match = re.search(self.hooks["electrons"], data[i])
        if match:
            mLogger.info("number of occupied MOs (a, b)", extra={"Parsed":V.n_occ})
            return tuple(map(int, match.groups()))

    @var_tag(V.nuc_repulsion)
    def nuc_rep(self, i, data):
        """ Get nuclear repulsion energy [a.u.] """
        mLogger.info("nuclear repulsion energy",
                     extra={"Parsed":V.nuc_repulsion})
        return float(data[i].split()[-2])

    @var_tag(V.n_bas)
    def n_basis(self, i, data):
        """ Parse number of basis functions """
        match = re.search(self.hooks["n_basis"], data[i])
        if match:
            mLogger.info("number of basis functions",
                         extra={"Parsed":V.n_bas})
            return int(match.group("bsf"))

    @var_tag(V.mulliken)
    def mulliken(self, i, data):
        """ Parse ground-state mulliken charges """
        chg, n = [], 0
        while True:
            if "-------" in data[i+4+n]:
                break
            else:
                chg.append(data[i+4+n].split()[1:])
                n += 1
        chg = [[x[0]]+[float(x[1])] for x in chg]
        mLogger.info("Ground-State Mulliken charges",
                     extra={"Parsed":V.mulliken})
        return chg

    @var_tag(V.chelpg)
    def chelpg(self, i, data):
        """ Parse ground-state ChElPG charges """
        chg, n = [], 0
        while True:
            if "-------" in data[i+4+n]:
                break
            else:
                chg.append(data[i+4+n].split()[1:])
                n += 1
        chg = [[x[0]]+[float(x[1])] for x in chg]
        mLogger.info("Ground-State ChElPG charges", extra={"Parsed":V.chelpg})
        return chg

    @var_tag(V.has_finished)
    def has_finished(self, i, data):
        """ Parse final statement that indicates if Q-Chem finished
        without errors. """
        mLogger.info("whether Q-Chem has finished successfully",
                     extra={"Parsed": V.has_finished})
        return True

    @var_tag(V.basis_name)
    def basis_name(self, i, data):
        """ Parse name of basis set. """
        mLogger.info("basis set name", extra={"Parsed": V.basis_name})
        return data[i].split()[-1]


class SCF(QCMethod):
    """ Parse scfman related quantities """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        # hooks as {function_name : hook_string}
        self.hooks = {"scf_energy" : "SCF   energy in the final basis set",
                      # "mo_energies": "Orbital Energies (a.u.)",
                      "mo_energies": r"\s(Alpha|Beta) MOs(, Restricted)?",
                      "overlap_matrix" : " Overlap Matrix",
                      "orthonorm_matrix" : " Orthonormalization Matrix",
                      "alpha_density_matrix" : r"\s(Alpha Density Matrix|Final Alpha density matrix\.)",
                      "mo_coefficients_r" : "RESTRICTED (RHF) MOLECULAR ORBITAL COEFFICIENTS",
                      "multipole_op" : "Multipole Matrix "}

    @var_tag(V.scf_energy)
    def scf_energy(self, i, data):
        """ Get Hartree-Fock energy [a.u.] from scfman """
        mLogger.info("SCF energy in [a.u.]", extra={"Parsed":V.scf_energy})
        return float(data[i].split()[-1])

    @var_tag(V.mo_energies)
    def mo_energies(self, i, data):
        """ Parse Hartree-Fock molecular orbital energies """
        n = 0
        v_start, v_end = 0, 0
        while True:
            if "-----------" in data[i+1+n] or\
            len(data[i+1+n]) < 3:
                v_end = n
                break
            if "-- Virtual --" in data[i+1+n]:
                v_start = n
            n += 1
        occ_s = "".join(data[i+1:i+1+v_start])
        vrt_s = "".join(data[i+1+v_start:i+1+v_end])
        a_occ = re.findall(r"-?\d+\.\d+", occ_s, re.M)
        a_vrt = re.findall(r"-?\d+\.\d+", vrt_s, re.M)
        mo_obj = MolecularOrbitals(a_occ, a_vrt)
        mLogger.info("molecular orbital energies",
                     extra={"Parsed" : V.mo_energies})
        return mo_obj

    @var_tag(V.overlap_matrix)
    def overlap_matrix(self, i, data):
        """ Parse overlap matrix S """
        mLogger.info("overlap matrix", extra={"Parsed":V.overlap_matrix})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.orthonorm_matrix)
    def orthonorm_matrix(self, i, data):
        """ Parse orthonormalization matrix X """
        mLogger.info("orthonormalization matrix",
                     extra={"Parsed":V.orthonorm_matrix})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.alpha_dens_mat)
    def alpha_density_matrix(self, i, data):
        """ Parse alpha density matrix P_alpha """
        mLogger.info("SCF alpha density matrix",
                     extra={"Parsed":V.alpha_dens_mat})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.mo_coefficients)
    def mo_coefficients_r(self, i, data):
        """ Parse MO coefficients C for restricted SCF"""
        C = parse_AO_matrix(data, i)
        mLogger.info("{0:} molecular orbital coefficients".format(C.shape[0]),
                     extra={"Parsed":V.mo_coefficients})
        return C

    @var_tag(V.multipole_operator)
    def multipole_op(self, i, data):
        """ Parse Multipole matrix (x,x,x)."""
        M = parse_symmetric_matrix(i, data)
        mLogger.info("Multipole matrix", extra={"Parsed": V.multipole_operator})
        return M

class GENSCF(QCMethod):
    """ Parse quantitites from gen_scfman """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"density_matrix" : "Density Matrix\n"}

    @var_tag(V.dens_mat)
    def density_matrix(self, i, data):
        """ Parse density matrix
        TODO: is that the total density?"""
        mLogger.info("SCF density matrix", extra={"Parsed": V.dens_mat})
        return parse_symmetric_matrix(i, data)

class ADC(QCMethod):
    """ Parse adcman related quantities """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"scf_energy": "HF Summary",
                      "mp_energy": r"(RI-)?MP\([2-3]\) Summary",
                      "mp_correction": (r"(MP e|E)nergy contribution:\s+"
                                        r"(?P<E>-\d+\.\d+) a.u."),
                      "exc_energies": "Excitation energy:",
                      "osc_strength": "Osc. strength:",
                      "has_converged" : r"Excited state\s+\d+\s+\(.*?\)\s+\[(.*?)\]",
                      "amplitudes": "Important amplitudes:",
                      "total_dipole" : "Total dipole [Debye]",
                      "dipole_moment" : (r"Dip\. moment \[a\.u\.\]:\s+\[\s+(?P<x>-?\d+\.\d+),"
                      r"\s+(?P<y>-?\d+\.\d+),\s+(?P<z>-?\d+\.\d+)\]"),
                      "mulliken_adc": "Mulliken Population Analysis",
                      "diff_detach_mean": "Exciton analysis of the difference density matrix",
                      "diff_attach_mean": "Exciton analysis of the difference density matrix",
                      "diff_da_dist": "Exciton analysis of the difference density matrix",
                      "diff_detach_size": "Exciton analysis of the difference density matrix",
                      "diff_attach_size": "Exciton analysis of the difference density matrix",
                      "diff_detach_mom": "Exciton analysis of the difference density matrix",
                      "diff_attach_mom": "Exciton analysis of the difference density matrix",
                      "trans_hole_mean": "Exciton analysis of the transition density matrix",
                      "trans_elec_mean": "Exciton analysis of the transition density matrix",
                      "trans_eh_dist": "Exciton analysis of the transition density matrix",
                      "trans_hole_size": "Exciton analysis of the transition density matrix",
                      "trans_elec_size": "Exciton analysis of the transition density matrix",
                      "trans_eh_sep": "Exciton analysis of the transition density matrix",
                      "trans_hole_mom": "Exciton analysis of the transition density matrix",
                      "trans_elec_mom": "Exciton analysis of the transition density matrix",
                      "trans_eh_sep_mom": "Exciton analysis of the transition density matrix",
                      "trans_eh_cov": "Exciton analysis of the transition density matrix",
                      "trans_eh_corr": "Exciton analysis of the transition density matrix"}

    @var_tag(V.exc_energy_rel)
    def exc_energies(self, i, data):
        """ Parse excitation energies [eV] """
        mLogger.info("relative ADC(x) excitation energy/-ies [eV]",
                     extra={"Parsed":V.exc_energy_rel})
        return float(data[i].split()[-2])

    @var_tag(V.osc_str)
    def osc_strength(self, i, data):
        """ Parse oscillator strengths """
        mLogger.info("ADC oscillator strength/s",
                     extra={"Parsed":V.osc_str})
        return float(data[i].split()[-1])

    @var_tag(V.scf_energy)
    def scf_energy(self, i, data):
        """ Parse SCF energy [a.u.] from adcman """
        mLogger.info("SCF energy in [a.u.] (adcman)",
                     extra={"Parsed":V.scf_energy})
        return float(data[i+2].split()[-2])

    @var_tag(V.mp_energy)
    def mp_energy(self, i, data):
        """ Parse MP(x) reference state energy [a.u.] from adcman """
        match = re.search(self.hooks["mp_energy"], data[i])
        if match:
            mLogger.info("MP(x) energy in [a.u.] (adcman)",
                         extra={"Parsed":V.mp_energy})
            return float(data[i+3].split()[2])

    @var_tag(V.mp_correction)
    def mp_correction(self, i, data):
        """Parse MP(x) energy contribution in [a.u.] from adcman."""
        match = re.search(self.hooks["mp_correction"], data[i])
        if match:
            mLogger.info("MP(x) correction in [a.u.] (adcman)",
                             extra={"Parsed": V.mp_correction})
            return float(match.group("E"))

    @var_tag(V.has_converged)
    def has_converged(self, i, data):
        """ Parse if state has converged. """
        match = re.search(self.hooks["has_converged"], data[i])
        if match:
            mLogger.info("if states converged (adcman)",
                         extra={"Parsed":V.has_converged})
            return match.group(1) == "converged"

    @var_tag(V.amplitudes)
    def amplitudes(self, i, data):
        """ Parse occ -> virt amplitudes """
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

        idx = i+3
        amplist = []# Format: [[occ_i, occ_j,..., virt_a, virt_b,..., ampl], ...]
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
        mLogger.info("ADC amplitudes", extra={"Parsed":V.amplitudes})
        return Amplitudes.from_list(amplist, factor=2.0)

    @var_tag(V.total_dipole)
    def total_dipole(self, i, data):
        """ Parse total dipole moment [Debye] for HF, MP2 and ADC states

        Returns
        -------
        float
            The total dipole moment [Debye]"""
        mLogger.info("Total dipole moment [Debye]",
                     extra={"Parsed":V.total_dipole})
        return float(data[i].split()[-1])

    @var_tag(V.dipole_moment)
    def dipole_moment(self, i, data):
        """ Parse dipole moment components in [a.u.]

        Returns
        -------
        numpy.ndarray
            Dipole moment in [a.u.]
        """
        match = re.search(self.hooks["dipole_moment"], data[i])
        if match:
            mLogger.info("Dipole moment [a.u.]",
                         extra={"Parsed":V.dipole_moment})
            g = list(map(float,match.groups()))
            return np.asarray(g)

    @var_tag(V.diff_detach_mean)
    def diff_detach_mean(self, i, data):
        """ Parse mean position of detachment density [Ang] """
        mLogger.info("mean position of detachment density [Ang]",
                extra={"Parsed": V.diff_detach_mean})
        return parse_libwfa_vec("<r_h> [Ang]:", data, i)

    @var_tag(V.diff_attach_mean)
    def diff_attach_mean(self, i, data):
        """ Parse mean position of attachment density [Ang] """
        mLogger.info("mean position of attachment density [Ang]",
                extra={"Parsed": V.diff_attach_mean})
        return parse_libwfa_vec("<r_e> [Ang]:", data, i)

    @var_tag(V.diff_da_dist)
    def diff_da_dist(self, i, data):
        """ Parse linear D/A distance [Ang] """
        mLogger.info("linear D/A distance [Ang]",
                extra={"Parsed": V.diff_da_dist})
        return parse_libwfa_float("|<r_e - r_h>| [Ang]:", data, i)

    @var_tag(V.diff_detach_size)
    def diff_detach_size(self, i, data):
        """ Parse detachment size [Ang] """
        mLogger.info("RMS size of detachment density [Ang]",
                extra={"Parsed": V.diff_detach_size})
        return parse_libwfa_float("Hole size [Ang]:", data, i)

    @var_tag(V.diff_attach_size)
    def diff_attach_size(self, i, data):
        """ Parse attachment size [Ang] """
        mLogger.info("RMS size of attachment density [Ang]",
                extra={"Parsed": V.diff_attach_size})
        return parse_libwfa_float("Electron size [Ang]:", data, i)

    @var_tag(V.diff_detach_mom)
    def diff_detach_mom(self, i, data):
        """ Parse cartesian components of detachment size [Ang] """
        mLogger.info("cartesian components of detachment size [Ang]",
                extra={"Parsed": V.diff_detach_mom})
        j = 1
        vec = []
        while True:
            ls = data[i+j].split()
            if len(ls) == 0:
                break
            elif "Hole size [Ang]:" in data[i+j] and \
            "Cartesian components [Ang]:" in data[i+j+1]:
                vec = extract_floats(data[i+j+1])
                break
            j += 1
        return vec

    @var_tag(V.diff_attach_mom)
    def diff_attach_mom(self, i, data):
        """ Parse cartesian components of attachment size [Ang] """
        mLogger.info("cartesian components of attachment size [Ang]",
                extra={"Parsed": V.diff_attach_mom})
        j = 1
        vec = []
        while True:
            ls = data[i+j].split()
            if len(ls) == 0:
                break
            elif "Electron size [Ang]:" in data[i+j] and \
            "Cartesian components [Ang]:" in data[i+j+1]:
                vec = extract_floats(data[i+j+1])
                break
            j += 1
        return vec

    @var_tag(V.mulliken)
    def mulliken_adc(self, i, data):
        """ Parse MP(x) and ADC(x) mulliken charges """
        chg, n = [], 0
        while True:
            if "-------" in data[i+3+n]:
                break
            else:
                chg.append(data[i+3+n].split()[1:])
                n += 1
        chg = [[x[0]]+[float(x[1])] for x in chg]
        mLogger.info("MP(x)/ADC(x) Mulliken charges",
                     extra={"Parsed":V.mulliken})
        return chg

    @var_tag(V.trans_hole_mean)
    def trans_hole_mean(self, i, data):
        """ Parse mean position of hole [Ang] """
        mLogger.info("mean position of hole [Ang]",
                extra={"Parsed": V.trans_hole_mean})
        return parse_libwfa_vec("<r_h> [Ang]:", data, i)

    @var_tag(V.trans_elec_mean)
    def trans_elec_mean(self, i, data):
        """ Parse mean position of electron [Ang] """
        mLogger.info("mean position of electron [Ang]",
                extra={"Parsed": V.trans_elec_mean})
        return parse_libwfa_vec("<r_e> [Ang]:", data, i)

    @var_tag(V.trans_eh_dist)
    def trans_eh_dist(self, i, data):
        """ Parse linear e/h distance [Ang] """
        mLogger.info("linear e/h distance [Ang]",
                extra={"Parsed": V.trans_eh_dist})
        return parse_libwfa_float("|<r_e - r_h>| [Ang]:", data, i)

    @var_tag(V.trans_hole_size)
    def trans_hole_size(self, i, data):
        """ Parse RMS hole size [Ang] """
        mLogger.info("RMS hole size [Ang]",
                extra={"Parsed": V.trans_hole_size})
        return parse_libwfa_float("Hole size [Ang]:", data, i)

    @var_tag(V.trans_elec_size)
    def trans_elec_size(self, i, data):
        """ Parse RMS electron size [Ang] """
        mLogger.info("RMS electron size [Ang]",
                extra={"Parsed": V.trans_elec_size})
        return parse_libwfa_float("Electron size [Ang]:", data, i)

    @var_tag(V.trans_eh_sep)
    def trans_eh_sep(self, i, data):
        """ Parse RMS electron-hole separation [Ang] """
        mLogger.info("RMS electron size [Ang]",
                extra={"Parsed": V.trans_eh_sep})
        return parse_libwfa_float("RMS electron-hole separation", data, i)

    @var_tag(V.trans_hole_mom)
    def trans_hole_mom(self, i, data):
        """ Parse cartesian components of hole size [Ang] """
        mLogger.info("cartesian components of hole size [Ang]",
                extra={"Parsed": V.trans_hole_mom})
        j = 1
        vec = []
        while True:
            ls = data[i+j].split()
            if len(ls) == 0:
                break
            elif "Hole size [Ang]:" in data[i+j] and \
            "Cartesian components [Ang]:" in data[i+j+1]:
                vec = extract_floats(data[i+j+1])
                break
            j += 1
        return vec

    @var_tag(V.trans_elec_mom)
    def trans_elec_mom(self, i, data):
        """ Parse cartesian components of electron size [Ang] """
        mLogger.info("cartesian components of electron size [Ang]",
                extra={"Parsed": V.trans_elec_mom})
        j = 1
        vec = []
        while True:
            ls = data[i+j].split()
            if len(ls) == 0:
                break
            elif "Electron size [Ang]:" in data[i+j] and \
            "Cartesian components [Ang]:" in data[i+j+1]:
                vec = extract_floats(data[i+j+1])
                break
            j += 1
        return vec

    @var_tag(V.trans_eh_sep_mom)
    def trans_eh_sep_mom(self, i, data):
        """ Parse cartesian components of RMS
        electron-hole separation [Ang] """
        mLogger.info("RMS electron-hole separation [Ang]",
                extra={"Parsed": V.trans_eh_sep_mom})
        j = 1
        vec = []
        while True:
            ls = data[i+j].split()
            if len(ls) == 0:
                break
            elif "RMS electron-hole separation" in data[i+j] and \
            "Cartesian components [Ang]:" in data[i+j+1]:
                vec = extract_floats(data[i+j+1])
                break
            j += 1
        return vec

    @var_tag(V.trans_eh_cov)
    def trans_eh_cov(self, i, data):
        """ Parse electron-hole covariance [Ang^2] """
        mLogger.info("electron-hole covariance [Ang^2]",
                extra={"Parsed": V.trans_eh_cov})
        return parse_libwfa_float("Covariance(r_h, r_e)", data, i)

    @var_tag(V.trans_eh_corr)
    def trans_eh_corr(self, i, data):
        """ Parse electron-hole correlation coefficient """
        mLogger.info("electron-hole correlation coefficient",
                extra={"Parsed": V.trans_eh_corr})
        return parse_libwfa_float("Correlation coefficient:", data, i)

class FDE_ADC(QCMethod):
    """ Parsing related to FDE-ADC implementation in Q-Chem """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"omega_ref": "FDE control parameter",
                      "omega_I": "Omega(FDE)",
                      "trust": "lambda(FDE)",
                      "delta_lin": "Delta_Lin:",
                      "timing": "FDE timings",
                      "expansion": "FDE-Expansion",
                      "method_Aref": "Method to generate rhoA_ref:",
                      "method_B": "Method to generate rhoB:",
                      "method_B_legacy": "Environment method",
                      "tfunc": "Embedding Pot. T functional",
                      "xcfunc": "Embedding Pot. XC functional",
                      "xfunc" : "Embedding Pot. X functional",
                      "cfunc" : "Embedding Pot. C functional",
                      "nonel_ref" : "LinFDET terms involving",
                      "subsysA"  : "Embedded System (A):",
                      "subsysB"  : "Environment (B):",
                      "E_xc_nad" : "non-additive E_xc:",
                      "Ts_nad"   : "non-additive T_s:",
                      "V_xc_nad_ref" : "# using -REFERENCE- potential",
                      "V_xc_nad"     : "# using -UPDATED- potential",
                      "V_t_nad_ref"  : "# using -REFERENCE- potential",
                      "V_t_nad"      : "# using -UPDATED- potential",
                      "J_int"        : "rho_A <-> rho_B:",
                      "AnucB"        : "rho_A <-> Nuc_B:",
                      "BnucA"        : "rho_B <-> Nuc_A:",
                      "V_NA_NB"      : "Nuc_A <-> Nuc_B:"}

    #TODO: implement new parsing for omega_ref, keep old one for legacy purpose
    @var_tag(V.fde_omega_ref)
    def omega_ref(self, i, data):
        """ Parse FDE trust parameter (before construction of embedding potential) [ppm] """
        mLogger.info("a priori FDE overlap parameter Omega_ref",
                     extra={"Parsed":V.fde_omega_ref})
        return float(data[i].split()[5])

    @var_tag(V.fde_omega_I)
    def omega_I(self, i, data):
        """ Parse state-specific trust parameter Omega_I [ppm]
            from final FDE output """
        mLogger.info("state specific FDE overlap parameter Omega",
                     extra={"Parsed":V.fde_omega_I})
        return float(data[i].split()[1])

    @var_tag(V.fde_trust)
    def trust(self, i, data):
        """ Parse FDE trust parameter [ppm] from final FDE output """
        mLogger.info("state specific FDE overlap parameter Lambda",
                     extra={"Parsed":V.fde_trust})
        return float(data[i].split()[1])

    @var_tag(V.fde_delta_lin)
    def delta_lin(self, i, data):
        """ Parse 1st order term [eV] of LinFDET approximation """
        mLogger.info("1st order term of linearized FDET",
                     extra={"Parsed":V.fde_delta_lin})
        return float(data[i].split()[1])

    @var_tag(V.fde_timing)
    def timing(self, i, data):
        """ Parses FDE timings from the FDE summary.
        The general order of appearance is FDE-method, RhoA_ref, RhoB, v_emb,
        although RhoA_ref or RhoB might not be present due to import.

        Returns
        -------
        timings : dict
            A dictionary mapping tuples (CPU, wall) containing the times in seconds
            according to their label.
        """
        timings = {}
        pattern = r"(?P<label>\b.+\S+)\s+(?P<cpu>\d+\.\d+)\s+\(.+\)\s+(?P<wall>\d+\.\d+)"
        n = 0
        while True:
            if "------" in data[i+4+n]:
                break
            match = re.search(pattern, data[i+4+n])
            if match:
                lbl = match.group("label")
                timings[lbl] = (float(match.group("cpu")),
                                float(match.group("wall")))
            elif match == None:
                break
            n += 1

        mLogger.info("final FDE timings", extra={"Parsed":V.fde_timing})
        return timings

    @var_tag(V.fde_expansion)
    def expansion(self, i, data):
        """ Parses FDE-Expansion type [ME, SE, RADM]

        Returns
        -------
        string
            Expansion type
        """
        mLogger.info("FDET Basis set expansion",
                     extra={"Parsed":V.fde_expansion})
        return data[i].split()[1]

    @var_tag(V.fde_method_rhoA)
    def method_Aref(self, i, data):
        """ Parses which method is used to generate rhoAref (the density used
        to construct the initial embedding potential)

        So far the choices are "HF", "DFT / Func_Name" or "imported"

        Returns
        -------
        string
            RhoAref method
        """
        mLogger.info("FDET Method for rhoAref",
                     extra={"Parsed":V.fde_method_rhoA})
        s = data[i].split()
        if s[-1] == "'Densmat_A.txt')":
            return "imported"
        else:
            return s[-1]

    @var_tag(V.fde_method_rhoB)
    def method_B(self, i, data):
        """ Parses method type for rhoB.

        So far the choices are "HF", "DFT / Func_Name" or "imported"

        Returns
        -------
        string
            RhoB method
        """
        mLogger.info("FDET Method for rhoB",
                     extra={"Parsed":V.fde_method_rhoB})
        s = data[i].split()
        if s[-1] == "'Densmat_B.txt')":
            return "imported"
        else:
            return s[-1]

    @var_tag(V.fde_method_rhoB)
    def method_B_legacy(self, i, data):
        """ Parses method type for rhoB.

        So far the choices are "HF", "DFT / Func_Name" or "imported"

        Returns
        -------
        string
            RhoB method
        """
        mLogger.info("FDET Method for rhoB", extra={"Parsed":V.fde_method_rhoB})
        return data[i].split()[-1]

    @var_tag(V.fde_Tfunc)
    def tfunc(self, i, data):
        """ Determine what kinetic energy functional is used to construct
        the embedding potential.

        Returns
        -------
        string
            Abbreviation of kinetic energy functional
        """
        mLogger.info("kinetic energy functional for v_emb",
                     extra={"Parsed":V.fde_Tfunc})
        return data[i].split()[-1]

    @var_tag(V.fde_XCfunc)
    def xcfunc(self, i, data):
        """ Determine what exchange-correlation functional is used to construct
        the embedding potential.

        Returns
        -------
        string
            Abbreviation of exchange-correlation energy functional
        """
        mLogger.info("exchange-correlation functional for v_emb",
                     extra={"Parsed":V.fde_XCfunc})
        return data[i].split()[-1]

    @var_tag(V.fde_Xfunc)
    def xfunc(self, i, data):
        """ Determine what exchange functional is used to construct
        the embedding potential.

        Returns
        -------
        string
            Abbreviation of exchange energy functional
        """
        mLogger.info("exchange functional for v_emb",
                     extra={"Parsed":V.fde_XCfunc})
        return data[i].split()[-1]

    @var_tag(V.fde_Cfunc)
    def cfunc(self, i, data):
        """ Determine what correlation functional is used to construct
        the embedding potential.

        Returns
        -------
        string
            Abbreviation of correlation energy functional
        """
        mLogger.info("correlation functional for v_emb",
                     extra={"Parsed":V.fde_XCfunc})
        return data[i].split()[-1]

    @var_tag(V.fde_non_elstat_ref)
    def nonel_ref(self, i, data):
        """Parses non additive properties of the reference density

        Returns
        -------
        list
            List of [non-additive E_xc, non-additive T_s, integrated v_xc nad, integrated v_T nad]
        """
        l=[]
        mLogger.info("reference density non additive properties",
                     extra={"Parsed":V.fde_non_elstat_ref})
        l.append(float(data[i+2].split()[-2]))#non-additive E_xc
        l.append(float(data[i+3].split()[-2]))#non-additive T_s
        l.append(float(data[i+5].split()[-2]))#integrated v_xc nad
        l.append(float(data[i+6].split()[-2]))#integrated v_T nad
        return l

    @var_tag(V.fde_sysA)
    def subsysA(self, i, data):
        """Parses total energy of subsystem A WITHOUT embedding potential.

        Returns
        -------
        float
            Energy of subsystem A in [a.u.]
        """
        mLogger.info("non-additive exchange-correlation bi-functional",
                     extra={"Parsed":V.fde_sysA})
        return float(data[i].split()[-2])

    @var_tag(V.fde_sysB)
    def subsysB(self, i, data):
        """Parses total energy of subsystem B if calculated within fdeman.

        Returns
        -------
        float
            Energy of subsystem B in [a.u.]
        """
        mLogger.info("non-additive exchange-correlation bi-functional",
                     extra={"Parsed":V.fde_sysB})
        return float(data[i].split()[-2])

    @var_tag(V.fde_Exc_nad)
    def E_xc_nad(self, i, data):
        """Parses values of the non-additive exchange-correlation bi-functional

        Returns
        -------
        float
        """
        mLogger.info("non-additive exchange-correlation bi-functional",
                     extra={"Parsed":V.fde_Exc_nad})
        return float(data[i].split()[-2])

    @var_tag(V.fde_Ts_nad)
    def Ts_nad(self, i, data):
        """Parses values of the non-additive kinetic bi-functional

        Returns
        -------
        float
        """
        mLogger.info("non-additive kinetic bi-functional",
                     extra={"Parsed":V.fde_Ts_nad})
        return float(data[i].split()[-2])

    @var_tag(V.fde_int_xc_nad_ref)
    def V_xc_nad_ref(self, i, data):
        """Parses integral of some rhoA with the non-additive
        exchange-correlation potential obtained with rhoA_ref:

        :math:`\\int \\rho_{A}\\cdot v_{xc}^{nad}[\\rho_{A}^{ref},
        \\rho_{B}](\\mathbf{r})`


        Returns
        -------
        float
        """
        mLogger.info(("expectation value of non-additive XC REFERENCE "
                      "potential"),
                     extra={"Parsed":V.fde_int_xc_nad_ref})
        return float(data[i+1].split()[-2])

    @var_tag(V.fde_int_xc_nad)
    def V_xc_nad(self, i, data):
        """Parses integral of some rhoA with the non-additive
        exchange-correlation potential obtained with rhoA_emb (HF):

        :math:`\\int \\rho_{A}\\cdot v_{xc}^{nad}[\\rho_{A}^{emb,HF},
        \\rho_{B}](\\mathbf{r})`


        Returns
        -------
        float
        """
        mLogger.info(("expectation value of non-additive XC UPDATED "
                      "potential"),
                     extra={"Parsed":V.fde_int_xc_nad})
        return float(data[i+1].split()[-2])

    @var_tag(V.fde_int_Ts_nad_ref)
    def V_t_nad_ref(self, i, data):
        """Parses integral of some rhoA with the non-additive kinetic
        potential obtained with rhoA_ref:

        :math:`\\int \\rho_{A}\\cdot v_{xc}^{nad}[\\rho_{A}^{ref},
        \\rho_{B}](\\mathbf{r})`


        Returns
        -------
        float
        """
        mLogger.info(("expectation value of non-additive kinetic REFERENCE "
                      "potential"),
                     extra={"Parsed":V.fde_int_Ts_nad_ref})
        return float(data[i+2].split()[-2])

    @var_tag(V.fde_int_Ts_nad)
    def V_t_nad(self, i, data):
        """Parses integral of some rhoA with the non-additive kinetic
        potential obtained with rhoA_emb (HF):

        :math:`\\int \\rho_{A}\\cdot v_{xc}^{nad}[\\rho_{A}^{emb,HF},
        \\rho_{B}](\\mathbf{r})`


        Returns
        -------
        float
        """
        mLogger.info(("expectation value of non-additive kinetic UPDATED "
                      "potential"),
                     extra={"Parsed":V.fde_int_Ts_nad})
        return float(data[i+2].split()[-2])

    @var_tag(V.fde_J)
    def J_int(self, i, data):
        """Parses Coulomb repulsion of rhoA and rhoB.

         :math:`\\int\\int\\frac{\\rho_{A}(\\mathbf{r})\\rho_{B}(\\mathbf{r}')}
         {|\\mathbf{r}-\\mathbf{r}'|}\\mathrm{d}\\mathbf{r}'\\mathrm{d}
         \\mathbf{r}`

        Returns
        -------
        float
        """
        mLogger.info("Coulomb repulsion of rhoA-rhoB",
                     extra={"Parsed":V.fde_J})
        return float(data[i].split()[-2])

    @var_tag(V.fde_AnucB)
    def AnucB(self, i, data):
        """Parses Coulomb attraction of rhoA and NucB.

        :math:`\\int \\rho_{A}\\cdot v_B\\mathrm{d}
        \\mathbf{r}`

        """
        mLogger.info("Coulomb attraction of rhoA-nucB",
                     extra={"Parsed":V.fde_AnucB})
        return float(data[i].split()[-2])

    @var_tag(V.fde_BnucA)
    def BnucA(self, i, data):
        """Parses Coulomb attraction of rhoB and NucA.

        :math:`\\int \\rho_{B}\\cdot v_A\\mathrm{d}
        \\mathbf{r}`

        """
        mLogger.info("Coulomb attraction of rhoB-nucA",
                     extra={"Parsed":V.fde_BnucA})
        return float(data[i].split()[-2])

    @var_tag(V.fde_VNN)
    def V_NA_NB(self, i, data):
        """Parses Coulomb repulsion of nuclei NucA and NucB.

        :math:`\\sum_{A}\\sum_{B} \\frac{Z_{A}Z_{B}}{|\\mathbf{R}_{A} -
        \\mathbf{R}_{B}|}`

        """
        mLogger.info("Coulomb repulsion of nucA-nucB",
                     extra={"Parsed":V.fde_VNN})
        return float(data[i].split()[-2])

class CIS(QCMethod):
    """ Parsing related to CI Singles implementation in Q-Chem """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"exc_energies": (r"Excited state\s+(?P<state>\d+):\s*"
            r"excitation energy \(eV\)\s*=\s*(?P<energy>[-+]?\d+\.\d+)"),
            "adiabatic_center": (r"ADIABATIC CENTER\d+\s*LOCATIONS:.*"
            r"Occ\sx\s+([-]?\d+\.\d+)\sy\s+([-]?\d+\.\d+)\sz\s+([-]?\d+\.\d+) "
            r"Virt\sx\s+([-]?\d+\.\d+)\sy\s+([-]?\d+\.\d+)\sz\s+([-]?\d+\.\d+)"),
            "adiabatic_ss_dist": (r"ADIABATIC:.*sum of square of distances.*"
            r"Occ\s*(?P<occ>\d+\.\d+)\s*Virt\s*(?P<virt>\d+\.\d+)"),
            "adiabatic_hamilt": "Printing H in adiabatic representation",
            "diabatic_hamilt" : "Getting H in diabatic representation",
            "diabatic_center": (r"^DIABATIC CENTER\d+\s*LOCATIONS:.*"
            r"Occ\sx\s+([-]?\d+\.\d+)\sy\s+([-]?\d+\.\d+)\sz\s+([-]?\d+\.\d+) "
            r"Virt\sx\s+([-]?\d+\.\d+)\sy\s+([-]?\d+\.\d+)\sz\s+([-]?\d+\.\d+)"),
            "diabatic_ss_dist": (r"^DIABATIC:.*sum of square of distances.*"
            r"Occ\s*(?P<occ>\d+\.\d+)\s*Virt\s*(?P<virt>\d+\.\d+)")}

    @var_tag(V.exc_energy_rel)
    def exc_energies(self, i, data):
        """ Parse excitation energies [eV] """
        mLogger.info("relative CIS excitation energy/-ies [eV]",
                     extra={"Parsed": V.exc_energy_rel})
        match = re.search(self.hooks["exc_energies"], data[i])
        if match:
            return float(match.group("energy"))
    @var_tag(V.adia_center)
    def adiabatic_center(self, i, data):
        """ Parse adiabatic center """
        mLogger.info("adiabatic center locations",
                     extra={"Parsed": V.adia_center})
        match = re.search(self.hooks["adiabatic_center"], data[i])
        if match:
            all_coord = [float(c) for c in match.groups()]
            return [all_coord[:3], all_coord[3:]]

    @var_tag(V.adia_ss_dist)
    def adiabatic_ss_dist(self, i, data):
        """ Parse sum of square distances to other states """
        mLogger.info("adiabatic sum of squares of distance",
                     extra={"Parsed": V.adia_ss_dist})
        match = re.search(self.hooks["adiabatic_ss_dist"], data[i])
        if match:
            return [float(match.group("occ")), float(match.group("virt"))]

    @var_tag(V.adia_hamilt)
    def adiabatic_hamilt(self, i, data):
        """ Parse the adiabatic representation of the Hamiltonian """
        mLogger.info("adiabatic representation of the Hamiltonian",
                extra={"Parsed": V.adia_hamilt})
        j = 1
        H = []
        while True:
            ls = data[i+j].split()
            if not "showmatrix" in data[i+j]:
                break
            elif not "=" in data[i+j]:
                break
            elif len(ls) == 0:
                break
            H.append(float(ls[-1]))
            j += 1
        H = np.asarray(H)
        if is_square(len(H)):
            n = int(np.sqrt(len(H)))
            H = H.reshape((n,n))
        return H

    @var_tag(V.dia_hamilt)
    def diabatic_hamilt(self, i, data):
        """ Parse the diabatic representation of the Hamiltonian """
        mLogger.info("diabatic representation of the Hamiltonian",
                extra={"Parsed": V.dia_hamilt})
        j = 1
        H = []
        while True:
            ls = data[i+j].split()
            if not "showmatrix" in data[i+j]:
                break
            elif not "=" in data[i+j]:
                break
            elif len(ls) == 0:
                break
            H.append(float(ls[-1]))
            j += 1
        H = np.asarray(H)
        if is_square(len(H)):
            n = int(np.sqrt(len(H)))
            H = H.reshape((n,n))
        return H

    @var_tag(V.dia_center)
    def diabatic_center(self, i, data):
        """ Parse diabatic center """
        mLogger.info("adiabatic center locations",
                     extra={"Parsed": V.dia_center})
        match = re.search(self.hooks["diabatic_center"], data[i])
        if match:
            all_coord = [float(c) for c in match.groups()]
            return [all_coord[:3], all_coord[3:]]

    @var_tag(V.dia_ss_dist)
    def diabatic_ss_dist(self, i, data):
        """ Parse sum of square distances to other states """
        mLogger.info("diabatic sum of squares of distance",
                     extra={"Parsed": V.dia_ss_dist})
        match = re.search(self.hooks["diabatic_ss_dist"], data[i])
        if match:
            return [float(match.group("occ")), float(match.group("virt"))]


class CDFTCI(QCMethod):
    """ Parsing related to CDFT-CI implementation in Q-Chem """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"overlap_matrix": "CDFT-CI overlap matrix",
                "non_orthogonal_hamilt": "CDFT-CI Hamiltonian matrix in non-orthogonal basis",
                "orthogonal_hamilt": "CDFT-CI Hamiltonian matrix in orthogonalized basis",
                "nonorthogonal_dip": r"dipole [xyz] component in nonorthogonal basis",
                "diabatic_dip": r"dipole [xyz] component in diabatic basis",
                "adiabatic_dip": r"dipole [xyz] component in adiabatic basis",
                "adiabatic_energy": r"CDFT-CI Energy state \d+"}

    @var_tag(V.wf_overlap)
    def overlap_matrix(self, i, data):
        """ Parse CDFT-CI overlap matrix. """
        mLogger.info("CDFT-CI overlap matrix", extra={"Parsed": V.wf_overlap})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.nonorthogonal_H)
    def non_orthogonal_hamilt(self, i, data):
        """ Parse CDFT-CI Hamiltonian matrix in non-orthogonal basis."""
        mLogger.info("CDFT-CI H in non-orth. basis", extra={"Parsed": V.nonorthogonal_H})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.orthogonal_H)
    def orthogonal_hamilt(self, i, data):
        """ Parse CDFT-CI Hamiltonian matrix in orthogonalized basis."""
        mLogger.info("CDFT-CI H in orthogonal basis", extra={"Parsed": V.orthogonal_H})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.nonorthogonal_dip)
    def nonorthogonal_dip(self, i, data):
        """ Parse x,y, or z component of dipole tensor in nonorthogonal basis"""
        mLogger.info("dipole matrix in non-orth. basis",
                extra={"Parsed": V.nonorthogonal_dip})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.dia_dip)
    def diabatic_dip(self, i, data):
        """ Parse x,y, or z component of dipole tensor in diabatic basis"""
        mLogger.info("dipole matrix in diabatic basis",
                extra={"Parsed": V.dia_dip})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.adia_dip)
    def adiabatic_dip(self, i, data):
        """ Parse x,y, or z component of dipole tensor in adiabatic basis"""
        mLogger.info("dipole matrix in adiabatic basis",
                extra={"Parsed": V.adia_dip})
        return parse_symmetric_matrix(i, data)

    @var_tag(V.adia_energy)
    def adiabatic_energy(self, i, data):
        """ Parse energy of CDFT-CI adiabatic states """
        mLogger.info("CDFT-CI energy of adiabatic state",
                extra={"Parsed": V.adia_energy})
        return float(data[i].split()[-1])

class TDDFT(QCMethod):
    """ Parsing related to TDDFT implementation in Q-Chem """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {"dip_mom_ES": "Electron Dipole Moments of Singlet Excited State",
                "trans_dip": r"GMH Couplings Between( Ground and)? Singlet Excited States",
                "coupling": r"GMH Couplings Between( Ground and)? Singlet Excited States",
                "oscillator_strength": r"Transition Moments Between( Ground and)? Singlet Excited States"}

    # TODO: not the same format as in other calculations... this gives
    # blocks of dipole moments instead of individual ones.
    @var_tag(V.sts_dip_mom)
    def dip_mom_ES(self, i, data):
        """ Parse dipole moments [a.u.] of singlet excited states"""
        mLogger.info("singlet ES dipole moment",
                extra={"Parsed": V.sts_dip_mom})
        return parse_STS_table(i, data, cols=[1,2,3], fmt=[float, float, float])

    @var_tag(V.sts_trans_dip)
    def trans_dip(self, i, data):
        """ Parse transition dipole moments [a.u.] between states"""
        mLogger.info("transition dipole moments [a.u.]",
                extra={"Parsed":V.sts_trans_dip})
        return parse_STS_table(i, data, cols=[0,1,2,3,4],
                fmt=[int, int, float, float, float])

    @var_tag(V.sts_coupling)
    def coupling(self, i, data):
        """ Parse coupling between states [eV]"""
        mLogger.info("coupling between states [eV]",
                extra={"Parsed": V.sts_coupling})
        return parse_STS_table(i, data, cols=[0,1,5], fmt=[int, int, float])

    # TODO: is this normal TDDFT output?
    @var_tag(V.sts_osc_str)
    def oscillator_strength(self, i, data):
        """ Parse oscillator strenghts"""
        mLogger.info("transition moment strength [a.u.]",
                extra={"Parsed": V.sts_osc_str})
        return parse_STS_table(i, data, cols=[0,1,5], fmt=[int, int, float])

class ALMO(QCMethod):
    """ Parse general ALMO output """
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {}
        # self.hooks = dict.fromkeys(["SCF_frz_simple", "SCF_pol_simple",
        #     "SCF_ct_simple", "SCF_tot_simple"],
        #         "SCF only Simplified EDA Summary (kJ/mol)")
        self.hooks.update(dict.fromkeys(["MP2_frz_corr", "MP2_pol_corr",
            "MP2_disp_corr", "MP2_ct_corr", "MP2_tot_corr"],
            "MP2 only EDA results (Hartree)"))
        self.hooks.update(dict.fromkeys(["E_frz", "E_pol", "E_disp", "E_ct",
            "E_tot"], "Total EDA results (Hartree)"))

        self.hooks["Escf_frz"] = "Energy of the unrelaxed supersystem initial determinant"
        # self.hooks["Eorth_decomp"] = "Orthogonal Decomposition of the Initial Supersystem Wavefunction"
        self.hooks["Escf_pol"]  = "for the determination of the polarized wavefunction"
        self.hooks["Escf_full"] = "for the determination of the CT-Allowed wavefunction"

        self.hooks["E_cls_elec"]  = "E_cls_elec  (CLS ELEC)  (kJ/mol)"
        self.hooks["E_cls_pauli"] = "E_cls_pauli (CLS PAULI) (kJ/mol)"
        self.hooks["SCF_disp"] = r"E_disp\s+\(DISP\)\s+\(kJ/mol\) = ([-]?\d+\.\d+)"
        self.hooks["SCF_frz"] = "E_frz (kJ/mol) ="
        self.hooks["SCF_pol"] = "E_pol (kJ/mol) ="
        self.hooks["SCF_ct"]  = r"E_vct \(kJ/mol\) = ([-]?\d+\.\d+) \(CP-corrected\)"
        self.hooks["SCF_tot"] = r"E_int \(kJ/mol\) = ([-]?\d+\.\d+) \(CP-corrected\)"
        # simplified EDA summary
        # self.hooks["prep_simple"] = r"PREPARATION\s+(-?\d+\.\d+)"
        self.hooks["frz_simple"]  = r"FROZEN\s+(-?\d+\.\d+)"
        self.hooks["disp_simple"] = r"DISPERSION\s+(-?\d+\.\d+)"
        self.hooks["pol_simple"]  = r"POLARIZATION\s+(-?\d+\.\d+)"
        self.hooks["ct_simple"]   = r"CHARGE TRANSFER\s+(-?\d+\.\d+)"
        self.hooks["tot_simple"]  = r"TOTAL\s+(-?\d+\.\d+)\s+\(.+\)"
        self.hooks["frag_energy"] = "Fragment Energies (Ha):"

    @var_tag(V.almo_frz_tot)
    def Escf_frz(self, i, data):
        """ Parse total SCF Frozen energy"""
        mLogger.info("total frozen energy E_frz [a.u.]",
                extra={"Parsed" : V.almo_frz_tot})
        return float(data[i].split()[-1])

    @var_tag(V.almo_pol_tot)
    def Escf_pol(self, i, data):
        """ Parse total SCF polarization energy"""
        mLogger.info("total polarization energy E_pol [a.u.]",
                extra={"Parsed" : V.almo_pol_tot})
        return parse_converged_genscfman(i, data, offset=5)

    @var_tag(V.almo_ful_tot)
    def Escf_full(self, i, data):
        """ Parse total (CT-allowed) SCF energy"""
        mLogger.info("total supersystem energy E_full [a.u.]",
                extra={"Parsed" : V.almo_ful_tot})
        return parse_converged_genscfman(i, data, offset=5)

    @var_tag(V.almo_cls_elec)
    def E_cls_elec(self, i, data):
        """ Parse classical electrostatics component of frozen energy (SCF)"""
        mLogger.info("classical electrostatic component of E_frz [kJ/mol]",
                extra={"Parsed" : V.almo_cls_elec})
        return float(data[i].split()[-1])

    @var_tag(V.almo_cls_pauli)
    def E_cls_pauli(self, i, data):
        """ Parse exchange component of frozen energy (SCF)"""
        mLogger.info("exchange component of E_frz [kJ/mol]",
                extra={"Parsed" : V.almo_cls_pauli})
        return float(data[i].split()[-1])

    @var_tag(V.almo_disp)
    def SCF_disp(self, i, data):
        """ Parse dispersion energy from Decomposition of frozen interaction energy [kJ/mol] """
        mLogger.info("ALMO-SCF dispersion energy [kJ/mol]",
                extra={"Parsed" : V.almo_disp})
        return float(data[i].split()[-1])

    @var_tag(V.almo_frz)
    def SCF_frz(self, i, data):
        """ Parse frozen energy from SCF EDA2 summary [kJ/mol]"""
        mLogger.info("ALMO-SCF frozen energy [kJ/mol]",
                extra={"Parsed" : V.almo_frz})
        return float(data[i].split()[-1])

    @var_tag(V.almo_pol)
    def SCF_pol(self, i, data):
        """ Parse frozen energy from SCF EDA2 summary [kJ/mol] """
        mLogger.info("ALMO-SCF polarization energy [kJ/mol]",
                extra={"Parsed" : V.almo_pol})
        return float(data[i].split()[-1])

    @var_tag(V.almo_ct)
    def SCF_ct(self, i, data):
        """ Parse charge-transfer energy from SCF EDA2 summary [kJ/mol] """
        mLogger.info("ALMO-SCF charge-transfer energy (incl. BSSE corr.) [kJ/mol]",
                extra={"Parsed" : V.almo_ct})
        return float(data[i].split()[-2])

    @var_tag(V.almo_tot)
    def SCF_tot(self, i, data):
        """ Parse interaction energy from simplified SCF EDA summary [kJ/mol] """
        mLogger.info("ALMO-SCF total interaction energy [kJ/mol]",
                extra={"Parsed" : V.almo_tot})
        return float(data[i].split()[-2])

    @var_tag(V.almo_frz)
    def frz_simple(self, i, data):
        """ Parse frozen energy from simplified SCF EDA summary [kJ/mol] """
        mLogger.info("ALMO-SCF frozen energy [kJ/mol]",
                extra={"Parsed" : V.almo_frz})
        # return float(data[i+3].split()[1])
        match = re.search(self.hooks["frz_simple"], data[i])
        if match:
            return float(match.groups()[0])

    @var_tag(V.almo_disp)
    def disp_simple(self, i, data):
        """ Parse dispersion energy from simplified EDA summary [kJ/mol] """
        mLogger.info("ALMO-SCF dispersion energy [kJ/mol]",
                extra={"Parsed" : V.almo_disp})
        match = re.search(self.hooks["disp_simple"], data[i])
        if match:
            return float(match.groups()[0])

    @var_tag(V.almo_pol)
    def pol_simple(self, i, data):
        """ Parse polarization energy from simplified EDA summary [kJ/mol] """
        mLogger.info("ALMO-SCF polarization energy [kJ/mol]",
                extra={"Parsed" : V.almo_pol})
        # return float(data[i+4].split()[1])
        match = re.search(self.hooks["pol_simple"], data[i])
        if match:
            return float(match.groups()[0])

    @var_tag(V.almo_ct)
    def ct_simple(self, i, data):
        """ Parse charge-transfer energy from simplified EDA summary [kJ/mol] """
        mLogger.info("ALMO-SCF charge-transfer energy (incl. BSSE corr.) [kJ/mol]",
                extra={"Parsed" : V.almo_ct})
        # return float(data[i+5].split()[2])
        match = re.search(self.hooks["ct_simple"], data[i])
        if match:
            return float(match.groups()[0])

    @var_tag(V.almo_tot)
    def tot_simple(self, i, data):
        """ Parse interaction energy from simplified EDA summary [kJ/mol] """
        mLogger.info("ALMO-SCF total interaction energy [kJ/mol]",
                extra={"Parsed" : V.almo_tot})
        # return float(data[i+6].split()[1])
        match = re.search(self.hooks["tot_simple"], data[i])
        if match:
            return float(match.groups()[0])

    @var_tag(V.almo_frz)
    def MP2_frz_corr(self, i, data):
        """ Parse correlation contribution to frozen energy (MP2) [Hartree] """
        mLogger.info("ALMO correlation contribution (MP2) to E_frz [a.u.]",
                extra={"Parsed" : V.almo_frz})
        return float(data[i+1].split()[1])

    @var_tag(V.almo_pol)
    def MP2_pol_corr(self, i, data):
        """ Parse correlation contribution to polarization energy (MP2) [Hartree] """
        mLogger.info("ALMO correlation contribution (MP2) to E_pol [a.u.]",
                extra={"Parsed" : V.almo_pol})
        return float(data[i+2].split()[1])

    @var_tag(V.almo_disp)
    def MP2_disp_corr(self, i, data):
        """ Parse correlation contribution to dispersion energy (MP2) [Hartree] """
        mLogger.info("ALMO correlation contribution (MP2) to E_disp [a.u.]",
                extra={"Parsed" : V.almo_disp})
        return float(data[i+3].split()[1])

    @var_tag(V.almo_ct)
    def MP2_ct_corr(self, i, data):
        """ Parse correlation contribution to charge-transfer energy (MP2) [Hartree] """
        mLogger.info("ALMO correlation contribution (MP2) to E_ct [a.u.]",
                extra={"Parsed" : V.almo_ct})
        return float(data[i+4].split()[1])

    @var_tag(V.almo_tot)
    def MP2_tot_corr(self, i, data):
        """ Parse correlation contribution to total interaction energy (MP2) [Hartree] """
        mLogger.info("ALMO correlation contribution (MP2) to E_int [a.u.]",
                extra={"Parsed" : V.almo_tot})
        return float(data[i+5].split()[2])

    @var_tag(V.almo_frz)
    def E_frz(self, i, data):
        """ Parse ALMO-EDA frozen energy [Hartree] """
        mLogger.info("ALMO-EDA frozen energy E_frz [a.u.]",
                extra={"Parsed" : V.almo_frz})
        return float(data[i+1].split()[1])

    @var_tag(V.almo_pol)
    def E_pol(self, i, data):
        """ Parse ALMO-EDA polarization energy [Hartree] """
        mLogger.info("ALMO-EDA polarization energy E_pol [a.u.]",
                extra={"Parsed" : V.almo_pol})
        return float(data[i+2].split()[1])

    @var_tag(V.almo_disp)
    def E_disp(self, i, data):
        """ Parse ALMO-EDA dispersion energy [Hartree] """
        mLogger.info("ALMO-EDA dispersion energy E_disp [a.u.]",
                extra={"Parsed" : V.almo_disp})
        return float(data[i+3].split()[1])

    @var_tag(V.almo_ct)
    def E_ct(self, i, data):
        """ Parse ALMO-EDA charge-transfer energy (incl. CP-corr) [Hartree] """
        mLogger.info("ALMO-EDA charge-transfer energy E_ct (incl. BSSE corr.) [a.u.]",
                extra={"Parsed" : V.almo_ct})
        return float(data[i+4].split()[1])

    # @var_tag(V.almo_ct)
    # def E_ct_nobsse(self, i, data):
    #     """ Parse ALMO-EDA charge-transfer energy (excl. CP-corr) [Hartree] """
    #     mLogger.info("ALMO-EDA charge-transfer energy E_ct (excl. BSSE corr.) [a.u.]",
    #             extra={"Parsed" : V.almo_frz})
    #     return float(data[i+4].split()[1])

    @var_tag(V.almo_tot)
    def E_tot(self, i, data):
        """ Parse total interaction energy [Hartree] """
        mLogger.info("ALMO-EDA total interaction energy [a.u.]",
                extra={"Parsed" : V.almo_tot})
        return float(data[i+5].split()[2])

    @var_tag(V.almo_frg_ene)
    def frag_energy(self, i, data):
        """ Parse total energies [Hartree] of isolated fragments """
        mLogger.info("fragment energies (isolated) [a.u.]",
                extra={"Parsed" : V.almo_frg_ene})
        k = 0
        e_frag = []
        while True:
            cwline = data[i+1+k]
            if "-----" in cwline:
                break
            e_frag.append(cwline.split()[-1])
            k += 1
        return list(map(float, e_frag))

class RIMP2(QCMethod):
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {
                "mp2_aaaa": "aaaa    correlation energy =",
                "mp2_abab": "abab    correlation energy =",
                "mp2_bbbb": "bbbb    correlation energy =",
                "mp2_nbs": "non-Brillouin singles      =",
                "mp2_tot_ss": "total same-spin energy     =",
                "mp2_tot_os": "total opposite-spin energy =",
                "mp2_corr": "Total  RIMP2   correlation energy =",
                "mp2_tot": "RIMP2         total energy =",
                "n_frz_occ": "# of frozen core orbitals:"}

    @var_tag(V.mp_aaaa)
    def mp2_aaaa(self, i, data):
        """ Parse MP2 same spin component (aaaa) [Hartree] """
        mLogger.info("MP2 same spin aaaa [a.u.]",
                extra={"Parsed" : V.mp_aaaa})
        return float(data[i].split()[-2])

    @var_tag(V.mp_bbbb)
    def mp2_bbbb(self, i, data):
        """ Parse MP2 same spin component (bbbb) [Hartree] """
        mLogger.info("MP2 same spin bbbb [a.u.]",
                extra={"Parsed" : V.mp_bbbb})
        return float(data[i].split()[-2])

    @var_tag(V.mp_abab)
    def mp2_abab(self, i, data):
        """ Parse MP2 same spin component (abab) [Hartree] """
        mLogger.info("MP2 same spin abab [a.u.]",
                extra={"Parsed" : V.mp_abab})
        return float(data[i].split()[-2])

    @var_tag(V.mp_nonBrill)
    def mp2_nbs(self, i, data):
        """ Parse MP2 non-Brillouin singles [Hartree] """
        mLogger.info("MP2 non-Brillouin singles [a.u.]",
                extra={"Parsed" : V.mp_nonBrill})
        return float(data[i].split()[-2])

    @var_tag(V.mp_ss)
    def mp2_tot_ss(self, i, data):
        """ Parse MP2 total same spin [Hartree] """
        mLogger.info("MP2 total same spin energy [a.u.]",
                extra={"Parsed" : V.mp_ss})
        return float(data[i].split()[-2])

    @var_tag(V.mp_os)
    def mp2_tot_os(self, i, data):
        """ Parse MP2 total opposite spin [Hartree] """
        mLogger.info("MP2 total opposite spin energy [a.u.]",
                extra={"Parsed" : V.mp_os})
        return float(data[i].split()[-2])

    @var_tag(V.mp_correction)
    def mp2_corr(self, i, data):
        """ Parse MP2 correlation energy [Hartree] """
        mLogger.info("MP2 correlation energy [a.u.]",
                extra={"Parsed" : V.mp_correction})
        return float(data[i].split()[-2])

    @var_tag(V.mp_energy)
    def mp2_tot(self, i, data):
        """ Parse MP2 total energy [Hartree] """
        mLogger.info("MP2 total energy [a.u.]",
                extra={"Parsed" : V.mp_energy})
        return float(data[i].split()[-2])

    @var_tag(V.n_frz_occ)
    def n_frz_occ(self, i, data):
        """ Number of frozen occupied orbitals """
        mLogger.info("number of frozen occupied MOs",
                extra={"Parsed" : V.n_frz_occ})
        return int(data[i].split()[-1])

class GeometryOpt(QCMethod):
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {
            "final_energy" : r"Final energy is\s+([-]?\d+\.\d+)",
            "final_xyz" : "OPTIMIZATION CONVERGED",
        }

    @var_tag(V.final_energy)
    def final_energy(self, i, data):
        """ Parse total energy of converged molecular structure [Hartree]. """
        mLogger.info("final energy (GeomOpt) [Hartree]",
                extra={"Parsed" : V.final_energy})
        return float(data[i].split()[-1])

    @var_tag(V.final_xyz)
    def final_xyz(self, i, data):
        """ Parse converged molecular structure in XYZ format. """
        mLogger.info("final XYZ (GeomOpt)", extra={"Parsed" : V.final_xyz})
        n = 0
        xyz = []
        while True:
            line = data[i+5+n]
            if "Z-matrix Print:" in line:
                break
            if len(line.split()) > 0:
                xyz.append(line.split()[1:])
            n +=1
        xyz = [[x[0]]+list(map(float,x[1:])) for x in xyz]
        return xyz

class Frequencies(QCMethod):
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {
            'frequencies' : 'VIBRATIONAL ANALYSIS',
            'ir_intensity': 'VIBRATIONAL ANALYSIS',
        }

    @var_tag(V.vib_freq)
    def frequencies(self, i, data):
        """ Parse vibrational frequencies in [cm^-1] """
        mLogger.info("vibrational frequencies [cm^-1]",
            extra={"Parsed" : V.vib_freq})
        n = 9
        freq = []
        while True:
            line = data[i+n]
            if 'STANDARD THERMODYNAMIC QUANTITIES AT' in line:
                break
            if "Frequency:" in line:
                freq += extract_floats(line)
            n += 1
        return freq

    @var_tag(V.vib_intensity)
    def ir_intensity(self, i, data):
        """ Parse IR intensities in [KM/mol] """
        mLogger.info("vibrational frequencies.",
            extra={"Parsed" : V.vib_intensity})
        n = 9
        intensities = []
        while True:
            line = data[i+n]
            if 'STANDARD THERMODYNAMIC QUANTITIES AT' in line:
                break
            if "IR Intens:" in line:
                intensities += extract_floats(line)
            n += 1
        return intensities

class CoupledCluster(QCMethod):
    def __init__(self):
        super().__init__()# necessary for every derived class of QCMethod
        self.hooks = {
            "libpt_mp2": r"MP2 energy\s+=\s+([-]?\d+\.\d+)",
            "libpt_ecorr": r"(CC[SDT\(\)]+) correlation energy\s+=\s+([-]?\d+\.\d+)",
            "libpt_etot": r"(CC[SDT\(\)]+) total energy\s+=\s+([-]?\d+\.\d+)",
        }

    @var_tag(V.mp_energy)
    def libpt_mp2(self, i, data):
        """ Parse MP2 energy [Hartree] from libpt module"""
        mLogger.info("total MP2 energy [a.u.]",
                extra={"Parsed" : V.mp_energy})
        return float(data[i].split()[-1])

    @var_tag(V.cc_correlation)
    def libpt_ecorr(self, i, data):
        """ Parse Coupled Cluster correlation energy [Hartree] in libpt module"""
        mLogger.info("Coupled Cluster correlation energy [a.u.]",
                extra={"Parsed" : V.cc_correlation})
        return float(data[i].split()[-1])

    @var_tag(V.cc_energy)
    def libpt_etot(self, i, data):
        """ Parse Coupled Cluster total energy [Hartree] in libpt module"""
        mLogger.info("Coupled Cluster total energy [a.u.]",
                extra={"Parsed" : V.cc_energy})
        return float(data[i].split()[-1])
