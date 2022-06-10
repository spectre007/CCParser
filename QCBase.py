#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import logging
"""
Created on Thu Jun 15 04:48:19 2017

@author: Alexander Zech
"""

class VarNames(object):
    """ General mapping of member variable names to parsing function names.
    After adding a new parsing function, add its result name here. It will
    then be callable as for instance ccp.Parser.results.my_new_var[0] or
    if it was parsed several times as ccp.Parser.results.my_new_var.get_last()"""
    # -- Energies --
    scf_energy      = "scf_energy"
    mp_energy       = "mp_energy"        # total MP energy
    mp_aaaa         = "mp_aaaa"          # mp2 same spin (alpha) component
    mp_bbbb         = "mp_bbbb"          # mp2 same spin (beta) component
    mp_abab         = "mp_abab"          # mp2 opposite spin component
    mp_nonBrill     = "mp_nonBrill"      # mp2 non-brillouin singles
    mp_ss           = "mp_ss"            # mp2 total same spin corr
    mp_os           = "mp_os"            # mp2 total opposite spin corr
    mp_correction   = "mp_correction"    # MP correlation energy
    cc_energy       = "cc_energy"        # total coupled cluster energy
    cc_correlation  = "cc_correlation"   # CC correlation energy
    exc_energy_abs  = "exc_energy_abs"
    exc_energy_rel  = "exc_energy_rel"
    mo_energies     = "mo_energies"
    nuc_repulsion   = "nuc_rep"

    # -- Excited States --
    osc_str           = "osc_str"
    amplitudes        = "ampl"
    total_dipole      = "tot_dip"
    dipole_moment     = "dip_moment"
    transition_dipole = "trans_dip"

    # -- Geometry Optimization --
    final_energy = "final_energy"
    final_xyz    = "final_xyz"

    # -- Vibrational Frequencies --
    vib_freq      = "freq"
    vib_intensity = "IR_int"

    # -- Orbital spaces --
    n_bas     = "nbas"
    n_occ     = "n_occ"
    # n_occb    = "n_occb"
    n_frz_occ = "n_frz_occ" #may need to be reworked for USCF

    # -- AO matrices --
    overlap_matrix     = "S"
    orthonorm_matrix   = "X"
    mo_coefficients    = "C"
    alpha_dens_mat     = "P_alpha"
    beta_dens_mat      = "P_beta"
    dens_mat           = "P_tot"
    multipole_operator = "multipole_op"

    # -- Charges --
    mulliken     = "mulliken"
    chelpg       = "chelpg"

    # -- General --
    xyz_coordinates    = "xyz"
    pt_order           = "pt_order"
    has_finished = "has_finished"
    basis_name   = "bas_name"
    state_label  = "state_label"
    version      = "version"

    # -- Issues & Warnings --
    has_converged = "state_cnvgd"

    # -- Frozen-Density Embedding --
    fde_omega_ref = "fde_omega_ref"
    fde_omega_I   = "fde_omega_I"
    fde_trust     = "fde_trust"
    fde_delta_lin = "fde_delta_lin"
    fde_electrostatic     = "fde_elstat"
    fde_non_electrostatic = "fde_non_elstat"
    fde_timing   = "fde_times"
    fde_scf_vemb = "fde_scf_vemb"
    fde_scf_vemb_components = "fde_scf_comp"
    fde_expansion    = "fde_expansion"
    fde_method_rhoA  = "fde_method_A"
    fde_method_rhoB  = "fde_method_B"
    fde_isA_imported = "fde_import_A"
    fde_isB_imported = "fde_import_B"
    fde_Tfunc  = "fde_Tfunc"
    fde_XCfunc = "fde_XCfunc"
    fde_Xfunc  = "fde_Xfunc"
    fde_Cfunc  = "fde_Cfunc"
    fde_non_elstat_ref = "fde_non_elstat_ref"
    # subsystem
    fde_sysA = "E_A"
    fde_sysB = "E_B"
    # nonelectrostatic
    fde_Exc_nad        = "Exc_nad"
    fde_Ts_nad         = "Ts_nad"
    fde_int_xc_nad_ref = "Vxc_nad_ref"
    fde_int_Ts_nad_ref = "Vt_nad_ref"
    fde_int_xc_nad     = "Vxc_nad"
    fde_int_Ts_nad     = "Vt_nad"
    # electrostatic
    fde_J     = "J_int"
    fde_AnucB = "AnucB"#rhoA * vB
    fde_BnucA = "BnucA"#rhoB * vA
    fde_VNN   = "V_AB"

    # -- electronic coupling --
    adia_center  = "adiabatic_center"
    adia_ss_dist = "adiabatic_ssd"
    adia_hamilt  = "adiabatic_H"
    dia_center   = "diabatic_center"
    dia_ss_dist  = "diabatic_ssd"
    dia_hamilt   = "diabatic_H"
    # from CDFT
    wf_overlap        = "overlap"
    nonorthogonal_H   = "nonorth_H"
    orthogonal_H      = "orth_H"
    nonorthogonal_dip = "nonorth_dip"
    dia_dip           = "diabatic_dip"
    adia_dip          = "adiabatic_dip"
    adia_energy       = "adiabatic_energy"
    # State-to-state
    sts_dip_mom   = "sts_dip_mom"
    sts_trans_dip = "sts_trans_dip"
    sts_coupling  = "sts_coupling"
    sts_osc_str   = "sts_osc_str"

    # -- wave function analysis --
    diff_attach_mean  = "diff_attach_mean"
    diff_detach_mean  = "diff_detach_mean"
    diff_da_dist      = "diff_da_dist"
    diff_attach_size  = "diff_attach_size"
    diff_detach_size  = "diff_detach_size"
    diff_attach_mom   = "diff_attach_mom"
    diff_detach_mom   = "diff_detach_mom"
    trans_elec_mean   = "trans_elec_mean"
    trans_hole_mean   = "trans_hole_mean"
    trans_eh_dist     = "trans_eh_dist"
    trans_elec_size   = "trans_elec_size"
    trans_hole_size   = "trans_hole_size"
    trans_elec_mom    = "trans_elec_mom"
    trans_hole_mom    = "trans_hole_mom"
    trans_eh_sep      = "trans_eh_sep"
    trans_eh_sep_mom  = "trans_eh_sep_mom"
    trans_eh_cov      = "trans_eh_cov"
    trans_eh_corr     = "trans_eh_corr"

    # -- Energy decomposition analysis --
    almo_frz_tot = "almo_frz_tot"
    almo_pol_tot = "almo_pol_tot"
    almo_ful_tot = "almo_full_tot"
    almo_frz  = "almo_frz"
    almo_pol  = "almo_pol"
    almo_ct   = "almo_ct"
    almo_disp = "almo_disp"
    almo_tot  = "almo_tot"
    almo_cls_elec  = "almo_cls_elec"
    almo_cls_pauli = "almo_cls_pauli"
    almo_mod_pauli = "almo_mod_pauli"
    almo_frg_ene   = "almo_frg_energy"
    # almo_components = "almo_comp"
    # almo_total      = "almo_tot"
    sapt0_total     = "sapt0_tot"
    sapt_components = "sapt_comp"

    # -- Input variables --
    params   = "input"
    frag_xyz = "frag_xyz"
    elconf = "elconf"


def var_tag(var_name):
    def var_decorator(func):
        def func_wrapper(self, i, data):
            self.map[func.__name__] = var_name
            return func(self, i, data)
        return func_wrapper
    return var_decorator

class QCMethod(object):
    """ Base class for Quantum Chemistry methods """
    def __init__(self):
        self.map = {}#map of function names to variable names

    def add_variable(self, func_name, var_name):
        """ Registers variable to the `out.results` ParseContainer"""
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

class AtomicBasis(object):
    def __init__(self, xyz, expo, coef, name=""):
        self.name = name
        self.center = xyz
        self.exponents = expo
        self.coefficients = coef

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

class PeriodicTable(object):
    """ Minimal implementation of a periodic table of elements. """
    PTE = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
           "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
           "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
           "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
           "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
           "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
           "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au",
           "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
           "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
           "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
           "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

    def get_atomic_num(self, symbol):
        """ Convert atomic symbol to atomic number.

        Parameters
        ----------
        symbol : string
            Atomic symbol.

        Returns
        -------
         : int
            Atomic number.
        """
        return PeriodicTable.PTE.index(symbol) + 1

    def get_atomic_sym(self, number):
        """ Convert atomic number to atomic symbol.

        Parameters
        ----------
        number : int
            Atomic number.

        Returns
        -------
         : string
            Atomic symbol.
        """
        return PeriodicTable.PTE[number-1]
