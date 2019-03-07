import pandas as pd
import numpy as np
from .QCBase import VarNames

class Exporter(object):
    """ Export class which writes parsed data to a certain format"""
    valid_formats = ["pdf", "xlsx", "txt", "csv", "dataframe"]
    
    def __init__(self, data=None):
        self.data = data
    
    # for later: add pandas independent functions to export arrays to file
        
    def arrays_to_dframe(self, **kwargs):
        """ Using keyworded arguments, expects arrays  """
        try:
            df = pd.DataFrame(kwargs)
        except ValueError: #if arrays do not have the same length
            d = {}
            for key, value in kwargs.items():
                d[key] = pd.Series(value)
            df = pd.DataFrame(d)
        return df
    
    def ExcitedStateSummary(self, results, fname="es_smry", fmt="csv",
                            ground_state=False):
        """ Exports energy related excited state quantities to file 
        
        Parameters
        ----------
        results : CCParser.ParseContainer
            Parsing container that holds parsed values.
        fname : string
            Filename prefix.
        fmt : string
            Output format ('csv', 'xlsx'/'xls' or 'df' for pandas.DataFrame).
        ground_state : bool
            Whether to include an empty line in the table for the ground state.
        """
        if fmt not in Exporter.valid_formats:
            raise ValueError("File format '{0:}' not recognized or supported!".format(fmt))
        if False in getattr(results, VarNames.has_converged).data:
            raise ValueError("Not converged state detected!")
        d = {}
        # (1) Excitation energies (default minimum)
        #if hasattr(results, VarNames.exc_energy_rel):
        d[VarNames.exc_energy_rel] = getattr(results, VarNames.exc_energy_rel).data
        n_states = len(d[VarNames.exc_energy_rel])
        # (2) Oscillator strengths
        if hasattr(results, VarNames.osc_str):
            d[VarNames.osc_str] = getattr(results, VarNames.osc_str).data
        # (3) Amplitudes
        if hasattr(results, VarNames.amplitudes):
            ampl = getattr(results, VarNames.amplitudes)
            pieces = [a.to_dataframe() for a in ampl]
            key = [x for x in range(1,len(pieces)+1)]
            amp_df = pd.concat(pieces, keys=key, names=["State", "Row ID"])
            
        # prepare MultiIndex (there has to be a better way to do that...)
        arrays = [[x for x in range(1, n_states+1)],
                  [0 for x in range(n_states)]]
        tuples = list(zip(*arrays))# asterisk unpacks
        df1 = pd.DataFrame(d)
        df1.index = pd.MultiIndex.from_tuples(tuples, names=["State", "Row ID"])
        
        df = pd.concat([df1, amp_df], axis=1)
        # add row to MultiIndex, see https://stackoverflow.com/q/24917700
        if ground_state:
            df.ix[(0,0),:] = np.nan
            df.sort_index(level=0, inplace=True)
        
        # EXPORT TO FILE or dataframe
        fout = fname + "." + fmt
        if fmt == "csv":
            df.to_csv(fout, encoding="utf-8")
        elif fmt == ("xlsx" or "xls"):
            writer = pd.ExcelWriter(fout)
            df.to_excel(writer, "Sheet1")
            writer.save()
        elif fmt.lower() == ("dataframe" or "df"):
            return df
        
    def ReducedWeights(self, results, nbsfA, extern=None, fmt="print",
                       fname="AmplAnl", silent=False):
        """ Calculate reduced weights based on fragment information.
        
        The reduced weight for a single excitation :math:`i \\rightarrow a` is defined as       
        :math:`v_{i}^{a} = 0.5\\cdot(c_{i,A}^{2} + c_{a,A}^{2})\\cdot w_{i}^{a}`, with
        c and w being the molecular orbital coefficient and transition weight,
        respectively.
        The MO coefficients from the output first have to be transformed to an
        orthonormal basis.
        
        Parameters
        ----------
        results : CCParser.ParseContainer
            Container object which contains excited state amplitudes
        nbsfA : int
            Number of basis functions on System A (assumes system A comes first!)
        extern : CCParser.ParseContainer
            Optional second container which contains orthonormalisation matrix and/or MO coefficients
        fmt : string
            Output format. Available are "print", "dataframe", "xlsx" or "csv"
        fname : string
            Output file name (basename only).
        silent : bool
            Whether to ignore lengthy printouts.

        """
        # consistency
        has_extern = True if extern != None else False
        if False in getattr(results, VarNames.has_converged).data:
            raise ValueError("Not converged state detected!")
        if not has_extern and not hasattr(results, VarNames.orthonorm_matrix):
            raise AttributeError("Could not find orthonormalization matrix! Was it parsed?")
        elif has_extern and not hasattr(extern, VarNames.orthonorm_matrix):
            raise AttributeError("Could not find orthonormalization matrix! Was it parsed?")
        elif not has_extern and not hasattr(results, VarNames.mo_coefficients):
            raise AttributeError("Could not find MO coefficients! Were they parsed?")
        elif has_extern and not hasattr(extern, VarNames.mo_coefficients):
            raise AttributeError("Could not find MO coefficients! Were they parsed?")
        elif not hasattr(results, VarNames.amplitudes):
            raise AttributeError("Could not find amplitudes! Were they parsed?")
        elif not hasattr(results, VarNames.n_bas):
            raise AttributeError("Could not find number of basis functions! Was it parsed?")
        else:
            # (1) Orthonormalization matrix, hardcoded last
            X = getattr(results, VarNames.orthonorm_matrix).get_last() if not \
                has_extern else getattr(extern, VarNames.orthonorm_matrix).get_last()
            X_inv = np.linalg.inv(X)
            # (2) MO coeffiecients, hardcoded last
            C = getattr(results, VarNames.mo_coefficients).get_last() if not \
                has_extern else getattr(extern, VarNames.mo_coefficients).get_last()
            C_prime = C * X_inv # Szabo, Ostlund, page 142
            max_mo = C.shape[0]
            # (3) Amplitudes
            ampl = getattr(results, VarNames.amplitudes)
            n_states = len(ampl)
            # (4) Number of basis functions
            nbsf = getattr(results, VarNames.n_bas).get_last()
            # (4) Output variables
            sum_weights = [0 for i in range(n_states)]
            sum_redweights = [0 for i in range(n_states)]
            # --------------
            sos_A = [0 for a in range(C_prime.shape[0])]
            sos_B = [0 for a in range(C_prime.shape[0])]
            for c, vect in enumerate(C_prime):
                for n in range(nbsf):
                    if n < nbsfA:
                        sos_A[c] += vect[0,n]**2
                    else:
                        sos_B[c] += vect[0,n]**2
            for i,a in enumerate(ampl):#state
                for t in range(len(a.occ)):#transition
                    if max(a.virt[t]) > max_mo:
                        if not silent:
                            print("State {0:>2d}: Omitting transition with weight \
{1:.1%} due to missing MO coefficients.".format(i+1, a.weights[t]))
                        continue
                    if len(a.occ[t]) == 1:#single amplitudes
                        rw = 0.5*(sos_A[a.occ[t][0]-1] + sos_A[a.virt[t][0]-1]) * a.weights[t]
                    elif len(a.occ[t]) == 2:#double amplitudes
                        rw = 0.25*(sos_A[a.occ[t][0]-1] + sos_A[a.occ[t][1]-1] +
                                  sos_A[a.virt[t][0]-1] + sos_A[a.virt[t][1]-1]
                                  )*a.weights[t]
                    else:
                        raise IndexError("Currently no more than double \
amplitudes are supported!")
                    sum_weights[i] += a.weights[t]
                    sum_redweights[i] += rw
            #----------------
            # Export as
            fout = fname + "." + fmt
            d = {"State": [i+1 for i in range(n_states)],
                 "sum_weight" : sum_weights,
                 "sum_red_weight" : sum_redweights}
            df = pd.DataFrame(d)
            df = df.assign(diff=df["sum_weight"]-df["sum_red_weight"],
                           ratio=df["sum_red_weight"]/df["sum_weight"])
            if fmt == "print":
                print("State | Sum(W) | Sum(P) | Sum(W) - Sum(P) | ratio P/W |\n",50*"-")
                for i in range(n_states):
                    print("  S{0:>2d} |  {1:.3f} |  {2:.3f} | {3:15.3f} | {4:.1%}".format(
                            i+1, sum_weights[i], sum_redweights[i], sum_weights[i] -
                            sum_redweights[i], sum_redweights[i]/sum_weights[i]))
            elif fmt == "dataframe":
                return df
            elif fmt == "csv":
                df.to_csv(fout, encoding="utf-8")
            elif fmt == "xlsx" or fmt == "xls":
                writer = pd.ExcelWriter(fout)
                df.to_excel(writer, "Sheet1")
                writer.save()
            else:
                raise ValueError("Output format not supported!")

    def MO_Molden(self, results, atom_basis, fname="molecular_orbitals",
                  tmp_5d=True):
        """ Writes molecular orbitals to a molden file.
        
        Expects molecular geometry in Angstrom.
        More information on the molden format at
        http://www.cmbi.ru.nl/molden/molden_format.html
        
        Parameters
        ----------
        results : CCParser.ParseContainer
            Container object which holds MO coefficients.
        exponents : dict
            Dictionary mapping GTO exponents/coefficients to atoms. Expected
            format of dictionary entry is list of strings.
        fname : string
            Output file name.
        """
        from .QCBase import PeriodicTable
        import re
        
        C = results.C.get_last()
        xyz = results.xyz.get_last()
        en = results.mo_energies.get_last()
        PeTa = PeriodicTable()
        #TODO: Permutator needed in case of different formats (Molcas, Gaussian)
        
        with open(fname+".molden", "w") as out:
            out.write("[Molden Format]\n")
            # write XYZ
            out.write("[Atoms] (Angs)\n")
            for i,atom in enumerate(xyz):
                num = PeTa.get_atomic_num(atom[0])
                out.write("{0:>3}{1:7d}{2:5d}".format(atom[0], i+1, num))
                out.write("".join("{0:16.8f}".format(c) for c in atom[1:])+"\n")
            # write basis exponents
            out.write("[GTO]\n")
            for n in range(len(xyz)):
                # atom sequence number, 0
                out.write("{0:d}{1:5d}\n".format(n+1, 0))
                symb = xyz[n][0].upper()
                #a = atom.upper()
                basis = atom_basis[symb]
                for coeff in basis:
                    # shell label, number of primitives, 1.00
                    if re.search(r"[SDPF]", coeff[0]):
                        out.write("{0:}{1:6d}{2:12.6f}\n".format(
                                coeff[0], int(coeff[1]), float(coeff[2])))
                    # exponent, contraction coefficient
                    else:
                        out.write("{0:18.8e}{1:18.8e}\n".format(
                                float(coeff[0]), float(coeff[1])))
                out.write("\n")
            for imo in range(C.shape[0]):#assumes counting from MO 1 !!
                out.write("[MO]\nSym=X\n")
                if imo < en.n_occ:#occupied
                    out.write("Ene={0:12.6f}\n".format(en.occ[imo]))
                    out.write("Spin=alpha\n")
                    out.write("Occup=1\n")
                else:#virtual
                    out.write("Ene={0:12.6f}\n".format(en.virt[imo]))
                    out.write("Spin=alpha\n")
                    out.write("Occup=0\n")
                for i in range(C.shape[1]):
                    out.write("{0:6d}{1: 22.12e}\n".format(i+1,C[imo, i]))
            if tmp_5d:
                out.write("[5D]\n")
            print("MOs written to Molden file.")
        