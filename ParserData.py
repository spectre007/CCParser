from . import constants as c
import numpy as np

class Struct(object):
    """ Struct-like container object """
    def __init__(self, **kwds): # keyword args define attribute names and values
        self.__dict__.update(**kwds)
        
class ParseContainer(object):
    """ Generic container object which keeps track of the parsed quantities.
    It allows to parse the same quantity several times.
    There should be one instance/parsed quantity. """
    def __init__(self):
        self.nversion = 0
        self.data = []
        self.lines = []
    
    def add(self, hook_line, new_obj):
        #self.data[hook_line] = new_pvalue
        self.data.append(new_obj)
        self.lines.append(hook_line)
        self.nversion += 1
    
    def get_first(self):
        idx = self.lines.index(min(self.lines))#not needed if assuming ordered parsing (line by line)
        #return self.data[0]
        return self.data[idx]
    
    def get_last(self):
        idx = self.lines.index(max(self.lines))#not needed if assuming ordered parsing (line by line)
        # return self.data[-1]
        return self.data[idx]
    
    def get_data(self):
        return self.data
    
    def __len__(self):
        assert len(self.data) == len(self.lines)
        return len(self.data)
    
    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return self.data[idx.start : idx.stop : idx.step]
        else:
            if idx >= len(self.data) or abs(idx) > len(self.data):
                raise IndexError("ParseContainer: Index out of range")
            return self.data[idx]
        
    def __setitem__(self, idx, value):
        """ Setter method which expects value tuple (line, parsed_obj) """
        self.lines[idx] = value[0]
        self.data[idx] = value[1]
        
    def __delitem__(self, idx):
        self.data.remove(idx)
        self.lines.remove(idx)
        
    def __iter__(self):
#        self.n = 0
        return iter(self.data)
        
#    def __next__(self):
#        if self.n <= self.nversion:
#            return self.data[self.n]
#        else:
#            raise StopIteration
    
    def __contains__(self, line):
        if type(line) == str:
            line = int(line)
        return True if line in self.lines else False
    
    def __str__(self):
        s = "\n"
        s+= "Line" + 3*" " + "Parsed Value\n"
        for i in range(self.nversion):
            s+= str(self.lines[i]) + 3*" " + str(self.data[i]) + "\n"
        return s
    
    
        
        
class MolecularOrbitals(object):# TODO: change name? "MolecularOrbitalEnergies"
    """ General molecular orbital class, which has more functionality than simple arrays """
    def __init__(self, o, v):
        self.occ = list(map(float,o))
        self.virt = list(map(float,v))
        self.n_occ = len(o)
        self.n_virt = len(v)
        self.homo = max(self.occ)
        self.lumo = min(self.virt)
        self.n_orb = self.n_occ + self.n_virt
        
    
    def __str__(self):
        n1 = [i for i in range(1,self.n_occ+1)]
        n2 = [i for i in range(self.n_occ +1, self.n_orb+1)]
        s = "\n"
        for i in range(0, len(self.occ), 10):
            s += 4*" " + " ".join("{:>8}".format(j) for j in n1[i:i+10])+"\n"
            if i == 0:
                s += " occ: " +' '.join("{:8.3f}".format(j) for j in self.occ[i:i+10])+"\n"
            else:
                s += 6*" "+' '.join("{:8.3f}".format(j) for j in self.occ[i:i+10])+"\n"
        s += 7*" "+88*"-"+"\n"
        for i in range(0, len(self.virt), 10):
            s += 4*" " + " ".join("{:>8}".format(j) for j in n2[i:i+10])+"\n"
            if i == 0:
                s += " virt:" +' '.join("{:8.3f}".format(j) for j in self.virt[i:i+10])+"\n"
            else:
                s += 6*" "+' '.join("{:8.3f}".format(j) for j in self.virt[i:i+10])+"\n"
        return s
    
    def RVS(self, gap):
        """ Determine amount of virtual orbitals to freeze based on energy gap (in eV) """
        if gap <= 0:
            raise ValueError("Negative or Zero energy gap not allowed for restriction of virtual space.")
        else:
            thr = gap/c.Hartree2eV + self.homo
#            print("THR: ",thr)
#            print("N_VIRT: ", self.n_virt)
            idx = min(range(len(self.virt)), key=lambda i: abs(self.virt[i]-thr))
            freeze = self.n_virt - (idx +1)
            part_of_v = float(freeze)/float(self.n_virt)
            s = "Index: {0:3d}, Number of frozen virtuals: {1:3d}, ({2:.1%})".format(idx, freeze, part_of_v)
            print(s)
            
class Amplitudes(object):
    """ General container for amplitudes of one state for easier access to and export of amplitude data """
    def __init__(self, occ, virt, v, factor=1.0):
        self.occ = occ # list of lists, even if only single int in sublist
        self.virt= virt
        self.v = v
        self.factor = factor
        self.weights = list(map(lambda x: self.factor * x**2, self.v))
    
    def __str__(self):
        s = "Amplitudes: Weights > 5%\n"
        for i in range(len(self.occ)):
            if self.weights[i] > 0.05:
                if len(self.occ[i]) == 1:
                    s += "{0:>4} -> {1:>4} : {2:.1f}\n".format(self.occ[i][0],
                          self.virt[i][0], 100*self.weights[i])
                elif len(self.occ[i]) == 2:
                    s += "{0:>4}, {1:>4} -> {2:>4}, {3:>4} : {4:.1f}\n".format(
                            self.occ[i][0], self.occ[i][1], self.virt[i][0],
                            self.virt[i][1], 100*self.weights[i])
        return s
        
    @classmethod
    def from_list(cls, allinone, factor=1.0):
        """ Alternative constructor which expects single list.
        Format: [[occ_i, occ_j,..., virt_a, virt_b,..., ampl], ...] """
        occ, virt, v = [], [], []
        for transition in allinone:
            assert(len(transition) % 2 != 0)
            n_mo = int((len(transition)-1)/2)
            occ.append(transition[0:n_mo])# slices yield list, even if only one element
            virt.append(transition[n_mo:-1])
            v.append(transition[-1])# index yields float
        return cls(occ, virt, v, factor)
        
    def to_dataframe(self, thresh=0.05):
        """ Converts the amplitude data to handy pandas.DataFrame object """
        try:
            import pandas as pd
        except ImportError:
            raise ImportError("Module 'pandas' needed for 'Amplitudes.to_dataframe()' ")
        # TODO: improve this clunky part
        max_exc = max(list(map(len,self.occ)))
        occ, virt = [], []
        for i in range(len(self.occ)):
            occ.append(self.occ[i] + [0]*(max_exc - len(self.occ[i])))
            virt.append(self.virt[i] + [0]*(max_exc - len(self.virt[i])))
        idx_o = list(map(lambda x: "occ_"+str(x), [n for n in range(1,max_exc+1)]))
        idx_v = list(map(lambda x: "virt_"+str(x), [n for n in range(1,max_exc+1)]))
        df = pd.concat([pd.DataFrame(occ, columns=idx_o),
                        pd.DataFrame(virt, columns=idx_v),
                        pd.Series(self.weights, name="weight")], axis=1)
        return df[(df["weight"] > thresh)] # trim DataFrame based on awesome function
    
    def get_single_list(self):
        """ Return single list of amplitude data in the format:
            [[occ_i, occ_j,..., virt_a, virt_b,..., ampl], ...] """
        return []

