from . import constants as c

class Struct(object):
    """ Generic container object """
    def __init__(self, **kwds): # keyword args define attribute names and values
        self.__dict__.update(**kwds)
        
class MolecularOrbitals(object):
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
            
        
#if __name__ == "__main__":
#    test = Struct(keyy="value1")