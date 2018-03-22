from .ParserData import Struct
from .ParserData import ParseContainer
#from .QChem import methods as m
import importlib as il
import os.path, inspect
import re

class Parser(object):
    def __init__(self, output, software=None):
        self.f_output = output
        self.software = software
        #self.output_folder 
        self.output_basename = os.path.basename(output)
        self.read_output()# read output to memory
        self.load_methods()# software dependent import
        
        self.results = Struct()# Set up container
#        self.methods = m.SCF()
#        print(self.methods.hooks.items())
        for i,line in enumerate(self.rawData):
            for mthd in self.methods:
                match, key = self.canParse(line, mthd)
                if match:
                    q = self.get_quantity(i, key, mthd)
                    #print(key, mthd.map)
                    if hasattr(self.results, mthd.map[key]):
                        obj = getattr(self.results, mthd.map[key])
#                        obj.add(i, self.get_quantity(i, key, mthd))
                        obj.add(i, q)
                    else:
                        obj = ParseContainer()
#                        obj.add(i, self.get_quantity(i, key, mthd))
                        obj.add(i, q)
                        setattr(self.results, mthd.map[key], obj)
    #                setattr(self.results, key, obj)
                    #print("Line "+str(i), obj.get_last())
                    
        
#        self.pData = Struct()

    def read_output(self):
        """ Read in output file """
        with open(self.f_output, "r") as f:
            self.rawData = f.readlines()
    
    def read_input(self, f_input):
        """ (Optional) Read input file """
        with open(f_input) as n:
            self.rawData.insert(0,n.readlines())
            
    def canParse(self, line, mthd):
        """ Check if line is parsable """
        found = False
#        for mthd in self.methods:
#        for key,value in self.methods.hooks.items():
        for key,value in mthd.hooks.items():
            if value in line:
                found = True
                return found, key
            else:
                match = re.search(value, line)
                if match:
                    found = True
                    return found, key
        if found == False:
            return found, None

            
    def get_quantity(self, i, key, mthd):
        """ Call function of method class. This is the actual parsing. """
        method_func = getattr(mthd, key)# needs to be method not list of methods
        result = method_func(i, self.rawData)
        return result
    
    def load_methods(self):
        """ Load correct module which contains parsing information
        based on which software was specified. """
        tmp = re.sub('[^A-Za-z0-9]+', '', self.software.lower())
        if tmp == "qchem":
            m_package = ".QChem"
        elif tmp == "gaussian":
            m_package = ".Gaussian"
        elif tmp == "molcas":
            m_package = ".Molcas"
        elif tmp == "turbomole":
            m_package = ".Turbomole"
        elif tmp == "psi":
            m_package = ".Psi4"
        else:
            raise Exception("The specified software is misspelled or not implemented yet!")
        global m
#        m = il.import_module(m_package+".methods",package="CCParser")
        m = il.import_module(m_package,package="CCParser")
        self.method_names = [k[0] for k in inspect.getmembers(m, inspect.isclass) if k[1].__module__ == "CCParser"+m_package]
        self.methods = [getattr(m,mname)() for mname in self.method_names]
        #print(self.methods)
        
        
    