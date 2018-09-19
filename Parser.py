from .ParserData import Struct
from .ParserData import ParseContainer
from .QCBase import GenFormatter, VarNames as V
import importlib as il
import os.path, inspect
import re
import logging

class Parser(object):
    def __init__(self, output, *, software=None, toConsole=True,
                 toFile=False, logname="CCParser.log"):#cf. PEP-3102
        self.f_output = output
        self.logger = logging.getLogger("CCParser")
        self.toConsole = toConsole
        self.toFile = toFile
        self.logname = logname
        self.setupLogger()
        
        if software != None:
            self.software = software
        else:
            raise ValueError("No software specified!")

        self.output_basename = os.path.basename(output)
        self.read_output()# read output to memory
        self.load_methods()# software dependent import
        
        self.results = Struct()# Set up container
        self.logger.warning("CCParser starts...")
        for i,line in enumerate(self.rawData):
            for mthd in self.methods:
#                match, key = self.canParse(line, mthd)
                match, keys = self.canParse(line, mthd)
                if match:
                    for key in keys:# if not 1-to-1 mapping
                        q = self.get_quantity(i, key, mthd)
                        if hasattr(self.results, mthd.map[key]):
                            obj = getattr(self.results, mthd.map[key])
                            obj.add(i, q)
                        else:
                            obj = ParseContainer()
                            obj.add(i, q)
                            setattr(self.results, mthd.map[key], obj)
        if not hasattr(self.results, V.has_finished):
            container = ParseContainer()
            container.add(0, False)
            setattr(self.results, V.has_finished, container)
            self.logger.warning("Output indicates abnormal exit. Added "+
                                "[results.has_finished] = False")
        self.logger.warning("CCParser has finished.")
        self.loggerCleanUp()


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
        keys = []#for cases where there's no 1-to-1 mapping
        for key, value in mthd.hooks.items():
            if value in line:
                found = True
                keys.append(key)
#                return found, key
            else:
                match = re.search(value, line)
                if match:
                    found = True
                    keys.append(key)
#                    return found, key
        if found == False:
            return found, None
        else:
            return found, keys

            
    def get_quantity(self, i, key, mthd):
        """ Call function of method class. This is the actual parsing. """
        method_func = getattr(mthd, key)# needs to be method not list of methods
        result = method_func(i, self.rawData)
        return result
    
    def load_methods(self):
        """ Load correct module which contains parsing information
        based on which software was specified. """
        tmp = re.sub('[^A-Za-z]+', '', self.software.lower())
        if tmp == "qchem":
            m_package = ".QChem"
        elif tmp == "gaussian":
            raise NotImplementedError("Gaussian parsingnot implemented yet!")
            m_package = ".Gaussian"
        elif tmp == "molcas":
            raise NotImplementedError("Molcas parsing not implemented yet!")
            m_package = ".Molcas"
        elif tmp == "turbomole":
            raise NotImplementedError("Turbomole parsing not implemented yet!")
            m_package = ".Turbomole"
        elif tmp == "psi":
            m_package = ".Psi4"
        else:
            raise ValueError("The specified software is misspelled or not implemented yet!")
        global m
#        m = il.import_module(m_package+".methods",package="CCParser")
        m = il.import_module(m_package,package="CCParser")
        self.method_names = [k[0] for k in inspect.getmembers(m,
                             inspect.isclass) if k[1].__module__ == "CCParser"+m_package]
        self.methods = [getattr(m,mname)() for mname in self.method_names]#this also instantiates!!

    def setupLogger(self):
        # Set main logger's minimum output level
        self.logger.setLevel(logging.INFO)
        # Set up Formatter
#        p_fmt = logging.Formatter("[results.%(Parsed)s] Parsed %(message)s")
        # TODO: change number of loggers
        # This is abusing the Formatter class a bit, but I wanted to avoid
        # one Logger for every format, maybe I'll change this in the future.
        p_fmt = GenFormatter(
                {logging.INFO: "[results.%(Parsed)s] Parsed %(message)s",
                 logging.WARNING: "==[%(asctime)s]== %(message)s",
                 logging.ERROR: "%(message)s"})
        # Set up Handlers
        if self.toFile:
            fh = logging.FileHandler(self.logname)
            fh.setLevel(logging.INFO)
            fh.setFormatter(p_fmt)
            self.logger.addHandler(fh)
        if self.toConsole:
            ch = logging.StreamHandler()
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(p_fmt)
            self.logger.addHandler(ch)
        # No output in case both booleans are False
        if not any([self.toConsole, self.toFile]):
            self.logger.setLevel(logging.CRITICAL)
            
    def loggerCleanUp(self):
        """In order to avoid multiplying handlers. """
        for i in range(len(self.logger.handlers)):
            self.logger.handlers.pop()
        
        
        
    
