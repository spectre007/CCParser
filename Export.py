import pandas as pd
from .QCBase import VarNames

class Exporter(object):
    """ Export class which writes parsed data to a certain format"""
    valid_formats = ["pdf", "xlsx", "txt", "csv"]
    
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
    
    def ExcitedStateSummary(self, results, fname="es_summary", fmt="csv"):
        """ Exports energy related excited state quantities to file """
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
        arrays = [[x for x in range(1,n_states+1)],[0 for x in range(n_states)]]
        tuples = list(zip(*arrays))# asterisk unpacks
        df1 = pd.DataFrame(d)
        df1.index = pd.MultiIndex.from_tuples(tuples, names=["State", "Row ID"])
        
        df = pd.concat([df1, amp_df], axis=1)
        
        # EXPORT TO FILE
        # TODO: write correct pd export
        fout = fname + "." + fmt
        if fmt == "csv":
            df.to_csv(fout, encoding="utf-8")
        elif fmt == ("xlsx" or "xls"):
            writer = pd.ExcelWriter(fout)
            df.to_excel(writer, "Sheet1")
            writer.save()
        
            
            
        
        
        
        
