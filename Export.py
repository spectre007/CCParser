class Export:
    """ Export class which writes parsed data to a certain format"""
    valid_formats = ["screen", "pdf", "excel", "txt", "csv"]
    
    def __init__(self, data, outformat=None, basename="CCP_Export"):
        if outformat == None:
            print(" [WARNING] No export format defined! Using default: screen print.")
            outformat = "screen"
        self.data = data

if __name__ == "__main__":
    e = Export([0,1], outformat="screen")