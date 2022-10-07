from Sweep import Sweep

class Volume:

    def __init__(self, sweeps, ranges):
        
        self.sweeps = sweeps
        self.ranges = ranges
    
    def getMetadata(self):
        self.dimensions = self.ncfile.dimensions.keys() # list
        self.variableNames = self.ncfile.variables.keys()   # list
        self.globalAttList = self.ncfile.ncattrs() #dir(self.ncfile)
