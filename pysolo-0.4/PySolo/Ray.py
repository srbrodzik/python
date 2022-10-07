class Ray:
    """ store a Ray of radar/lidar data """
    def __init__(self, az, el, time, varDict):
        self.az = az
        self.el = el
        self.time = time
        self.varDict = varDict
