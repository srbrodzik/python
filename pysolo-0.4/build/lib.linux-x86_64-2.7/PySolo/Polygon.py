import matplotlib.patches as patches 
from matplotlib.pyplot import figure, show, rc, grid
import matplotlib.pyplot as plt

class Polygon:

    def __init__(self, verts, color):
        self.verts = verts
        self.color = color
        
    def draw(self,ax):
        poly = patches.Polygon(self.verts,facecolor=self.color,closed=True,fill=True,edgecolor=self.color) 
        ax.add_patch(poly)

