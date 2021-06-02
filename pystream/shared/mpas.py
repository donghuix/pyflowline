from abc import ABCMeta, abstractmethod
from osgeo import gdal, osr, ogr
from pystream.shared.vertex import pyvertex
from pystream.shared.edge import pyedge
from pystream.shared.cell import pycell
class pympas(pycell):
    dLength=0.0
    dArea=0.0
    aEdge=None
    aVertex=None
    aFlowline=None

    def __init__(self, aEdge):    
        pass
    
   
    def calculate_cell_area(self, aEdge):

        pass