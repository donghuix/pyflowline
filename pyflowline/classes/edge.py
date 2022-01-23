from abc import ABCMeta, abstractmethod
import json
import numpy as np
from osgeo import gdal, osr, ogr
from pyflowline.algorithms.auxiliary.gdal_functions import  calculate_angle_betwen_vertex
from pyflowline.classes.vertex import pyvertex

class pyedge(object):
    __metaclass__ = ABCMeta 
    lEdgeID=-1
    lIndex=-1
    pVertex_start = None
    pVertex_end = None
    dLength=0.0
    
    lIndex_upstream=-1
    lIndex_downstream=-1

    def __init__(self, pVertex_start_in, pVertex_end_in):
        

        if pVertex_start_in == pVertex_end_in:
            print('The two vertices are the same')    
        else:
            self.pVertex_start = pVertex_start_in
            self.pVertex_end = pVertex_end_in

            #dict.__init__(self, pVertex_start_in=pVertex_start_in, pVertex_end_in=pVertex_end_in)
            self.dLength = self.calculate_length()

        return

    def calculate_length(self):
        dLength =0.0
        dLength = self.pVertex_start.calculate_distance( self.pVertex_end)
        self.dLength= dLength

        return dLength

    def check_shared_vertex(self, other):
        '''
        check whether two edges are sharing the same vertex
        '''
        iFlag_shared =-1

        v0 = self.pVertex_start
        v1 = self.pVertex_end

        v2 = other.pVertex_start
        v3 = other.pVertex_end

        if v0 == v2 or v0 == v2 or v1==v2 or v1==v3:
            iFlag_shared =1
        else:
            iFlag_shared=0

        return iFlag_shared
    
    def check_upstream(self, other):
        iFlag_upstream =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end

        v2 = other.pVertex_start
        v3 = other.pVertex_end

        if v0 == v3:
            iFlag_upstream =1
        else:
            iFlag_upstream=0

        return iFlag_upstream

    def check_downstream(self, other):
        iFlag_downstream =-1
        v0 = self.pVertex_start
        v1 = self.pVertex_end

        v2 = other.pVertex_start
        v3 = other.pVertex_end

        if v1 == v2:
            iFlag_downstream =1
        else:
            iFlag_downstream=0

        return iFlag_downstream

    def reverse(self):
        v0 = self.pVertex_start
        v1 = self.pVertex_end

        self.pVertex_start = v1
        self.pVertex_end = v0
    
    def is_overlap(self, pEdge_in):
        iFlag_overlap = 0
        pVertex_start1 = self.pVertex_start
        pVertex_end1 = self.pVertex_end

        pVertex_start2 = pEdge_in.pVertex_start
        pVertex_end2 = pEdge_in.pVertex_end

        if pVertex_start1 == pVertex_start2 and pVertex_end1 == pVertex_end2:
            iFlag_overlap = 1
        else:
            if  pVertex_start1 == pVertex_end2 and pVertex_end1 == pVertex_start2:
                iFlag_overlap = 1
            else:
                iFlag_overlap = 0

        return iFlag_overlap

    def check_vertex_on_edge(self, pVertex_in):
        iFlag =0 

        pVertex_start = self.pVertex_start
        pVertex_end = self.pVertex_end
        self.dLength = pVertex_start.calculate_distance(pVertex_end)
        if pVertex_in != pVertex_start and pVertex_in!=pVertex_end:
            d1 = pVertex_start.calculate_distance(pVertex_in)
            d2 = pVertex_end.calculate_distance(pVertex_in)
            d3 = d1 +d2 -self.dLength

            angle3deg = calculate_angle_betwen_vertex(\
                 pVertex_start.dLongitude_degree, pVertex_start.dLatitude_degree,\
                 pVertex_in.dLongitude_degree, pVertex_in.dLatitude_degree,\
                 pVertex_end.dLongitude_degree,pVertex_end.dLatitude_degree)

        
            if  angle3deg > 175 and d3 < 0.1:
                iFlag = 1
            else:
                iFlag = 0
            #if (d1<self.dLength and d2<self.dLength and d3 < 0.1 ):#care
                
                
            #else:
            #    if (d3 < 0.1):
            #        print('debug')
            #    iFlag = 0 
        else:
                iFlag = 0 

        return iFlag 
    
    def __eq__(self, other):                
        iFlag_overlap = self.is_overlap(other)  
        return iFlag_overlap

    def __ne__(self, other):
        return not self.__eq__(other)

