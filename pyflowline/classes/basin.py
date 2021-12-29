import os
from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
from pathlib import Path
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.wkt import loads
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib import cm

import cartopy.crs as ccrs


from pyflowline.algorithm.auxiliary.text_reader_string import text_reader_string

from pyflowline.formats.convert_shapefile_to_json import convert_shapefile_to_json
from pyflowline.formats.export_flowline import export_flowline_to_json

desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)
desired_proj = ccrs.PlateCarree()

class BasinClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
        return JSONEncoder.default(self, obj)


class pybasin(object):
    lBasinID =1 
    sBasin=''
    lCellID_outlet=-1
    iFlag_disconnected =0
    iFlag_dam=0
    dLongitude_outlet_degree = -9999.
    dLatitude_outlet_degree = -9999.
    dAccumulation_threshold= 100000.0
    dThreshold_small_river = 10000
    dLength_flowline_total = 0.0
    sWorkspace_output_basin=''
    sFilename_flowline_raw=''    
    sFilename_flowline_filter=''
    sFilename_flowline_filter_json=''
    sFilename_dam=''
    sFilename_flowline_topo=''
    #before intersect
    sFilename_flowline_segment_order_before_intersect=''
    sFilename_flowline_segment_index_before_intersect=''
    sFilename_flowline_final=''

    aFlowline_basin=None
    def __init__(self, aParameter):

        if 'lBasinID' in aParameter:            
            self.lBasinID             = int(aParameter['lBasinID'])
        else:
            self.lBasinID   = 1
        
        self.sBasin =  "{:03d}".format(self.lBasinID )
        
        if 'lCellID_outlet' in aParameter:            
            self.lCellID_outlet             = int(aParameter['lCellID_outlet'])
        else:
            self.lCellID_outlet   = -1

        if 'iFlag_disconnected' in aParameter:            
            self.iFlag_disconnected             = int(aParameter['iFlag_disconnected'])
        else:
            self.iFlag_disconnected   = 0
        
        if 'iFlag_dam' in aParameter:            
            self.iFlag_dam             = int(aParameter['iFlag_dam'])
        else:
            self.iFlag_dam   = 0
        
        if 'dLongitude_outlet_degree' in aParameter:            
            self.dLongitude_outlet_degree             = float(aParameter['dLongitude_outlet_degree'])
        else:
            self.dLongitude_outlet_degree   = -9999.
        
        if 'dLatitude_outlet_degree' in aParameter:            
            self.dLatitude_outlet_degree             = float(aParameter['dLatitude_outlet_degree'])
        else:
            self.dLatitude_outlet_degree   = -9999.
        
        if 'dThreshold_small_river' in aParameter:            
            self.dThreshold_small_river             = float(aParameter['dThreshold_small_river'])
        else:
            self.dThreshold_small_river   =10000.0

        if 'dAccumulation_threshold' in aParameter:            
            self.dAccumulation_threshold             = float(aParameter['dAccumulation_threshold'])
        else:
            self.dAccumulation_threshold = 100000.0

        if 'sFilename_flowline_raw' in aParameter:
            self.sFilename_flowline_raw = aParameter['sFilename_flowline_raw']
        else:
            self.sFilename_flowline_raw   =''
       
        if 'sFilename_flowline_filter' in aParameter:
            self.sFilename_flowline_filter = aParameter['sFilename_flowline_filter']
        else:
            self.sFilename_flowline_filter   = ''

        if 'sWorkspace_output_basin' in aParameter:
            self.sWorkspace_output_basin = aParameter['sWorkspace_output_basin']
        else:
            self.sWorkspace_output_basin   = ''
            print('The basin output path is not specified!')


        self.sFilename_flowline_filter_json = os.path.join(str(self.sWorkspace_output_basin ), "flowline_filter.json"  )

        if 'sFilename_dam' in aParameter:
            self.sFilename_dam = aParameter['sFilename_dam']
        else:
            self.sFilename_dam   = ''

        if 'sFilename_flowline_topo' in aParameter:
            self.sFilename_flowline_topo = aParameter['sFilename_flowline_topo']
        else:
            self.sFilename_flowline_topo   =''

        sBasinID = "{:03d}".format(self.lBasinID)

        self.sFilename_flowline_segment_index_before_intersect = 'flowline_segment_index_before_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_segment_order_before_intersect = 'flowline_segment_order_before_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_intersect  = 'flowline_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_final = 'flowline_final_' + sBasinID + '.json'
        
        return
    
    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)
        return sJson
    
    def convert_flowline_to_json(self):
        sFilename_raw = self.sFilename_flowline_filter            
        sFilename_out = self.sFilename_flowline_filter_json
        convert_shapefile_to_json(1, sFilename_raw, sFilename_out)
    
    def export_flowline(self, aFlowline_in, sFilename_json_in,iFlag_projected_in = None,  pSpatial_reference_in = None):

        export_flowline_to_json(aFlowline_in, iFlag_projected_in= iFlag_projected_in, \
            pSpatial_reference_in = pSpatial_reference_in)



    def calculate_flowline_length(self, aFlowline_in):

        dLength = 0.0

        nflowline = len(aFlowline_in)

        for i in range(nflowline):

            pFlowline= aFlowline_in[i]

            pFlowline.calculate_length()

            dLength = dLength + pFlowline.dLength

        self.dLength_flowline_total = dLength

        return dLength
    
    def plot(self, sVariable_in=None):

        if sVariable_in is not None:
            if sVariable_in == 'flowline_filter_json':
                sFilename_json = self.sFilename_flowline_filter_json
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_out = self.sFilename_flowline_segment_index_before_intersect
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                else:
                    sFilename_out = self.sFilename_flowline_final
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                pass
        else:
            #default 
            sFilename_json = self.sFilename_flowline_filter_json
        
        fig = plt.figure( dpi=150 )
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=desired_proj  )
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
    
        lID = 0
        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180
        n_colors = pLayer.GetFeatureCount()
        
        colours = cm.rainbow(np.linspace(0, 1, n_colors))
        for pFeature_shapefile in pLayer:
            pGeometry_in = pFeature_shapefile.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)
                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon
                    
                    if dLon < dLon_min:
                        dLon_min = dLon
                    
                    if dLat > dLat_max:
                        dLat_max = dLat
    
                    if dLat < dLat_min:
                        dLat_min = dLat
    
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= colours[lID])
                lID = lID + 1
                
    
        pDataset = pLayer = pFeature  = None      
    
        ax.set_extent([dLon_min  , dLon_max , dLat_min , dLat_max ])
        
        sDirname = os.path.dirname(sFilename_json)
        sFilename  = Path(sFilename_json).stem + '.png'
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')
        plt.show()
    
        return

    def preprocess_flowline(self):

        try:
            sFilename_flowline_filter = self.sFilename_flowline_filter
            sFilename_flowline_filter_json = self.sFilename_flowline_filter_json
            aFlowline_basin, pSpatial_reference = read_flowline_geojson( sFilename_flowline_filter_json )   
            sWorkspace_output_basin = self.sWorkspace_output_basin
            if self.iFlag_dam ==1:
                sFilename_dam = self.sFilename_dam
                aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1,cDelimiter_in=',' )
                sFilename_flowline_topo = self.sFilename_flowline_topo
                aData_flowline_topo = text_reader_string(sFilename_flowline_topo, iSkipline_in =1,cDelimiter_in=',' )
                aFromFlowline = aData_flowline_topo[:,1].astype(int).ravel()
                aToFlowline = aData_flowline_topo[:,2].astype(int).ravel()
                sFilename_flowline_raw = self.sFilename_flowline_raw
                aNHDPlusID_filter = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_filter)
                aNHDPlusID_raw = read_nhdplus_flowline_shapefile_attribute(sFilename_flowline_raw)
                ndam = len(aData_dam)
                aNHDPlusID_dams_headwater = list()
                aNHDPlusID_dams_nonheadwater = list()
                for j in range(0, ndam):
                    #print(j)
                    dLon = float(aData_dam[j][1])
                    dLat = float(aData_dam[j][0])
                    sDam = aData_dam[j][4]            
                    lNHDPlusID = int(aData_dam[j][5])
                    aNHDPlusID_dams_headwater.append(lNHDPlusID)

                    if lNHDPlusID in aNHDPlusID_filter:
                        #remove by id
                        for k in range(len(aFlowline_basin)):
                            if aFlowline_basin[k].lNHDPlusID == lNHDPlusID:
                                aFlowline_basin.pop(k)
                                break
                        pass
                    else:                                
                        aNHDPlusID_dam_nonheadwater = track_nhdplus_flowline(aNHDPlusID_filter, aFromFlowline, aToFlowline, lNHDPlusID)
                        aNHDPlusID_filter = aNHDPlusID_filter + aNHDPlusID_dams_headwater+ aNHDPlusID_dam_nonheadwater  
                        aNHDPlusID_dams_nonheadwater = aNHDPlusID_dams_nonheadwater + aNHDPlusID_dam_nonheadwater

                aFlowline_dams_headwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_headwater )
                for i in range(len(aFlowline_dams_headwater)):
                    aFlowline_dams_headwater[i].iFlag_dam = 1

                aFlowline_dams_nonheadwater = extract_nhdplus_flowline_shapefile_by_attribute(sFilename_flowline_raw, aNHDPlusID_dams_nonheadwater )
                aFlowline_basin = aFlowline_basin + aFlowline_dams_headwater + aFlowline_dams_nonheadwater
            else:
                pass

            if iFlag_disconnected == 1:                
                #aThreshold = np.full(2, 300.0, dtype=float)
                #aFlowline_basin = connect_disconnect_flowline(aFlowline_basin, aVertex, aThreshold)
                #sFilename_out = 'flowline_connect.json'
                #sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)    
                #export_flowline_to_json(iFlag_projected, aFlowline_basin,pSpatial_reference_gcs, sFilename_out)
                pass
            else:
                pass

            sFilename_out = 'flowline_before_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)            
            export_flowline(aFlowline_basin, sFilename_out)

            #calculate length
            calculate_flowline_length(aFlowline_basin)
            print(self.dLength_flowline_total)

            aVertex = find_flowline_vertex(aFlowline_basin)
            sFilename_out = 'flowline_vertex_without_confluence_before_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_vertex_to_json(iFlag_projected, aVertex,pSpatial_reference_gcs, sFilename_out)

            aFlowline_basin = split_flowline(aFlowline_basin, aVertex)
            sFilename_out = 'flowline_split_by_point_before_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, aFlowline_basin,pSpatial_reference_gcs, sFilename_out)

            #ues location to find outlet

            point= dict()   
            point['lon'] = pBasin.dLongitude_outlet_degree
            point['lat'] = pBasin.dLatitude_outlet_degree
            pVertex_outlet=pyvertex(point)

            aFlowline_basin = correct_flowline_direction(aFlowline_basin,  pVertex_outlet )

            pVertex_outlet = aFlowline_basin[0].pVertex_end

            sFilename_out = 'flowline_direction_before_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatial_reference_gcs, sFilename_out)

            #step 4: remove loops

            aFlowline_basin = remove_flowline_loop(aFlowline_basin)    
            sFilename_out = 'flowline_loop_before_intersect_' + sBasin + '.json'
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, aFlowline_basin,pSpatial_reference_gcs, sFilename_out)

            #using loop to remove small river, here we use 5 steps

            for i in range(3):
                sStep = "{:02d}".format(i+1)
                aFlowline_basin = remove_small_river(aFlowline_basin, dThreshold)
                sFilename_out = 'flowline_large_'+ sStep +'_before_intersect_' + sBasin + '.json'
                sFilename_out =os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatial_reference_gcs, sFilename_out)


                aVertex, lIndex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence, aConnectivity = find_flowline_confluence(aFlowline_basin,  pVertex_outlet)
                sFilename_out = 'flowline_vertex_with_confluence_'+ sStep +'_before_intersect_' + sBasin + '.json'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_vertex_to_json(iFlag_projected, aVertex, pSpatial_reference_gcs, sFilename_out, aAttribute_data=aConnectivity)

                aFlowline_basin = merge_flowline( aFlowline_basin,aVertex, pVertex_outlet, aIndex_headwater,aIndex_middle, aIndex_confluence  )  
                sFilename_out = 'flowline_merge_'+ sStep +'_before_intersect_' + sBasin + '.json'
                sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
                export_flowline_to_json(iFlag_projected, aFlowline_basin, pSpatial_reference_gcs, sFilename_out)

                if len(aFlowline_basin) ==1:
                    break

            
            dLength_total_new = calculate_flowline_length(aFlowline_basin)
            print(dLength_total_new)

            #build segment index
            aFlowline_basin, aStream_segment = define_stream_segment_index(aFlowline_basin)
            sFilename_out = pBasin.sFilename_flowline_segment_index_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, \
                aFlowline_basin, pSpatial_reference_gcs, sFilename_out, \
                aAttribute_data=[aStream_segment], aAttribute_field=['iseg'], aAttribute_dtype=['int'])

            #build stream order 
            aFlowline_basin, aStream_order = define_stream_order(aFlowline_basin)
            sFilename_out = pBasin.sFilename_flowline_segment_order_before_intersect
            sFilename_out = os.path.join(sWorkspace_output_basin, sFilename_out)
            export_flowline_to_json(iFlag_projected, \
                aFlowline_basin, pSpatial_reference_gcs, sFilename_out, \
                aAttribute_data=[aStream_segment, aStream_order], aAttribute_field=['iseg','iord'], aAttribute_dtype=['int','int'])
        except:
            print('Flowline preprocess failed')

        return aFlowline_basin