import sys
import numpy as np
import osgeo
import math
from math import radians, cos, sin, asin, sqrt
from osgeo import ogr, osr, gdal, gdalconst
from numpy import arctan2, cos, sin, sqrt, pi, power, append, diff

#https://stackoverflow.com/questions/8204998/how-to-check-if-a-pointlonc-latc-lie-on-a-great-circle-running-from-lona-lata
def calculate_distance_to_plane(x1, y1, x2, y2, x3, y3):
    
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) )) #this is the middle one
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)
    #The formula is x+b*y+c*z=0 
    x1,y1,z1 = a3
    x2,y2,z2 = b3
    x3,y3,z3 = c3
    c = (-x1*y3 + x3* y1)/( z1*y3 - z3*y1 )
    b = (  -x1*z3 + x3 * z1 ) / (y1 * z3 - y3*z1)

    distance = abs(  x2 + b * y2 + c * z2 )

    return distance

def calculate_angle_betwen_vertex_normal(x1, y1, x2, y2, x3, y3):
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) ))
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)

    a3vec = a3 - b3
    c3vec = c3 - b3 

    dot=np.dot(a3vec, c3vec)
    g = np.cross(a3vec, c3vec)
    det = np.dot(b3, g)
    angle = np.arctan2(det, dot)
    f = np.degrees(angle) 

    if f < 0:
        f = 360 + f
    
    return f

def calculate_angle_betwen_vertex(x1, y1, x2, y2, x3, y3):
    #all in degree
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) ))
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)
    # Vectors in 3D space
    a3vec = a3 - b3
    c3vec = c3 - b3        
    angle3deg = angle_between_vectors_degrees(a3vec, c3vec)

    return  angle3deg
    
def longlat_to_3d(lonr, latr):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""
    a = math.cos(latr) * math.cos(lonr)
    b = math.cos(latr) * math.sin(lonr)
    c = math.sin(latr)
    return np.array((a,b,c))

def angle_between_vectors_degrees(u, v):
    """Return the angle between two vectors in any dimension space,
    in degrees."""

    a = np.dot(u, v)
    b = np.linalg.norm(u)
    c = np.linalg.norm(v)
    d = a / (b* c)
    if d > 1:
        d = 1
    if d < -1:
        d = -1
    e = math.acos(d)
    f = np.degrees(e)
    #if d < 0:
    #    f = f + 180.0    

    return f
    #return np.degrees(
    #    math.acos(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))))

def convert_360_to_180(dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """
 
    a = int(dLongitude_in /180)
    dLongitude_out = dLongitude_in - a*360.0

    return dLongitude_out



def convert_180_to_360(dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """

    dLongitude_out = (dLongitude_in + 360.0) % 360.0

    return dLongitude_out

def retrieve_geotiff_metadata(sFilename_geotiff_in):
    """[retrieve the metadata of a geotiff file]

    Args:
        sFilename_geotiff_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    pDriver = gdal.GetDriverByName('GTiff')
   
    pDataset = gdal.Open(sFilename_geotiff_in, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_geotiff_in)
        sys.exit("Try again!")
    else: 
        pProjection = pDataset.GetProjection()
        pSpatial_reference = osr.SpatialReference(wkt=pProjection)
    
    
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        #nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]       
        
        print( dPixelWidth, dOriginX, dOriginY, nrow, ncolumn)
        return dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatial_reference, pProjection, pGeotransform


def degree_to_meter(dLatitude_in, dResolution_degree_in):
    """[Conver a degree-based resolution to meter-based resolution]

    Args:
        dLatitude_in ([type]): [description]
        dResolution_degree_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    dRadius = 6378137.0  #unit: m earth radius
    dRadius2 = dRadius * np.cos( dLatitude_in / 180.0 * np.pi)
    dResolution_meter = dResolution_degree_in / 360.0 * (2*np.pi * dRadius2)

    return dResolution_meter

def meter_to_degree(dResolution_meter_in, dLatitude_mean_in):
    """[Convert a meter-based resolution to degree-based resolution]

    Args:
        dResolution_meter_in ([type]): [description]
        dLatitude_mean_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    dLatitude_mean = abs(dLatitude_mean_in)

    dRadius = 6378137.0 # #unit: m earth radius
    dRadius2 = dRadius * np.cos( dLatitude_mean / 180.0 * np.pi)

    ##dResolution_meter = dResolution_degree / 360.0 * 2*np.pi * dRadius2

    dResolution_degree= dResolution_meter_in/(2*np.pi * dRadius2) * 360.0

    return dResolution_degree



def reproject_coordinates(x_in, y_in, spatial_reference_source, spatial_reference_target=None):
    """[Reproject coordinates from one reference to another. By default to WGS84.]

    Args:
        x_in ([type]): [description]
        y_in ([type]): [description]
        spatial_reference_source ([type]): [description]
        spatial_reference_target ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """    
  

    if spatial_reference_target is not None:

        pass
    else:
        spatial_reference_target = osr.SpatialReference()
        spatial_reference_target.ImportFromEPSG(4326)
        
        pass

    
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    
        spatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    
    pTransform = osr.CoordinateTransformation( spatial_reference_source, spatial_reference_target)
   
    x_new,y_new, z = pTransform.TransformPoint( x_in,y_in)
    
    return x_new,y_new

def reproject_coordinates_batch(aX_in, aY_in, spatial_reference_source, spatial_reference_target=None):
    """[Reproject coordinates from one reference to another in batch mode. By default to WGS84.]

    Args:
        aX_in ([type]): [description]
        aY_in ([type]): [description]
        spatial_reference_source ([type]): [description]
        spatial_reference_target ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    #Reproject a list of x,y coordinates. 

    if spatial_reference_target is not None:

        pass
    else:
        spatial_reference_target = osr.SpatialReference()
        spatial_reference_target.ImportFromEPSG(4326)
        
        pass

    
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    
        spatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    
    pTransform = osr.CoordinateTransformation( spatial_reference_source, spatial_reference_target)

    npoint = len(aX_in)
    x_new=list()
    y_new=list()
    for i in range(npoint):
        x0 = aX_in[i]
        y0 = aY_in[i]
   
        x1,y1, z = pTransform.TransformPoint( x0,y0)

        x_new.append(x1)
        y_new.append(y1)
    
    return x_new,y_new


def calculate_distance_based_on_lon_lat(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6378137.0
    return c * r

def calculate_polygon_area(lons, lats,  algorithm = 0, radius = 6378137.0):
    """
    Computes area of spherical polygon, assuming spherical Earth. 
    Returns result in ratio of the sphere's area if the radius is specified. Otherwise, in the units of provided radius.
    lats and lons are in degrees.
    """
    #TODO: take into account geodesy (i.e. convert latitude to authalic sphere, use radius of authalic sphere instead of mean radius of spherical earth)
    lats = np.deg2rad(lats)
    lons = np.deg2rad(lons)

    if algorithm==0:
        # Line integral based on Green's Theorem, assumes spherical Earth
        

        #close polygon
        if lats[0]!=lats[-1]:
            lats = append(lats, lats[0])
            lons = append(lons, lons[0])

        # Get colatitude (a measure of surface distance as an angle)
        a = sin(lats/2)**2 + cos(lats)* sin(lons/2)**2
        colat = 2*arctan2( sqrt(a), sqrt(1-a) )

        #azimuth of each point in segment from the arbitrary origin
        az = arctan2(cos(lats) * sin(lons), sin(lats)) % (2*pi)

        # Calculate step sizes
        # daz = diff(az) % (2*pi)
        daz = diff(az)
        daz = (daz + pi) % (2 * pi) - pi

        # Determine average surface distance for each step
        deltas=diff(colat)/2
        colat=colat[0:-1]+deltas

        # Integral over azimuth is 1-cos(colatitudes)
        integrands = (1-cos(colat)) * daz

        # Integrate and save the answer as a fraction of the unit sphere.
        # Note that the sum of the integrands will include a factor of 4pi.
        area = abs(sum(integrands))/(4*pi) # Could be area of inside or outside

        area = min(area,1-area)
        if radius is not None: #return in units of radius
            return area * 4*pi*radius**2
        else: #return in ratio of sphere total area
            return area
    elif algorithm==2:
        #L'Huilier Theorem, assumes spherical earth
        #see:
        # https://mathworld.wolfram.com/SphericalPolygon.html
        # https://web.archive.org/web/20160324191929/http://forum.worldwindcentral.com/showthread.php?20724-A-method-to-compute-the-area-of-a-spherical-polygon
        # https://github.com/spacetelescope/spherical_geometry/blob/master/spherical_geometry/polygon.py
        # https://github.com/tylerjereddy/spherical-SA-docker-demo/blob/master/docker_build/demonstration.py
        #TODO
        pass
    elif algorithm==3:
        #https://trs.jpl.nasa.gov/handle/2014/41271
        #TODO
        pass