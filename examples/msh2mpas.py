import numpy as np
from scipy.spatial import Voronoi
import fiona
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import pyproj

def msh2mpas(mesh, bnd_file, debug):

    radius_earth_m = 6371220

    coord = mesh.point["coord"]
    tria3 = mesh.tria3["index"]
    print('Generating Voronoi...')
    voronoi = Voronoi(coord)

    c = fiona.open(bnd_file)
    boundary = [np.array(poly['geometry']['coordinates']) for poly in c.values()]
    bndPoly = Polygon(boundary[0][0])

    if debug:
        import matplotlib.pyplot as plt
        plt.plot(boundary[0][0][:,0],boundary[0][0][:,1],'k-',linewidth=1)
        plt.triplot(coord[:,0],coord[:,1],tria3,linewidth=0.1)

    proj_wgs84 = pyproj.CRS('EPSG:4326')
    proj_utm10n = pyproj.CRS('EPSG:32610')
    project = pyproj.Transformer.from_crs(proj_wgs84, proj_utm10n, always_xy=True).transform
    bndProj = Polygon(project(x, y) for x, y in boundary[0][0])

    xc = []; yc = []; xv = []; yv = []; frac = []; area = []
    for idx, point_region in enumerate(voronoi.point_region):
        region = np.array(voronoi.regions[point_region])
        if np.all(region> -1) and len(region) > 0:
            pt = Point(voronoi.points[idx,0],voronoi.points[idx,1])
            if np.min(bndPoly.distance(pt)) < 1e-5:
                region = np.append(region,region[0])
                cell = voronoi.vertices[np.array(region),:]
                cell =  np.vstack([cell, cell[0]])
                cell_utm = Polygon(project(x, y) for x, y in cell)
                cell_int = cell_utm.intersection(bndProj)
                area_int = cell_int.area
                aream2   = cell_utm.area
                area.append(aream2/(radius_earth_m**2))
                ff = area_int/aream2
                if abs(ff - 1) < 1e-5:
                    ff = 1.0
                frac.append(ff)
                if debug:
                    print("frac is " + str(frac[-1]) + ", area is " + str(aream2) + " [m^2], " + str(area[-1]) + " [radians]...")
                    plt.plot(voronoi.vertices[np.array(region),0], voronoi.vertices[np.array(region),1],'r-',linewidth=0.5)
                    plt.plot(voronoi.points[idx,0],voronoi.points[idx,1],'g.',markersize=0.5)
                xc.append(voronoi.points[idx,0])
                yc.append(voronoi.points[idx,0])

                a = np.empty((9,))
                a[:] = np.nan
                a[:len(region)] = voronoi.vertices[np.array(region),0]
                xv.append(a)
                a[:len(region)] = voronoi.vertices[np.array(region),1]
                yv.append(a)

    if debug:
        plt.savefig("debug.pdf") 
    
    return xc,yc,xv,yv,area,frac

