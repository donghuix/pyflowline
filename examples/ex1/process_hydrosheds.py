import sys
sys.path.insert(0, '/Users/xudo627/Developments/reach')

def download_file(args):    
    import requests
    import zipfile

    filename = args.shp_file.split("/")[1] + ".zip"
    url = "https://data.hydrosheds.org/file/HydroRIVERS/" + filename
    
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                f.write(chunk)
    with zipfile.ZipFile(filename, 'r') as zip_ref:
        zip_ref.extractall("./dat")

    return filename

def clip_hydrosheds(args):
    import geopandas as gpd

    bnd = gpd.read_file(args.bnd_file)
    shp = gpd.read_file(args.shp_file)
    shp = shp.to_crs(bnd.crs)
    clp = gpd.overlay(shp, bnd, how="intersection")
    clp.to_file(args.clp_file)
    
    # Plot the clipped polylines
    #fig, ax = plt.subplots(figsize=(8, 6))
    #clp.plot(ax=ax, color='green', label='Clipped Polylines')
    #bnd.boundary.plot(ax=ax, color='red', linewidth=1, label='Polygon Boundary')
#
    #ax.set_title('Clipped Polylines')
    #ax.legend()
    #plt.show()

    return

def output_jgsw(args, rnet, hfun, tags, rsph=1.):
    """
    Output filtered reach network for mesh-gen. using jigsaw.

    """
    import numpy as np
    import jigsawpy  # here, so jigsaw is not a "strict" dep.
    import pyproj
    import fiona

    geom = jigsawpy.jigsaw_msh_t()
    
    geom.mshID = "euclidean-mesh"
    #geom.radii = \
    #    rsph * np.ones(3, dtype=geom.REALS_t)

    nset = []; eset = []; last = 0

    # add boundary
    if len(args.bnd_file) > 0:
        print("Adding boundary in the geometry")
        c = fiona.open(args.bnd_file)
        coords = [np.array(poly['geometry']['coordinates']) for poly in c.values()]

        nbox = 4
        temp = jigsawpy.jigsaw_msh_t()
        temp.vert2 = np.zeros(
            (nbox + 0), dtype=temp.VERT2_t)
        temp.edge2 = np.zeros(
            (nbox + 0), dtype=temp.EDGE2_t)
        
        temp.edge2["IDtag"][:] = 1
        
        indx = np.arange(0, nbox - 0) + last
        last = last + nbox

        temp.vert2["coord"][:, 0] = np.array([args.box_xmin, args.box_xmax, args.box_xmax, args.box_xmin])
        temp.vert2["coord"][:, 1] = np.array([args.box_ymin, args.box_ymin, args.box_ymax, args.box_ymax])

        temp.edge2["index"][:, 0] = indx + 0
        temp.edge2["index"][:, 1] = np.append(np.arange(1,nbox), 0)

        nset.append(temp.vert2)
        eset.append(temp.edge2)

        nbnd = len(coords[0][0])

        temp = jigsawpy.jigsaw_msh_t()
        temp.vert2 = np.zeros(
            (nbnd + 0), dtype=temp.VERT2_t)
        temp.edge2 = np.zeros(
            (nbnd + 0), dtype=temp.EDGE2_t)
        
        temp.edge2["IDtag"][:] = 2

        indx = np.arange(0, nbnd - 0) + last

        temp.vert2["coord"][:, 0] = coords[0][0][:,0]
        temp.vert2["coord"][:, 1] = coords[0][0][:,1]

        temp.edge2["index"][:, 0] = indx + 0
        temp.edge2["index"][:, 1] = np.append(np.arange(1,nbnd), 0) + last

        last = last + nbnd

        nset.append(temp.vert2)
        eset.append(temp.edge2)

    # add river network
    for rdat in rnet:
        if (rdat.flag >= +1):

            npts = len(rdat.vert)

            temp = jigsawpy.jigsaw_msh_t()
            temp.vert2 = np.zeros(
                (npts + 0), dtype=temp.VERT2_t)
            temp.edge2 = np.zeros(
                (npts - 1), dtype=temp.EDGE2_t)

            temp.edge2["IDtag"][:] = rdat.hpos

            indx = np.arange(0, npts - 1) + last

            last = last + npts

            temp.vert2["coord"][:, 0] = \
                [x.xlon for x in rdat.vert]

            temp.vert2["coord"][:, 1] = \
                [x.ylat for x in rdat.vert]

            temp.edge2["index"][:, 0] = indx + 0
            temp.edge2["index"][:, 1] = indx + 1

            nset.append(temp.vert2)
            eset.append(temp.edge2)
    
    geom.vert2 = np.concatenate(nset, axis=+0)
    geom.edge2 = np.concatenate(eset, axis=+0)

    #geom.vert2["coord"] *= np.pi / 180.

    __, fmap, rmap = np.unique(
        geom.vert2["coord"], 
        return_index=True, return_inverse=True, axis=0)

    geom.vert2 = geom.vert2[fmap]
    geom.edge2["index"] = rmap[geom.edge2["index"]]

    # Create boundary for the geometry
    # Outer boundary
    nbox = len(np.where(geom.edge2["IDtag"] == 1)[0])
    nbnd = len(np.where(geom.edge2["IDtag"] == 2)[0])
    print("nbox = " + str(nbox))
    print("nbnd = " + str(nbnd))
    geom.bound = np.zeros((nbox+nbnd), dtype=geom.BOUND_t)
    et = jigsawpy.\
        jigsaw_def_t.JIGSAW_EDGE2_TAG
    ibound = np.concatenate((np.ones(nbox),np.ones(nbnd)*2),axis=-1)
    geom.bound["IDtag"] = ibound
    geom.bound["index"] = np.arange(0,nbox+nbnd)
    geom.bound["cells"] = et

    jigsawpy.savemsh(tags + "_geom.msh", geom)
    jigsawpy.savevtk(tags + "_geo2.vtk", geom)

    init = jigsawpy.jigsaw_msh_t()

    init.mshID = "euclidean-mesh"
    init.radii = \
        rsph * np.ones(3, dtype=init.REALS_t)

    nset = []; eset = []; last = 0
    for rdat in rnet:
        if (rdat.flag >= 1 
                and rdat.vert[+0].seed >= 1):

            npts = +1
            temp = jigsawpy.jigsaw_msh_t()
            temp.vert2 = np.zeros(
                (npts + 0), dtype=temp.VERT2_t)

            temp.vert2["coord"][0, 0] = \
                rdat.vert[+0].xlon
            temp.vert2["coord"][0, 1] = \
                rdat.vert[+0].ylat

            nset.append(temp.vert2)

        if (rdat.flag >= 1 
                and rdat.vert[-1].seed >= 1):

            npts = +1
            temp = jigsawpy.jigsaw_msh_t()
            temp.vert2 = np.zeros(
                (npts + 0), dtype=temp.VERT2_t)

            temp.vert2["coord"][0, 0] = \
                rdat.vert[-1].xlon
            temp.vert2["coord"][0, 1] = \
                rdat.vert[-1].ylat

            nset.append(temp.vert2)
            
    init.vert2 = np.concatenate(nset, axis=+0)

    #init.vert2["coord"] *= np.pi / 180.
    
    jigsawpy.savemsh(tags + "_init.msh", init)


    spac = jigsawpy.jigsaw_msh_t()

    spac.mshID = "euclidean-grid"
    spac.radii = \
        rsph * np.ones(3, dtype=init.REALS_t)

    spac.xgrid = hfun.xpos #* np.pi / 180.
    spac.ygrid = hfun.ypos #* np.pi / 180.
    spac.value = hfun.vals

    jigsawpy.savemsh(tags + "_spac.msh", spac)

    return

def filter_hydrosheds(args):
    import os
    import time

    import netCDF4 as nc
    import numpy as np
    import fiona
    import argparse

    from reach.filter import filter_init, filter_eval, \
                            filter_core, reach_dat, reach_xyz
    from reach.interp import interp_2dim, interp_grid

    print("Loading filter spacing...")
    ttic = time.time()

    hfun = interp_grid()
    with nc.Dataset(args.spc_file) as data:
        hfun.xpos = \
            np.asarray(data.variables["xlon"][:])
        hfun.ypos = \
            np.asarray(data.variables["ylat"][:])
        hfun.vals = \
            np.asarray(data.variables["vals"][:])*100000.0

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")

    print("Loading stream network...")

    ttic = time.time()

    rnet = []; rpos = 0
    for feat in fiona.open(args.clp_file, "r"):   
        if (args.flt_endo and 
                feat["properties"]["ENDORHEIC"]):
            continue

        rdat = reach_dat()
        rdat.flag = + 1
        rdat.rpos = rpos
        rdat.rank = feat["properties"]["UPLAND_SKM"]
        rdat.hpos = feat["properties"]["HYRIV_ID"]
        rdat.dpos = feat["properties"]["NEXT_DOWN"]

        xmin, xmax = +180., -180.
        ymin, ymax = +90.0, -90.0
        for vert in feat["geometry"]["coordinates"]:
            xmin = min(xmin, vert[0])
            xmax = max(xmax, vert[0])
            ymin = min(ymin, vert[1])
            ymax = max(ymax, vert[1])

            rdat.vert.append(
                reach_xyz(vert[0], vert[1]))

        rpos = rpos + 1

        rnet.append (rdat)
        
    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")

    print("Form stream filtration...")

    ttic = time.time()

    filter_init(rnet, args.sph_size)
    filter_eval(rnet, hfun)
    filter_core(rnet, args.sph_size)
    hfun.vals = hfun.vals / 100000.0
    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")

    print("Writing stream network...")

    ttic = time.time()

    dir_ = os.path.dirname(os.path.abspath(args.out_file))

    if (not os.path.exists(dir_)):
        os.makedirs(dir_)

    src_ = fiona.open(args.clp_file)
    dst_ = fiona.open(
        args.out_file, "w", crs=src_.crs, schema=src_.schema, 
        driver="ESRI Shapefile")

    hmap = set()
    for rdat in rnet:
        if (rdat.flag >= +1):
            for rpos in rdat.rsub:
                hmap.add(rnet[rpos].hpos)

    for feat in src_:
        if (feat["properties"]["HYRIV_ID"] in hmap):
            dst_.write(feat)

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")

    if (args.msh_tags == ""): return

    print("Writing jigsaw outputs...")

    ttic = time.time()

    dir_ = os.path.dirname(os.path.abspath(args.msh_tags))

    if (not os.path.exists(dir_)):
        os.makedirs(dir_)

    output_jgsw(args, rnet, hfun, args.msh_tags, args.sph_size)

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")