import os
import numpy as np
import netCDF4 as nc
import sys


sys.path.insert(0, '/Users/xudo627/Developments/jigsaw-python')
import jigsawpy
import argparse
import matplotlib.pyplot as plt
import math
from msh2mpas import msh2mpas
from process_hydrosheds import download_file, clip_hydrosheds, filter_hydrosheds
import geopandas as gpd
class base: pass

name = "susquehanna"
tag_ = "filtered_5k"

args = base()

args.bnd_file = os.path.join( # location of Hy-SHEDS
    "dat",
    'wbd',
    "Susquehanna.shp")

args.shp_file = os.path.join( # location of Hy-SHEDS
    "dat", 
    "HydroRIVERS_v10_na_shp",
    "HydroRIVERS_v10_na.shp")

args.clp_file = os.path.join( # location of filtered
    "out",  
    tag_,
    "clipped_rivers.shp")

args.out_file = os.path.join( # location of filtered
    "out",  
    tag_,
    "filtered_rivers.shp")

args.spc_file = os.path.join( # location for spacing
    "out",
    tag_,
    "resolution.nc")

args.msh_tags = os.path.join( # location for meshes
    "msh",
    tag_,
    name+"_msh")

if (not os.path.exists("./dat")): os.makedirs("./dat")
dir_ = os.path.dirname(os.path.abspath(args.spc_file))
if (not os.path.exists(dir_)): os.makedirs(dir_)

if (not os.path.isfile(args.shp_file)):
    print('Downloading HydroRIVERS...')
    download_file(args)

if (not os.path.isfile(args.clp_file)):
    print('Clipping HydroRIVERS...')
    clip_hydrosheds(args)

nlon = 100; nlat = 100
bnd = gpd.read_file(args.bnd_file)

# Calculate bounds
args.box_xmin, args.box_ymin, args.box_xmax, args.box_ymax = bnd.total_bounds
args.box_xmin = math.floor(args.box_xmin) - 0.5
args.box_ymin = math.floor(args.box_ymin) - 0.5
args.box_xmax = math.ceil(args.box_xmax)  + 0.5
args.box_ymax = math.ceil(args.box_ymax)  + 0.5

data = nc.Dataset(args.spc_file, "w")
data.createDimension("nlon", nlon)
data.createVariable("xlon", "f8", ("nlon"))
data["xlon"][:] = np.linspace(args.box_xmin, args.box_xmax, nlon)
data.createDimension("nlat", nlat)
data.createVariable("ylat", "f8", ("nlat"))
data["ylat"][:] = np.linspace(args.box_ymin, args.box_ymax, nlat)
data.createVariable("vals", "f4", ("nlat", "nlon"))

xvec = np.asarray(data["xlon"][:]) #* np.pi / 180.
yvec = np.asarray(data["ylat"][:]) #* np.pi / 180.

xmat, ymat = np.meshgrid(xvec, yvec, sparse=True)

# 25km globally, 5km on Mississippi basin
filt = np.ones((nlat, nlon)) * 0.05

data["vals"][:, :] = filt  # in [m]
data.close()

args.sph_size = 6371220.  # earth radius
args.flt_endo = True  # strip no-outlet rivers??   

filter_hydrosheds(args)

opts = jigsawpy.jigsaw_jig_t()  
mesh = jigsawpy.jigsaw_msh_t()

opts.geom_file = args.msh_tags + "_geom.msh"
opts.hfun_file = args.msh_tags + "_spac.msh"
#opts.init_file = args.msh_tags + "_init.msh"
opts.mesh_file = args.msh_tags + "_mesh.msh"

opts.jcfg_file = args.msh_tags + "_opts.jig"

opts.verbosity = +1

opts.hfun_scal = "absolute"
opts.hfun_hmax = float("inf")       # null HFUN limits
#opts.hfun_hmin = float(+0.00)
opts.optm_qlim = +.95
opts.mesh_dims = +2
opts.mesh_rad2 = +1.20              # relax edge-ratio
#opts.mesh_eps1 = +1.00
#opts.mesh_top1 = True
    
#opts.optm_iter = +64                # tight optim. tol
#opts.optm_qtol = +5.E-05
#opts.optm_cost = "skew-cos"
#opts.optm_dual = True

jigsawpy.cmd.jigsaw(opts, mesh)

xc, yc, xv, yv, area, frac = msh2mpas(mesh, args.bnd_file, True)








