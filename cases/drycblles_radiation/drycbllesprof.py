import numpy as np
import netCDF4 as nc

# Settings
dthetadz = 0.003
site = 0
expt = 0
float_type = "f8"

# Get number of vertical levels and size from .ini file
with open('drycblles.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z   = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl = np.zeros(np.size(z))
qt  = np.zeros(np.size(z))

# linearly stratified profile
for k in range(kmax):
    thl[k] = 300. + dthetadz*z[k]

# Save all the input data to NetCDF
nc_file = nc.Dataset("drycblles.nc", mode="w", datamodel="NETCDF4", clobber=True)
nc_file_rfmip = nc.Dataset("rfmip.nc", mode="r", datamodel="NETCDF4", clobber=False)

nc_file.createDimension("z", kmax)
nc_z = nc_file.createVariable("z", float_type, ("z"))
nc_z[:] = z[:]

# Create a group called "init" for the initial profiles.
nc_group_init = nc_file.createGroup("init")

nc_thl = nc_group_init.createVariable("thl", float_type, ("z"))
nc_qt  = nc_group_init.createVariable("qt" , float_type, ("z"))
nc_thl[:] = thl[:]
nc_qt [:] = qt [:]

nc_group_radiation = nc_file.createGroup("radiation")
nc_group_radiation.createDimension("level", nc_file_rfmip.dimensions["level"].size)
nc_group_radiation.createDimension("layer", nc_file_rfmip.dimensions["layer"].size)
nc_group_radiation.createDimension("col", 1)

nc_pres_level = nc_group_radiation.createVariable("pres_level", float_type, ("level"))
nc_pres_layer = nc_group_radiation.createVariable("pres_layer", float_type, ("layer"))
nc_temp_level = nc_group_radiation.createVariable("temp_level", float_type, ("level"))
nc_temp_layer = nc_group_radiation.createVariable("temp_layer", float_type, ("layer"))

nc_pres_level[:] = nc_file_rfmip.variables["pres_level"][site,:]
nc_pres_layer[:] = nc_file_rfmip.variables["pres_layer"][site,:]
nc_temp_level[:] = nc_file_rfmip.variables["temp_level"][expt,site,:]
nc_temp_layer[:] = nc_file_rfmip.variables["temp_layer"][expt,site,:]

nc_surface_emissivity = nc_group_radiation.createVariable("surface_emissivity", float_type, ("col"))
nc_surface_temperature = nc_group_radiation.createVariable("surface_temperature", float_type, ("col"))

nc_surface_emissivity[:] = nc_file_rfmip.variables["surface_emissivity"][site]
nc_surface_temperature[:] = nc_file_rfmip.variables["surface_temperature"][expt,site]

nc_file_rfmip.close()
nc_file.close()
