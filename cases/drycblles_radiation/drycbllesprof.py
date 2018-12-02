import numpy as np
import netCDF4 as nc

# Settings
float_type = "f8"

dthetadz = 0.003

site = 0
expt = 0


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

# Create a group for the radiation and set up the values.
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

nc_surface_emissivity = nc_group_radiation.createVariable("surface_emissivity", float_type)
nc_surface_temperature = nc_group_radiation.createVariable("surface_temperature", float_type)

nc_surface_emissivity[:] = nc_file_rfmip.variables["surface_emissivity"][site]
nc_surface_temperature[:] = nc_file_rfmip.variables["surface_temperature"][expt,site]

nc_h2o = nc_group_radiation.createVariable("h2o", float_type, ("layer"))
nc_o3  = nc_group_radiation.createVariable("o3" , float_type, ("layer"))
nc_co2 = nc_group_radiation.createVariable("co2", float_type)
nc_n2o = nc_group_radiation.createVariable("n2o", float_type)
nc_co  = nc_group_radiation.createVariable("co" , float_type)
nc_ch4 = nc_group_radiation.createVariable("ch4", float_type)
nc_o2  = nc_group_radiation.createVariable("o2" , float_type)
nc_n2  = nc_group_radiation.createVariable("n2" , float_type)

nc_h2o[:] = nc_file_rfmip.variables["water_vapor"][expt,site,:] * float(nc_file_rfmip.variables["water_vapor"].units)
nc_o3 [:] = nc_file_rfmip.variables["ozone"][expt,site,:]       * float(nc_file_rfmip.variables["ozone"].units)
nc_co2[:] = nc_file_rfmip.variables["carbon_dioxide_GM"][expt]  * float(nc_file_rfmip.variables["carbon_dioxide_GM"].units)
nc_n2o[:] = nc_file_rfmip.variables["nitrous_oxide_GM"][expt]   * float(nc_file_rfmip.variables["nitrous_oxide_GM"].units)
nc_co [:] = nc_file_rfmip.variables["carbon_monoxide_GM"][expt] * float(nc_file_rfmip.variables["carbon_monoxide_GM"].units)
nc_ch4[:] = nc_file_rfmip.variables["methane_GM"][expt]         * float(nc_file_rfmip.variables["methane_GM"].units)
nc_o2 [:] = nc_file_rfmip.variables["oxygen_GM"][expt]          * float(nc_file_rfmip.variables["oxygen_GM"].units)
nc_n2 [:] = nc_file_rfmip.variables["nitrogen_GM"][expt]        * float(nc_file_rfmip.variables["nitrogen_GM"].units)

nc_file_rfmip.close()
nc_file.close()
