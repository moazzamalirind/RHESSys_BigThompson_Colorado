#Source: https://github.com/Fire-and-Dryland-Ecosystems-Lab/EcohydrologyScripts/blob/main/Rscripts/CARB_Workflow/clim_2_ncdf_makebasestation.R

# clim_ncdf_makebasestation
#
# Create a RHESSys basestation file based on Gridmet netcdf inputs
#
# Requires
#   - createbaseinfo_netcdf.c compiled binary
#
# Inputs:
#   - Processed/edited/subset aggregated (through time), netcdf climate data, via: #http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
#     (easies to use/download spatial subset to begin with, via the NetcdfSubset download option)
#     Has already been processed and modified via clim_ncdf_processgridmet
#   - Basin and DEM raster maps

library(ncdf4)
library(terra)

plots = T


# ------------------------------ INPUTS ------------------------------
# uncompiled source of the createbaseinfo_netcdf.c code
ncbase_src = "Data/createbaseinfo_netcdf.c"

# USE SAME BASIN RASTER AS PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
basin = rast("Data/Basin")      #use same which you using for RHESSys.
DEM = rast("Data/DEM")

# location for maps to be used with base station creation
map_dest = "clim/netcdfmaps"

# where to output a new zone map based on Netcdf grid
zone_dest = "clim"

# specify climate files individualy or via pattern, USE PROCESSED CLIMATE INPUTS FROM PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
clim_files = list.files(path = "clim",pattern = "crop_agg_met_", full.names = T)

# CHECK THE EDITS TO THE BASE STATION FILE AT THE END, FILE NAMES MAY NEED TO BE CORRECTED

# ------------------------------ OUTPUT NETCDF AS TIF FOR USE IN CREATING GRID ------------------------------
# we know extent should be since we just cropped using it
basin_vect = as.polygons(basin)
basin_vect_unproj = project(basin_vect, "+proj=longlat +datum=WGS84 +no_defs ")

clim_tmp = rast(clim_files[1])
#ext(clim_tmp) = ext(basin_vect_unproj)  #This line is creating problem. It is changing resolution
#the way around is generate data for area larger then watershed and then mask extra data in Grass gis.

# output geotif raster using first day of data from ncdf
# writeRaster(clim_tmp[[1]], "preprocessing/spatial_source/GRIDSOURCE_crop_agg_met_pr_1979_CurrentYear_CONUS.tif", overwrite=T)
nc_grid = clim_tmp[[1]]

# ------------------------------ GENERATE GRID MAP IN PROJECTED AND UNPROJECTED CRS  ------------------------------
#nc_grid = rast("preprocessing/spatial_source/GRIDSOURCE_crop_agg_met_pr_1979_CurrentYear_CONUS.tif")
id_grid = rast(nc_grid)
values(id_grid) = seq_along(values(nc_grid))
names(id_grid) = "ID"
id_grid_proj = project(id_grid, basin) # back to original projection
writeRaster(id_grid_proj, "Data/id_map.tif", overwrite=T)


if (plots) {
  # projected maps
  plot(id_grid_proj)
  plot(basin_vect, add=T)
  # unprojected maps
  plot(id_grid)
  plot(basin_vect_unproj, add=T)
}

# ------------------------------ GENERATE INPUTS FOR NETCDF BASE STATION CREATION  ------------------------------
# Maps need to be based on the netcdf resolution + extent
dem_nc = resample(DEM, nc_grid)
# LAI
#I need to figure out how to use data from real LAI map.. replicate file from ming.
#lai_nc = rast ("Data/LAI.tif")
#lai_nc = mask (lai_nc, DEM, filename="MASK_LAI.tif")

lai_nc = dem_nc

values(lai_nc) = 3.0
names(lai_nc) = "lai"

if (!file.exists(map_dest)) {
  dir.create(map_dest)
}

writeRaster(dem_nc, filename = file.path(map_dest,"dem.asc"), filetype="AAIGrid", gdal = c("FORCE_CELLSIZE=TRUE"), NAflag=-9999, overwrite=T)
writeRaster(lai_nc, filename = file.path(map_dest,"lai.asc"), filetype="AAIGrid", gdal = c("FORCE_CELLSIZE=TRUE"), NAflag=-9999, overwrite=T)
writeRaster(id_grid, filename = file.path(map_dest,"cellid.asc"), filetype="AAIGrid", datatype =  "INT4S", gdal = c("FORCE_CELLSIZE=TRUE"), NAflag=-9999, overwrite=T)
file.remove(list.files(path = map_dest, pattern = ".prj|.aux.xml", full.names = T))

# format: ID ID Y X X Y
xyid_loc = cbind(1:nrow(crds(id_grid)), 1:nrow(crds(id_grid)), crds(id_grid)[,c("y","x")], crds(id_grid))
write.table(xyid_loc, file.path(map_dest,"xyid_loc.txt"), row.names = F, col.names = F)

# output the zone grid for use in preprocessing
zone = mask(x = id_grid_proj, mask = basin)
writeRaster(zone, file.path(zone_dest, "zone.tif"))

# ------------------------------ RUN C BIN TO GET BASE STATION  ------------------------------
# COMPILE
system(noquote(paste("bash -c \"", paste0("gcc ", ncbase_src, " -o ", file.path(map_dest,"create_netcdfbase")), "\"", sep = "")))

# RUN -- YOU MAY NEED TO EDIT THIS, command should look like: ./create_netcdfbase cellid.asc lai.asc dem.asc xyid_loc.txt clim/ netcdf.base
cmd = paste0(file.path(map_dest,"create_netcdfbase"), " ", file.path(map_dest,"cellid.asc"), " ", file.path(map_dest,"lai.asc"), " ",
             file.path(map_dest,"dem.asc"), " ", file.path(map_dest,"xyid_loc.txt"), " ", dirname(clim_files[1]), "/ ", file.path(map_dest,"netcdf.base"))
system(noquote(paste("bash -c \"",cmd ,"\"", sep = "")))

# ------------------------------ EDIT/FIX BASE STATION  ------------------------------

varnames = rep(NA, length(clim_files))
for (i in seq_along(clim_files)) {
  tmp = nc_open(clim_files[i])
  varnames[i] = names(tmp$var)
}

ncbase = read.table(file.path(map_dest,"netcdf.base"))

ncbase[ncbase[,2]=="year_start_index", 1] = "1900"

ncbase[ncbase[,2]=="netcdf_tmax_filename", 1] = clim_files[grepl("tmmx", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_tmax", 1] = varnames[grepl("tmmx", clim_files)]

ncbase[ncbase[,2]=="netcdf_tmin_filename", 1] = clim_files[grepl("tmmn", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_tmin", 1] = varnames[grepl("tmmn", clim_files)]

ncbase[ncbase[,2]=="netcdf_rain_filename", 1] = clim_files[grepl("pr", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_rain", 1] = varnames[grepl("pr", clim_files)]

ncbase[ncbase[,2]=="netcdf_huss_filename", 1] = clim_files[grepl("daily_mean_specific_humidity", varnames)]
ncbase[ncbase[,2]=="netcdf_var_huss", 1] = varnames[grepl("daily_mean_specific_humidity", varnames)]

ncbase[ncbase[,2]=="netcdf_rmax_filename", 1] = clim_files[grepl("rmax", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_rmax", 1] = varnames[grepl("rmax", clim_files)]

ncbase[ncbase[,2]=="netcdf_rmin_filename", 1] = clim_files[grepl("rmin", clim_files)]
ncbase[ncbase[,2]=="netcdf_var_rmin", 1] = varnames[grepl("rmin", clim_files)]

ncbase[ncbase[,2]=="netcdf_rsds_filename", 1] = clim_files[grepl("daily_mean_shortwave_radiation_at_surface", varnames)]
ncbase[ncbase[,2]=="netcdf_var_rsds", 1] = varnames[grepl("daily_mean_shortwave_radiation_at_surface", varnames)]

ncbase[ncbase[,2]=="netcdf_was_filename", 1] = clim_files[grepl("daily_mean_wind_speed", varnames)]
ncbase[ncbase[,2]=="netcdf_var_was", 1] = varnames[grepl("daily_mean_wind_speed", varnames)]

write.table(ncbase, file.path(map_dest,"netcdf.base"), row.names = F, col.names = F, quote = F)




