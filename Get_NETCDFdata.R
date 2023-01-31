#part 1
# clim_ndf_processgridmet
#
# This workflow takes aggregated (through time) GridMet netcdf climate data as input, and
# subsets and processes the netcdf data to be usable in RHESSys
#
# Requires
#   - NCO (installed on linux with apt with 'sudo apt install NCO')
#
# Inputs:
#   - Aggregated (through time), netcdf climate data, via: #http://thredds.northwestknowledge.net:8080/thredds/reacch_climate_MET_aggregated_catalog.html
#     (easies to use/download spatial subset to begin with, via the NetcdfSubset download option)
#   - Basin raster

#install.packages("ncdf4")
#install.packages("terra")
library(ncdf4)
library(terra)
library(raster)
library(sp)
library(dplyr)
plot = T


#download data ...I did manually.

# source netcdfs (these are default file names from THREADS aggregate download)
# pr = "../data/gridmet/agg_met_pr_1979_CurrentYear_CONUS.nc"
# tmin = "../data/gridmet/agg_met_tmmn_1979_CurrentYear_CONUS.nc"
# tmax = "../data/gridmet/agg_met_tmmx_1979_CurrentYear_CONUS.nc"
# rmax = "../data/gridmet/agg_met_rmax_1979_CurrentYear_CONUS.nc"
# rmin = "../data/gridmet/agg_met_rmin_1979_CurrentYear_CONUS.nc"
# vpd = "../data/gridmet/agg_met_vpd_1979_CurrentYear_CONUS.nc"
# sph = "../data/gridmet/agg_met_sph_1979_CurrentYear_CONUS.nc"
# srad = "../data/gridmet/agg_met_srad_1979_CurrentYear_CONUS.nc"
# wdir = "../data/gridmet/agg_met_th_1979_CurrentYear_CONUS.nc"
# wspd = "../data/gridmet/agg_met_vs_1979_CurrentYear_CONUS.nc"
# clim_files = c(pr, tmin, tmax, rmax, rmin, vpd, sph, srad, wdir, wspd)

# OR
setwd("E:/Fire_Research/RHESSys_Bigthompson/Climate_data")

Basin = rast("Data/Basin_Raster.tif")

clim_files = list.files(path = "Data",pattern = "agg_met_", full.names = T)

# destination folder for cropped + edited ncdf clim files
clim_dest=  "Clim"

# ------------------------------ CROP + PROCESS NETCDF CLIMATE DATA USING BASIN EXTENT AND NCO ------------------------------
# FOR NETCDF TO WORK WITH RHESSYS:
#   Must use 'time' as temporal dimension, below we change every instance of 'day' to 'time' to be safe
#   Must have variable as floating point/double data type, NOT INTEGER


basin_vect = as.polygons(Basin)
# writeVector(basin_vect, "preprocessing/spatial_source/basin_vect.shp", overwrite=T)

# projection for netcdf data is: +proj=longlat +datum=WGS84 +no_defs
basin_vect_unproj = terra::project(basin_vect, "+proj=longlat +datum=WGS84 +no_defs ")
# writeVector(basin_vect_unproj, "preprocessing/spatial_source/basin_vect_unproj.shp", overwrite=T)

# this should be in lon lat wgs84 to match ncdf - see xy axes
if (plot) {plot(basin_vect_unproj)}

basin_extent = ext(basin_vect_unproj)

#check whether data is uploaded correctly
#infile= clim_files[[1]]

# iterate through each climate file
for (infile in clim_files) {
  # get the variable name for when modifying data type
  nctmp = nc_open(infile)

  # TO CHECK THAT LAT AND LON IS -180 TO 180 (IT SEEMS TO BE - this is different from MACA data)
  # dims = unlist(strsplit(ncatt_get(nctmp, "precipitation_amount", "dimensions")$value, " "))
  # lonvals = nctmp$var$precipitation_amount$dim[[which(dims=="lon")]]$vals

  varname = names(nctmp$var)
  nc_close(nctmp)
  outfile = paste0(file.path(clim_dest, "crop_"), basename(infile))


  crop_cmd = paste0("ncks -O -d lon,",basin_extent[1],",",basin_extent[2]," -d lat,",basin_extent[3],",",basin_extent[4]," ",infile," ",outfile,"\n",
                    "ncrename -d day,time -v day,time -O ", outfile, " ", outfile, "\n",
                    "ncatted -O -a coordinates,,m,c,'time lat lon' ", outfile, "\n",
                    "ncap2 -O -s '",varname,"=double(",varname,")' ", outfile," ", outfile)
  # ASSUMES EITHER UNIX OR IF ON WINDOWS USING WSL
  if (.Platform$OS.type == "windows") {
    tmp = noquote(paste("bash -c \"", crop_cmd, "\"", sep = ""))
  }

  system(tmp)
}



#########################################################
######################## Part 2######################################################################

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
ncbase_src =  file.path("Data","createbaseinfo_netcdf.c")


# USE SAME BASIN RASTER AS PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
basin = rast("Data/Basin_Raster.tif")
DEM = rast("Data/Basin_DEM.tif")

# location for maps to be used with base station creation
map_dest = "Clim/netcdfmaps"

# where to output a new zone map based on Netcdf grid
zone_dest = "Clim"

# specify climate files individual or via pattern, USE PROCESSED CLIMATE INPUTS FROM PREVIOUS SCRIPT (clim_1_ncdf_processgridmet.R)
clim_files = list.files(path = "Clim",pattern = "crop_agg_met_", full.names = T)

# CHECK THE EDITS TO THE BASE STATION FILE AT THE END, FILE NAMES MAY NEED TO BE CORRECTED

# ------------------------------ OUTPUT NETCDF AS TIF FOR USE IN CREATING GRID ------------------------------
# we know extent should be since we just cropped using it
basin_vect = as.polygons(basin)
basin_vect_unproj = terra::project(basin_vect, "+proj=longlat +datum=WGS84 +no_defs ")

clim_tmp = rast(clim_files[1])
ext(clim_tmp) = ext(basin_vect_unproj)
# output geotif raster using first day of data from ncdf
# writeRaster(clim_tmp[[1]], "preprocessing/spatial_source/GRIDSOURCE_crop_agg_met_pr_1979_CurrentYear_CONUS.tif", overwrite=T)
nc_grid = clim_tmp[[1]]

# ------------------------------ GENERATE GRID MAP IN PROJECTED AND UNPROJECTED CRS  ------------------------------
#nc_grid = rast("preprocessing/spatial_source/GRIDSOURCE_crop_agg_met_pr_1979_CurrentYear_CONUS.tif")
id_grid = rast(nc_grid)
values(id_grid) = seq_along(values(nc_grid))
names(id_grid) = "ID"
id_grid_proj = terra::project(id_grid, basin) # back to original projection

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
writeRaster(zone, file.path(zone_dest, "zone.tif"), overwrite=TRUE)

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

