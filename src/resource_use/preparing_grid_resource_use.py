# -*- coding: utf-8 -*-
"""
Spyder Editor

# script to prepare the grid with the resource use information
# in the output file we will have a
# a 1 at the first position if pixel has to be discarded
# a 2 at the first position if the pixel is considered
# a 1 at the second position for IDN
# a 2 at the second position for MYS 
# a 0 - 5 at the 3rd position depending on the resource use category
# the grid_id from the 4rth to the 10th position (6 positions)
# 
# we will now clip with populations and then also do this in the 
# other resource use table, let's see how this looks like
"""


import numpy as np
import MacPyver as mp
import os
 
 
"""
# absence 
os.system("gdal_rasterize \
-l absence_diss_exp \
-burn 1 \
-ot Float32 \
-te -1751798.359 1267454.241  -629798.359 2560454.241 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep_expanded_22_08_17.shp \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep_expanded_22_08_17.tif")

os.system("gdal_rasterize \
-l absence_shape_rep_expanded_22_08_17 \
-burn 1 \
-ot Float32 \
 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep_expanded_22_08_17.shp \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep_expanded_22_08_17.tif")

os.system("gdalwarp \
-t_srs '+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' \
-ot Float64 -r bilinear \
-te -1751798.359 1267454.241 -629798.359 2560454.241 \
-ts 1122 1293 \
-overwrite  \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep_expanded_22_08_17.tif \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_expanded_22_08_17_repro_res.tif"
)

os.system("gdalwarp \
-t_srs '+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' \
-ot Float64 -r near \
-te -1751798.359 1267454.241 -629798.359 2560454.241 \
-ts 1122 1293 \
-overwrite  \
/homes/mv39zilo/work/Borneo/data/future/grid.tif \
/homes/mv39zilo/work/Borneo/data/future/grid_repro_res.tif"
)


os.system("gdalwarp \
-t_srs '+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' \
-ot Float64 -r near \
-te -1751798.359 1267454.241 -629798.359 2560454.241 \
-ts 1122 1293 \
-overwrite  \
/homes/mv39zilo/work/Borneo/data/predictors/cleaned_data/cleaned_predictors/dem.tif \
/homes/mv39zilo/work/Borneo/data/predictors/cleaned_data/cleaned_predictors/dem_repro_res.tif"
)

"""





absence_path = "/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_expanded_22_08_17_repro_res.tif"
absence = mp.tiff.read_tif(absence_path, 1)

grid_layer_path = "/homes/mv39zilo/work/Borneo/data/future/grid_repro_res.tif"
grid = mp.tiff.read_tif(grid_layer_path, 1)


grid = np.where(absence == 0, grid, 0)    
np.unique(grid)


# now assigning number for each resource use
# 0 - absence       
# 1 - plantation
# 2 - deforestation
# 3 - logging
# 4 - primary forest < 750m
# 5 - primary forest > 750
# 6 - regrowth
# 7 - plantations before 2000
# 8 - other 


gaveau_defor_path = "/homes/mv39zilo/work/Borneo/data/predictors/raw_data/Remaining_forest_in_2015_and_deforestation_1973-2015/REGIONBorneo_FCDefDeg_1973to2015_CIFOR_repro_res.tif"
gaveau_defor_raw = mp.tiff.read_tif(gaveau_defor_path, 1)
gaveau_deforested = np.where((gaveau_defor_raw > 6 ) & 
                              (gaveau_defor_raw < 10), 1, 0)
gaveau_logged = np.where(gaveau_defor_raw == 2, 1, 0)
gaveau_primary_forest = np.where(gaveau_defor_raw == 1, 1, 0)
gaveau_regrowth = np.where(gaveau_defor_raw == 4, 1, 0)

#plantation                       
plantation_path = "/homes/mv39zilo/work/Borneo/data/predictors/raw_data/Land_cover_trajectory_before_oil-palm_and_pulpwood_establishment/raster_plantation.tif"
plantations = mp.tiff.read_tif(plantation_path, 1)                              

# old plantations
old_plantations_path = "/homes/mv39zilo/work/Borneo/data/predictors/raw_data/Land_cover_trajectory_before_oil-palm_and_pulpwood_establishment/raster_plantation_old_repro_res.tif"
old_plantations = mp.tiff.read_tif(old_plantations_path, 1)



dem_path = "/homes/mv39zilo/work/Borneo/data/predictors/cleaned_data/cleaned_predictors/dem_repro_res.tif"
dem = mp.tiff.read_tif(dem_path, 1)



# elevation
# --> how much of primary forest in high elevation

# we can save some computing time by only considering
# the grid and not everything, because abundance will only be
# in rid
grid = np.where((plantations == 1) & 
                 (grid == 1), 1 , grid)
grid = np.where((plantations == 0) & 
                 (gaveau_deforested  == 1)&
                 (grid == 1), 2 , grid)               
                 
# 3 logged
grid = np.where((plantations == 0) & 
                 (gaveau_deforested == 0)&
                 (gaveau_logged == 1) & 
                 (grid == 1), 3 , grid)
# 4 primary forest
grid = np.where((plantations == 0) & 
                 (gaveau_deforested == 0)&
                 (gaveau_logged == 0) & 
                 (gaveau_primary_forest == 1) &
                 (grid == 1), 4 , grid)   
                 
# add here primary forest at high altitudes
# because everywhere where we have primary forest 
# we no longer have a 1 in grid but a six we do it differently   
# 5 primary montane
grid  = np.where((grid == 4) & 
                (dem > 750), 5, grid)    


# 6 regrowth
grid = np.where((plantations == 0) & 
                 (gaveau_deforested == 0)&
                 (gaveau_logged == 0) & 
                 (gaveau_primary_forest == 0) &
                 (gaveau_regrowth ==1) &
                 (grid == 1), 6 , grid)  
                 
# 7 is old plantations
            
grid = np.where((plantations == 0) & 
                 (gaveau_deforested == 0)&
                 (gaveau_logged == 0) & 
                 (gaveau_primary_forest == 0) &
                 (gaveau_regrowth == 0) &
                 (old_plantations == 1) & 
                 (grid == 1), 7 , grid)              
                 
                 

# 8 is other
grid = np.where((plantations == 0) & 
               (gaveau_deforested == 0)&
               (gaveau_logged == 0) & 
               (gaveau_primary_forest == 0) &
               (gaveau_regrowth ==0) &
               (grid == 1), 8, grid) 
                 

np.unique(grid)

mp.tiff.write_tif(file_with_srid = grid_layer_path, 
                   full_output_name = "/homes/mv39zilo/work/Borneo/data/resource_use/resource_use_grid.tif",
                   data =  grid, 
                   dtype = 0)  
       

# grid map with four categories
# 0 - absence                 - 0
# 1 - plantation              - 1
# 2 - deforestation           - 1
# 3 - logging                 - 1
# 4 - primary forest < 750m   - 2
# 5 - primary forest > 750    - 2
# 6 - regrowth                - 2
# 7 - plantations before 2000 - 3
# 8 - other                   - 3

grid_map = np.where((grid <= 3) & 
                    (grid > 0), 1, 0)  
grid_map = np.where((grid <= 6) & 
                    (grid > 3), 2, grid_map)                      
grid_map = np.where((grid > 6), 3, grid_map)    


# need to clip this with the area in which we have orangutans in 1999
year = 1999
abundance_path = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/abundMod/testing_ae_and_absence/pipeline_results/ppln_ae75m_50-2017-02-28T18-00-52/prediction_map_" + str(year) + "_2017-02-28_repro_res.tif"
abundance = mp.tiff.read_tif(abundance_path, 1)
  
  
grid_map = np.where(abundance > 0.1, grid_map, 0)
np.unique(grid_map)

mp.tiff.write_tif(file_with_srid = grid_layer_path, 
                   full_output_name = "/homes/mv39zilo/work/Borneo/data/resource_use/resource_use_grid_map.tif",
                   data =  grid_map, 
                   dtype = 0)  
"""
# test
indir = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use"
abundance_layer_path = indir + "/abundance_pred_2015_no_absence.tif"
# read the data
abundance_data = mp.tiff.read_tif(abundance_layer_path, 1)

test = np.where((grid == 0) &
                (abundance_data > 0), abundance_data, 0)
np.sum(test)
mp.tiff.write_tif(file_with_srid = grid_layer_path, 
                   full_output_name = "/homes/mv39zilo/work/Borneo/data/bootstrap/test.tif",
                   data =  test, 
                   dtype = 4)  
"""                
# rasterize borneo with info on country or province
                   
"""
os.system("gdal_rasterize \
-a ID_1 \
-l Borneo \
-ot Float32 \
-te -1751798.359 1267454.241  -629798.359 2560454.241 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo.shp \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_province.tif")

os.system("gdal_rasterize \
-a ID_1 \
-l Borneo \
-ot Float32 \
 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo.shp \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_province.tif")

os.system("gdalwarp \
-r near -t_srs '+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' \
-te -1751798.359 1267454.241  -629798.359 2560454.241 -ts 1122 1293 \
-overwrite \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_province.tif \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo__province_repro_res.tif")

os.system("gdal_rasterize \
-a ID_0 \
-l Borneo \
-ot Float32 \
-te -1751798.359 1267454.241  -629798.359 2560454.241 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo.shp \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_country.tif")

os.system("gdal_rasterize \
-a ID_0 \
-l Borneo \
-ot Float32 \
 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo.shp \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_country.tif")

os.system("gdalwarp \
-r near -t_srs '+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' \
-te -1751798.359 1267454.241  -629798.359 2560454.241 -ts 1122 1293 \
-overwrite \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_country.tif \
/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_country_repro_res.tif")

"""
