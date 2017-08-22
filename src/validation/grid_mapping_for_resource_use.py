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
# now assigning number for each resource use
# now assigning number for each resource use
# 0 - absence       
# 1 - plantation
# 2 - deforestation
# 3 - landcover change
# 4 - logging
# 5 - primary forest < 750m
# 6 - primary forest > 750
# 7 - regrowth
# 8 - plantations before 2000
# 9 - other 
"""


import numpy as np
import MacPyver as mp
import os
 
 
"""
# absence 
os.system("gdal_rasterize \
-l absence_shape_rep \
-burn 1 \
-ot Float32 \
-te -1751798.359 1267454.241  -629798.359 2560454.241 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep.shp \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep.tif")

os.system("gdal_rasterize \
-l absence_shape_rep \
-burn 1 \
-ot Float32 \
 -ts 1122 1293 \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep.shp \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep.tif")

os.system("gdalwarp \
-t_srs '+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs' \
-ot Float64 -r near \
-te -1751798.359 1267454.241 -629798.359 2560454.241 \
-ts 1122 1293 \
-overwrite  \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_rep.tif \
/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_repro_res.tif"
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



"""

absence_path = "/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_expanded_22_08_17_repro_res.tif"
absence = mp.tiff.read_tif(absence_path, 1)

# input grid with grid_ids as values

grid = mp.tiff.read_tif("/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/future/grid_with_id_repro_res.tif", 1)
np.unique(grid)

resource_use = mp.tiff.read_tif("/homes/mv39zilo/work/Borneo/data/resource_use/resource_use_grid.tif", 1)
np.unique(resource_use)  

country_layer_path = "/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_country_repro_res.tif"
borneo = mp.tiff.read_tif(country_layer_path, 1)
"""
# format
1        0    0     0   0   0 000 000
(1 / 2)  0 (3 / 6)  0 (1-9) 0 700 000
out/in,    MYS/IDN,  category,grid_id
"""

# pixel is out / out
resource_grid_1 = np.where(absence == 1, 100000000000, 200000000000)

# country
resource_grid_2 = np.where(borneo == 136, (resource_grid_1 + 3 * 1000000000), 
                resource_grid_1)
resource_grid_3 = np.where(borneo == 106, (resource_grid_2 + 6 * 1000000000), 
                resource_grid_2)           
                
                
resource_grid_4 = np.where(resource_use > 0, (resource_grid_3 + resource_use * 10000000), 
                resource_grid_3)
np.unique(resource_grid_4 )




 # now assigning number for each resource use

resource_grid_5 = np.where(grid > 0, (resource_grid_4 + grid), 
                resource_grid_4 )
np.unique(resource_grid_5)


mp.tiff.write_tif(file_with_srid = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/future/grid_with_id_repro_res.tif", 
                   full_output_name = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/resource_grid_absence_country_category_id.tif",
                   data = resource_grid_5, 
                   dtype = 5)
     
# use this in R script "prepare_boot_grid.R"
                   

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
