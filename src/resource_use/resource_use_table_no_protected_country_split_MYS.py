# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 12:51:35 2017

@author: mv39zilo
script to make the table for the graph and the table
# 0 - absence       
# 1 - plantation
# 2 - deforestation
# 3 - logging
# 4 - primary forest < 750m
# 5 - primary forest > 750
# 6 - regrowth
# 7 - plantations before 2000
# 8 - other 
"""

import numpy as np
import MacPyver as mp
import os
 
year =  1999
 
grid_layer_path = "/homes/mv39zilo/work/Borneo/data/resource_use/resource_use_grid.tif"
resource_use = mp.tiff.read_tif(grid_layer_path, 1)
np.unique(resource_use)

country_layer_path = "/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_country_repro_res.tif"
borneo = mp.tiff.read_tif(country_layer_path, 1)
np.unique(borneo)

province_layer_path = "/homes/mv39zilo/work/Borneo/data/auxiliary_additional_data/Borneo_shape/cleaned_data/Borneo_province_repro_res.tif"
borneo_province = mp.tiff.read_tif(province_layer_path, 1)
np.unique(borneo_province)


populations = mp.tiff.read_tif("/homes/mv39zilo/work/Borneo/data/populations_phva/meta_kalsarsab_20002015_diss_no_reintro_repro_res.tif", 1)

#

absence_path = "/homes/mv39zilo/work/Borneo/data/response/cleaned_data/absence_shape/absence_shape_expanded_22_08_17_repro_res.tif"
absence = mp.tiff.read_tif(absence_path, 1)
if year==1999:
    abundance_path = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/abundMod/testing_ae_and_absence/pipeline_results/ppln_ae75m_50-2017-02-28T18-00-52/prediction_map_" + str(year) + "_2017-02-28_repro_res.tif"
if year==2015:
    abundance_path = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/abundMod/testing_ae_and_absence/pipeline_results/ppln_ae75m_50-2017-02-28T18-00-52/prediction_map_" + str(year) + "_2017-03-01_repro_res.tif"
  
print abundance_path
# read the data
abundance_data = mp.tiff.read_tif(abundance_path, 1)


# we want to know mean density of ou in plantation
# so all pixels where no ou is 0, not na, also for absence

abundance_data = np.where(populations != 0, abundance_data, 0)


resource_use_MYS = np.where(borneo == 136, resource_use, 0)
np.unique(resource_use_MYS)



resource_use_SAB = np.where(borneo_province == 10, resource_use, 0)
np.unique(resource_use_SAB)

resource_use_SAW = np.where(borneo_province == 11, resource_use, 0)
np.unique(resource_use_SAW)


#IDN 
resource_use_IDN = np.where(borneo == 106, resource_use, 0)




"""
mp.tiff.write_tif(file_with_srid = grid_layer_path, 
                   full_output_name = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/test.tif",
                   data =  resource_use_SAW , 
                   dtype = 4)

"""




#MYS

file_path_table = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/abundance_density_resource_use_" + str(year) + "no_absence.csv"
resource_use = open(file_path_table, 'w')

resource_use.write("type, country, abundance, percent_total, area, density \n")


abundance_data_MYS = np.where(borneo == 136, abundance_data, 0 )
count_abundance_MYS = np.count_nonzero(abundance_data_MYS)
abundance_sum_MYS = np.sum(abundance_data_MYS)
mean_abundance_MYS = abundance_sum_MYS/count_abundance_MYS 

resource_use.write( "total abundance" + "," + "MYS" + "," + str(abundance_sum_MYS) + "," + 
  str(abundance_sum_MYS/abundance_sum_MYS * 100 ) + "," + 
  str(count_abundance_MYS) +"," + 
  str(mean_abundance_MYS) + "\n")


# planted
planted_abundance_MYS = np.where(resource_use_MYS == 1, abundance_data_MYS, 0)
abundance_count_MYS_planted = np.count_nonzero(planted_abundance_MYS)
abundance_sum_MYS_planted = np.sum(planted_abundance_MYS)
abundance_mean_MYS_planted = abundance_sum_MYS_planted / abundance_count_MYS_planted


resource_use.write( "plantation" +  "," + "MYS" + "," + str(abundance_sum_MYS_planted) + "," +
 str(abundance_sum_MYS_planted/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_planted) +"," + 
 str(abundance_mean_MYS_planted) + "\n")



# deforested
deforested_abundance_MYS = np.where(resource_use_MYS == 2, abundance_data_MYS, 0)
abundance_count_MYS_deforested = np.count_nonzero(deforested_abundance_MYS)
abundance_sum_MYS_deforested = np.sum(deforested_abundance_MYS)
abundance_mean_MYS_deforested = abundance_sum_MYS_deforested/abundance_count_MYS_deforested

resource_use.write( "deforestation" +  "," + "MYS" + "," + str(abundance_sum_MYS_deforested) + "," +
 str(abundance_sum_MYS_deforested/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_deforested) +"," + 
 str(abundance_mean_MYS_deforested) + "\n")
 
 
 

 # logged

logged_abundance_MYS = np.where(resource_use_MYS == 3, abundance_data_MYS, 0)
abundance_count_MYS_logged = np.count_nonzero(logged_abundance_MYS)
abundance_sum_MYS_logged = np.sum(logged_abundance_MYS)
abundance_mean_MYS_logged = abundance_sum_MYS_logged/abundance_count_MYS_logged

resource_use.write( "logging" +  "," + "MYS" + "," + str(abundance_sum_MYS_logged) + "," +
 str(abundance_sum_MYS_logged/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_logged) +"," + 
 str(abundance_mean_MYS_logged) + "\n")
 
 # primary
primary_abundance_MYS = np.where(resource_use_MYS == 4, abundance_data_MYS, 0)
abundance_count_MYS_primary = np.count_nonzero(primary_abundance_MYS)
abundance_sum_MYS_primary = np.sum(primary_abundance_MYS)
abundance_mean_MYS_primary = abundance_sum_MYS_primary/abundance_count_MYS_primary


resource_use.write( "primary forest" +  "," + "MYS" + "," + str(abundance_sum_MYS_primary) + "," +
 str(abundance_sum_MYS_primary/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_primary) +"," + 
 str(abundance_mean_MYS_primary) + "\n")


# primary montane
primary_montane_abundance_MYS = np.where(resource_use_MYS == 5, abundance_data_MYS, 0)
abundance_count_MYS_primary_montane = np.count_nonzero(primary_montane_abundance_MYS)
abundance_sum_MYS_primary_montane = np.sum(primary_montane_abundance_MYS)
abundance_mean_MYS_primary_montane = abundance_sum_MYS_primary_montane/abundance_count_MYS_primary_montane


resource_use.write( "primary montane forest" +  "," + "MYS" + "," + str(abundance_sum_MYS_primary_montane) + "," +
 str(abundance_sum_MYS_primary_montane/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_primary_montane) +"," + 
 str(abundance_mean_MYS_primary_montane) + "\n")


 # regrowth
regrowth_abundance_MYS = np.where(resource_use_MYS == 6, abundance_data_MYS, 0)
abundance_count_MYS_regrowth = np.count_nonzero(regrowth_abundance_MYS)
abundance_sum_MYS_regrowth = np.sum(regrowth_abundance_MYS)
abundance_mean_MYS_regrowth = abundance_sum_MYS_regrowth/abundance_count_MYS_regrowth


resource_use.write( "regrowth forest" +  "," + "MYS" + "," + str(abundance_sum_MYS_regrowth) + "," +
 str(abundance_sum_MYS_regrowth/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_regrowth) +"," + 
 str(abundance_mean_MYS_regrowth) + "\n")


# old plantations
old_plantations_abundance_MYS = np.where(resource_use_MYS == 7, abundance_data_MYS, 0)
abundance_count_MYS_old_plantations = np.count_nonzero(old_plantations_abundance_MYS)
abundance_sum_MYS_old_plantations = np.sum(old_plantations_abundance_MYS)
abundance_mean_MYS_old_plantations = abundance_sum_MYS_old_plantations/abundance_count_MYS_old_plantations


resource_use.write( "old plantations" +  "," + "MYS" + "," + str(abundance_sum_MYS_old_plantations) + "," +
 str(abundance_sum_MYS_old_plantations/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_old_plantations) +"," + 
 str(abundance_mean_MYS_old_plantations) + "\n")


# nothing

other_abundance_MYS = np.where(resource_use_MYS == 8, abundance_data_MYS, 0)
abundance_count_MYS_other = np.count_nonzero(other_abundance_MYS)
abundance_sum_MYS_other = np.sum(other_abundance_MYS)
abundance_mean_MYS_other = abundance_sum_MYS_other/abundance_count_MYS_other


resource_use.write( "other" +  "," + "MYS" + "," + str(abundance_sum_MYS_other) + "," +
 str(abundance_sum_MYS_other/abundance_sum_MYS * 100 ) + "," + 
 str(abundance_count_MYS_other) +"," + 
 str(abundance_mean_MYS_other) + "\n")
 
 
#######
# SAB #
#######

#SAB


abundance_data_SAB = np.where(borneo_province == 10, abundance_data, 0 )
count_abundance_SAB = np.count_nonzero(abundance_data_SAB)
abundance_sum_SAB = np.sum(abundance_data_SAB)
mean_abundance_SAB = abundance_sum_SAB/count_abundance_SAB 

resource_use.write( "total abundance" + "," + "SAB" + "," + str(abundance_sum_SAB) + "," + 
  str(abundance_sum_SAB/abundance_sum_SAB * 100 ) + "," + 
  str(count_abundance_SAB) +"," + 
  str(mean_abundance_SAB) + "\n")





# planted
planted_abundance_SAB = np.where(resource_use_SAB == 1, abundance_data_SAB, 0)
abundance_count_SAB_planted = np.count_nonzero(planted_abundance_SAB)
abundance_sum_SAB_planted = np.sum(planted_abundance_SAB)
abundance_mean_SAB_planted = abundance_sum_SAB_planted / abundance_count_SAB_planted


resource_use.write( "plantation" +  "," + "SAB" + "," + str(abundance_sum_SAB_planted) + "," +
 str(abundance_sum_SAB_planted/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_planted) +"," + 
 str(abundance_mean_SAB_planted) + "\n")



# deforested
deforested_abundance_SAB = np.where(resource_use_SAB == 2, abundance_data_SAB, 0)
abundance_count_SAB_deforested = np.count_nonzero(deforested_abundance_SAB)
abundance_sum_SAB_deforested = np.sum(deforested_abundance_SAB)
abundance_mean_SAB_deforested = abundance_sum_SAB_deforested/abundance_count_SAB_deforested

resource_use.write( "deforestation" +  "," + "SAB" + "," + str(abundance_sum_SAB_deforested) + "," +
 str(abundance_sum_SAB_deforested/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_deforested) +"," + 
 str(abundance_mean_SAB_deforested) + "\n")
 
 
 

 # logged

logged_abundance_SAB = np.where(resource_use_SAB == 3, abundance_data_SAB, 0)
abundance_count_SAB_logged = np.count_nonzero(logged_abundance_SAB)
abundance_sum_SAB_logged = np.sum(logged_abundance_SAB)
abundance_mean_SAB_logged = abundance_sum_SAB_logged/abundance_count_SAB_logged

resource_use.write( "logging" +  "," + "SAB" + "," + str(abundance_sum_SAB_logged) + "," +
 str(abundance_sum_SAB_logged/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_logged) +"," + 
 str(abundance_mean_SAB_logged) + "\n")
 
 # primary
primary_abundance_SAB = np.where(resource_use_SAB == 4, abundance_data_SAB, 0)
abundance_count_SAB_primary = np.count_nonzero(primary_abundance_SAB)
abundance_sum_SAB_primary = np.sum(primary_abundance_SAB)
abundance_mean_SAB_primary = abundance_sum_SAB_primary/abundance_count_SAB_primary


resource_use.write( "primary forest" +  "," + "SAB" + "," + str(abundance_sum_SAB_primary) + "," +
 str(abundance_sum_SAB_primary/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_primary) +"," + 
 str(abundance_mean_SAB_primary) + "\n")


# primary montane
primary_montane_abundance_SAB = np.where(resource_use_SAB == 5, abundance_data_SAB, 0)
abundance_count_SAB_primary_montane = np.count_nonzero(primary_montane_abundance_SAB)
abundance_sum_SAB_primary_montane = np.sum(primary_montane_abundance_SAB)
abundance_mean_SAB_primary_montane = abundance_sum_SAB_primary_montane/abundance_count_SAB_primary_montane


resource_use.write( "primary montane forest" +  "," + "SAB" + "," + str(abundance_sum_SAB_primary_montane) + "," +
 str(abundance_sum_SAB_primary_montane/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_primary_montane) +"," + 
 str(abundance_mean_SAB_primary_montane) + "\n")


 # regrowth
regrowth_abundance_SAB = np.where(resource_use_SAB == 6, abundance_data_SAB, 0)
abundance_count_SAB_regrowth = np.count_nonzero(regrowth_abundance_SAB)
abundance_sum_SAB_regrowth = np.sum(regrowth_abundance_SAB)
abundance_mean_SAB_regrowth = abundance_sum_SAB_regrowth/abundance_count_SAB_regrowth


resource_use.write( "regrowth forest" +  "," + "SAB" + "," + str(abundance_sum_SAB_regrowth) + "," +
 str(abundance_sum_SAB_regrowth/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_regrowth) +"," + 
 str(abundance_mean_SAB_regrowth) + "\n")


# old plantations
old_plantations_abundance_SAB = np.where(resource_use_SAB == 7, abundance_data_SAB, 0)
abundance_count_SAB_old_plantations = np.count_nonzero(old_plantations_abundance_SAB)
abundance_sum_SAB_old_plantations = np.sum(old_plantations_abundance_SAB)
abundance_mean_SAB_old_plantations = abundance_sum_SAB_old_plantations/abundance_count_SAB_old_plantations


resource_use.write( "old plantations" +  "," + "SAB" + "," + str(abundance_sum_SAB_old_plantations) + "," +
 str(abundance_sum_SAB_old_plantations/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_old_plantations) +"," + 
 str(abundance_mean_SAB_old_plantations) + "\n")


# nothing

other_abundance_SAB = np.where(resource_use_SAB == 8, abundance_data_SAB, 0)
abundance_count_SAB_other = np.count_nonzero(other_abundance_SAB)
abundance_sum_SAB_other = np.sum(other_abundance_SAB)
abundance_mean_SAB_other = abundance_sum_SAB_other/abundance_count_SAB_other


resource_use.write( "other" +  "," + "SAB" + "," + str(abundance_sum_SAB_other) + "," +
 str(abundance_sum_SAB_other/abundance_sum_SAB * 100 ) + "," + 
 str(abundance_count_SAB_other) +"," + 
 str(abundance_mean_SAB_other) + "\n") 
 
#######
# SAW #
#######


abundance_data_SAW = np.where(borneo_province == 11, abundance_data, 0 )
count_abundance_SAW = np.count_nonzero(abundance_data_SAW)
abundance_sum_SAW = np.sum(abundance_data_SAW)
mean_abundance_SAW = abundance_sum_SAW/count_abundance_SAW 

resource_use.write( "total abundance" + "," + "SAW" + "," + str(abundance_sum_SAW) + "," + 
  str(abundance_sum_SAW/abundance_sum_SAW * 100 ) + "," + 
  str(count_abundance_SAW) +"," + 
  str(mean_abundance_SAW) + "\n")





# planted
planted_abundance_SAW = np.where(resource_use_SAW == 1, abundance_data_SAW, 0)
abundance_count_SAW_planted = np.count_nonzero(planted_abundance_SAW)
abundance_sum_SAW_planted = np.sum(planted_abundance_SAW)
abundance_mean_SAW_planted = abundance_sum_SAW_planted / abundance_count_SAW_planted


resource_use.write( "plantation" +  "," + "SAW" + "," + str(abundance_sum_SAW_planted) + "," +
 str(abundance_sum_SAW_planted/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_planted) +"," + 
 str(abundance_mean_SAW_planted) + "\n")



# deforested
deforested_abundance_SAW = np.where(resource_use_SAW == 2, abundance_data_SAW, 0)
abundance_count_SAW_deforested = np.count_nonzero(deforested_abundance_SAW)
abundance_sum_SAW_deforested = np.sum(deforested_abundance_SAW)
abundance_mean_SAW_deforested = abundance_sum_SAW_deforested/abundance_count_SAW_deforested

resource_use.write( "deforestation" +  "," + "SAW" + "," + str(abundance_sum_SAW_deforested) + "," +
 str(abundance_sum_SAW_deforested/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_deforested) +"," + 
 str(abundance_mean_SAW_deforested) + "\n")
 
 
 

 # logged

logged_abundance_SAW = np.where(resource_use_SAW == 3, abundance_data_SAW, 0)
abundance_count_SAW_logged = np.count_nonzero(logged_abundance_SAW)
abundance_sum_SAW_logged = np.sum(logged_abundance_SAW)
abundance_mean_SAW_logged = abundance_sum_SAW_logged/abundance_count_SAW_logged

resource_use.write( "logging" +  "," + "SAW" + "," + str(abundance_sum_SAW_logged) + "," +
 str(abundance_sum_SAW_logged/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_logged) +"," + 
 str(abundance_mean_SAW_logged) + "\n")
 
 # primary
primary_abundance_SAW = np.where(resource_use_SAW == 4, abundance_data_SAW, 0)
abundance_count_SAW_primary = np.count_nonzero(primary_abundance_SAW)
abundance_sum_SAW_primary = np.sum(primary_abundance_SAW)
abundance_mean_SAW_primary = abundance_sum_SAW_primary/abundance_count_SAW_primary


resource_use.write( "primary forest" +  "," + "SAW" + "," + str(abundance_sum_SAW_primary) + "," +
 str(abundance_sum_SAW_primary/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_primary) +"," + 
 str(abundance_mean_SAW_primary) + "\n")


# primary montane
primary_montane_abundance_SAW = np.where(resource_use_SAW == 5, abundance_data_SAW, 0)
abundance_count_SAW_primary_montane = np.count_nonzero(primary_montane_abundance_SAW)
abundance_sum_SAW_primary_montane = np.sum(primary_montane_abundance_SAW)
abundance_mean_SAW_primary_montane = abundance_sum_SAW_primary_montane/abundance_count_SAW_primary_montane


resource_use.write( "primary montane forest" +  "," + "SAW" + "," + str(abundance_sum_SAW_primary_montane) + "," +
 str(abundance_sum_SAW_primary_montane/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_primary_montane) +"," + 
 str(abundance_mean_SAW_primary_montane) + "\n")


 # regrowth
regrowth_abundance_SAW = np.where(resource_use_SAW == 6, abundance_data_SAW, 0)
abundance_count_SAW_regrowth = np.count_nonzero(regrowth_abundance_SAW)
abundance_sum_SAW_regrowth = np.sum(regrowth_abundance_SAW)
abundance_mean_SAW_regrowth = abundance_sum_SAW_regrowth/abundance_count_SAW_regrowth


resource_use.write( "regrowth forest" +  "," + "SAW" + "," + str(abundance_sum_SAW_regrowth) + "," +
 str(abundance_sum_SAW_regrowth/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_regrowth) +"," + 
 str(abundance_mean_SAW_regrowth) + "\n")


# old plantations
old_plantations_abundance_SAW = np.where(resource_use_SAW == 7, abundance_data_SAW, 0)
abundance_count_SAW_old_plantations = np.count_nonzero(old_plantations_abundance_SAW)
abundance_sum_SAW_old_plantations = np.sum(old_plantations_abundance_SAW)
abundance_mean_SAW_old_plantations = abundance_sum_SAW_old_plantations/abundance_count_SAW_old_plantations


resource_use.write( "old plantations" +  "," + "SAW" + "," + str(abundance_sum_SAW_old_plantations) + "," +
 str(abundance_sum_SAW_old_plantations/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_old_plantations) +"," + 
 str(abundance_mean_SAW_old_plantations) + "\n")


# nothing

other_abundance_SAW = np.where(resource_use_SAW == 8, abundance_data_SAW, 0)
abundance_count_SAW_other = np.count_nonzero(other_abundance_SAW)
abundance_sum_SAW_other = np.sum(other_abundance_SAW)
abundance_mean_SAW_other = abundance_sum_SAW_other/abundance_count_SAW_other


resource_use.write( "other" +  "," + "SAW" + "," + str(abundance_sum_SAW_other) + "," +
 str(abundance_sum_SAW_other/abundance_sum_SAW * 100 ) + "," + 
 str(abundance_count_SAW_other) +"," + 
 str(abundance_mean_SAW_other) + "\n")
 
 
 
########
# IDN #
#######


abundance_data_IDN = np.where(borneo == 106, abundance_data, 0 )
count_abundance_IDN = np.count_nonzero(abundance_data_IDN)
abundance_data_IDN = np.where(borneo == 106, abundance_data, 0)
abundance_sum_IDN = np.sum(abundance_data_IDN)
mean_abundance_IDN = np.mean(abundance_data_IDN)

count_abundance_IDN = np.count_nonzero(abundance_data_IDN)
abundance_sum_IDN = np.sum(abundance_data_IDN)
mean_abundance_IDN = abundance_sum_IDN/count_abundance_IDN 

resource_use.write( "total abundance" + "," + "IDN" + "," + str(abundance_sum_IDN) + "," + 
  str(abundance_sum_IDN/abundance_sum_IDN * 100 ) + "," + 
  str(count_abundance_IDN) +"," + 
  str(mean_abundance_IDN) + "\n")





# planted
planted_abundance_IDN = np.where(resource_use_IDN == 1, abundance_data_IDN, 0)
abundance_count_IDN_planted = np.count_nonzero(planted_abundance_IDN)
abundance_sum_IDN_planted = np.sum(planted_abundance_IDN)
abundance_mean_IDN_planted = abundance_sum_IDN_planted / abundance_count_IDN_planted


resource_use.write( "plantation" +  "," + "IDN" + "," + str(abundance_sum_IDN_planted) + "," +
 str(abundance_sum_IDN_planted/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_planted) +"," + 
 str(abundance_mean_IDN_planted) + "\n")



# deforested
deforested_abundance_IDN = np.where(resource_use_IDN == 2, abundance_data_IDN, 0)
abundance_count_IDN_deforested = np.count_nonzero(deforested_abundance_IDN)
abundance_sum_IDN_deforested = np.sum(deforested_abundance_IDN)
abundance_mean_IDN_deforested = abundance_sum_IDN_deforested/abundance_count_IDN_deforested

resource_use.write( "deforestation" +  "," + "IDN" + "," + str(abundance_sum_IDN_deforested) + "," +
 str(abundance_sum_IDN_deforested/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_deforested) +"," + 
 str(abundance_mean_IDN_deforested) + "\n")
 
 
 

 # logged

logged_abundance_IDN = np.where(resource_use_IDN == 3, abundance_data_IDN, 0)
abundance_count_IDN_logged = np.count_nonzero(logged_abundance_IDN)
abundance_sum_IDN_logged = np.sum(logged_abundance_IDN)
abundance_mean_IDN_logged = abundance_sum_IDN_logged/abundance_count_IDN_logged

resource_use.write( "logging" +  "," + "IDN" + "," + str(abundance_sum_IDN_logged) + "," +
 str(abundance_sum_IDN_logged/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_logged) +"," + 
 str(abundance_mean_IDN_logged) + "\n")
 
 # primary
primary_abundance_IDN = np.where(resource_use_IDN == 4, abundance_data_IDN, 0)
abundance_count_IDN_primary = np.count_nonzero(primary_abundance_IDN)
abundance_sum_IDN_primary = np.sum(primary_abundance_IDN)
abundance_mean_IDN_primary = abundance_sum_IDN_primary/abundance_count_IDN_primary


resource_use.write( "primary forest" +  "," + "IDN" + "," + str(abundance_sum_IDN_primary) + "," +
 str(abundance_sum_IDN_primary/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_primary) +"," + 
 str(abundance_mean_IDN_primary) + "\n")


# primary montane
primary_montane_abundance_IDN = np.where(resource_use_IDN == 5, abundance_data_IDN, 0)
abundance_count_IDN_primary_montane = np.count_nonzero(primary_montane_abundance_IDN)
abundance_sum_IDN_primary_montane = np.sum(primary_montane_abundance_IDN)
abundance_mean_IDN_primary_montane = abundance_sum_IDN_primary_montane/abundance_count_IDN_primary_montane


resource_use.write( "primary montane forest" +  "," + "IDN" + "," + str(abundance_sum_IDN_primary_montane) + "," +
 str(abundance_sum_IDN_primary_montane/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_primary_montane) +"," + 
 str(abundance_mean_IDN_primary_montane) + "\n")


 # regrowth
regrowth_abundance_IDN = np.where(resource_use_IDN == 6, abundance_data_IDN, 0)
abundance_count_IDN_regrowth = np.count_nonzero(regrowth_abundance_IDN)
abundance_sum_IDN_regrowth = np.sum(regrowth_abundance_IDN)
abundance_mean_IDN_regrowth = abundance_sum_IDN_regrowth/abundance_count_IDN_regrowth


resource_use.write( "regrowth forest" +  "," + "IDN" + "," + str(abundance_sum_IDN_regrowth) + "," +
 str(abundance_sum_IDN_regrowth/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_regrowth) +"," + 
 str(abundance_mean_IDN_regrowth) + "\n")


# old plantations
old_plantations_abundance_IDN = np.where(resource_use_IDN == 7, abundance_data_IDN, 0)
abundance_count_IDN_old_plantations = np.count_nonzero(old_plantations_abundance_IDN)
abundance_sum_IDN_old_plantations = np.sum(old_plantations_abundance_IDN)
abundance_mean_IDN_old_plantations = abundance_sum_IDN_old_plantations/abundance_count_IDN_old_plantations


resource_use.write( "old plantations" +  "," + "IDN" + "," + str(abundance_sum_IDN_old_plantations) + "," +
 str(abundance_sum_IDN_old_plantations/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_old_plantations) +"," + 
 str(abundance_mean_IDN_old_plantations) + "\n")


# nothing

other_abundance_IDN = np.where(resource_use_IDN == 8, abundance_data_IDN, 0)
abundance_count_IDN_other = np.count_nonzero(other_abundance_IDN)
abundance_sum_IDN_other = np.sum(other_abundance_IDN)
abundance_mean_IDN_other = abundance_sum_IDN_other/abundance_count_IDN_other


resource_use.write( "other" +  "," + "IDN" + "," + str(abundance_sum_IDN_other) + "," +
 str(abundance_sum_IDN_other/abundance_sum_IDN * 100 ) + "," + 
 str(abundance_count_IDN_other) +"," + 
 str(abundance_mean_IDN_other) + "\n")
 
"""
mp.tiff.write_tif(file_with_srid = grid_layer_path, 
                   full_output_name = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/test.tif",
                   data =  abundance_data_MYS , 
                   dtype = 4)
np.max(planted_abundance)
"""

resource_use.close()
# CONTINUE WITH R SCRIPT calculating_overlap_densities_country.R

"""
mp.tiff.write_tif(file_with_srid = grid_layer_path, 
                   full_output_name = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/resource_use_MYS.tif",
                   data =  resource_use_MYS , 
                   dtype = 4)
                   
"""

year =  1999
 

if year==1999:
    abundance_path = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/abundMod/testing_ae_and_absence/pipeline_results/ppln_ae75m_50-2017-02-28T18-00-52/prediction_map_" + str(year) + "_2017-02-28_repro_res.tif"

print abundance_path
# read the data
abundance_data = mp.tiff.read_tif(abundance_path, 1)  

resource_use = mp.tiff.read_tif(grid_layer_path, 1)

abundance_data = np.where(absence == 0, abundance_data, 0)
#abundance_data = np.where(abundance_data > 0.1, abundance_data, 0)

# there are some pixels in the outside shape of borneo,
# where we have abundance but no grid, this we will clip
abundance_data = np.where((resource_use == 0) &
                           (abundance_data > 0), 0, abundance_data)
                                               
abundance_data_MYS = np.where(borneo == 136, abundance_data, 0 )
abundance_data_SAB = np.where(borneo_province == 10, abundance_data, 0 )
abundance_data_SAW = np.where(borneo_province == 11, abundance_data, 0 )
abundance_data_IDN = np.where(borneo == 106, abundance_data, 0 )

                   
file_path_table_affected = "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/OU_area_affected_" + str(year) + ".csv"
affected = open(file_path_table_affected, 'w')

affected.write("type, country, ou_affected, perc_ou_affected, area_affected, perc_area_affected \n")


# IT DOESNT MAKE A DIFFERENCE IF INCLUDING ONLY LARGER THAN 0.1

# how many OU have been affected

gone_ou_all = np.where(#(abundance_data > 0.1) & 
                ((resource_use == 1) | 
                (resource_use == 2) |
                (resource_use == 3)), abundance_data, 0)
                
 # # how much "area" was affected
total_ae = np.where(abundance_data > 0, 1, 0)
gone_ae_all = np.where(
              #  (abundance_data > 0.1) & 
                ((resource_use == 1) | 
                (resource_use == 2) |
                (resource_use == 3)), 1, 0)

affected.write( "all" +  "," + "borneo" + "," +  
str(np.sum(gone_ou_all) ) + "," + 
str(np.sum(gone_ou_all)*100/np.sum(abundance_data)) + "," + 
str(np.sum(gone_ae_all) * 0.0001) + "," + str(np.sum(gone_ae_all)*100/ np.sum(total_ae)) + "\n")
 
# how many OU have been affected
gone_ou_defor = np.where(#(abundance_data > 0.1) & 
                ((resource_use == 1) | 
                (resource_use == 2) ), abundance_data, 0)
                
 # # how much "area" was affected
gone_ae_defor = np.where(
               # (abundance_data > 0.1) & 
               ((resource_use == 1) | 
                (resource_use == 2) ), 1, 0)

                            
affected.write( "cover_change" +  "," + "borneo" + "," + 
 str(np.sum(gone_ou_defor) ) + "," +
 str(np.sum(gone_ou_defor)*100/np.sum(abundance_data)) + "," + 
 str(np.sum(gone_ae_defor)* 0.0001) + "," +
 str(np.sum(gone_ae_defor)*100/ np.sum(total_ae)) + "\n")
 
 
                
gone_ou_logging = np.where(#(abundance_data > 0.1) & 
                               (resource_use == 3), abundance_data, 0)
gone_ae_logging = np.where(#(abundance_data > 0.1) &
                               (resource_use == 3), 1, 0)


affected.write( "logging" +  "," + "borneo" + "," + 
 str(np.sum(gone_ou_logging) ) + "," +
 str(np.sum(gone_ou_logging)*100/np.sum(abundance_data)) + "," + 
 str(np.sum(gone_ae_logging)* 0.0001) + "," +
 str(np.sum(gone_ae_logging)*100/ np.sum(total_ae))+ "\n")
 
 
gone_ou_IDN_all = np.where(#(abundance_data_IDN > 0.1) & 
                ((resource_use_IDN == 1) | 
                (resource_use_IDN == 2) |
                (resource_use_IDN == 3)), abundance_data_IDN, 0)
                
 # # how much "area" was affected
total_ae_IDN = np.where(abundance_data_IDN > 0, 1, 0)
gone_ae_IDN_all = np.where(
              #  (abundance_data_IDN > 0.1) & 
                ((resource_use_IDN == 1) | 
                (resource_use_IDN == 2) |
                (resource_use_IDN == 3)), 1, 0)

affected.write( "all" +  "," + "IDN" + "," +  
str(np.sum(gone_ou_IDN_all) ) + "," + 
str(np.sum(gone_ou_IDN_all)*100/np.sum(abundance_data_IDN)) + "," + 
str(np.sum(gone_ae_IDN_all) * 0.0001) + "," + str(np.sum(gone_ae_IDN_all)*100/ np.sum(total_ae_IDN)) + "\n")
 
# how many OU have been affected
gone_ou_IDN_defor = np.where(#(abundance_data_IDN > 0.1) & 
                ((resource_use_IDN == 1) | 
                (resource_use_IDN == 2) ), abundance_data_IDN, 0)
                
 # # how much "area" was affected
gone_ae_IDN_defor = np.where(
               # (abundance_data_IDN > 0.1) & 
               ((resource_use_IDN == 1) | 
                (resource_use_IDN == 2) ), 1, 0)

                            
affected.write( "cover_change" +  "," + "IDN" + "," + 
 str(np.sum(gone_ou_IDN_defor) ) + "," +
 str(np.sum(gone_ou_IDN_defor)*100/np.sum(abundance_data_IDN)) + "," + 
 str(np.sum(gone_ae_IDN_defor)* 0.0001) + "," +
 str(np.sum(gone_ae_IDN_defor)*100/ np.sum(total_ae_IDN)) + "\n")
 
 
                
gone_ou_IDN_logging = np.where(#(abundance_data_IDN > 0.1) & 
                               (resource_use_IDN == 3), abundance_data_IDN, 0)
gone_ae_IDN_logging = np.where(#(abundance_data_IDN > 0.1) &
                               (resource_use_IDN == 3), 1, 0)


affected.write( "logging" +  "," + "IDN" + "," + 
 str(np.sum(gone_ou_IDN_logging) ) + "," +
 str(np.sum(gone_ou_IDN_logging)*100/np.sum(abundance_data_IDN)) + "," + 
 str(np.sum(gone_ae_IDN_logging)* 0.0001) + "," +
 str(np.sum(gone_ae_IDN_logging)*100/ np.sum(total_ae_IDN))+ "\n")
 
                
gone_ou_MYS_all = np.where(#(abundance_data_MYS > 0.1) & 
                ((resource_use_MYS == 1) | 
                (resource_use_MYS == 2) |
                (resource_use_MYS == 3)), abundance_data_MYS, 0)
                
 # # how much "area" was affected
total_ae_MYS = np.where(abundance_data_MYS > 0, 1, 0)
gone_ae_MYS_all = np.where(
               # (abundance_data_MYS > 0.1) & 
                ((resource_use_MYS == 1) | 
                (resource_use_MYS == 2) |
                (resource_use_MYS == 3)), 1, 0)

affected.write( "all" +  "," + "MYS" + "," + 
 str(np.sum(gone_ou_MYS_all) ) + "," +
 str(np.sum(gone_ou_MYS_all)*100/np.sum(abundance_data_MYS)) + "," + 
 str(np.sum(gone_ae_MYS_all) * 0.0001) + "," +
 str(np.sum(gone_ae_MYS_all)*100/ np.sum(total_ae_MYS)) + "\n")
 
# how many OU have been affected
gone_ou_MYS_defor = np.where(#(abundance_data_MYS > 0.1) & 
                ((resource_use_MYS == 1) | 
                (resource_use_MYS == 2)), abundance_data_MYS, 0)
                
 # # how much "area" was affected
gone_ae_MYS_defor = np.where(
               # (abundance_data_MYS > 0.1) & 
               ((resource_use_MYS == 1) | 
                (resource_use_MYS == 2)), 1, 0)

                            
affected.write( "cover_change" +  "," + "MYS" + "," + 
 str(np.sum(gone_ou_MYS_defor) ) + "," +
 str(np.sum(gone_ou_MYS_defor)*100/np.sum(abundance_data_MYS)) + "," + 
 str(np.sum(gone_ae_MYS_defor) * 0.0001) + "," +
 str(np.sum(gone_ae_MYS_defor)*100/ np.sum(total_ae_MYS)) + "\n")
 
 
                
gone_ou_MYS_logging = np.where(#(abundance_data_MYS > 0.1) & 
                               (resource_use_MYS == 3), abundance_data_MYS, 0)
gone_ae_MYS_logging = np.where(#(abundance_data_MYS > 0.1) &
                               (resource_use_MYS == 3), 1, 0)


affected.write( "logging" +  "," + "MYS" + "," + 
 str(np.sum(gone_ou_MYS_logging) ) + "," +
 str(np.sum(gone_ou_MYS_logging)*100/np.sum(abundance_data_MYS)) + "," + 
 str(np.sum(gone_ae_MYS_logging) * 0.0001) + "," +
 str(np.sum(gone_ae_MYS_logging)*100/ np.sum(total_ae_MYS)) + "\n")    

# SAB

gone_ou_SAB_all = np.where(#(abundance_data_SAB > 0.1) & 
                ((resource_use_SAB == 1) | 
                (resource_use_SAB == 2) |
                (resource_use_SAB == 3)), abundance_data_SAB, 0)
                
 # # how much "area" was affected
total_ae_SAB = np.where(abundance_data_SAB > 0, 1, 0)
gone_ae_SAB_all = np.where(
               # (abundance_data_SAB > 0.1) & 
                ((resource_use_SAB == 1) | 
                (resource_use_SAB == 2) |
                (resource_use_SAB == 3)), 1, 0)

affected.write( "all" +  "," + "SAB" + "," + 
 str(np.sum(gone_ou_SAB_all) ) + "," +
 str(np.sum(gone_ou_SAB_all)*100/np.sum(abundance_data_SAB)) + "," + 
 str(np.sum(gone_ae_SAB_all) * 0.0001) + "," +
 str(np.sum(gone_ae_SAB_all)*100/ np.sum(total_ae_SAB)) + "\n")
 
# how many OU have been affected
gone_ou_SAB_defor = np.where(#(abundance_data_SAB > 0.1) & 
                ((resource_use_SAB == 1) | 
                (resource_use_SAB == 2)), abundance_data_SAB, 0)
                
 # # how much "area" was affected
gone_ae_SAB_defor = np.where(
               # (abundance_data_SAB > 0.1) & 
               ((resource_use_SAB == 1) | 
                (resource_use_SAB == 2)), 1, 0)

                            
affected.write( "cover_change" +  "," + "SAB" + "," + 
 str(np.sum(gone_ou_SAB_defor) ) + "," +
 str(np.sum(gone_ou_SAB_defor)*100/np.sum(abundance_data_SAB)) + "," + 
 str(np.sum(gone_ae_SAB_defor) * 0.0001) + "," +
 str(np.sum(gone_ae_SAB_defor)*100/ np.sum(total_ae_SAB)) + "\n")
 
 
                
gone_ou_SAB_logging = np.where(#(abundance_data_SAB > 0.1) & 
                               (resource_use_SAB == 3), abundance_data_SAB, 0)
gone_ae_SAB_logging = np.where(#(abundance_data_SAB > 0.1) &
                               (resource_use_SAB == 3), 1, 0)


affected.write( "logging" +  "," + "SAB" + "," + 
 str(np.sum(gone_ou_SAB_logging) ) + "," +
 str(np.sum(gone_ou_SAB_logging)*100/np.sum(abundance_data_SAB)) + "," + 
 str(np.sum(gone_ae_SAB_logging) * 0.0001) + "," +
 str(np.sum(gone_ae_SAB_logging)*100/ np.sum(total_ae_SAB)) + "\n")    

gone_ou_SAW_all = np.where(#(abundance_data_SAW > 0.1) & 
                ((resource_use_SAW == 1) | 
                (resource_use_SAW == 2) |
                (resource_use_SAW == 3)), abundance_data_SAW, 0)
                
 # # how much "area" was affected
total_ae_SAW = np.where(abundance_data_SAW > 0, 1, 0)
gone_ae_SAW_all = np.where(
               # (abundance_data_SAW > 0.1) & 
                ((resource_use_SAW == 1) | 
                (resource_use_SAW == 2) |
                (resource_use_SAW == 3)), 1, 0)

affected.write( "all" +  "," + "SAW" + "," + 
 str(np.sum(gone_ou_SAW_all) ) + "," +
 str(np.sum(gone_ou_SAW_all)*100/np.sum(abundance_data_SAW)) + "," + 
 str(np.sum(gone_ae_SAW_all) * 0.0001) + "," +
 str(np.sum(gone_ae_SAW_all)*100/ np.sum(total_ae_SAW)) + "\n")
 
# how many OU have been affected
gone_ou_SAW_defor = np.where(#(abundance_data_SAW > 0.1) & 
                ((resource_use_SAW == 1) | 
                (resource_use_SAW == 2)), abundance_data_SAW, 0)
                
 # # how much "area" was affected
gone_ae_SAW_defor = np.where(
               # (abundance_data_SAW > 0.1) & 
               ((resource_use_SAW == 1) | 
                (resource_use_SAW == 2)), 1, 0)

                            
affected.write( "cover_change" +  "," + "SAW" + "," + 
 str(np.sum(gone_ou_SAW_defor) ) + "," +
 str(np.sum(gone_ou_SAW_defor)*100/np.sum(abundance_data_SAW)) + "," + 
 str(np.sum(gone_ae_SAW_defor) * 0.0001) + "," +
 str(np.sum(gone_ae_SAW_defor)*100/ np.sum(total_ae_SAW)) + "\n")
 
 
                
gone_ou_SAW_logging = np.where(#(abundance_data_SAW > 0.1) & 
                               (resource_use_SAW == 3), abundance_data_SAW, 0)
gone_ae_SAW_logging = np.where(#(abundance_data_SAW > 0.1) &
                               (resource_use_SAW == 3), 1, 0)


affected.write( "logging" +  "," + "SAW" + "," + 
 str(np.sum(gone_ou_SAW_logging) ) + "," +
 str(np.sum(gone_ou_SAW_logging)*100/np.sum(abundance_data_SAW)) + "," + 
 str(np.sum(gone_ae_SAW_logging) * 0.0001) + "," +
 str(np.sum(gone_ae_SAW_logging)*100/ np.sum(total_ae_SAW)) + "\n")    




affected.close()

