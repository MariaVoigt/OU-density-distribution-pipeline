rm(list = ls())
source("/homes/mv39zilo/projects/orangutan_density_distribution/src/preload/preload.R")
source("/homes/mv39zilo/projects/orangutan_density_distribution/src/functions/generic/path.to.current.R")
source("/homes/mv39zilo/projects/orangutan_density_distribution/src/functions/project_functions/values.to.tif.R")
outdir <- "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/future/"
# for the future we need a raster with the grid_ids
geography_path <- "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/abundMod/testing_ae_and_absence/pipeline_results/ppln_ae75m_50-2017-02-28T18-00-52/geography_2015_2017-02-28.rds"
geography_2015 <- readRDS(geography_path)%>%
  mutate(grid_id = id, id = 1:nrow(.))%>%
  dplyr::select(id, x_start, y_start, grid_id)


#--------------#
# resource use #
#--------------#
outdir <- "/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/"

resource_use <- raster("/homes/mv39zilo/work/Borneo/analysis/model_prep_and_running/results/resource_use/resource_grid_absence_country_category_id.tif")
resource_use_grid <- as.data.frame(resource_use)
names(resource_use_grid) <- "category"
resource_use_grid$category <- as.character(resource_use_grid$category)
resource_use_grid$population <- as.numeric(substr(resource_use_grid$category, start= 1, stop = 1))
resource_use_grid$country <- as.numeric(substr(resource_use_grid$category, start= 2, stop = 2))
resource_use_grid$province <- as.numeric(substr(resource_use_grid$category, start= 3, stop = 4))
resource_use_grid$resource_use <- as.numeric(substr(resource_use_grid$category, start= 5, stop =5))
resource_use_grid$grid_id <- as.numeric(substr(resource_use_grid$category, start= 6, stop = nchar(resource_use_grid$category[1])))

summary(resource_use_grid)

resource_use_grid[resource_use_grid$country == 3, "country"] <- "MYS"
resource_use_grid[resource_use_grid$country == 6, "country"] <- "IDN"

# for now not distinguishing Kalimantan
resource_use_grid[resource_use_grid$country == 6, "province"] <- "KAL"
resource_use_grid[resource_use_grid$province == 10, "province"] <- "SAB"
resource_use_grid[resource_use_grid$province == 11, "province"] <- "SAW"

summary(resource_use_grid)
plot(table(resource_use_grid$resource_use))

# # testing duplicates
# test <- resource_use_grid
# test <- dplyr::filter(test, grid_id != 0)
# testtest <- test[duplicated(test$grid_id) & test$grid_id != 0, ]

# we want a category with the pop number for all pixels with
# forest
# and the grid_id in order
resource_use_grid[resource_use_grid$population == 2, "resource_use"] <- NA

resource_use_grid[resource_use_grid$grid_id == 0 &
                    !is.na(resource_use_grid$resource_use), "resource_use"] <- NA

resource_use_grid[resource_use_grid$resource_use == 0 &
                    !is.na(resource_use_grid$resource_use), "resource_use"] <- NA

resource_use_grid <- resource_use_grid %>%
  dplyr::select(resource_use, country,  grid_id)

# SORT WITH GRID ID
resource_use_grid <- resource_use_grid %>%
  dplyr::select( resource_use, country, grid_id) %>%
  arrange(grid_id) %>%
  filter(grid_id != 0)

grid_ids <- resource_use_grid[resource_use_grid$grid_id != 0, "grid_id"]

# DEAL WITH NON UNIQUE GRID-IDS, MAYBE PROBLEM WITH REPROJECTION OR SOMETHING
# not_unique <- grid_ids[duplicated(grid_ids)]
# resource_use_grid_not_unique <- resource_use_grid[resource_use_grid$grid_id %in% not_unique,]
grid_ids_missing <- geography_2015$grid_id[!geography_2015$grid_id %in% grid_ids]
grid_ids_missing <- as.data.frame(cbind(rep(NA, times = length(grid_ids_missing)),
                                        rep(NA, times = length(grid_ids_missing)),
                                        grid_ids_missing))
names(grid_ids_missing) <- c("resource_use", "country", "grid_id")

resource_use_grid <- resource_use_grid %>%
  rbind(grid_ids_missing) %>%
  arrange(grid_id)

resource_use_grid[!is.na(resource_use_grid$resource_use) &
                    resource_use_grid$country != 0, "category"] <-
   paste0(resource_use_grid[!is.na(resource_use_grid$resource_use)&
                              resource_use_grid$country != 0, "resource_use"],
          "_",
          resource_use_grid[!is.na(resource_use_grid$resource_use)&
                              resource_use_grid$country != 0, "country"])
resource_use_grid$category <- as.factor(resource_use_grid$category)
summary(resource_use_grid)
unique(resource_use_grid$category)
resource_use_grid <- dplyr::select(resource_use_grid, category,
                                   resource_use, country, grid_id)

write.csv(resource_use_grid, file.path(outdir, "resource_use_grid.csv" ), row.names = F)
