##### Iteration 3.5 #####
#### Extract Data and Run Through "phenofit" ####

library(terra)
library(sf)
library(tmap)
library(tidyverse)
library(furrr)
library(data.table)
library(phenofit)


# ##### Single Point #####
# 
# # Set the directory where your HDF files are located
# hdf_folder <- "Project/Data/MCD15A3H/modis/hdf/subset"
# 
# # List all HDF files in the directory
# hdf_files <- list.files(hdf_folder, pattern = "\\.hdf$", full.names = TRUE)
# 
# # Specify the year or set of months to subset
# # Example: subset by year
# subset_year <- ""
# 
# # Load the spatial vector point 'rbp'
# rbp <- st_read("../../../Red Bluff/GIS/Layers/Fences/Pastures Current/rb_pastures_current.shp") %>%
#   dplyr::filter(Name%in%"Rye") %>%
#   st_centroid() %>%
#   # vect() %>%
#   st_transform("PROJCRS[\"unnamed\",\n    BASEGEOGCRS[\"Unknown datum based upon the custom spheroid\",\n        DATUM[\"Not specified (based on custom spheroid)\",\n            ELLIPSOID[\"Custom spheroid\",6371007.181,0,\n                LENGTHUNIT[\"metre\",1,\n                    ID[\"EPSG\",9001]]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Sinusoidal\",\n        METHOD[\"Sinusoidal\"],\n        PARAMETER[\"Longitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"Meter\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"Meter\",1]]]")
# rbp
# 
# # Initialize an empty list to store dataframe objects
# dataframe_list <- list()
# 
# # Initialize a vector to store the dates
# dates <- character()
# 
# # Loop through each HDF file
# for (file in hdf_files) {
#   # Extract the year and julian day information from the file name
#   file_parts <- unlist(strsplit(basename(file), "\\."))
#   year_jday <- substr(file_parts[2], 2, 8)
#   year <- substr(year_jday, 1, 4)
#   # Check if the file matches the specified year or set of months
#   if (grepl(subset_year, year)) {
#     # Load the raster from HDF file
#     raster <- rast(file)
#     # Find the cell containing the spatial vector point 'rbp'
#     cell <- cellFromXY(raster, st_coordinates(rbp))
#     # Extract values from layers 2 and 4 at the spatial vector point 'rbp'
#     values_layer2 <- terra::extract(raster[[2]], cell)
#     values_layer4 <- terra::extract(raster[[4]], cell)
#     # Get the date in yyyy-mm-dd format
#     date <- as.Date(year_jday, format = "%Y%j")
#     # Combine values and date into a dataframe
#     df <- data.frame(date = rep(date, length(values_layer2)),
#                      y = values_layer2,
#                      qc = values_layer4)
#     # Append the dataframe to the dataframe list
#     dataframe_list[[length(dataframe_list) + 1]] <- df
#     # Keep track of the dates
#     dates <- c(dates, date)
#   }
# }
# 
# # Combine all dataframes into a single dataframe
# combined_df <- do.call(rbind, dataframe_list)
# 
# # Now 'combined_df' is your dataframe where each row represents a unique observation
# # with columns 'date', 'y' (layer 2), and 'qc' (layer 4)
# 
# # head(combined_df)
# 
# # ggplot(combined_df,aes(x=date,y=Lai_500m))+
# #   geom_line()+
# #   scale_x_date(breaks = "1 year")
# 
# # write_csv(combined_df, "Project/Data/Output/Sample/rb_vd_lai.csv")
# 
# read_csv("Project/Data/Output/Sample/rb_vd_lai.csv")
# 
# # Arrange data and assign weights
# d <- combined_df  %>% setDT() %>%
#   cbind(combined_df[, as.list(qc_FparLai(FparExtra_QC))]) %>%
#   .[, .(date, y = Lai_500m, QC_flag, w)] %>%
#   mutate(t = date) %>%
#   filter(!year(date)==2002&!year(date)==2024)
# d
# 
# 
# # Parameters - need to change?
# # lambda         <- 8
# nptperyear     <- 92
# # minExtendMonth <- 0.5
# # maxExtendMonth <- 1
# minPercValid   <- 0
# wFUN           <- wTSM # wBisquare
# wmin           <- 0.2
# methods_fine <- c("AG", "Zhang", "Beck", "Elmore", "Gu")
# 
# 
# INPUT <- check_input(d$t, d$y, d$w,
#                      QC_flag = d$QC_flag,
#                      nptperyear = nptperyear,
#                      maxgap = nptperyear / 4, wmin = 0.2
# )
# 
# brks <- season_mov(INPUT,
#                    list(FUN = "smooth_wWHIT", wFUN = wFUN,
#                         # maxExtendMonth = 3,
#                         wmin = wmin, r_min = 0.1
#                    ))
# # plot_season(INPUT, brks)
# 
# ## 2.4 Curve fitting
# fit <- curvefits(INPUT, brks,
#                  list(
#                    methods = methods_fine, # ,"klos",, 'Gu'
#                    wFUN = wFUN,
#                    iters = 2,
#                    wmin = wmin,
#                    # constrain = FALSE,
#                    nextend = 2,
#                    # maxExtendMonth = maxExtendMonth, minExtendMonth = minExtendMonth,
#                    minPercValid = minPercValid
#                  ))
# 
# ## check the curve fitting parameters
# # l_param <- get_param(fit)
# # print(l_param$Beck)
# # 
# # dfit <- get_fitting(fit)
# # print(dfit)
# 
# ## 2.5 Extract phenology
# TRS <- c(0.1, 0.2, 0.5)
# # l_pheno <- get_pheno(fit, TRS = TRS, IsPlot = T) # %>% map(~melt_list(., "meth"))
# # print(l_pheno$doy$Beck)
# # 
# # pheno <- l_pheno$doy %>% melt_list("meth")
# # 
# # 
# # ggplot(pheno,aes(group=meth,col=meth))+
# #   geom_point(aes(x=Maturity,y=DER.pos))+
# #   geom_abline(intercept = 0, slope = 1)
# 
# ## AG seems like the best overall curvefit method
# 
# #### Visualization
# 
# # growing season dividing
# # plot_season(INPUT, brks, ylab = "LAI")
# 
# 
# # Ipaper::write_fig({  }, "Figure4_seasons.pdf", 9, 4)
# 
# # fine curvefitting
# # g <- plot_curvefits(dfit, brks, title = NULL, cex = 1.5, ylab = "EVI", angle = 0)
# # grid::grid.newpage()
# # grid::grid.draw(g)
# 
# # Ipaper::write_fig(g, "Figure5_curvefitting.pdf", 8, 6, show = TRUE)
# 
# 
# # extract phenology metrics, only the first 3 year showed at here
# # write_fig({
#   l_pheno <- get_pheno(fit, method = "AG",
#                        TRS = TRS, IsPlot = F, show.title = F)
# # }, "Project/Data/Output/phenology_metrics.png", 10, 8, show = TRUE)
# l_pheno$doy$AG
# 

##### Polygon #####

# Set the directory where your HDF files are located
hdf_folder <- "C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/MCD15A3H/modis/hdf"

# List all HDF files in the directory

# Year Subset
# year <- 2003
# hdf_files <- list.files(hdf_folder, pattern = paste0("MCD15A3H.A", year, ".*\\.hdf$"), full.names = TRUE)

# Load all
hdf_files <- list.files(hdf_folder, pattern = paste0("MCD15A3H.A", ".*\\.hdf$"), full.names = TRUE)


# Retrieve Data from File Name
file_info <- lapply(hdf_files, function(x) {
  file_parts <- unlist(strsplit(basename(x), "\\."))
  year_jday <- substr(file_parts[2], 2, 8)
  year <- substr(year_jday, 1, 4)
  date <- as.Date(year_jday, format = "%Y%j")
  tibble(year,date)
  }) %>%
  bind_rows(., .id = "file")


# Load the AOI spatial polygon 'rbp'
rbp <- st_read("C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/mt_madison_n.shp") %>%
  # dplyr::filter(NAME%in%"DEER LODGE") %>%
  # st_centroid() %>%
  # vect() %>%
  # st_union() %>%
  # st_convex_hull() %>%
  st_transform("PROJCRS[\"unnamed\",\n    BASEGEOGCRS[\"Unknown datum based upon the custom spheroid\",\n        DATUM[\"Not specified (based on custom spheroid)\",\n            ELLIPSOID[\"Custom spheroid\",6371007.181,0,\n                LENGTHUNIT[\"metre\",1,\n                    ID[\"EPSG\",9001]]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433,\n                ID[\"EPSG\",9122]]]],\n    CONVERSION[\"Sinusoidal\",\n        METHOD[\"Sinusoidal\"],\n        PARAMETER[\"Longitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"False easting\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"Meter\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"Meter\",1]]]") %>%
  vect()
rbp


# Generate raster to use as a structure
template_raster <- rast(hdf_files[1]) %>%
  terra::crop(ext(rbp))


## Method 2
# Unknown error when steps combined.
# Further debugging needed.
# "Error in h(simpleError(msg, call)) :      
#   error in evaluating the argument 'x' in selecting a method for function 'as.list': [crop] file exists. You can use 'overwrite=TRUE' to overwrite it"

# AOI Extent
extent <- ext(rbp)

# Load all files
r_sds <- sds(hdf_files) %>%
  as.list()

# Crop and Subset files to variables of interest
process_raster <- function(x) {
  terra::crop(x, extent) %>%
    terra::subset(.,subset = c("Lai_500m","FparExtra_QC")) %>%
    as.data.frame(cells=F)
}

# Each list item is one file
r_list <- map(r_sds, process_raster)

# terra::crop(r_list, extent) %>%
# as.list() %>%
# lapply(., function(x) terra::subset(x,subset = c("Lai_500m","FparExtra_QC")) %>%
#          as.data.frame(cells=F))

seq_indices <- seq_along(r_list)

# Collapse list into dataframe with cell ID annotations
r_df <- map_df(.x = seq_indices,
               .f = function(x) {
                 bind_rows(r_list[x]) %>%
                   mutate(cell = row_number(),
                          file = as.character(x))
               }) %>%
  arrange(cell)

# Join file information with cell data
d <- full_join(r_df, file_info, by = "file") %>%
  filter(date > mdy("01-01-2003") & date < mdy("01-01-2024")) %>%
  setDT() %>%
  # Extract quality controld data from "FparExtra_QC"
  cbind(.[, as.list(qc_FparLai(FparExtra_QC))]) %>%
  .[, .(t = date, y = Lai_500m, QC_flag, w, cell)] %>%
  split(by = "cell")


# Parameters - need to change?
# lambda         <- 8
nptperyear     <- 92
# minExtendMonth <- 0.5
# maxExtendMonth <- 1
minPercValid   <- 0
wFUN           <- wTSM # wBisquare
wmin           <- 0.2
methods_fine <- c("AG", "Zhang", "Beck", "Elmore") # Methods to calculate phenological metrics


#### Using furrr: purrr and future ####

# Set up parallel processing

# This seems to be the upper limit to use most of the total CPU capacity without crashing it
plan(multisession, workers = 18)

# Define the function to be applied
# t = date
# y = vegetation index
# w = weight
process_data <- function(df) {
  
  # Replace NA values in df$w with 1
  replaced_indices <- which(is.na(df$w))
  if (length(replaced_indices) > 0) {
    # Output warning message indicating replaced NAs and their locations
    message("Warning: NAs in df$w were replaced with 0.5 at indices ", replaced_indices) # This line needs work
    df$w[replaced_indices] <- 1
  }
  
  # Check input data
  INPUT <- check_input(df$t, df$y, df$w,
                       QC_flag = df$QC_flag,
                       nptperyear = nptperyear,
                       maxgap = nptperyear / 4, wmin = 0.2
  )
  
  # Define growing season breaks
  brks <- season_mov(INPUT,
                     list(FUN = "smooth_wWHIT", wFUN = wFUN,
                          # maxExtendMonth = 3,
                          wmin = wmin, r_min = 0.1
                     ))
  
  # Fit curves to data
  fit <- curvefits(INPUT, brks,
                   list(
                     methods = methods_fine,
                     wFUN = wFUN,
                     iters = 2,
                     wmin = wmin,
                     nextend = 2,
                     minPercValid = minPercValid
                   ))
  
  # Define threshold extraction points
  TRS <- c(0.1, 0.2, 0.5)
  
  l_pheno <- get_pheno(fit, method = "AG",
                       TRS = TRS, IsPlot = FALSE, show.title = FALSE)
  l_pheno$doy$AG
}

# Apply the function in parallel - this is the most time consuming part of the process
pheno_list <- future_map(d, process_data, .options = furrr_options(seed = TRUE), .progress = TRUE)



## Summarise output data

pheno_summ <- lapply(pheno_list, function(dat) {
  dplyr::select(dat,3:19) %>%
    mutate_all(~ifelse(. <= 0, NA, .)) %>%
    summarise_all(list(mean=mean,sd=sd), na.rm = T)
})

# Combine data frames into a single data frame
combined_df <- rbindlist(pheno_summ, idcol = NULL)


## Create Raster of summarized data

# Create an empty raster with the same extent and resolution as the template raster
empty_raster <- rast(ext(template_raster), resolution=res(template_raster),
                     crs = crs(template_raster),
                     nlyrs = length(pheno_summ[[1]]))
# Get the column names from the first element of pheno_list
layer_names <- names(pheno_summ[[1]])

# Set the layer names of the empty raster
names(empty_raster) <- layer_names

# Fill with data
filled_raster <- setValues(empty_raster, combined_df)

# Project the raster
proj_ras <- filled_raster %>% project("EPSG:26912")

# # Save output!!
# terra::writeRaster(filled_raster, "C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/Output/madison_valley_output.tif")
# 
# terra::writeRaster(proj_ras, "C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/Output/madison_valley_output_projected.tif")

# Test if saved
# test_read <- rast("C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/Output/madison_valley_output.tif")


## Visualize Data

# Load supporting layers
mad_pop <- st_read("C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/madison_select_pop_places.shp") %>%
  filter(Name %in% "Ennis")

county <- st_read("C:/Users/noahg/OneDrive - Montana State University (1)/Noah_Davis/Research/GIS Layers/Administrative Boundaries/MontanaCounties_shp/County.shp") %>%
  st_crop()

rb_past <- st_read("C:/Users/noahg/OneDrive - Montana State University (1)/Red Bluff/GIS/Layers/Fences/Pastures Current/rb_pastures_current.shp")

ras_mean <- proj_ras[[1:17]] %>%
  setNames(c("Threshold Start of Season 0.1", "Threshold End of Season 0.1", "Threshold Start of Season 0.2", "Threshold End of Season 0.2", "Threshold Start of Season 0.5", "Threshold End of Season 0.5", "Derivative Start of Season", "Derivative Peak of Season", "Derivative End of Season", "Upturn Date", "Senescence Date", "Downturn Date", "Regression Date", "Greenup", "Maturity", "Senescence", "Dormancy"))

ras_sd <- proj_ras[[18:34]] %>%
  setNames(c("Threshold Start of Season 0.1", "Threshold End of Season 0.1", "Threshold Start of Season 0.2", "Threshold End of Season 0.2", "Threshold Start of Season 0.5", "Threshold End of Season 0.5", "Derivative Start of Season", "Derivative Peak of Season", "Derivative End of Season", "Upturn Date", "Senescence Date", "Downturn Date", "Regression Date", "Greenup", "Maturity", "Senescence", "Dormancy"))


# Generate figures
tmap_mean <- tm_shape(ras_mean[[c(1,3,5,6,4,2,7,8,9,10,11,12,13,14,15,16,17)]])+tm_raster(palette = "-Spectral",
                                                                                          breaks = c(0,31,59,90,120,151,181,212,243,273,304,334,365),
                                                                                          labels = month.abb,
                                                                                          title = "Mean Month")+
  tm_shape(mad_pop)+tm_dots(size = 0.35)+
  tm_shape(county)+tm_borders()

tmap_sd <- tm_shape(ras_sd[[c(1,3,5,6,4,2,7,8,9,10,11,12,13,14,15,16,17)]])+tm_raster(palette = "-Spectral", title = "Standard Deviation, days")+
  tm_shape(mad_pop)+tm_dots(size = 0.35)+
  tm_shape(county)+tm_borders()

tmap_save(tmap_mean, "C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/Output/madison_mean.jpg",
          height = 4, width = 6, units = "in")

tmap_save(tmap_sd, "C:/Users/noahg/OneDrive - Montana State University (1)/Classes/Spring 2024/LRES_525/Project/Data/Output/madison_sd.jpg",
          height = 4, width = 6, units = "in")

# tm_shape(proj_ras[[c(1,3,5, 6,4,2)]])+tm_raster(palette = "-Spectral", breaks = c(0,31,59,90,120,151,181,212,243,273,304,334,365), labels = month.abb)+
#   # tm_shape(mad_pop)+tm_dots(size = 0.5, labels = c("Norris","Ennis","Twin Bridges"))+
#   tm_shape(county)+tm_borders()+
#   tm_facets(nrow = 2)
# 
# tm_shape(proj_ras[[c(18,20,22, 23,21,19)]])+tm_raster(palette = "-Reds")+
#   tm_shape(mad_pop)+tm_dots(labels = c("Norris","Ennis","Twin Bridges"))+
#   tm_facets(nrow = 2)
# 
# tm_shape(proj_ras[[c(1,3,5,6,4,2)]])+tm_raster(palette = "Greens", breaks = c(0,31,59,90,120,151,181,212,243,273,304,334,365), labels = month.abb)+
#   tm_shape(mad_pop)+tm_dots(labels = c("Norris","Ennis","Twin Bridges"))+
#   tm_facets(nrow = 2)
# 
# tm_shape(proj_ras[[18:34]])+tm_raster(palette = "Spectral", style = "cont")+
#   tm_shape(rb_past)+tm_borders()