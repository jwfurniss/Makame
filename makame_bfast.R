
packages<-(c("sf", "gdalcubes", "rstac", "terra",
             "torch", "dplyr", "dismo", "exactextractr", 
             "ape", "reshape2","ggplot2", "parallel"))

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

rm(packages)
rm(installed_packages)
## General ## 
setwd("D://Makame") # Set working directory to whatever folder you would like where data can be written and retrieved

gdalcubes_options(parallel = detectCores()-2) # Option to enable parallel processing for gdalcubes

output.dir<-paste("./Bfast") # create a directory to store files 
if (!file.exists(output.dir)) {
  dir.create(output.dir)
}

#sf_pa<-st_read("ForestBenchmark_2022_LeakageBelt.shp")  #Define Project Area
sf_pa<-st_read("D://TerraCarbon_DropBox//TerraCarbon Dropbox//John Furniss//Makame REDD Project Folder//Makame Internal Folder//Makame GIS data//Monitoring_2023//MK_2022_PALB_F.shp")

## Image Prep Input ##

## Short Dry
#start_date = "2023-01-01"
#end_date = "2023-03-31"

## Long Dry
start_date = "2020-01-01"
end_date = "2023-08-03"  

max_cc = 30
band_list=c("B04","B08","SCL")


sf_pa<-st_transform(sf_pa$geometry, 4326) #transform PA to lat/long for search purposes
sf_pa_buffer<-st_buffer(st_zm(sf_pa), 10, 10) #### 10m buffer around LB

roi<-st_bbox(sf_pa_buffer) #get Bounding Box of PA Buffer area
roi_list<-c(roi[[1]], roi[[2]], roi[[3]], roi[[4]]) # transform Bounding Box to list of coordinates


s_obj <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") #connect to Microsoft Planetary Computer API
it_obj <- s_obj %>%
  stac_search(collections = "sentinel-2-l2a",
              bbox = roi_list,
              datetime = paste0(as.character(start_date),"/",as.character(end_date), sep = ""),
              limit = 1000) # sets limit to number of images called (default is 500)

max_cc<-as.numeric(max_cc)
it_obj$params$query = paste0("{\"eo:cloud_cover\": {\"lt\": ",max_cc,"}}",sep="")#limits images to under X% cloud cover
stac <- NULL
attempt <- 0
while( is.null(stac) && attempt <= 20 ) {
  attempt <- attempt + 1 #this has never actually timed out, but I figured an escape hatch would be helpful if something weird happens (like MPC going down)
  try(
    stac <- it_obj %>% get_request() %>%
      items_sign(sign_fn = sign_planetary_computer()), #actual call to get images
    silent = TRUE)
} 
#Create Image Collection from stac
col<-stac_image_collection(stac$features, 
                           out_file = tempfile(fileext = ".sqlite"),
                           asset_names = band_list,
                           property_filter = function(x) {x[["datetime"]] !="2022-01-28T07:40:49.024000Z"}) #convert stac to image collection, keep qa_pixel band for cloud masking, remove any images from Landsat-7


col.prj<<-stac$features[[1]]$properties$`proj:epsg` #function to call projection from first tile to create data cube
#get image collection extent and projection information
col.epsg<<-paste0("EPSG:",stac$features[[1]]$properties$`proj:epsg`, sep = "") #function to call projection from first tile to create data cube
col.ext<<-gdalcubes::extent(col, srs = col.epsg) #get image collection spatial and temporal extent
#Project Project Area buffer to match image collection for Extent
sf_pa_prj<-sf::st_transform(sf_pa, crs = col.prj)
sf_pa_buffer_prj<-sf::st_transform(sf_pa_buffer, crs = col.prj)
roi_prj<<-sf::st_bbox(sf_pa_buffer_prj) #get Bounding Box of Projected Donor Area
roi_ext<<-list(left = roi_prj[[1]], 
               right = roi_prj[[3]], 
               top = roi_prj[[4]], 
               bottom = roi_prj[[2]], 
               t0=substr(col.ext$t0,0,10), 
               t1=substr(col.ext$t1,0,10))


roi_ext$t0 = start_date
roi_ext$t1 = end_date
v.overview = cube_view(extent=roi_ext, dx = 30, dy = 30, dt="P1M", srs = col.epsg,
                       aggregation = "mean", resampling = "bilinear")


cloud_mask = image_mask("SCL", values =c(3,8,9)) # Mask Sentinel clouds  # Mask Landsat clouds 

#Function to create data cube and write .tif files to output directory
#Need final pipe function as write_tif to make X into list of file paths, which we can export into raster stack

start.time.stac<-Sys.time()
print(start.time.stac)

x = raster_cube(col, v.overview, mask = cloud_mask)|>
  select_bands(c("B04","B08"))|>
  filter_geom(sf_pa, "EPSG:4326")|>
  reduce_time(names = c("change_date", "change_magnitude"), FUN = function(x) {
    knr <- exp(-((x["B08",]/10000)-(x["B04",]/10000))^2/(2))
    kndvi <- (1-knr) / (1+knr)   
    if (all(is.na(kndvi))) {
      return(c(NA,NA))
    }
    kndvi_ts = ts(kndvi, start = c(2020, 1), frequency = 12)
    library(bfast)
    tryCatch({
      result = bfastmonitor(kndvi_ts, start = c(2022,7), level = 0.01)
      return(c(result$breakpoint, result$magnitude))
    }, error = function(x) {
      return(c(NA,NA))
    })
  }) |>
  write_tif(dir = output.dir,prefix = "2023_Makame_Bfast_kndvi_30m_30cc_")

end.time.stac<-Sys.time()
print(end.time.stac-start.time.stac)
plotRGB(rast(x), r=4, g=3, b = 2)
