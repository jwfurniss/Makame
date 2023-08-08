
packages<-(c("sf", "gdalcubes", "rstac", "terra",
             "torch", "dplyr", "dismo", "exactextractr", 
             "reshape2","ggplot2", "parallel",))

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

rm(packages)
rm(installed_packages)
## General ## 
setwd("C://GIS//Tanzania//Makame") # Set working directory to whatever folder you would like where data can be written and retrieved

gdalcubes_options(parallel = detectCores()-2) # Option to enable parallel processing for gdalcubes

output.dir<-paste("./Sent2_bfast") # create a directory to store files 
if (!file.exists(output.dir)) {
  dir.create(output.dir)
}

sf_pa<-st_read("MK_LB_BBox.shp")  #Define Project Area

## Image Prep Input ##

## Short Dry
#start_date = "2023-01-01"
#end_date = "2023-03-31"

## Long Dry
start_date = "2020-01-01"
end_date = "2020-12-31"  

max_cc = 30
#band_list=c("B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B11","B12","SCL")
band_list=c("B01","B03","B05","B09","B12","SCL")

  sf_pa<-st_transform(sf_pa$geometry, 4326) #transform PA to lat/long for search purposes
  #sf_pa_buffer<-st_buffer(st_zm(sf_pa), 10000, 10) #### 10k buffer around LB
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
                                 property_filter = function(x) {x[["platform"]] !='landsat-7'}) #convert stac to image collection, keep qa_pixel band for cloud masking, remove any images from Landsat-7
      
      
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
      #v.overview = cube_view(extent=roi_ext, dx = 10, dy = 10, dt="P15D", srs = col.epsg,
      #                       aggregation = "mean", resampling = "bilinear")
      v.overview = cube_view(extent=roi_ext, dx = 10, dy = 10, dt="P15D", srs = col.epsg,
                             aggregation = "mean", resampling = "bilinear")
      
      cloud_mask = image_mask("SCL", values =c(3,8,9)) # Mask Sentinel clouds  # Mask Landsat clouds 
      
      #Function to create data cube and write .tif files to output directory
      #Need final pipe function as write_tif to make X into list of file paths, which we can export into raster stack
      
      #start.time.stac<-Sys.time()
      #print(start.time.stac)
      
      #x = raster_cube(col, v.overview, mask = cloud_mask)|>
      # select_bands(c("B01","B02","B03","B04","B05","B06","B07","B08","B8A","B09","B11","B12"))|>
      #  reduce_time(c("median(B01)","median(B02)","median(B03)","median(B04)","median(B05)","median(B06)",
      #                "median(B07)","median(B08)","median(B8A)","median(B09)","median(B11)","median(B12)"))|>
      #  write_tif(dir = output.dir, prefix = paste0("2023_c",max_cc,"_P15D_median_LDry_all","_", sep=""))
      #end.time.stac<-Sys.time()
      #print(end.time.stac-start.time.stac)
      #plotRGB(rast(x), r=4, g=3, b = 2)
      
      start.time.stac<-Sys.time()
      print(start.time.stac)
      
      x = raster_cube(col, v.overview, mask = cloud_mask)|>
       select_bands(c("B01","B03","B05","B09","B12"))|>
       reduce_time(c("median(B01)","median(B03)","median(B05)","median(B09)","median(B12)"))|>
       write_tif(dir = output.dir, prefix = paste0("2020_c",max_cc,"_P1Y_pca_in_median","_", sep=""))
      end.time.stac<-Sys.time()
      print(end.time.stac-start.time.stac)

      
      
      ##### NDFI ######
      img<-raster::brick(x)
      
      em<-(data=matrix(c(119,1514,1799,4031,
                         475,1597,2479,8714,
                         169,1421,3158,7900,
                         6250,3053,5437,8989,
                         2399,7707,7707,7002,
                         675,1975,6646,6607), nrow=4))
      
      
      colnames(em)<-c("B02_mean","B03_mean","B04_mean","B08_mean","B11_mean","B12_mean")
      rownames(em)<-c("gv","npv","soil","cloud")
      print("Spectral Unmixing Begun")
      time.start<-Sys.time()
      probs<-RStoolbox::mesma(img, em, method = "NNLS")
      time.end<-Sys.time()
      print("Spectral Unmixing Finished")
      print(time.end-time.start)
      
      summed<-(probs$gv+probs$npv+probs$soil+probs$cloud)*100
      gvs<-((probs$gv*100)/summed)*100
      npvSoil<- (probs$npv*100) + (probs$soil*100)
      ndfi <- (gvs-npvSoil)/(gvs+npvSoil)
      writeRaster(ndfi,paste0(output.dir,"//c",max_cc,"_MK_",start_date,"_SDry_ndfi.tif", sep=""))
      writeRaster(probs$gv,paste0(output.dir,"//c",max_cc,"_MK_",start_date,"_SDry_gv.tif", sep=""))
      writeRaster(probs$npv,paste0(output.dir,"//c",max_cc,"_MK_",start_date,"_SDry_npv.tif", sep=""))
      print("Finished Writing NDFI .tif to output")
      
      

img<-rast("Sent2//2023_c10_P15D_median_Ldry_all_2023-07-01.tif")
source<-"MK_LDry_S2"
###### Sentinel ######
    BLUE<- img$B02_median/10000
    GREEN<-img$B03_median/10000
    RED<- img$B04_median/10000
    RE1<-img$B05_median/10000
    RE2<-img$B06_median/10000
    RE3<-img$B07_median/10000
    NIR<- img$B08_median/10000
    RE4<-img$B8A_median/10000
    SWIR1<- img$B11_median/10000
    SWIR2<-img$B12_median/10000      
      
#### KNDVI ####
knr <- exp(-((NIR)-(RED))^2/(2))
kNDVI <- (1-knr) / (1+knr)
plot(kNDVI, main = "kNDVI")
writeRaster(kNDVI,paste0(output.dir,"//c",max_cc,"_",source,"_kndvi.tif", sep=""))

##### EVI #####
EVI<-2.5*((NIR-RED)/(1+NIR+6*RED-7.5*BLUE))
plot(EVI, main = "EVI")
writeRaster(EVI,paste0(output.dir,"//c",max_cc,"_",source,"_EVI.tif", sep=""))

##### MSAVI #####
MSAVI<-NIR + 0.5 - (0.5 * sqrt((2 * NIR + 1)^2 - 8 * (NIR - (2 * RED))))
plot(MSAVI, main = "MSAVI")
writeRaster(MSAVI,paste0(output.dir,"//c",max_cc,"_",source,"_MSAVI.tif", sep=""))


###### NDWI #####
NDWI<-(NIR - SWIR1)/(NIR + SWIR1)
plot(NDWI, main = "NDWI")
writeRaster(NDWI, paste0(output.dir,"//c",max_cc,"_",source,"_NDWI.tif"))

##### GNDVI #####

GNDVI<-(NIR-GREEN)/(NIR+GREEN)
plot(GNDVI, main = "GNDVI")
writeRaster(GNDVI, paste0(output.dir,"//c",max_cc,"_",source,"_GNDVI.tif"))

##### Bare Soil Index (BSI) #####
BSI<-((RED+SWIR1)-(NIR+BLUE))/((RED+SWIR1)+(NIR+BLUE))
plot(BSI, main = "BSI")
writeRaster(BSI, paste0(output.dir,"//c",max_cc,"_",source,"_BSI.tif"))

##### Atmospherically Resistant Vegetation Index #####
ARVI<-(NIR-(2*RED)+BLUE)/(NIR+(2*RED)+BLUE)
writeRaster(ARVI, paste0(output.dir,"//c",max_cc,"_",source,"_ARVI.tif"))

##### Green Coverage Index (GCI) #####
GCI<-(NIR)/(GREEN)-1
writeRaster(GCI, paste0(output.dir,"//c",max_cc,"_",source,"_GCI.tif"))


##### Structure Insensitive Pigment Index (SIPI) #####
SIPI<-(NIR-BLUE)/(NIR-RED)
plot(SIPI)
writeRaster(SIPI, paste0(output.dir,"//c",max_cc,"_",source,"_SIPI.tif"))


##### Normalized Burned Ratio Index (NBRI) #####
NBR<- (NIR-SWIR2)/(NIR+SWIR2)
plot(NBR)
writeRaster(NBR, paste0(output.dir,"//c",max_cc,"_",source,"_NBR.tif"))


##### Green Leaf Index (GLI) #####
GLI<-((GREEN-RED)+(GREEN-BLUE))/((2*GREEN)+RED+BLUE)
writeRaster(GLI, paste0(output.dir,"//c",max_cc,"_",source,"_GLI.tif"))


################################################################################
############################## Sentinel-2 Only #################################
################################################################################

##### Normalized Difference Red Edge #####
NDRE<-(NIR-RE1)/(NIR+RE1)
writeRaster(NDRE, paste0(output.dir,"//c",max_cc,"_",source,"_NDRE.tif"))

##### Red-edge Inflection Point (REIP) #####
REIP<-700+40*((((RED+RE3)/2)-RE1)/(RE2-RE1))
writeRaster(REIP, paste0(output.dir,"//c",max_cc,"_",source,"_REIP.tif"))

plot(NDRE, range = c(0,1))







#### BFAST ####
start_date = "2020-01-01"
end_date = "2023-07-25"  

max_cc = 20
band_list=c("B02","B03","B04","B06","B07","B08","B11","B12","SCL")


sf_pa<-st_transform(sf_pa$geometry, 4326) #transform PA to lat/long for search purposes
sf_pa_buffer<-st_buffer(st_zm(sf_pa), 10000, 10) #### 10k buffer around LB

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
v.overview = cube_view(extent=roi_ext, dx = 10, dy = 10, dt="P1M", srs = col.epsg,
                       aggregation = "mean", resampling = "bilinear")


cloud_mask = image_mask("SCL", values =c(3,8,9)) # Mask Sentinel clouds  # Mask Landsat clouds 


x = raster_cube(col, v.overview, mask = cloud_mask)|>
  select_bands(c("B04","B08"))|>
  reduce_time(names = c("change_date", "change_magnitude"), FUN = function(x) {
    knr <- exp(-((x["B08",]/10000)-(x["B04",]/10000))^2/(2))
    kndvi <- (1-knr) / (1+knr)   
    if (all(is.na(kndvi))) {
      return(c(NA,NA))
    }
    kndvi_ts = ts(kndvi, start = c(2020, 1), frequency = 12)
    library(bfast)
    tryCatch({
      result = bfastmonitor(kndvi_ts, start = c(2023,1), level = 0.01)
      return(c(result$breakpoint, result$magnitude))
    }, error = function(x) {
      return(c(NA,NA))
    })
  }) |>
  write_tif("makame_bfast_2023_10km_v1.tif")
