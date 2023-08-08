setwd("C://GIS//Tanzania//Makame") # Set working directory to whatever folder you would like where data can be written and retrieved

packages<-(c("sf", "gdalcubes", "rstac", "terra",
             "torch", "dplyr", "dismo", "exactextractr", 
             "reshape2","ggplot2", "parallel",
             "caret","randomForest","VSURF","RSNNS"))

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

rm(packages)
rm(installed_packages)


output.dir<-paste("./0802WF") # create a directory to store files 
if (!file.exists(output.dir)) {
  dir.create(output.dir)
}


opt<-rast("./Sent2/2023_c5_P15D_median_SDry_all_2023-01-01.tif")
s.1<-rast("./0731WF//MK_2023_SDry_SAR.tif")
si<-rast("./0731WF//3MK_2023_SDry_SI_stack.tif")
stack<-c(opt,si,s.1)
stack<-stack[[-28:-29]]

poly<-vect("MK_2023_training.shp")
table(poly$MK_Class)

###########Set seed to make sure the same random sample is selected next time########
set.seed(2) 
###Split the dataset into test and training dataset######

polysplit <- createDataPartition(poly$MK_Class,
                                  p = 0.5, 
                                  list = FALSE,
                                  times = 1)


training_poly <- poly[ polysplit,]
testing_poly  <- poly[-polysplit,]

sample_rast<-rasterize(training_poly,stack[[1]],"MK_Class")

training_points<-spatSample(sample_rast, size = 1000, method = "stratified", as.points = TRUE, xy=TRUE)

table(training_points$MK_Class)

training_sample<-st_as_sf(training_points) %>% group_by(MK_Class) %>% slice_sample(n = min(table(training_points$MK_Class)))

training_buffer<-st_buffer(training_sample, dist = 15)



training_data<-exact_extract(stack, training_buffer, 'mean')


training_data$Class<-as.factor(training_sample$MK_Class)

levels(training_data$Class)<-c("Ag","Forest","Grass","L.Forest","Water")

trainDF<-Filter(function(x)!all(is.na(x)), training_data)


testDF <- exact_extract(stack, st_as_sf(testing_poly), 'mean')
testDF$Class<-as.factor(testing_poly$MK_Class)


Fitcontrol <- trainControl("repeatedcv", 
                           number=10, 
                           repeats=10,
                           classProbs = TRUE)

####### Un-Reduced Classifications #####
    
  #### RF ####

rf<-caret::train(Class ~.,
                  data = trainDF, 
                  method= "rf",
                  trControl = Fitcontrol,
                  preProcess = c("center", "scale"),
                  importance = TRUE)

print(rf)
rf$finalModel
rfVarImp2<-importance(rf$finalModel, type=1, scale=TRUE)
plot(rfVarImp2)
importance(rf$finalModel)

rfplsClasses <- predict(rf, newdata = testDF)

caret::confusionMatrix(data = rfplsClasses, testDF$Class)


    #### SVM #####
# svm<-caret::train(Class ~.,
#                   data = trainDF, 
#                   method= "svmRadial",
#                   trControl = Fitcontrol,
#                   preProcess = c("center", "scale"),
#                   tunelength = 10)

# print(svm)
# svm$finalModel
# summary(svm)

# svmplsClasses <- predict(svm, newdata = testDF)
# caret::confusionMatrix(data = svmplsClasses, testDF$Class)


####### Variable Reduction ######


### VSURF ###
v.rf<-VSURF(x=subset(trainDF,select=-Class), y=trainDF$Class)

summary(v.rf)
print(v.rf)
v.rf$mean.perf

head(trainDF[c(v.rf$varselect.pred)])


# Normalization of VSURF selected variables for further classification and prediction
vsurf.stack<-stack[[c(v.rf$varselect.pred)]]

for(i in seq(1:nlyr(vsurf.stack))){
  band<-vsurf.stack[[i]]
  band<-setMinMax(band)
  band_min<-minmax(band)[1]
  band_max<-minmax(band)[2]
  band_norm<-(band-band_min)/(band_max-band_min)
  if(i==1){
    vsurf_norm<-band_norm
  } else {
    vsurf_norm<-c(vsurf_norm, band_norm)
  }
  gc()
}

# Create new training dataset from normalized data
vrf.training_data<-exact_extract(vsurf_norm, training_buffer, 'mean')
vrf.training_data$Class<-as.factor(training_buffer$MK_Class)

levels(vrf.training_data$Class)<-c("Ag","Forest","Grass","L.Forest","Water")

vrf.trainDF<-vrf.training_data

# Create testing data from retained polygons
vrf.testDF <- exact_extract(vsurf_norm, st_as_sf(testing_poly), 'mean')
vrf.testDF$Class<-as.factor(testing_poly$MK_Class)
levels(vrf.testDF$Class)<-c("Ag","Forest","Grass","L.Forest","Water")



#### RF using VSURF reduced Variables #####
vrf<-caret::train(Class ~.,
           data = vrf.trainDF, 
           method= "rf",
           trControl = Fitcontrol,
           preProcess = c("center", "scale"),
           importance = TRUE)

print(vrf)
vrf$finalModel
vrfVarImp2<-importance(vrf$finalModel, type=1, scale=TRUE)
plot(vrfVarImp2)
importance(vrf$finalModel)

vrfplsClasses <- predict(vrf, newdata = vrf.testDF)

levels(vrfplsClasses)<-c("Ag","Forest","Grass","L.Forest","Water")

caret::confusionMatrix(data = vrfplsClasses, as.factor(vrf.testDF$Class))


#### SVM  using VSURF reduced variables #####
vsvm<-caret::train(Class ~.,
            data = vrf.trainDF, 
            method= "svmPoly",
            trControl = Fitcontrol,
            preProcess = c("center", "scale"),
            tunelength = 10)

print(vsvm)
vsvm$finalModel
summary(vsvm)

vsvmplsClasses <- predict(vsvm, newdata = vrf.testDF)

levels(vsvmplsClasses)<-c("Ag","Forest","Grass","L.Forest","Water")

caret::confusionMatrix(data = vsvmplsClasses, vrf.testDF$Class)


#### MLP using VSURF reduced Variables #####

train_x<-vrf.trainDF[,1:nlyr(vsurf_norm)]
train_y<-decodeClassLabels(vrf.trainDF$Class)

vmlp<-mlp(x = train_x,
              y = train_y,
              size = c(10,10),
              maxit = 10000)


vmlpplsClasses<-as.factor(encodeClassLabels(predict(vmlp,vrf.testDF[,1:nlyr(vsurf_norm)])))

levels(vmlpplsClasses)<-c("Ag","Forest","Grass","L.Forest","Water")

caret::confusionMatrix(vmlpplsClasses, vrf.testDF$Class)

###### Predict ####### 

names(vsurf_norm)<-names(subset(vrf.trainDF, select= -Class))

# Classification Raster
vrf.class<-terra::predict(vsurf_norm, vrf, na.rm = TRUE)

writeRaster(vrf.class, "C://GIS//Tanzania//Makame//0807WF//MK_vrf_vsurf_class_0807.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))

# Class Probibility Raster
vrf.class.probs<-terra::predict(vsurf_norm, vrf ,type = 'prob', na.rm = TRUE)

writeRaster(vrf.class.probs, "C://GIS//Tanzania//Makame//0807WF//MK_vrf_vsurf_probs_0807.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))
