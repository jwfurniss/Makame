setwd("C://GIS//Tanzania//Makame") # Set working directory to whatever folder you would like where data can be written and retrieved

packages<-(c("sf", "gdalcubes", "rstac", "terra",
             "torch", "dplyr", "dismo", "exactextractr", 
             "ape", "reshape2","ggplot2", "parallel","RStoolbox",
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


opt<-rast("./Sent2/2023_c10_P15D_median_LDry_all_2023-07-01.tif")
s.1<-rast("./0731WF//MK_2023_VVVH_LDry_prj_rnm.tif")
si<-rast("./0731WF//2MK_2023_LDry_SI_stack.tif")
stack<-c(opt,si,s.1)
stack<-stack[[-28:-29]]

poly<-vect("MK_2023_training.shp")
table(poly$MK_Class)

###########Set seed to make sure the same random sample is selected next time########
set.seed(25) 
###Split the dataset into test and training dataset######

polysplit <- createDataPartition(poly$MK_Class,
                                  p = 0.5, 
                                  list = FALSE,
                                  times = 1)


training_poly <- poly[ polysplit,]
testing_poly  <- poly[-polysplit,]

training_points<-spatSample(training_poly, size = 70000, method = "random", strata = "MK_Class")
table(training_points$MK_Class)

training_sample<-st_as_sf(training_points) %>% group_by(MK_Class) %>% slice_sample(n = min(table(training_points$MK_Class)))

training_buffer<-st_buffer(training_sample, dist = 15)



training_data<-exact_extract(stack, training_buffer, 'mean')


training_data$Class<-as.factor(training_sample$MK_Class)
trainDF<-Filter(function(x)!all(is.na(x)), training_data)


testDF <- exact_extract(stack, st_as_sf(testing_poly), 'mean')
testDF$Class<-as.factor(testing_poly$MK_Class)


Fitcontrol <- trainControl("repeatedcv", 
                           number=10, 
                           repeats=10,
                           classProbs = TRUE)
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
svm<-caret::train(Class ~.,
                   data = trainDF, 
                   method= "svmRadial",
                   trControl = Fitcontrol,
                   preProcess = c("center", "scale"),
                   tunelength = 10)

print(svm)
svm$finalModel
summary(svm)

svmplsClasses <- predict(svm, newdata = testDF)
caret::confusionMatrix(data = svmplsClasses, testDF$Class)



### VSURF ###
v.rf<-VSURF(x=subset(trainDF,select=-Class), y=trainDF$Class)

summary(v.rf)
print(v.rf)
v.rf$mean.perf

head(trainDF[c(v.rf$varselect.pred)])

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

vrf.training_data<-exact_extract(vsurf_norm, training_buffer, 'mean')

vrf.training_data$Class<-training_buffer$MK_Class

vrf.trainDF<-vrf.training_data

vrf.testDF <- exact_extract(vsurf_norm, st_as_sf(testing_poly), 'mean')
vrf.testDF$Class<-as.factor(testing_poly$MK_Class)



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

rfplsClasses <- predict(vrf, newdata = vrf.testDF)

caret::confusionMatrix(data = rfplsClasses, as.factor(vrf.testDF$Class))


#### SVM #####
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
caret::confusionMatrix(data = vsvmplsClasses, vrf.testDF$Class)


#### MLP #####

train_x<-vrf.trainDF[,1:nlyr(vsurf_norm)]
train_y<-decodeClassLabels(vrf.trainDF$Class)

mlp.nnet<-mlp(x = train_x,
              y = train_y,
              size = c(50,50),
              maxit = 10000)


mlpplsClasses<-as.factor(encodeClassLabels(predict(mlp.nnet,vrf.testDF[,1:nlyr(vsurf_norm)])))

levels(mlpplsClasses)<-c("Ag","Forest","Grass","L.Forest","Water")

caret::confusionMatrix(mlpplsClasses, testDF$Class)

###### Predict ####### 

names(vsurf_norm)<-names(subset(vrf.trainDF, select= -Class))

vrf.class<-terra::predict(vsurf_norm, vrf, na.rm = TRUE)

writeRaster(vrf.class, "C://GIS//Tanzania//Makame//0802WF//MK_vrf_vsurf_class_0802_2.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))

vrf.class.probs<-terra::predict(vsurf_norm, vrf ,type = 'prob', na.rm = TRUE)
writeRaster(vrf.class.probs, "C://GIS//Tanzania//Makame//0802WF//MK_vrf_vsurf_probs_0802_2.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))
###
vsvm.class<-terra::predict(vsurf_norm, vsvm, na.rm = TRUE)

writeRaster(vsvm.class, "C://GIS//Tanzania//Makame//0802WF//MK_vsvm_vsurf_class_0802.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))

vsvm.class.probs<-terra::predict(vsurf_norm, vsvm ,type = 'prob', na.rm = TRUE)
writeRaster(vsvm.class.probs, "C://GIS//Tanzania//Makame//0802WF//MK_vsvm_vsurf_probs_0802.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))


vmlp.class<-terra::predict(vsurf_norm, mlp.nnet, na.rm = TRUE)
vmlp.class.encode<-encodeClassLabels(vlmp.class)

writeRaster(vmlp.class, "C://GIS//Tanzania//Makame//0802WF//MK_vmlp_vsurf_class_0802.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))

vmlp.class.probs<-terra::predict(vsurf_norm, vlmp ,type = 'prob', na.rm = TRUE)
writeRaster(vrf.class.probs, "C://GIS//Tanzania//Makame//0802WF//MK_vmlp_vsurf_probs_0802.tif", gdal=c("COMPRESS=DEFLATE","PREDICTOR=2"))
