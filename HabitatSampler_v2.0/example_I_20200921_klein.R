# HaSa - HabitatSampler
# for documentation see "HabitatSampler_Usage.txt"

##0)##
#####

##0.1##
setwd("E:/USER/")
#-> to source file location

##0.2##
inPath<-"schmehl/HaSa/HabitatSampler-master/HabitatSampler_v2.0/Funktionen/"
dataPath <- "drohnen_daten/20200921_TB/"
#outPath <- "schmehl/HaSa/HabitatSampler-master/HabitatSampler_v2.0/test/Results/"
outPath <- "E:/USER/schmehl/HaSa/klein1/" #absoluten Pfad angeben!!

##0.3##
dd = dir(paste0(inPath,"package_modified/"), pattern = ".r|R", full.names = T, all.files = T)
for(i in 1:length(dd)) {source(dd[i])}
source(paste(inPath,"install_packages.r",sep=""))

#zu ladende Pakete
usePack("rgdal","raster","maptools","spatialEco","randomForest","e1071","devtools","installr","velox","rgeos","leaflet","htmlwidgets")
rasterOptions(tmpdir="./RasterTmp/")

########################################################################################
###load data
a1<-brick(paste(dataPath,"20200921_I_filt1.tif",sep=""))
a1 = a1[[-4]] #Lösche "4. Kanal" weil keine Varianz drin (alle =255)
##1.a.2##
cut<-readOGR("E:/USER/schmehl/test_site_I_hasa_klein.shp")
a1<-clip(a1,cut)
proj4string(a1)<-proj4string(cut)


###texture
require(glcm)

texture<-glcm(a1[[2]], window = c(5, 5), shift = c(1, 1), statistics = c("mean", "variance", "homogeneity", "contrast", "dissimilarity", "entropy", "second_moment", "correlation"), na_opt="center", na_val=NA,scale_factor=10000, asinteger=TRUE)

###RGB Transformation to remove brightness diagonal (hsv color space)
brightCor<-function(x) {
  
  if (any ( is.na(x) )) {
    rep(NA, length(x))
  } else {
    
    flo1<-rgb2hsv(x[1],x[2],x[3],maxColorValue = 255)
    back_rgb<-col2rgb(hsv(h =flo1[1], s = flo1[2], v = 0.8, alpha=1), alpha = FALSE)
    #back_rgb<-col2rgb(test, alpha = FALSE)
    #flower<-rgb(back_rgb[1,],back_rgb[2,],back_rgb[3,],maxColorValue = 255)
    #return(flower)
  }
}

a1_br_cor<-calc(a1, brightCor, forceapply=T)

###indices

###Visible Atmospheric Resistant Index (VARI)
rgb = a1
rgb$vari = (rgb[[2]]-rgb[[1]])/(rgb[[2]]+rgb[[1]]-rgb[[3]])
###Triangular Greenness Index (TGI)
rgb$tgi<-rgb[[2]]-0.39*rgb[[1]]-0.61*rgb[[3]]
###Green Chromatic Value (gcv)
rgb$gcv<-rgb[[2]]/(rgb[[1]]+rgb[[2]]+rgb[[3]])
rgb$bcv<-rgb[[3]]/(rgb[[1]]+rgb[[2]]+rgb[[3]])
rgb$rcv<-rgb[[1]]/(rgb[[1]]+rgb[[2]]+rgb[[3]])
##others
rgb$ngrdi<-(rgb[[2]]-rgb[[1]])/(rgb[[2]]+rgb[[1]])
rgb$mgrvi<-(rgb[[2]]^2-rgb[[1]]^2)/(rgb[[2]]^2+rgb[[1]]^2)
rgb$rgbvi<-rgb[[2]]^2 - (rgb[[1]]*rgb[[3]])/rgb[[2]]^2 + (rgb[[1]]*rgb[[3]])
rgb$gli<-(2*rgb[[2]]-rgb[[1]]-rgb[[3]]) /(2*rgb[[2]]+ rgb[[1]]+ rgb[[3]])
rgb$exg<-2*rgb[[2]]-rgb[[1]]-rgb[[3]]
rgb$ng<-rgb[[2]]/((rgb[[1]]^0.667)*(rgb[[3]]^(1 - 0.667)))


###stack all together

a3 = brick(list(a1, texture[[c(1,2,3,5,6,7)]], a1_br_cor)) #not correlated features
a3 = brick(list(a1, texture[[c(1,2,3,5,6,7)]], rgb[[-c(1:3)]], a1_br_cor)) #not correlated features
proj4string(a3)<-proj4string(cut)

#a3<-stack(a2,texture[[c(1,2,3,5,6,7)]]) #not correlated features

#############################################################################

##get spectra of ref_points
##1.b.2##
shp<-readOGR("schmehl/HaSa/ref_shp.shp")
shp = shp[c(1,6,7,10),] #reduce to most important classes
ref<-as.data.frame(extract(a3,shp))

##1.c.2##
r=1; g=2; b=3;
plotRGB(a3,r=r,g=g,b=b,stretch="lin", axes=T)
plot(shp,pch=21,bg="red",col="yellow",cex=1.9,lwd=2.5,add=T)
##1.c.3##
col<-colorRampPalette(c("lightgrey","orange","yellow","limegreen","forestgreen"))
##1.c.4##
#classNames<-c("deciduous","coniferous","heather_young","heather_old","heather_shrub","bare_ground","xeric_grass")
classNames <- sub(pattern = "/", replacement = ".", x = shp$cover)

##2.a.1##
init.samples = 50
nb_models = 200 #kann zum testen auch <200 sein
nb_it = 10
buffer = 0.1 #in m -> muss an Bildauflösung angepassten werden. etwas über 1 bis 2 Pixel
mtry = 8  #muss kleiner als Anzahl kanäle sein, damit zufallsauswahl getroffen werden kann
init.seed = "sample"
n_classes = 4

multi_Class_Sampling(in.raster=a3, init.samples=init.samples, sample_type="regular",
                     nb_models=nb_models, nb_it=nb_it, buffer=buffer, reference=ref, model="rf",
                     mtry=mtry, last=F, seed=3, init.seed=init.seed, outPath=outPath, 
                     step=1, classNames=classNames, n_classes=n_classes, multiTest=1)

##2.b.1##
multi_Class_Sampling(in.raster=out.raster, init.samples=200, sample_type="regular",
                     nb_models=230, nb_it=nb_it, buffer=buffer, reference=out.reference, model="rf",
                     mtry=mtry, last=F, seed=3, init.seed=init.seed, outPath=outPath, 
                     step=3, #anpassen! bei letztem Schritt, last=T setzen; falls 
                     classNames=out.names, n_classes=n_classes, multiTest=1)

########################################################################################
##3)##
######

##3.a.1##
plot_Results(inPath=outPath)





#############################################################################

#############################################################################

#############################################################################

#############################################################################




##1)##
######


##1.a.1##


##1.b.1##
#ref<-read.table(paste(dataPath,"Example_Reference_table.txt", sep=""),header=T)
########################################################################################
##2)##
######

