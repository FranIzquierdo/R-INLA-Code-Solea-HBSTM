rm(list=ls()) # Clean workspace
library(INLA)

## ----load data----------------------------------------------------------------

data <- read_delim("9_modelo_tallas/datos/dataSize.csv", ";", escape_double = FALSE, trim_ws = TRUE)

data$year_id<-as.numeric(factor(data$year_id))# create year_id
library(dplyr)
data <- data %>%
  filter(bathy<198.5) # filter bathy obs


## ----COORDS UTM--------------------------------------------------------------- 

# Spatial polygons and coordinates already transformed to UTM
coords<-as.matrix(data[,6:7])# UTM COORDS

# Distance Matern ----------------------------------------------------------

# This function help us to guess the magnitude of the spatial range

D <- dist(coords)
par(mfrow = c(1,2), mar = c(5,5,2,2), cex.lab = 1.5)
hist(D, 
     freq = TRUE,
     main = "", 
     xlab = "Dist. between sites (km)",
     ylab = "Frequency")

plot(x = sort(D), 
     y = (1:length(D))/length(D), 
     type = "l",
     xlab = "Dist. sites (km)",
     ylab = "Cumulative proportion")

# The maximum distance between 2 points is 600 km
# We know that the majority of points will be separated at least by 100-150 km = prior spatial range


# mesh ------------------------------------------------------------------------------------------------

# This code help us to fit a good mesh. A good mesh implies that not too many points should be located inside of a single mesh triangle
# If the mesh is too finner, a good option to reduce computation is by using convex argument, to allow the dense inner part be mor focused over sampling points

RangeGuess<- 100 # rango sera mayor de 100 km 
MaxEdge<- RangeGuess/5
ConvHull <- inla.nonconvex.hull(coords,convex=-0.05)# mesh reducido
mesh <- inla.mesh.2d(boundary = ConvHull, max.edge = c(1,4) * MaxEdge, cutoff=MaxEdge/5)
dev.off()
plot(mesh, asp=1, main="n=740")
points(coords, col=1, pch=16, cex=0.7)
mesh$n

### Prepare datapred ----------------------------------------------------------

# We need to prepare stack.pred
# Transform raster of bathy to UTM, extract values and then apply inla.group
# Inla.group: we need to have the same bathy groups in the prediction and observation stack

library(raster)
bathy<-raster("9_modelo_tallas/datos/bathy.asc")# lat long without CRS

# CRS lat long
projection(bathy) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs") #project raster in WGS 84

# Project in UTM
sr <- "+init=epsg:4326 +proj=utm +zone=30 +units=km +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
bathy_utm<-projectRaster(bathy, crs=sr)
plot(bathy_utm)

# WRITE RASTER projected in UTM
dir<-"./9_modelo_tallas/datos"
writeRaster(bathy_utm, filename=file.path(dir, "bathy_utm"), bylayer=TRUE,format="GTiff", overwrite=TRUE)

# MASK bathy_utm 200 metros for SOLEA
bathy_utm_200 <- mask(x = bathy_utm , mask = ( bathy_utm[[1]]<(200) & bathy_utm[[1]]>(0) ), maskvalue = 0)
writeRaster(bathy_utm_200, filename=file.path(dir, "bathy_utm_200"), bylayer=TRUE,format="GTiff", overwrite=TRUE)

# EXTRACT BATHY PRED at MESH LOCS
plot(bathy_utm_200)
points(mesh$loc[,1:2])
meshvars<-extract(bathy_utm_200, as.matrix(mesh$loc[,1:2]))
summary(meshvars)# from 0 to 200 m

# Create bathy pred
bathy_pred<-rep(meshvars, times=19)# bathy values en mesh locs rep * 19 year

# Create YEAR PRED
year_pred<-rep(seq(from=1, to=19), times=1, each=mesh$n)

# UTM COORDINATES (mesh coords repeated each month and year)
lat_utm<-rep(mesh$loc[,2], times=19)#  utm zone 30 
lon_utm<-rep(mesh$loc[,1], times=19)

# CREATE DATAPRED
datapred<-cbind(bathy_pred,year_pred,lat_utm,lon_utm)
datapred<-as.data.frame(datapred)

write.table(datapred ,"./9_modelo_tallas/datos/9_datapred.txt" )

# INLA GROUP --------------------------------------------------------------

# We need to ensure that groups are exactly the same in data.obs and data.pred

data # obs
length(data$bathy)# 351
datapred # pred
length(datapred$bathy_pred) # 14060

all<-c(data$bathy,datapred$bathy_pred) #JOIN DATASETS 
length(all)# 14411
all<-as.data.frame(all)

idxgbath<-inla.group(all[,1], n=25, idx.only = TRUE) # ID
groupbath<-inla.group(all[,1], n=25) # bath groups
allin<-cbind(all, idxgbath, groupbath)

dataobs<-cbind(data, allin$groupbath[1:351],allin$idxgbath[1:351])
names(dataobs)<-c("year","lon","lat","bathy","size","lon_utm","lat_utm","year_id","groupbath","idxgbath")

datapred<-cbind(datapred, allin$groupbath[352:14411],allin$idxgbath[352:14411])
names(datapred)<-c("bathy_pred","year_pred","lat_utm", "lon_utm","groupbath_pred","idxgbath_pred")

write.table(dataobs, file="./9_modelo_tallas/datos/9_dataobs_IG.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)
write.table(datapred, file="./9_modelo_tallas/datos/9_datapred_IG.txt", append = FALSE, sep = "\t", dec = ".",row.names = TRUE, col.names = TRUE)


# LOAD datasets ------------------------------------------------------------------

data<-read.table(file="./9_modelo_tallas/datos/9_dataobs_IG.txt", dec=".", header=TRUE)
datapred<-read.table(file="./9_modelo_tallas/datos/9_datapred_IG.txt", dec=".", header=TRUE)

# SPDE ------------------------------------------------------------------------------------------------

Range0<-100 # we force inla to come up with a range larger than 100 
# P(Range < 100 km) = 0.001  and P(sigma > 0.5) = 0.5
# If we don't have a clue or priior information about sigma we can give the same probability of being higher/lower than 0.5

AlphaRange0<- 0.001
Sigma0<- 0.5 #
AlphaSigma0 <-0.05
spde<-inla.spde2.pcmatern(mesh, prior.range = c(Range0, AlphaRange0), 
                          prior.sigma = c(Sigma0,AlphaSigma0))

# A Matrix --------------------------------------------------------------------------------------------

est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$lon_utm,data$lat_utm),ncol=2),
                             #index=index.est
                             group=data$year_id,
                             n.group=length(unique(data$year_id)))

Astp <- inla.spde.make.A(mesh = mesh,
                         index=rep(1:mesh$n, times=19),
                         group=datapred$year_pred,
                         n.group = length(unique(datapred$year_pred)))

dim(est.temp)# A matrix est
dim(Astp)# A matrix for prediction

## ----idx-----------------------------------------------------------------
idx<- inla.spde.make.index("s", n.spde=spde$n.spde, n.group=max(data$year_id))


# PRIORS  -----------------------------------------------

pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))


# informative priors from previous studies/knowledge
# Prob (\sigma > U) = \alpha = 0.05

hyper.prec <- list(theta = list(prior="pc.prec", 
                                param = c(1, 0.05))) #allow smaller values prec iid
hyper.prec.bath <- list(theta = list(prior="pc.prec", 
                                     param = c(1, 0.05)))

# Stack-----------------------------------------------------------------------

# INLA does estimation and prediction at the same time
# estimation stack: it contains all the data information
# prediction stack: it contains data for all components but not for response variable

stk.e <- inla.stack(data=list(y=data$size, link=1), 
                    A=list(est.temp, 1),
                    effects=list(idx,list(a0=rep(1, length=length(data$size)), 
                                          bath=data$groupbath,year=data$year_id)), tag="est")

stk.p <- inla.stack(  data = list(y = NA, link=1),
                      A=list(Astp, 1),
                      effects=list(idx,list(a0=rep(1, length=length(datapred$lat)), 
                                            bath=datapred$groupbath_pred, year=datapred$year_pred)), tag="pred")# vessel ID con y sin
stk.full <- inla.stack(stk.e, stk.p)

# MODEL -----------------------------------------------------------------------

form <- y ~ 0 + a0 +
  f(s, model = spde, group = s.group, control.group = list(model = 'ar1', hyper = list(theta = pcrho))) +
  f(bath, model="rw2", hyper=hyper.prec.bath) + f(year, model="iid", hyper=hyper.prec)

abu.res <- inla(form, family = 'gaussian', 
  data = inla.stack.data(stk.full),
  control.predictor = list(A = inla.stack.A(stk.full),link=1),
  control.compute = list(waic=TRUE, cpo=TRUE, dic=TRUE, config=TRUE),
  control.inla = list(strategy = 'adaptive'), verbose=TRUE, num.threads=2)

saveRDS(abu.res, "./9_modelo_tallas/9_rds/9_modelo_tallas_gaussian_stkpred.rds")


# Results --------------------------------------------------------------

shared<-readRDS("./9_modelo_tallas/9_rds/9_modelo_tallas_gaussian_stkpred.rds")# complejo trend
summary(shared)

dev.off()
plot(shared, plot.fixed.effects = TRUE, plot.lincomb = FALSE, plot.random.effects = FALSE, plot.hyperparameters = FALSE,
     plot.predictor = FALSE, plot.q = FALSE, plot.cpo = FALSE, single = FALSE)

plot(shared, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=TRUE, plot.q=FALSE, plot.cpo=FALSE,
     single=FALSE)

plot(shared, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=TRUE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=FALSE,
     single=FALSE)

plot(shared, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE,
     plot.hyperparameters=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=TRUE,
     single=FALSE)


library(ggplot2)
suabinm <- shared$summary.random$year$mean
suabin2 <- shared$summary.random$year$`0.025quant`
suabin9 <-shared$summary.random$year$`0.975quant`
suabinID<-shared$summary.random$year$ID
year<-seq(from=2001, to=2019)
year<-as.factor(year)
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID, year)

ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
  geom_line(aes(x = suabinID, y = suabinm), color="black")+
  geom_ribbon(aes(x = suabinID, ymin = suabin2, ymax = suabin9), alpha = 0.3,fill="dodgerblue4")+
  ggtitle("")+
  xlab("Year")+
  ylab("Temporal shared effect")+
  theme_bw() + geom_hline(yintercept = 0, linetype="dashed", 
                          color = "gray39", size=0.7)+ theme(axis.text.x = element_text(angle =90 ))+
  scale_x_continuous(name="Year", limits=c(1, 19), breaks=suabinID, labels=year) 



## ----polygons--------------------------------------------------------------

library(sf)
library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)

# get shapefiles from http://www.naturalearthdata.com
shape_path <- "./data/oceandatafiles/50m_physical_separated_zip/"
ocean_shapefile <- paste(shape_path, "ne_50m_ocean/ne_50m_ocean.shp", sep="")
layer <- ogrListLayers(ocean_shapefile)
ogrInfo(ocean_shapefile, layer=layer)
ocean_poly <- readOGR(ocean_shapefile, layer=layer)# read the shape file
ocean <- ocean_poly
bbx <- readWKT("POLYGON((-9.65 44.20, -1 44.20, -1 41.83, -9.65 41.83, -9.65 44.20))") 
proj4string(bbx)<-proj4string(ocean) #  CRS from pen_rec (WGS84)
ocean.cut <- gIntersection(ocean, bbx)
pen.cut<-gDifference(bbx, ocean.cut)

# Polygons
CRS.new <- CRS("+init=epsg:4326 +proj=utm +zone=30 +units=km +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
ocean.cut.km <- spTransform(ocean.cut, CRS.new)
pen.cut.km <- spTransform(pen.cut, CRS.new)
bbx.km<- spTransform(bbx, CRS.new)
summary(ocean.cut.km)

## Projector MATRIX------------------------------------------------------------------------

library(splancs)
box<-bbox(ocean.cut.km)
Xrange<-c(box[1,1], box[1,2])
Yrange<-c(box[2,1], box[2,2])
nxy=c(120,120)

prj <- inla.mesh.projector(mesh, xlim=range(Xrange),
                           ylim=range(Yrange), dims=nxy) 


## ----PLOT spatial effect----------------------------------------------------------

# Repeat for SD

k<-max(unique(data$year_id))
m<-mesh$n
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$s$mean[1:m + (j - 1) * m])
  return(r) 
})


# create the breaks- and label vectors
ewbrks <- seq(-100,600,100)
nsbrks <- seq(4650,4950,100)


zlm <- range(unlist(m.prj), na.rm = TRUE)
library(RColorBrewer)
colors <- brewer.pal(5, "YlOrRd")
library(fields)

dev.off()
par(mfrow = c(5, 4), mar = c(0, 0, 1.2, 0))
library(fields)
for (j in 1:k) {
  jpeg(paste("./9_modelo_tallas/9_images/efesp/",j,"efesp.jpeg"), width=1700, height=1030, res=300)
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  bath<-raster("./9_modelo_tallas/datos/bathy_utm_200.tif")
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  
  image.plot(aa ,xaxt= "n", yaxt= "n",asp=1, xlab = '', ylab=' ',zlim = zlm, 
             axes = TRUE, col = colors, 
             breaks=c(-3.5,-2,-1.0,1,2,3.5),
             main = paste0("Year ", j), cex=0.7,
             axis.args=list(cex.axis=0.8))
  
  axis(1, at=ewbrks, cex.axis=0.8 )
  axis(2, at=nsbrks, cex.axis=0.8, las=1)
  
  #points(coords[igr == j, 1:2], pch = 19, cex=0.7)
  plot(pen.cut.km, col=c("grey","darkgrey"), add=TRUE)
  dev.off()
  
}


# PRED. SPATIAL -----------------------------------------------------------------

index.p <-inla.stack.index(stk.full, tag = "pred")$data # index of the pred. values
pred_mean <-shared$summary.linear.predictor[index.p, "0.5quant"] # fitted values MOMENTOS mas estables
pred_meanq2 <- shared$summary.fitted.values[index.p, "0.025quant"]
pred_meanq975 <- shared$summary.fitted.values[index.p, "0.975quant"]

# The next loop code allows us to run, cut and save the yearly prediction maps automatically

# MEDIAN ------------------------------------------------------------------

#Repeat this step for each quantile

m<-mesh$n
k<-which.max(unique(datapred$year_pred))

# create the breaks- and label vectors
ewbrks <- seq(-100,600,100)
nsbrks <- seq(4650,4950,100)

### EE BIN
m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         pred_mean[1:m + (j - 1) * m])# REPETIR CAMBIANDO EL pred_MEAN
  return(r) 
})

library(RColorBrewer)
library(fields)
colors <- brewer.pal(7, "YlOrRd")# 7
(zlm <- range(unlist(m.prj), na.rm = TRUE))# help for setting the breaks
dev.off()
par(mfrow = c(5, 4), mar = c(0, 1.3, 1.2, 0.8))
for (j in 1:k) {
  jpeg(paste("./9_modelo_tallas/9_images/prediction/median/",j,"pred_median.jpeg"), width=1700, height=1030, res=300)
  
  library(raster)
  sp.mean.raster1<-raster(list(x=prj$x,
                               y = prj$y,
                               z = m.prj[[j]]))
  bath<-raster("./9_modelo_tallas/datos/bathy_utm_200.tif")
  res<-resample(bath,sp.mean.raster1)
  aa<-mask(sp.mean.raster1,res)
  
  image.plot(aa ,xaxt= "n", yaxt= "n",asp=1, xlab = '', ylab=' ',zlim = zlm, axes = TRUE, col = colors, 
             breaks=c(26,28,30,32,34,36,38,40),
             main = paste0("Year ", j), cex=0.7,
             axis.args=list(cex.axis=0.8))
  
  axis(1, at=ewbrks, cex.axis=0.8 )
  axis(2, at=nsbrks, cex.axis=0.8, las=1)
  
  #points(coords[igr == j, 1:2], pch = 19, cex=0.7)
  plot(pen.cut.km, col=c("grey","darkgrey"), add=TRUE)
  dev.off()
}
