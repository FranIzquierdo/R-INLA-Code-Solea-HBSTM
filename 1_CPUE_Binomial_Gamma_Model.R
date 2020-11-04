rm(list=ls()) # Clean workspace
library(INLA)

# Load dataset ------------------------------------------------

data<-read.table(file="./data_trend/IG_14_data19.txt", dec=".", header=TRUE)

# Prepare polygons -------------------------------------------------------------

library(sp)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)
library(rasterVis)

# get shapefiles from http://www.naturalearthdata.com
shape_path <- "./2_oceandatafiles/50m_physical_separated_zip/"
ocean_shapefile <- paste(shape_path, "ne_50m_ocean/ne_50m_ocean.shp", sep="")
layer <- ogrListLayers(ocean_shapefile)
ogrInfo(ocean_shapefile, layer=layer)
ocean_poly <- readOGR(ocean_shapefile, layer=layer) # Read shapefile
ocean <- ocean_poly
bbx <- readWKT("POLYGON((-9.65 44.20, -1 44.20, -1 41.83, -9.65 41.83,
               -9.65 44.20))") 
proj4string(bbx)<-proj4string(ocean) # Set CRS from pen_rec (WGS84)
ocean.cut <- gIntersection(ocean, bbx) # Intersect 
pen.cut<-gDifference(bbx, ocean.cut) # Difference: Obtain study area

par(mfrow=c(2,2)) # PLOT files
plot(ocean_poly, col="lightblue", main="Ocean Map")
plot(bbx, main="bbx") # DOMAIN
plot(ocean.cut, main="ocean.cut") # STUDY AREA
points(data[,1:2], pch=20,cex=0.7, col="darkred")
plot(pen.cut, main="pen.cut") # PENINSULA

# 2UTM ---------- Transform to UTM: ZONE 30-------------------------------------

# Polygons
CRS.new <- CRS("+init=epsg:4326 +proj=utm +zone=30 +units=km +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
ocean.cut.km <- spTransform(ocean.cut, CRS.new)
pen.cut.km <- spTransform(pen.cut, CRS.new)
bbx.km<- spTransform(bbx, CRS.new)
summary(ocean.cut.km)

# Sampling locations
d <- data.frame(lon=data[,1], lat=data[,2])
coordinates(d) <- c("lon", "lat")
proj4string(d) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") # WGS 84
CRS.new <- CRS("+init=epsg:4326 +proj=utm +zone=30 +units=km +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
d.ch1903 <- spTransform(d, CRS.new)
coords4<-coordinates(d.ch1903)
coords4<-as.matrix(coords4)

par(mfrow=c(1,2)) # PLOT files UTM
plot(ocean.cut.km, main="ocean.cut", asp=1) # STUDY AREA
points(coords4[,1:2], pch=20,cex=0.7, col="darkred")
plot(pen.cut.km, main="pen.cut", asp=1) # PENINSULA
points(coords4[,1:2], pch=20,cex=0.7, col="darkred")

# add UTM to dataset ------------------------------------------------------

data1<-cbind(data,coords4)
colnames(data1)<-c("lon","lat","weight_kg","pres","bathy","year_id","year","haul","camp","ID","weight","gbathy","Igbathy","lon_utm","lat_utm")
data<-data1

# Distance Matern ----------------------------------------------------------

# This function help us to guess the magnitude of the spatial range

D <- dist(coords4)
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

RangeGuess<- 100 # try 100 and 150 km
MaxEdge<- RangeGuess/5
ConvHull <- inla.nonconvex.hull(coords4,convex=-0.05)# mesh reducido
mesh <- inla.mesh.2d(boundary = ConvHull, max.edge = c(1,4) * MaxEdge, cutoff=MaxEdge/5)
dev.off()
plot(mesh, asp=1, main=" ")
points(coords4, col=1, pch=16, cex=0.7)
plot(pen.cut.km, col="grey", add=TRUE)
mesh$n


# SPDE ------------------------------------------------------------------------------------------------

Range0<-100# we force inla to come up with a range larger than 100 
# P(Range < 100 km) = 0.001  and P(sigma > 0.5) = 0.5
# If we don't have a clue or priior information about sigma we can give the same probability of being higher/lower than 0.5

AlphaRange0<- 0.001
Sigma0<- 0.5 #
AlphaSigma0 <-0.05
spde<-inla.spde2.pcmatern(mesh, prior.range = c(Range0, AlphaRange0), 
                          prior.sigma = c(Sigma0,AlphaSigma0))
# A matrix --------------------------------------------------------------------------------------------

est.temp <- inla.spde.make.A(mesh, loc=matrix(c(data$lon_utm,data$lat_utm),
                             ncol=2),
                             #index=index.est
                             group=data$year_id,
                             n.group=length(unique(data$year_id)))

# index--------------------------------------------------------------------------------------

mesh.index.bin<- inla.spde.make.index("i.bin", n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))
mesh.index.con<- inla.spde.make.index("i.con", n.spde=spde$n.spde, 
                                      n.group=max(data$year_id))

# Stack ------------------------------------------------------------

est.bin<-inla.stack(data=list(y=cbind(data$pres,NA)),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.bin,
                                 list(bin.b0=1,
                                      bath.bin=data$gbathy, 
                                      year.bin=data$year_id)),tag='est.bin')


est.con<-inla.stack(data=list(y=cbind(NA,ifelse(data$pres>0,data$weight_kg,NA))),
                    A=list(est.temp, 1),
                    effects=list(mesh.index.con,
                                 list(con.b0=1,
                                      bath.con=data$gbathy, 
                                      year.con=data$year_id)), tag='est.con')

est.stack<-inla.stack(est.bin,est.con)

# Link -------------------------------------------------------------------------------------

link=c(rep(1,length(est.stack$data$data[,1])/2), 
       rep(2,length(est.stack$data$data[,1])/2))

# Priors --------------------------------------------------------------------------

pcrho <- list(prior = 'pccor1', param = c(0.7, 0.7))

# informative priors from previous studies/knowledge
# Prob (\sigma > U) = \alpha = 0.05

hyper.prec.bath <- list(theta = list(prior="pc.prec", param = c(0.5, 0.05)))
hyper.prec.year <- list(theta = list(prior="pc.prec", param = c(1, 0.05)))


# MODEL with bath shared, and not shared spatial effect neither year trend-----------

prog_est_B_S_T_ns<- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(bath.bin, model="rw2", hyper=hyper.prec.bath) + f(bath.con, copy="bath.bin", fixed=F) +
  f(year.bin, model="rw2", hyper=hyper.prec.year) + f(year.con, model="rw2",hyper=hyper.prec.year) +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con,model=spde,group = i.con.group,control.group = list(model="ar1"))

prog_est_B_S_T_ns<-inla(prog_est_B_S_T_ns,
                     family=c('binomial',"gamma"),
                     #control.fixed=list(expand.factor.strategy="inla"), 
                     data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T, config=TRUE),
                     control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
                     verbose=TRUE,control.inla = list(strategy = "gaussian"), num.threads = 2)

saveRDS(prog_est_B_S_T_ns, "8_cpue_trend_ns_utm.rds")


# MODEL with bath and trend year shared, not shared spatial effect ------------------------------------------------------------------

prog_est_B_S_T_s<- y ~ -1 + 
  f(bin.b0,model="linear", prec.linear=.1,mean.linear=-.5) + con.b0 +
  f(bath.bin, model="rw2", hyper=hyper.prec.bath) + f(bath.con, copy="bath.bin", fixed=F) +
  f(year.bin, model="rw2", hyper=hyper.prec.year) +f (year.con, copy="year.bin") +
  f(i.bin,model=spde,group = i.bin.group,control.group = list(model="ar1")) +
  f(i.con,model=spde,group = i.con.group,control.group = list(model="ar1"))

prog_est_B_S_T_s<-inla(prog_est_B_S_T_s,
                     family=c('binomial',"gamma"),
                     #control.fixed=list(expand.factor.strategy="inla"), 
                     data=inla.stack.data(est.stack), control.compute=list(dic=TRUE,cpo=TRUE,waic=T, config=TRUE),
                     control.predictor=list(A=inla.stack.A(est.stack), compute=TRUE,link=link),
                     verbose=TRUE,control.inla = list(strategy = "gaussian"), num.threads = 2)

saveRDS(prog_est_B_S_T_s, "8_cpue_trend_s_utm.rds")# 8 años en 46 minutos

samples_s=inla.posterior.sample(n=1000, prog_est_B_S_T_s)
#saveRDS(samples_s, "8_cpue_s_IPS_1000_utm.rds")# 8 años en 46 minutos

# Plot results model TREND SHARED ------------------------------------------------------------

shared<- readRDS("./rds_trend/8_cpue_trend_s_utm.rds") 
summary(shared)

plot(shared$summary.random$year.con[,1:2])

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


# Trend marginal ----------------------------------------------------------------------------------------------------------------------------------------------------------

library(ggplot2)
suabinm <- shared$summary.random$year.con$mean
suabin2 <- shared$summary.random$year.con$`0.025quant`
suabin9 <-shared$summary.random$year.con$`0.975quant`
suabinID<-shared$summary.random$year.con$ID
year<-seq(from=2001, to=2019)
#year<-as.factor(year)
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID, year)

ggplot(data = suabin, aes(x = suabinID, y = suabinm)) +
  geom_line(aes(x = suabinID, y = suabinm), color="dodgerblue4") +
  #geom_ribbon(aes(x = suabinID, ymin = suabin2, ymax = suabin9), alpha = 0.3,fill="dodgerblue4")+
  ggtitle("") +
  xlab("Year") +
  ylab("Temporal shared effect") +
  theme_bw() + 
  geom_hline(yintercept = 0, linetype="dashed", 
                          color = "grey", size=1) + 
  scale_x_discrete(name ="Year", 
                   breaks=suabinID, labels=year)



## Projector MATRIX------------------------------------------------------------------------

library(splancs)
box<-bbox(bbx.km)
Xrange<-c(box[1,1], box[1,2])
Yrange<-c(box[2,1], box[2,2])
nxy=c(120,120)
prj <- inla.mesh.projector(mesh, xlim=range(Xrange),
                           ylim=range(Yrange), dims=nxy) 

# Loops Spatial effect plot -------------------------------------------------------------------------------

m<-mesh$n
k<-which.max(unique(data$year_id))

### EE BIN mean

m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$i.bin$mean[1:m + (j - 1) * m])
  return(r) 
})
library(RColorBrewer)
library(fields)
colors <- brewer.pal(9, "YlGnBu")
zlm <- range(unlist(m.prj), na.rm = TRUE)
par(mfrow = c(3, 3), mar = c(0, 1.3, 1.2, 0.8))
for (j in 1:k) {
  image.plot(x = prj$x, y = prj$y, z = m.prj[[j]], asp = 1, 
             xlab = '', zlim = zlm, axes = FALSE, col = colors,
             main = paste0("Year ", j), cex=0.7, legend=TRUE,legend.cex = 5, legend.shrink = 0.38)
  #points(coords4[igr == j, 1:2], pch = 19, cex=0.7)
  plot(pen.cut.km, col=c("grey","darkgrey"), add=TRUE)
}

### EE GAMMA mean

m.prj <- lapply(1:k, function(j) {
  r <- inla.mesh.project(prj,
                         shared$summary.ran$i.con$mean[1:m + (j - 1) * m])
  return(r) 
})

colors <- brewer.pal(9, "YlOrRd")
zlm <- range(unlist(m.prj), na.rm = TRUE)
par(mfrow = c(3, 3), mar = c(0, 1.3, 1.2, 0.8))

for (j in 1:k) {
  image.plot(x = prj$x, y = prj$y, z = m.prj[[j]], asp = 1, 
             xlab = '', zlim = zlm, axes = FALSE, col = colors,
             main = paste0("Year ", j), cex=0.7,legend=TRUE,legend.shrink = 0.37)
  #points(coords4[igr == j, 1:2], pch = 19, cex=0.7)
  plot(pen.cut.km, col=c("grey","darkgrey"), add=TRUE)
}

### REPEAT THE SAME FOR SD
