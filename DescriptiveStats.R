# ------------------------------------------------------------------------------
# Title: Descriptive Statistics
#
# Author: Christoph Schönenberger
# Date: 30.10.17
# Version:
# Purpose:
# Instructions:
#
# ------------------------------------------------------------------------------
#
# Settings and Variables
rm(list=ls()) # Clean the environment
options(scipen=6) # Display digits, not the scientific version
dataFolder <- file.path("data") # Data folder
RFolder <- file.path("R") # RScript folder
figureFolder <- file.path("figures")
#
# ------------------------------------------------------------------------------
#
# Required R packages (add those on the top of your file)
# library()
require(RPostgreSQL)
require(rgeos)
require(rgdal)
require(raster)
require(maptools)
require(geosphere)
require(plotrix)
require(velox)
#
# install.packages("")
# install.packages("RPostgreSQL")
# install.packages("sf")
# install.packages("Plotrix")
# R Scripts
# source(file.path(RFolder,"myScript.r")) # load your scripts

############################################################################
# Database initialization --------------------------------------------------
############################################################################
pw <- {
  "m:9%sX-Kd8h"
}

# loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "whiterisk",
                 host = "sc12.geo.uzh.ch", port = 5432,
                 user = "cschoene", password = pw)
rm(pw) # removes the password

CRS <- "+init=epsg:4326"

############################################################################
# Data loading -------------------------------------------------------------
############################################################################
# Query routes from Switzerland
# The filtered routes don't contain routes with a length > 30km and < 1m, and none with a step length > 1km
routes <- dbGetQuery(con, "SELECT route_id, tour_id, st_astext(route_path) FROM route_filtered;")
# Load DHM25 and derivates
dhm25 <- raster("./data/DHM25_4326.tif")
slope25 <- raster("./data/slope_dhm25.tif")
aspect25 <- raster("./data/aspect_dhm25.tif")
curvature25 <- raster("./data/curv_dhm25.tif")
planCurv25 <- raster("./data/planCurv_dhm25.tif")
profCurv25 <- raster("./data/profCurv_dhm25.tif")
# Load borders of Switzerland
switzerland <- readOGR(dsn = "./data/borderCH_4326.shp")

############################################################################
# Data Conversion ----------------------------------------------------------
############################################################################
# Allocate space for lines
routesSpatialLines <- vector("list",length(routes[,3]))
# Convert routes from WKT to lines object
for(i in 1:length(routes[,3])){
  routesSpatialLines[[i]] <- readWKT(routes[i,3],id=i)@lines
}
rm(i)
# Convert to SpatialLines object
routesSpatialLines <- SpatialLines(unlist(routesSpatialLines),proj4string = CRS(CRS))
# Convert to SpatialLinesDataFrame
routesDF <- SpatialLinesDataFrame(routesSpatialLines,routes[,1:2])

############################################################################
# Calculations ------------------------------------------------------------
############################################################################

# Length ------------------------------------------------------------------

# Calculate length of routes
routesDF@data$"length" <- SpatialLinesLengths(routesDF, longlat = TRUE)*1000 #km to m

# Step Length --------------------------------------------------------------
# Average segment length
stepLength <- vector("numeric",length(routesDF))
for(i in 1:length(routesDF)){
  lengthI <- length(routesDF@lines[[i]]@Lines[[1]]@coords[,1])
  stepLength[i] <- routesDF@data$length[i]/lengthI
}
routesDF@data$stepLength <- stepLength
rm(lengthI,stepLength)
  
#====================================================================
# Subset for height calculations until more efficient method is found
#====================================================================
routesTest <- routesDF[sample(nrow(routesDF), 2958),]

# Height ------------------------------------------------------------------
# Allocate space for lists
mins <- vector("numeric",length(routesTest))
maxs <- vector("numeric",length(routesTest))
means <- vector("numeric",length(routesTest))

# Loop through chunks of whole dataset
for(i in seq(1,length(routesTest),986)){
  # buffer around lines
  buffer <- suppressWarnings(buffer(routesTest[i:(i+985),],0.001))
  # Crop digital elevation model to decrease size
  dhm25_cr <- crop(dhm25,buffer)
  # Create velox object
  dhm25_cr_vx <- velox(dhm25_cr)
  # Extract values
  vals <- dhm25_cr_vx$extract(routesTest[i:(i+985),])
  # Store aggregated values in allocated list
  for(j in 1:length(vals)){
    mins[i+j-1] <- min(unlist(vals[j]))
    maxs[i+j-1] <- max(unlist(vals[j]))
    means[i+j-1] <- mean(unlist(vals[j]))
  }
}
# Add values to dataframe
routesTest$minHeight <- mins
routesTest$maxHeight <- maxs
routesTest$meanHeight <- means
rm(i,j,mins,maxs,means,dhm25_cr,dhm25_cr_vx,vals)


# Slope -------------------------------------------------------------------
# Extract all slope values along lines from Slope (from DHM)
# Allocate space for lists
mins <- vector("numeric",length(routesTest))
maxs <- vector("numeric",length(routesTest))
means <- vector("numeric",length(routesTest))

# Loop through chunks of whole dataset
for(i in seq(1,length(routesTest),986)){
  # buffer around lines
  buffer <- suppressWarnings(buffer(routesTest[i:(i+985),],0.001))
  # Crop digital elevation model to decrease size
  slope25_cr <- crop(slope25,buffer)
  # Create velox object
  slope25_cr_vx <- velox(slope25_cr)
  # Extract values
  vals <- slope25_cr_vx$extract(routesTest[i:(i+985),])
  # Store aggregated values in allocated list
  for(j in 1:length(vals)){
    mins[i+j-1] <- min(unlist(vals[j]))
    maxs[i+j-1] <- max(unlist(vals[j]))
    means[i+j-1] <- mean(unlist(vals[j]))
  }
}
# Add values to dataframe
routesTest$minSlope <- mins
routesTest$maxSlope <- maxs
routesTest$meanSlope <- means
rm(i,j,mins,maxs,means,slope25_cr,slope25_cr_vx,vals)

# Aspect --------------------------------------------------------------
# Extract all aspect values along lines from aspect25 (from DHM)
# Allocate space
aspect <- rep(list(numeric()),length(routesTest))

for(i in seq(1,length(routesTest),986)){
  buffer <- suppressWarnings(buffer(routesTest[i:(i+985),],0.001))
  aspect25_cr <- crop(aspect25,buffer)
  aspect25_cr_vx <- velox(aspect25_cr)
  vals <- aspect25_cr_vx$extract(routesTest[i:(i+985),])
  for(j in 1:length(vals)){
    aspect[[i+j-1]] <- as.vector(unlist(vals[j]))
  }
}
routesTest@data$aspect <- aspect

aspectCount <- data.frame(N=numeric(),NE=numeric(),E=numeric(),SE=numeric(),S=numeric(),SW=numeric(),W=numeric(),NW=numeric())
for(i in 1:length(routesTest)){
  aspectCount[i,]$N <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]<=22.5 | routesTest@data$aspect[[i]]>337.5)])
  aspectCount[i,]$NE <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>22.5 & routesTest@data$aspect[[i]]<=67.5)])
  aspectCount[i,]$E <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>67.5&routesTest@data$aspect[[i]]<=112.5)])
  aspectCount[i,]$SE <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>112.5&routesTest@data$aspect[[i]]<=157.5)])
  aspectCount[i,]$S <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>157.5&routesTest@data$aspect[[i]]<=202.5)])
  aspectCount[i,]$SW <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>202.5&routesTest@data$aspect[[i]]<=247.5)])
  aspectCount[i,]$W <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>247.5&routesTest@data$aspect[[i]]<=292.5)])
  aspectCount[i,]$NW <- length(routesTest@data$aspect[[i]][which(routesTest@data$aspect[[i]]>292.5&routesTest@data$aspect[[i]]<=337.5)])
}

aspectSums <- colSums(aspectCount)
aspectSums <- aspectSums/max(aspectSums)

rm(i,j,aspect25_cr,aspect25_cr_vx)

# Plan Curvature ---------------------------------------------------------------
# Extract all plan curvature values along lines from curvature (from DHM)
# Allocate space for lists
mins <- vector("numeric",length(routesTest))
maxs <- vector("numeric",length(routesTest))
means <- vector("numeric",length(routesTest))

# Loop through chunks of whole dataset
for(i in seq(1,length(routesTest),986)){
  # buffer around lines
  buffer <- suppressWarnings(buffer(routesTest[i:(i+985),],0.001))
  # Crop digital elevation model to decrease size
  planCurv25_cr <- crop(planCurv25,buffer)
  # Create velox object
  planCurv25_cr_vx <- velox(planCurv25_cr)
  # Extract values
  vals <- planCurv25_cr_vx$extract(routesTest[i:(i+985),])
  # Store aggregated values in allocated list
  for(j in 1:length(vals)){
    mins[i+j-1] <- min(unlist(vals[j]))
    maxs[i+j-1] <- max(unlist(vals[j]))
    means[i+j-1] <- mean(unlist(vals[j]))
  }
}
# Add values to dataframe
routesTest$minPlanCurv <- mins
routesTest$maxPlanCurv <- maxs
routesTest$meanPlanCurv <- means
rm(i,j,mins,maxs,means,planCurv25_cr,planCurv25_cr_vx,vals)

# Profile Curvature ---------------------------------------------------------------
# Extract all profile curvature values along lines from curvature (from DHM)
# Allocate space for lists
mins <- vector("numeric",length(routesTest))
maxs <- vector("numeric",length(routesTest))
means <- vector("numeric",length(routesTest))

# Loop through chunks of whole dataset
for(i in seq(1,length(routesTest),986)){
  # buffer around lines
  buffer <- suppressWarnings(buffer(routesTest[i:(i+985),],0.001))
  # Crop digital elevation model to decrease size
  profCurv25_cr <- crop(profCurv25,buffer)
  # Create velox object
  profCurv25_cr_vx <- velox(profCurv25_cr)
  # Extract values
  vals <- profCurv25_cr_vx$extract(routesTest[i:(i+985),])
  # Store aggregated values in allocated list
  for(j in 1:length(vals)){
    mins[i+j-1] <- min(unlist(vals[j]))
    maxs[i+j-1] <- max(unlist(vals[j]))
    means[i+j-1] <- mean(unlist(vals[j]))
  }
}
# Add values to dataframe
routesTest$minProfCurv <- mins
routesTest$maxProfCurv <- maxs
routesTest$meanProfCurv <- means
rm(i,j,mins,maxs,means,profCurv25_cr,profCurv25_cr_vx,vals)


# Heading -----------------------------------------------------------------
# Allocate space 
headings <- vector("list",length(routesTest))
# For each line get length and allocate space
for(i in 1:length(routesTest)){
  lengthI <- length(routesTest@lines[[i]]@Lines[[1]]@coords[,1])-1
  heading <- vector("numeric",length(lengthI))
  # For each point get bearing to next point
  for(j in 1:lengthI){
    he <- bearing(routesTest@lines[[i]]@Lines[[1]]@coords[j,],routesTest@lines[[i]]@Lines[[1]]@coords[j+1,])
    if(he >= 0){
      heading[j] <- he
    } else {
      heading[j] <- 360-(180+he) 
    }
  }
  headings[[i]] <- heading
}
# Store in dataframe
routesTest@data$headings <- headings
rm(he,heading,headings,i,j,lengthI)

headCount <- data.frame(N=numeric(),NE=numeric(),E=numeric(),SE=numeric(),S=numeric(),SW=numeric(),W=numeric(),NW=numeric())
for(i in 1:length(routesTest)){
  headCount[i,]$N <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]<=22.5 | routesTest@data$headings[[i]]>337.5)])
  headCount[i,]$NE <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>22.5 & routesTest@data$headings[[i]]<=67.5)])
  headCount[i,]$E <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>67.5&routesTest@data$headings[[i]]<=112.5)])
  headCount[i,]$SE <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>112.5&routesTest@data$headings[[i]]<=157.5)])
  headCount[i,]$S <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>157.5&routesTest@data$headings[[i]]<=202.5)])
  headCount[i,]$SW <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>202.5&routesTest@data$headings[[i]]<=247.5)])
  headCount[i,]$W <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>247.5&routesTest@data$headings[[i]]<=292.5)])
  headCount[i,]$NW <- length(routesTest@data$headings[[i]][which(routesTest@data$headings[[i]]>292.5&routesTest@data$headings[[i]]<=337.5)])
}

headSums <- colSums(headCount)
headSums <- headSums/max(headSums)

# Straightness -------------------------------------------------------------
straightness <- vector("numeric",length(routesTest))
for(i in 1:length(routesTest)){
  lengthI <- length(routesTest@lines[[i]]@Lines[[1]]@coords[,1])
  # Diffusion distance (first to last coordinate)
  diffDist <- distGeo(routesTest@lines[[i]]@Lines[[1]]@coords[1,],routesTest@lines[[i]]@Lines[[1]]@coords[lengthI,])
  straightness[i] <- 1-(diffDist/routesTest@data$length[i])
}
routesTest@data$straightness <- straightness
rm(lengthI,diffDist,straightness)

############################################################################
# Visualization ------------------------------------------------------------
############################################################################

png("./img/boxplots/lengthAll.png")
boxplot(routesDF@data$length, main="Boxplot Route Length", ylab="Length of routes [m]") 
dev.off()
# png("./img/boxplots/lengthSubset.png")
# boxplot(routesDF@data$length[which(routesDF@data$length<30000)], main="Boxplot Route Length", ylab="Length of routes [m]") 
# dev.off()

png("./img/boxplots/stepLengthAll.png")
boxplot(routesDF@data$stepLength, main="Boxplot Route Step Width", ylab="Widths of steps [m]") 
dev.off()
# png("./img/boxplots/stepLengthsubsetLength.png")
# boxplot(routesDF@data$stepLength[which(routesDF@data$length<30000)], main="Boxplot Route Step Width", ylab="Widths of steps [m]") 
# dev.off()
# png("./img/boxplots/stepLengthsubsetStLength.png")
# boxplot(routesDF@data$stepLength[which(routesDF@data$stepLength<1000)], main="Boxplot Route Step Width", ylab="Widths of steps [m]") 
# dev.off()

png("./img/boxplots/heightMean.png")
boxplot(routesTest@data$meanHeight, main="Boxplot Mean Height", ylab="Height [m]") 
dev.off()
png("./img/boxplots/heightMin.png")
boxplot(routesTest@data$minHeight, main="Boxplot Minimum Height", ylab="Height [m]") 
dev.off()
png("./img/boxplots/heightMax.png")
boxplot(routesTest@data$maxHeight, main="Boxplot Maximum Height", ylab="Height [m]") 
dev.off()

png("./img/boxplots/slopeMean.png")
boxplot(routesTest@data$meanSlope, main="Boxplot Mean Slope", ylab="Slope [°]") 
dev.off()
png("./img/boxplots/slopeMin.png")
boxplot(routesTest@data$minSlope, main="Boxplot Minimum Slope", ylab="Slope [°]") 
dev.off()
png("./img/boxplots/slopeMax.png")
boxplot(routesTest@data$maxSlope, main="Boxplot Maximum Slope", ylab="Slope [°]") 
dev.off()

png("./img/boxplots/planCurvatureMean.png")
boxplot(routesTest@data$meanPlanCurv, main="Boxplot Mean Plan Curvature", ylab="Curvature [1/100 m]") 
dev.off()
png("./img/boxplots/planCurvatureMin.png")
boxplot(routesTest@data$minPlanCurv, main="Boxplot Minimum Plan Curvature", ylab="Curvature [1/100 m]") 
dev.off()
png("./img/boxplots/planCurvatureMax.png")
boxplot(routesTest@data$maxPlanCurv, main="Boxplot Maximum Plan Curvature", ylab="Curvature [1/100 m]") 
dev.off()

png("./img/boxplots/profileCurvatureMean.png")
boxplot(routesTest@data$meanProfCurv, main="Boxplot Mean profile Curvature", ylab="Curvature [1/100 m]") 
dev.off()
png("./img/boxplots/profileCurvatureMin.png")
boxplot(routesTest@data$minProfCurv, main="Boxplot Minimum profile Curvature", ylab="Curvature [1/100 m]") 
dev.off()
png("./img/boxplots/profileCurvatureMax.png")
boxplot(routesTest@data$maxProfCurv, main="Boxplot Maximum profile Curvature", ylab="Curvature [1/100 m]") 
dev.off()

png("./img/boxplots/straightness.png")
boxplot(routesTest@data$straightness, main="Boxplot Straightness", ylab="Straightness") 
dev.off()

# Aspect
library(RColorBrewer)
#cols <- colorRampPalette(brewer.pal(4,"Blues"))(length(aspectSums))
cols <- colorRampPalette(c(rgb(0.4,0.7,1,0.7),rgb(0,0.2,1,0.7),rgb(0.4,0.7,1,0.7)),alpha=T)(length(aspectSums)+1)
#cols <- c("#EFF3FF","#C9DEED","#93C2DE","#589ECD","#2171B5","#589ECD","#93C2DE","#C9DEED")
png("./img/boxplots/aspect.png")
radial.pie(aspectSums, labels=c('N','NW','W','SW','S','SE','E','NE'), clockwise = T,
           start = 0.625*pi, label.pos = seq(2.125,0.375,-0.25)*pi, radial.lim = c(0,1),
           boxed.radial = F, sector.colors = cols, title(main="Distribution Aspect"))
dev.off()

radial.plot(c(aspectSums,aspectSums[1]), labels=c('N','NW','W','SW','S','SE','E','NE'), clockwise = T,
           start = 0.5*pi, label.pos = seq(2,0.25,-0.25)*pi, radial.lim = c(0,1),
           boxed.radial = F, rp.type = "l", line.col = "blue", lwd = 2)


# Heading
png("./img/boxplots/heading.png")
radial.pie(headSums, labels=c('N','NW','W','SW','S','SE','E','NE'), clockwise = T,
           start = 0.625*pi, label.pos = seq(2.125,0.375,-0.25)*pi, radial.lim = c(0,1),
           boxed.radial = F, sector.colors = cols, main="Distribution Heading")
dev.off()
rm(cols)

############################################################################
# Test SHP ------------------------------------------------------------
############################################################################
# Select 10 random lines
#routesTester <- dbGetQuery(con, "SELECT route_id, tour_id, route_type, st_astext(route_path) FROM (SELECT * FROM route ORDER BY random() LIMIT 10) AS routeTest, border WHERE st_within(routeTest.route_path,border.geom);")
routesTester <- dbGetQuery(con, "SELECT route_id, tour_id, route_type, st_astext(route_path) FROM (SELECT * FROM route ORDER BY random() LIMIT 10000) AS routeTest;")
# Allocate space for lines
routesSpatialLines2 <- vector("list",length(routesTester[,4]))
# Convert routes from WKT to lines object
for(i in 1:length(routesTester[,4])){
  routesSpatialLines2[[i]] <- readWKT(routesTester[i,4],id=i)@lines
}
rm(i)
# Convert to SpatialLines object
routesSpatialLines2 <- SpatialLines(unlist(routesSpatialLines2),proj4string = CRS(CRS))
# Convert to SpatialLinesDataFrame
routesTester <- SpatialLinesDataFrame(routesSpatialLines2,routesTester[,1:3])
writeOGR(obj=routesTester, dsn=".", layer="routesTester2", driver="ESRI Shapefile")
rm(routesSpatialLines2,routesTester)

############################################################################
# Test Rasterize -----------------------------------------------------------
############################################################################
rasterList <- list()
for(i in 170:length(routesDF)){
  routesRas <- rasterize(routesDF[i,],dhm25)
  routesRas <- crop(routesRas,routesDF[i,])
  rasterList[[i]] <- routesRas
  print(i)
}
rm(i,routesRas)

# Overlay
test <- overlay(rasterList[[1]],crop(dhm25,rasterList[[1]]),fun=function(x,y){return(x*y)})
cellStats(test,min) #mean, min, max

mins <- vector("numeric",length(routesTest))
maxs <- vector("numeric",length(routesTest))
means <- vector("numeric",length(routesTest))
for(i in 1:1){
  profCurv <- overlay(rasterList[[i]],crop(profCurvNorm,rasterList[[i]]),fun=function(x,y){return(x*y)})
  planCurv <- overlay(rasterList[[i]],crop(planCurvNorm,rasterList[[i]]),fun=function(x,y){return(x*y)})
  mins[i] <- cellStats(profCurv,min)
  maxs[i] <- cellStats(profCurv,max)
  means[i] <- cellStats(profCurv,mean)
}

ptm <- proc.time()
mins <- vector("numeric",1000)#length(routesDF))
maxs <- vector("numeric",1000)#length(routesDF))
means <- vector("numeric",1000)#length(routesDF))
for(i in 1:1000){
  buffer <- suppressWarnings(buffer(routesDF[i,],0.001))
  dhm25_cr <- crop(dhm25,buffer)
  dhm25_cr_vx <- velox(dhm25_cr)
  test_cr_vx <- unlist(dhm25_cr_vx$extract(routesDF[i,]))
  mins[i] <- min(test_cr_vx)
  maxs[i] <- max(test_cr_vx)
  means[i] <- mean(test_cr_vx)
  if(i%%100 == 0){
    cat(i,"lines processed \n")
  }
}
proc.time() - ptm

routesTest$minHeight <- mins
routesTest$maxHeight <- maxs
routesTest$meanHeight <- means
