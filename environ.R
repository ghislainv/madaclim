##=====================================================
##
## Environmental data for Madagascar
##
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
##
## March 2016
##
##=====================================================

## gdal library is needed to run this script
## http://www.gdal.org/

## Load libraries
library(sp)
library(raster)
library(rgdal)

## Create some directories
dir.create("gisdata") ## To save GIS data
dir.create("environ") ## To save 1km environmental data
dir.create("temp") ## Temporary folder

## gdalwrap options
Extent <- "298000 7155000 1101000 8683000"
Res <- "1000"
nodat <- "-9999"
proj.s <- "EPSG:4326"
proj.t <- "EPSG:32738"

##==========================================================
##
## Simplified geological map by Kew
##
##==========================================================

## Reference: Kew (http://www.kew.org/gis/projects/madagascar/download.html)

## Create directory
dir.create("gisdata/geol")
## Download and unzip
url.geol <- "http://www.kew.org/gis/projects/madagascar/downloads/geolsimp.zip"
download.file(url=url.geol,destfile="gisdata/geol/geolsimp.zip",method="wget",quiet=TRUE)
unzip("gisdata/geol/geolsimp.zip",exdir="temp/geol",overwrite=TRUE)
## Reproject in UTM 38S
geol.latlong <- readOGR("temp/geol","geolsimp")
crs(geol.latlong) <- "+init=epsg:4326"
geol <- spTransform(geol.latlong,CRS("+init=epsg:32738"))
writeOGR(geol,dsn="temp/geol",layer="geolsimp_38S",driver="ESRI Shapefile")
## Rasterize with gdal
system("gdal_rasterize -ot Int16 -a_nodata -32768 -te 298000 7155000 1101000 8683000 \\
        -tr 1000 1000 -a RECLASS_ID -l geolsimp_38S \\
        temp/geol/geolsimp_38S.shp \\
        environ/geol_1km.tif")
## Attribute table
dgeol <- geol@data
tgeol <- as.data.frame(tapply(as.character(dgeol$RECLASS_NA),dgeol$RECLASS_ID,unique))
names(tgeol) <- "type"
sink("environ/geol_attribute.txt")
tgeol
sink()

##==========================================================
##
## Watersheds
##
##==========================================================

## Reference: Pearson, R. G., & Raxworthy, C. J. (2009). The evolution of local
## endemism in Madagascar: watershed versus climatic gradient hypotheses
## evaluated by null biogeographic models. Evolution, 63(4), 959-967.

## Create directory
dir.create("gisdata/watersheds")
## Download and unzip
url.wshed <- "http://onlinelibrary.wiley.com/store/10.1111/j.1558-5646.2008.00596.x/asset/supinfo/EVO_596_sm_AppendixS2.tif?v=1&s=62c6690369a2b83a224560f08be91462993d485d"
download.file(url=url.wshed,destfile="gisdata/watersheds/EVO_596_sm_AppendixS2.tif",method="wget",quiet=TRUE)
system("gdalwarp -overwrite -ot Int16 -srcnodata 255 -dstnodata -32768 -s_srs EPSG:4326 -t_srs EPSG:32738 \\
        -r near -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        gisdata/watersheds/EVO_596_sm_AppendixS2.tif \\
        environ/wshed_1km.tif")

##==========================================================
##
## Soil
##
##==========================================================

## Reference: Delenne M., Pelletier F., 1981. Carte du Potentiel des Unités
## Physiques de Madagascar, au 1:1 000 000e, ORSTOM, Bondy

unzip("gisdata/soil/morphopedo_mada.zip",exdir="temp/soil",overwrite=TRUE)
soil.latlong <- readOGR("temp/soil/shape/","geomorph_wgs")
crs(soil.latlong) <- "+init=epsg:4326"
soil <- spTransform(soil.latlong,CRS("+init=epsg:32738"))
soil$SOLDT_ID <- as.numeric(soil$SOLDT)
writeOGR(soil,dsn="temp/soil/shape/",layer="soil_38S",driver="ESRI Shapefile")
system("gdal_rasterize -ot Int16 -a_nodata -32768 -te 298000 7155000 1101000 8683000 \\
        -tr 1000 1000 -a SOLDT_ID -l soil_38S \\
        temp/soil/shape/soil_38S.shp \\
        environ/soil_1km.tif")

##==========================================================
##
## Vegetation map by Kew
##
##==========================================================

## Reference: The CEFP Madagascar Vegetation mapping project (http://vegmad.org)

## Create directory
dir.create("gisdata/vegmada_kew")
## Download and unzip
url.veggeol <- "http://www.vegmad.org/downloads/utm/veg_tif.zip"
download.file(url=url.veggeol,destfile="gisdata/vegmada_kew/vegmada.zip",method="wget",quiet=TRUE)
unzip("gisdata/vegmada_kew/vegmada.zip",exdir="temp/veg",overwrite=TRUE)
## VegMada
system("gdalwarp -overwrite -srcnodata 0 -dstnodata -32768 -ot Int16 -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r near -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        temp/veg/vegetation.tif \\
        environ/vegmada_1km.tif")

##==========================================================
##
## Forest: percentage of forest in 1km2
##
##==========================================================

## Set region on CI map
system("g.region rast=Forest_CI_905_30m -ap")
## Create Mada raster at 30m
system("r.mapcalc 'Mada=if(!isnull(Forest_CI_905_30m),1,null())'")

## Create grid with 1 x 1 km cells
system("g.region w=298000 e=1101000 s=7155000 n=8683000 res=1000 -ap")
## Use r.resamp.stat
system("r.resamp.stats input=for2010_30m output=forest_n method=sum")
system("r.resamp.stats input=Mada output=land_n method=sum")
## Mada 1km if land_n > 555 (more than 50% of the 1km cell is covered by land)
system("r.mapcalc 'mada_1km=Mada'")
## correction Mada if land_n > 555 (more than 50% of the 1km cell is covered by land)
system("r.mapcalc 'mada_1km=if(land_n>555 &&& !isnull(land_n),1,mada_1km)'")
## Percentage of forest
system("r.mask -o 'mada_1km'")
system("r.mapcalc 'percfor10=round(100*forest_n/land_n)'")
system("r.info percfor10")
system("r.mask -r")
system("r.mapcalc 'percfor10=if(!isnull(mada_1km) && isnull(forest_n),0,percfor10)'")

## Export
system("r.out.gdal input=percfor10 output=percfor2010.tif type=UInt16 createopt='compress=lzw,predictor=2'")

##==========================================================
##
## Altitude, Slope, Aspect and Solar radiation
##
##==========================================================

## ## ## location: frb.defor.mada
## system("g.region rast=altitude -ap")
## system("r.out.gdal input=altitude output=/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/altitude.tif \\
##         type=Int16 createopt='compress=lzw,predictor=2'")
## system("r.out.gdal input=slope output=/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/slope.tif \\
##         type=Int16 createopt='compress=lzw,predictor=2'")
## system("r.out.gdal input=aspect output=/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/aspect.tif \\
##         type=Int16 createopt='compress=lzw,predictor=2'")
## system("r.out.gdal input=global_rad output=/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/solar.tif \\
##         type=Int16 createopt='compress=lzw,predictor=2'")

## Resample to grid
## altitude
system("gdalwarp -overwrite -srcnodata -32768 -dstnodata -9999 -ot Int16 -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/altitude.tif \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/altitude_1km.tif")
## slope
system("gdalwarp -overwrite -srcnodata -32768 -dstnodata -9999 -ot Int16 -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/slope.tif \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/slope_1km.tif")
## aspect
system("gdalwarp -overwrite -srcnodata -32768 -dstnodata -9999 -ot Int16 -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/aspect.tif \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/aspect_1km.tif")
## solar radiation
system("gdalwarp -overwrite -srcnodata -32768 -dstnodata -9999 -ot Int16 -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/solar.tif \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/solar_1km.tif")

## Import
altitude <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/altitude_1km.tif")
slope <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/slope_1km.tif")
aspect <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/aspect_1km.tif")
## Caution: in GRASS, aspect is calculated counterclockwise from east in degrees.
## Transformation to clockwise from north in degrees
asp.north <- (360+(90-values(aspect)))%%360
values(aspect) <- asp.north
solar <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/solar_1km.tif")
geol <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/geol_1km.tif")
soil <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/soil_1km.tif")
veg <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/vegmada_1km.tif")
wshed <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/wshed_1km.tif")
percfor2010 <- raster("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/percfor2010.tif")

## Stack and output
s.env <- stack(altitude,slope,aspect,solar,geol,soil,veg,wshed,percfor2010)
writeRaster(s.env,filename="environ.tif",overwrite=TRUE,
            datatype="INT2S",format="GTiff",options=c("COMPRESS=LZW","PREDICTOR=2"))

## Plot
s <- stack("environ.tif")
names(s) <- c("alt","slop","asp","solar","geol","soil","veg","wshed","percfor2010")
## All layers
png("environ.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(5,4,3,2))
plot(s)
dev.off()
## Altitude
png("altitude.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(s$alt,col=terrain.colors(255))
dev.off()
## Watershed
png("watershed.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(s$wshed,col=terrain.colors(255))
dev.off()

##= END
