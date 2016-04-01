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

## GRASS GIS 7.x.x is also needed to run this script
## https://grass.osgeo.org/

## Read argument for download
## Set "down" to TRUE if you want to download the sources. Otherwise, the data already provided in the gisdata repository will be used.
arg <- commandArgs(trailingOnly=TRUE)
down <- FALSE
if (length(arg)>0) {
  down <- arg[1]
}

## Load libraries
library(sp)
library(raster)
library(rgdal)
library(rgrass7)

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

## Reference: 
## Kew: http://www.kew.org/gis/projects/madagascar/download.html

## Create directory
dir.create("gisdata/geol")
## Download and unzip
if (down) {
  url.geol <- "http://www.kew.org/gis/projects/madagascar/downloads/geolsimp.zip"
  download.file(url=url.geol,destfile="gisdata/geol/geolsimp.zip",method="wget",quiet=TRUE)
}
unzip("gisdata/geol/geolsimp.zip",exdir="temp/geol",overwrite=TRUE)
## Reproject in UTM 38S
geol.latlong <- readOGR("temp/geol","geolsimp")
crs(geol.latlong) <- "+init=epsg:4326"
geol <- spTransform(geol.latlong,CRS("+init=epsg:32738"))
writeOGR(geol,dsn="temp/geol",layer="geolsimp_38S",driver="ESRI Shapefile",overwrite_layer=TRUE)
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

## Reference: 
## Pearson, R. G., & Raxworthy, C. J. (2009). The evolution of local endemism in
## Madagascar: watershed versus climatic gradient hypotheses evaluated by null
## biogeographic models. Evolution, 63(4), 959-967.

## Create directory
dir.create("gisdata/watersheds")
## Download
if (down) {
  url.wshed <- "http://onlinelibrary.wiley.com/store/10.1111/j.1558-5646.2008.00596.x/asset/supinfo/EVO_596_sm_AppendixS2.tif?v=1&s=62c6690369a2b83a224560f08be91462993d485d"
  download.file(url=url.wshed,destfile="gisdata/watersheds/EVO_596_sm_AppendixS2.tif",method="wget",quiet=TRUE)
}
## Resample
system("gdalwarp -overwrite -ot Int16 -srcnodata 255 -dstnodata -32768 -s_srs EPSG:4326 -t_srs EPSG:32738 \\
        -r near -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        gisdata/watersheds/EVO_596_sm_AppendixS2.tif \\
        environ/wshed_1km.tif")

##==========================================================
##
## Soil
##
##==========================================================

## Reference: 
## Delenne M., Pelletier F., 1981. Carte du Potentiel des Unités Physiques de
## Madagascar, au 1:1 000 000e, ORSTOM, Bondy

unzip("gisdata/soil/morphopedo_mada.zip",exdir="temp/soil",overwrite=TRUE)
soil.latlong <- readOGR("temp/soil/shape/","geomorph_wgs")
crs(soil.latlong) <- "+init=epsg:4326"
soil <- spTransform(soil.latlong,CRS("+init=epsg:32738"))
soil$SOLDT_ID <- as.numeric(soil$SOLDT)
writeOGR(soil,dsn="temp/soil/shape/",layer="soil_38S",driver="ESRI Shapefile",overwrite_layer=TRUE)
system("gdal_rasterize -ot Int16 -a_nodata -32768 -te 298000 7155000 1101000 8683000 \\
        -tr 1000 1000 -a SOLDT_ID -l soil_38S \\
        temp/soil/shape/soil_38S.shp \\
        environ/soil_1km.tif")

##==========================================================
##
## Vegetation map by Kew
##
##==========================================================

## Reference: 
## The CEFP Madagascar Vegetation mapping project: http://vegmad.org

## Create directory
dir.create("gisdata/vegmada_kew")
## Download and unzip
if (down) {
  url.veggeol <- "http://www.vegmad.org/downloads/utm/veg_tif.zip"
  download.file(url=url.veggeol,destfile="gisdata/vegmada_kew/vegmada.zip",method="wget",quiet=TRUE)
}
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

## References:
## (1) Hansen, M. C., Potapov, P. V., Moore, R., Hancher, M., Turubanova, S. A.,
## Tyukavina, A., ... & Kommareddy, A. (2013). High-resolution global maps of
## 21st-century forest cover change. Science, 342(6160), 850-853.
## (2) Harper, G. J., Steininger, M. K., Tucker, C. J., Juhn, D., & Hawkins, F.
## (2007). Fifty years of deforestation and forest fragmentation in Madagascar.
## Environmental Conservation, 34(4), 325.
## (3) BioSceneMada project: http://bioscenemada.net

## Create directory
dir.create("gisdata/forest")
## Download
if (down) {
  url.for2010 <- "http://bioscenemada.net/FileTransfer/for2010.tif"
  url.forCI.905 <- "http://bioscenemada.net/FileTransfer/Forest_CI_905.tif"
  download.file(url=url.for2010,destfile="gisdata/forest/for2010.tif",method="wget",quiet=TRUE)
  download.file(url=url.forCI.905,destfile="gisdata/forest/Forest_CI_905.tif",method="wget",quiet=TRUE)
}
## Create new grass location in UTM 38S
dir.create("grassdata")
system("grass70 -c epsg:32738 grassdata/environ.mada")
## Connect R to grass location
initGRASS(gisBase="/usr/local/grass-7.0.1",home=tempdir(), 
          gisDbase="grassdata",
          location="environ.mada",mapset="PERMANENT",
          override=TRUE)
## Import rasters in grass
system("r.in.gdal input=gisdata/forest/Forest_CI_905.tif output=Forest_CI_905")
system("r.in.gdal input=gisdata/forest/for2010.tif output=for2010")
## Set region on CI map
system("g.region rast=Forest_CI_905 -ap")
## Create Mada raster at 30m
system("r.mapcalc 'Mada=if(!isnull(Forest_CI_905),1,null())'")

## Create grid with 1 x 1 km cells
system("g.region w=298000 e=1101000 s=7155000 n=8683000 res=1000 -ap")
## Use r.resamp.stat
system("r.resamp.stats input=for2010 output=forest_n method=sum")
system("r.resamp.stats input=Mada output=land_n method=sum")
## Mada 1km
system("r.mapcalc 'mada_1km = Mada'")
## Correction Mada if land_n > 555 (more than 50% of the 1km cell is covered by land)
## Indeed, 555*30*30 = 499500 km2 and 556*30*30 = 500400 km2
system("r.mapcalc --o 'mada_1km = if(land_n>555 &&& !isnull(land_n),1,mada_1km)'")
## Percentage of forest
system("r.mask 'mada_1km'")
system("r.mapcalc 'percfor10 = round(100*forest_n/land_n)'")
system("r.mask -r")
system("r.mapcalc --o 'percfor10 = if(!isnull(mada_1km) && isnull(forest_n),0,percfor10)'")

## Export
system("r.out.gdal --o input=percfor10 output=environ/percfor2010.tif type=UInt16 createopt='compress=lzw,predictor=2'")

##==========================================================
##
## Altitude, Slope, Aspect and Solar radiation
##
##==========================================================

## Download and unzip CGIAR-CSI 90m DEM data
tiles <- c("45_18","46_18","45_17","46_17","45_16","46_16","47_16","46_15","47_15")
for (i in 1:length(tiles)) {
  dst <- paste0("gisdata/altitude/srtm_",tiles[i],".zip")
  if (down) {
    url.tile <- paste0("http://srtm.csi.cgiar.org/SRT-ZIP/SRTM_V41/SRTM_Data_GeoTiff/srtm_",tiles[i],".zip")
    download.file(url=url.tile,destfile=dst,method="wget",quiet=TRUE)
  }
  unzip(dst,exdir="temp/altitude",overwrite=TRUE)
}

## Mosaic with gdalbuildvrt
system("gdalbuildvrt temp/altitude/altitude.vrt temp/altitude/*.tif")
## Resample (dstnodata need to be set to 32767 as we pass from Int16 (nodata=-32768) to INT2S (nodata=-32767) in R)
system("gdalwarp -overwrite -s_srs EPSG:4326 -t_srs EPSG:32738 -srcnodata -32768 -dstnodata -32767 \\
        -r bilinear -tr 90 90 -te 298000 7155000 1101000 8683000 -ot Int16 -of GTiff \\
        temp/altitude/altitude.vrt \\
        temp/altitude/altitude.tif")

## Import raster in grass
system("r.in.gdal --o input=temp/altitude/altitude.tif output=altitude")

## Compute slope, aspect and global radiation in grass
system("g.region rast=altitude -ap")
system("r.slope.aspect --o elevation=altitude slope=slope aspect=aspect format=degrees")
## Computing radiation with r.sun at 90m resolution (very long to run: ~ some hours, so we set res=1000)
system("g.region rast=altitude res=1000 -ap")
system("r.sun --o --verbose elevation=altitude aspect=aspect slope=slope day=79 glob_rad=global_rad")

## Export
system("g.region rast=altitude -ap")
system("r.out.gdal --o input=altitude output=temp/altitude/altitude.tif \\
        type=Int16 nodata=-32767 createopt='compress=lzw,predictor=2'")
system("r.out.gdal -f --o input=slope output=temp/altitude/slope.tif \\
        type=Int16 nodata=-32767 createopt='compress=lzw,predictor=2'")
system("r.out.gdal -f --o input=aspect output=temp/altitude/aspect.tif \\
        type=Int16 nodata=-32767 createopt='compress=lzw,predictor=2'")
system("r.out.gdal -f --o input=global_rad output=temp/altitude/solar.tif \\
        type=Int16 nodata=-32767 createopt='compress=lzw,predictor=2'")

## Resample to grid
## altitude
system("gdalwarp -overwrite -ot Int16 -srcnodata -32767 -dstnodata -32767 \\
        -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        temp/altitude/altitude.tif \\
        environ/altitude_1km.tif")
## slope
system("gdalwarp -overwrite -ot Int16 -srcnodata -32767 -dstnodata -32767 \\
        -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        temp/altitude/slope.tif \\
        environ/slope_1km.tif")
## aspect
system("gdalwarp -overwrite -ot Int16 -srcnodata -32767 -dstnodata -32767 \\
        -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        temp/altitude/aspect.tif \\
        environ/aspect_1km.tif")
## solar radiation
system("gdalwarp -overwrite -ot Int16 -srcnodata -32767 -dstnodata -32767 \\
        -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r bilinear -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        temp/altitude/solar.tif \\
        environ/solar_1km.tif")

## Import
altitude <- raster("environ/altitude_1km.tif")
slope <- raster("environ/slope_1km.tif")
aspect <- raster("environ/aspect_1km.tif")
## Caution: in GRASS, aspect is calculated counterclockwise from east in degrees.
## Transformation to clockwise from north in degrees
asp.north <- (360+(90-values(aspect)))%%360
values(aspect) <- asp.north
solar <- raster("environ/solar_1km.tif")
geol <- raster("environ/geol_1km.tif")
soil <- raster("environ/soil_1km.tif")
veg <- raster("environ/vegmada_1km.tif")
wshed <- raster("environ/wshed_1km.tif")
percfor2010 <- raster("environ/percfor2010.tif")

## Stack and export
s.env <- stack(altitude,slope,aspect,solar,geol,soil,veg,wshed,percfor2010)
## Remove data for Comoro Islands
bbCom <- extent(xmin(s.env),600000,8500000,ymax(s.env)) # bounding-box
cellsCom <- cellsFromExtent(s.env,bbCom)
values(s.env)[cellsCom,] <- NA
s.env <- stack(s.env) # Transform back from RasterBrick to RasterStack
## Export
writeRaster(s.env,filename="environ/environ.tif",overwrite=TRUE,
            datatype="INT2S",format="GTiff",options=c("COMPRESS=LZW","PREDICTOR=2"))

## Plot
s <- stack("environ/environ.tif")
names(s) <- c("alt","slop","asp","solar","geol","soil","veg","wshed","percfor2010")
## All layers
png("environ/environ.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(5,4,3,1))
plot(s,legend.mar=7)
dev.off()
## Altitude
png("environ/altitude.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(s$alt,col=terrain.colors(255))
dev.off()
## Watershed
png("environ/watershed.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(s$wshed,col=terrain.colors(255))
dev.off()

##==========================================================
##
## Clean directory if necessary
##
##========================================================== 

## unlink("grassdata",recursive=TRUE)
## unlink("temp",recursive=TRUE)

##========================================
## End of script
##========================================
