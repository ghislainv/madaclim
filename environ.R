##=====================================================
## Environmental data for Madagascar
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## July 2015
##=====================================================

##= libraries
library(sp)
library(raster)
library(rgdal)

##= gdalwrap options
Extent <- "298000 7155000 1101000 8683000"
Res <- "1000"
nodat <- "-9999"
proj.s <- "EPSG:4326"
proj.t <- "EPSG:32738"

##==========================================================
##
##= Vegetation map by Kew
##
##==========================================================

## ## VegGeol
## kew <- readOGR("/media/ghislain/BioSceneMada2/gisdata/mada/vegetation_kew/","VEGGEOL")
## projection(kew)
## kew.38S <- spTransform(kew,CRS("+init=epsg:32738"))
## writeOGR(kew.38S,"/media/ghislain/BioSceneMada2/gisdata/mada/vegetation_kew/",layer="kew_38S",driver="ESRI Shapefile")
## system("gdal_rasterize -ot Float -a_nodata -32768 -te 298000 7155000 1101000 8683000 \\
##         -tr 1000 1000 -a ID -l kew_38S \\
##         /media/ghislain/BioSceneMada2/gisdata/mada/vegetation_kew/kew_38S.shp \\
##         /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/kew_1km.tif")

## VegMada
system("gdalwarp -overwrite -srcnodata 0 -dstnodata -32768 -ot Int16 -s_srs EPSG:32738 -t_srs EPSG:32738 \\
        -r near -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        /media/ghislain/BioSceneMada2/gisdata/mada/vegmada_kew/vegetation.tif \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/vegmada_1km.tif")

##==========================================================
##
##= forest: percentage of forest in 1km2
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
##= altitude, slope, aspect and solar radiation
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

## geology
geol.latlong <- readOGR("/media/ghislain/BioSceneMada2/gisdata/mada/geolsimp/","geolsimp")
crs(geol.latlong) <- "+init=epsg:4326"
geol <- spTransform(geol.latlong,CRS("+init=epsg:32738"))
writeOGR(geol,dsn="/media/ghislain/BioSceneMada2/gisdata/mada/geolsimp/",layer="geolsimp_38S",driver="ESRI Shapefile")
system("gdal_rasterize -ot Int16 -a_nodata -32768 -te 298000 7155000 1101000 8683000 \\
        -tr 1000 1000 -a RECLASS_ID -l geolsimp_38S \\
        /media/ghislain/BioSceneMada2/gisdata/mada/geolsimp/geolsimp_38S.shp \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/geol_1km.tif")
# Attribute table
dgeol <- geol@data
tgeol <- as.data.frame(tapply(as.character(dgeol$RECLASS_NA),dgeol$RECLASS_ID,unique))
names(tgeol) <- "type"
sink("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/geol_attribute.txt")
tgeol
sink()

## soil
soil.latlong <- readOGR("/media/ghislain/BioSceneMada2/gisdata/mada/Soil_morphopedo/shape/","geomorph_wgs")
crs(soil.latlong) <- "+init=epsg:4326"
soil <- spTransform(soil.latlong,CRS("+init=epsg:32738"))
soil$SOLDT_ID <- as.numeric(soil$SOLDT)
writeOGR(soil,dsn="/media/ghislain/BioSceneMada2/gisdata/mada/Soil_morphopedo/shape/",layer="soil_38S",driver="ESRI Shapefile")
system("gdal_rasterize -ot Int16 -a_nodata -32768 -te 298000 7155000 1101000 8683000 \\
        -tr 1000 1000 -a SOLDT_ID -l soil_38S \\
        /media/ghislain/BioSceneMada2/gisdata/mada/Soil_morphopedo/shape/soil_38S.shp \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/soil_1km.tif")

## watershed
wshed <- raster("/media/ghislain/BioSceneMada2/gisdata/mada/Watersheds/EVO_596_sm_AppendixS2.tif")
system("gdalwarp -overwrite -ot Int16 -srcnodata 255 -dstnodata -32768 -s_srs EPSG:4326 -t_srs EPSG:32738 \\
        -r near -tr 1000 1000 -te 298000 7155000 1101000 8683000 -of GTiff \\
        /media/ghislain/BioSceneMada2/gisdata/mada/Watersheds/EVO_596_sm_AppendixS2.tif \\
        /home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/environ/wshed_1km.tif")

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
