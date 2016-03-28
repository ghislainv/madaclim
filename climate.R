##=====================================================
## Climate data for Madagascar
## Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
## July 2015
##=====================================================

##= Set working directory
setwd("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/climate/")

##= libraries
library(sp)
library(raster)
library(dismo) # for function biovars()
library(insol) # for function daylength

#== Thornthwaite functions
thorn.indices <- function(clim) {
    I <- rep(0,ncell(clim))
    for (i in 1:12) {
        Tmin <- values(subset(clim,i))/10
        Tmax <- values(subset(clim,i+12))/10
        Tmean <- (Tmin+Tmax)/2
        I <- I+(Tmean/5)^(1.514)
    }
    alpha <- (6.75e-7)*I^3-(7.71e-5)*I^2+(1.792e-2)*I+0.49239
    return(list(I=I,alpha=alpha))
}
thorn.f <- function(Tm,I,alpha,Jday,lat,long) {
    L <- daylength(lat,long,Jday,tmz=3)[,3]
    PET <- 1.6*(L/12)*(10*Tm/I)^alpha
    return(PET)
}

##= Variables
yr <- c("2050","2080") # For 2050, 2080
## mod.ccafs <- c("ipsl_cm5a_lr","miroc_miroc5","ncc_noresm1_m") # For global climate models (GCMs): GISS-E2-R, HadGEM2-ES, CCSM4
mod.ccafs <- c("csiro_access1_0","giss_e2_r","ipsl_cm5a_lr","miroc_miroc5","mohc_hadgem2_es","ncar_ccsm4","ncc_noresm1_m")
## mod <- c("ip","mc","no")
mod <- c("ac","gs","ip","mc","he","cc","no")
rcp.ccafs <- c("4_5","8_5") # For representative concentration pathways (RCPs): RCP 4.5, RCP 8.5
rcp <- c("45","85")
var <- c("tmin","tmax","prec")

##= gdalwrap options
Extent <- "298000 7155000 1101000 8683000"
Res <- "1000"
nodat <- "-9999"
proj.s <- "EPSG:4326"
proj.t <- "EPSG:32738"

##= function to compute PET, CWD and NDM
## PET: potential evapotranspiration (Thornthwaite equation,1948)
## CWD: climatic water deficit
## NDM: number of dry months
pet.cwd.ndm.f <- function(clim) {
    # get latitude in radians
    xy.utm <- SpatialPoints(coordinates(clim), proj4string=CRS("+init=epsg:32738"))
    xy <- spTransform(xy.utm,CRS("+init=epsg:4326"))
    long_deg <- coordinates(xy)[,1]
    lat_deg <- coordinates(xy)[,2]
    # initialize
    cwd <- rep(0,ncell(clim)) 
    ndm <- rep(0,ncell(clim))
    pet <- rep(0,ncell(clim))
    # thorn.index
    ind <- thorn.indices(clim)
    # loop on months
    for (i in 1:12) {
        cat(paste("Month: ",i,"\n",sep=""))
        evap.thorn <- clim[[1]] # Evap Thornthwaite
        Tmin <- values(subset(clim,i))/10
        Tmax <- values(subset(clim,i+12))/10
        Prec <- values(subset(clim,i+24))
        d <- data.frame(day=(30*i)-15,Tmin,Tmax,lat_deg,long_deg)
        d[is.na(d)] <- 0
        ## Thornthwaite
        pet.thorn <- thorn.f(Tm=(d$Tmax+d$Tmin)/2,lat=d$lat_deg,long=d$long_deg,
                             I=ind$I,alpha=ind$alpha,Jday=d$day)*10
        pet.thorn[is.na(Tmin)] <- NA # to correct for NA values
        values(evap.thorn) <- pet.thorn
        if (i==1) {
            PET12.thorn <- stack(evap.thorn)
        }
        if (i>1) {
            PET12.thorn <- addLayer(PET12.thorn, evap.thorn)
        }
        pet <- pet+pet.thorn # annual PET
        pe.diff <- Prec-pet.thorn
        cwd <- cwd+pmin(pe.diff,0.0) # climatic water deficit
        dm <- rep(0,ncell(clim)) # dry month
        dm[pe.diff<0] <- 1
        ndm <- ndm+dm
    }
    # make rasters
    PET <- CWD <- NDM <- clim[[1]]
    values(PET) <- pet
    values(CWD) <- -cwd
    values(NDM) <- ndm
    NDM[is.na(PET)] <- NA # to account for NA values
    return (list(PET12=PET12.thorn,PET=PET,CWD=CWD,NDM=NDM))
}

##==================================================================================
##
## Present
## - WorlClim data
## - http://www.worldclim.org/current
##
##==================================================================================

##= Source folder
## sf.pres <- "/media/ghislain/BioSceneMada2/gisdata/mada/worldclim/"
worldclim.folder <- "/home/ghislain/gisdata/mada/worldclim/"

## Loop on variables
for (k in 1:length(var)) {
    ## unzip
    nf.ccafs.zip <- paste(worldclim.folder,var[k],"_37_tif.zip",sep="")
    ffile <- unzip(nf.ccafs.zip,junkpaths=TRUE)
    ffile <- ffile[c(1,5:12,2:4)]
    ## gdalwrap
    for (f in 1:length(ffile)) {
        Input <- ffile[f]
        Output <- paste(var[k],"_",f,".tif",sep="")
        system(paste("gdalwarp -ot Int16 -dstnodata ",nodat," -s_srs ",proj.s," -t_srs ",proj.t,
                     " -te ",Extent," -tr ",Res," ",Res," -r bilinear -overwrite ",Input," ",Output,sep=""))
    }
    ## clean
    system("rm ./*_37.tif")               
}

## stack
tn <- paste("tmin_",1:12,".tif",sep="")
tx <- paste("tmax_",1:12,".tif",sep="")
pr <- paste("prec_",1:12,".tif",sep="")
s <- stack(c(tn,tx,pr))
## bioclim
b <- biovars(tmin=subset(s,1:12),tmax=subset(s,13:24),prec=subset(s,25:36))
## pet.cwd.ndm
pet.cwd.ndm <- pet.cwd.ndm.f(s,slope,aspect)
## comparing Thorthwaite and Priestley-Taylor

## output stack
os <- stack(s,b,pet.cwd.ndm$PET12,pet.cwd.ndm$PET,pet.cwd.ndm$CWD,pet.cwd.ndm$NDM)
writeRaster(os,filename="current.tif",overwrite=TRUE,
            datatype="INT2S",format="GTiff",options=c("COMPRESS=LZW","PREDICTOR=2"))
## clean
system("rm ./tmin_*.tif")
system("rm ./tmax_*.tif")
system("rm ./prec_*.tif")

##= Plot
## PET
os <- stack("current.tif")
names(os) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),paste("prec",1:12,sep=""),
               paste("bio",1:19,sep=""),paste("pet",1:12,sep=""),"pet","cwd","ndm")
png("pet.png",width=600,height=600,res=72,pointsize=16)
plot(subset(os,56:67),zlim=c(0,210))
dev.off()

## NDM
png("ndm.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(os$ndm,col=terrain.colors(255))
dev.off()

## CWD
png("cwd.png",width=600,height=600,res=72,pointsize=16)
par(mar=c(3,3,1,1))
plot(os$cwd,col=terrain.colors(255))
dev.off()

##==================================================================================
##
## Future
## - CGIAR CCAFS data
## - http://www.ccafs-climate.org/data/
##
##==================================================================================

##= Import data from CGIAR CCAFS website
## ccafs.folder <- "/media/ghislain/BioSceneMada2/gisdata/ccafs.climate.tntxpr/"
ccafs.folder <- "/home/ghislain/gisdata/ccafs.climate.tntxpr/"
setwd(ccafs.folder)

##= Loop on models, scenarios, years and variables
for (i in 1:length(mod)) {
    for (j in 1:length(rcp)) {
        for (l in 1:length(yr)) {
            for (k in 1:length(var)) {
                target <- paste("http://cgiardata.s3.amazonaws.com/ccafs/ccafs-climate/data/ipcc_5ar_ciat_tiled/rcp",
                                rcp.ccafs[j],"/",yr[l],"s/",mod.ccafs[i],"/30s/",mod.ccafs[i],"_rcp",rcp.ccafs[j],
                                "_",yr[l],"s_",var[k],"_30s_r1i1p1_c4_asc.zip",sep="")
                system(paste("wget ",target,sep=""))
            }
        }
    }
}

## Reset wd
setwd("/home/ghislain/Documents/Ghislain-CIRAD/FRB_Mada/madaclim/climate/")

##= Loop on models, scenarios and years
for (i in 1:length(mod)) {
    for (j in 1:length(rcp)) {
        for (l in 1:length(yr)) {
            for (k in 1:length(var)) {
                ## unzip
                nf.ccafs.zip <- paste(ccafs.folder,mod.ccafs[i],"_rcp",rcp.ccafs[j],"_",yr[l],
                                      "s_",var[k],"_30s_r1i1p1_c4_asc.zip",sep="")
                ffile <- unzip(nf.ccafs.zip,junkpaths=TRUE)
                if (length(grep(".prj",ffile))>0) { ffile <- ffile[-grep(".prj",ffile)] }
                ffile <- ffile[c(1,5:12,2:4)]
                ## gdalwrap
                for (f in 1:length(ffile)) {
                    Input <- ffile[f]
                    Output <- paste(var[k],"_",f,".tif",sep="")
                    system(paste("gdalwarp -ot Int16 -dstnodata ",nodat," -s_srs ",proj.s," -t_srs ",proj.t,
                                 " -te ",Extent," -tr ",Res," ",Res," -r bilinear -overwrite ",Input," ",Output,sep=""))
                }
                ## clean
                system("rm ./*.asc")
                if (length(grep(".prj",ffile))>0) { system("rm ./*.prj") } 
            }
            # stack
            tn <- paste("tmin_",1:12,".tif",sep="")
            tx <- paste("tmax_",1:12,".tif",sep="")
            pr <- paste("prec_",1:12,".tif",sep="")
            s <- stack(c(tn,tx,pr))
            # bioclim
            b <- biovars(tmin=subset(s,1:12),tmax=subset(s,13:24),prec=subset(s,25:36))
            # pet.cwd.ndm
            pet.cwd.ndm <- pet.cwd.ndm.f(s)
            # output stack
            os <- stack(s,b,pet.cwd.ndm$PET12,pet.cwd.ndm$PET,pet.cwd.ndm$CWD,pet.cwd.ndm$NDM)
            names(os) <- c(paste("tmin",1:12,sep=""),paste("tmax",1:12,sep=""),paste("prec",1:12,sep=""),
                           paste("bio",1:19,sep=""),paste("pet",1:12,sep=""),"pet","cwd","ndm")
            writeRaster(os,filename=paste(mod[i],"_",rcp[j],"_",yr[l],".tif",sep=""),overwrite=TRUE,
                        datatype="INT2S",format="GTiff",options=c("COMPRESS=LZW","PREDICTOR=2"))
            # clean
            system("rm ./tmin_*.tif")
            system("rm ./tmax_*.tif")
            system("rm ./prec_*.tif")
        }
    }
}

##= END
