rm(list=ls())
library(raster)
library(mgcv)

## simulated fishery observations
loc <- read.csv("supp_csp_sim.csv",stringsAsFactors = T)
## amon = layer position in the following raster object

## a list of seascape raster stacks
load(file="supp_csp_env.rda")
## each layer is one month
## bathy= bathymetry
## sst = sea surface temp.
## sst2 = quadratic sst
## fpi = frontal probability index
## ssh = sea surface height
## chl = chlorophyll-a index

library(leaflet)
library(mapview)
library(sp)
#library(sp708M)
source("visio_2.R")

dummy <- loc
coordinates(dummy) <- ~long+lat
m <- visio(obj=dummy,leg.col="Set1",field="group",
           legend.title = "Organizations")
mapshot(m,file="spot_raw.png")


#' Extract environmental data from rasters used in model output

## change to a list from each month 
track_list <- split(x=loc,f=loc$amon,drop=TRUE)
env_list <- vector("list",length(track_list))
for(i in 1:length(track_list)){
  ## locations
  xy <- as.matrix(track_list[[i]][,c("long","lat"),drop=FALSE])
  ## month
  imonth <- track_list[[i]][1,"amon"]
  cat("extracting data from ",names(cov__[[1]])[imonth],".\n",sep="")
  #' bathy
  l.bathy <- extract(cov__[[1]][[imonth]],xy)
  #' SST
  l.sst <- extract(cov__[[2]][[imonth]],xy)
  #' SSH
  l.ssh <- extract(cov__[[5]][[imonth]],xy)
  #' chla
  l.chla <- extract(cov__[[6]][[imonth]],xy)
  env_list[[i]] <- cbind(track_list[[i]],
                         BATHY_median_nn=l.bathy,
                         allSST_median_nn=l.sst,
                         allSSH_median_nn=l.ssh,
                         allCHLA_median_nn=l.chla)
}

##############################
##############################
#' Plot Environmental Maps
#' 
#' Melt list into dataframe for ggplot2 because it can't handle lists
tracks.melt = do.call(rbind,env_list)
coordinates(tracks.melt) <- ~long+lat
#' Map BATHY
tracks.melt$bathy2 <- abs(tracks.melt$BATHY_median_nn)
m <- visio(obj=subset(tracks.melt,!is.na(BATHY_median_nn)),
           leg.col=rainbow(8),field="bathy2",
           choropleth=FALSE,pt.col="black",legend.title = "Depth (x1,000m)")
m
mapshot(m,file="bathy_raw.png")


#' MAP SST
library(RColorBrewer)
m <- visio(obj=subset(tracks.melt,!is.na(allSST_median_nn)),
           leg.col=rev(brewer.pal(n=9,"Spectral")),
           field="allSST_median_nn",
           choropleth = FALSE,pt.col="black",legend.title="SST (DegC)")
mapshot(m,file="sst_raw.png")


#' MAP SSH
m <- visio(obj=subset(tracks.melt,!is.na(allSSH_median_nn)),
           leg.col="RdBu",field="allSSH_median_nn",
           choropleth=FALSE,pt.col="black",legend.title="SSH (meter)")
mapshot(m,file="ssh_raw.png")


#' MAP Chl-a
m <- visio(obj=subset(tracks.melt,!is.na(allCHLA_median_nn)),
           leg.col="BuGn",field="allCHLA_median_nn",
           choropleth = TRUE,bins = c(0,0.15,0.5,1,5,30),
           pt.col="black",legend.title="CHLA (mg/m3)")
mapshot(m,file="chla_raw.png")

cat("processing fishery observations.\n")

## load point process codes
source("cellFromXY_nbr.R")
source("extract_1.R")
source("ppm_st_glm_6f.R")
source("ppm_st_utils_2.R")

## study area raster
## missing = outside study area
## non-missing = inside study area
sp <- raster(cov__[["bathy"]],layer=1)

# generate pseudo-absence data 
pa1 <- ppm.st.glm(
  loc, window = sp, eps = 0.5,
  covariates = cov__,
  trend=~bathy+sst+sst2+fpi+ssh+chl,
  coordinate=~long+lat+amon)

## maximum likelihood estimation
ppm1 <- gam(z~bathy+s(sst,bs="cr",k=5)+fpi+ssh,data=pa1,
            weights=wt,
            family=quasi(link="log",variance="mu"))

## summary statistics
summary(ppm1)

## paritial effect plots
plot(ppm1)

## likelihood approximation at 0.5 resolution
ppm.st.lik(ppm1,na.omit(pa1))

# generate pseudo-absence data 
pa2 <- ppm.st.glm(
  loc, window = sp, eps = 1.0,
  covariates = cov__,
  trend=~bathy+sst+sst2+fpi+ssh+chl,
  coordinate=~long+lat+amon)

## maximum likelihood estimation
ppm2 <- gam(z~bathy+s(sst,bs="cr",k=5)+fpi+ssh,data=pa2,
            weights=wt,
            family=quasi(link="log",variance="mu"))

## likelihood approximation at 0.5 resolution
ppm.st.lik(ppm2,na.omit(pa2))

## Akaike Information Criterion
qaic(ppm1)

## prediction
cat("preparing prediction locations...\n")
newd <- ppm.st.glm.new(
  window=sp,
  newdata=cov__,
  coord=~long+lat+amon)

pred <- predict(ppm1,newd,type="response",se=T)
newd$fit <- pred$fit
newd$se <- pred$se.fit

## make raster objects for raw prediction
predMonth <- ppm.st.glm.fill(
  window = sp,
  newdata = newd,
  formula = fit~long+lat+amon,
  verbose = T,
)

## load scripts
library(RColorBrewer)
source("choropleth_raster_2.R")
pal <- brewer.pal(6,"RdPu")


## quantiles across months
k_ <- 1 #  per cell
brks <- k_*quantile(
  c(as.array(predMonth)),
  probs=seq(0,1,len=length(pal)+1),na.rm=T)
brks[c(1,length(brks))] <- range(brks)+c(-100,100)


png(file="monthPredictFigLegend.png",width = 3000,
    height = 3000,res=500)
legend <- choropleth.raster(
  predMonth[[1]],pal,expression(paste("Intensity per cell ",10^-3,"")),
  legend=TRUE,
  leg.label=c("1.7 to 6.0","6.0 to 7.4","7.4 to 8.6",
              "8.7 to 10.3","10.3 to 12.9","12.9 to 62.2"),
  style="fixed",
  fixedBreaks=brks)
dev.off()

r <- 3/1.369619/1.5
pdf(file="PredMonth.pdf",height=10,width=10*r)
op <- par(mar=c(0.1,0.1,2.1,0.1),oma=c(0,0,0,0),mfrow=c(2,5))
for(i in 1:dim(predMonth)[3]){
  cat(i,"\n")
  line <- choropleth.raster(
    k_*predMonth[[i]],pal,"test",legend=FALSE,
    style="fixed",fixedBreaks=brks)
  plot(line,col=pal,main=names(predMonth)[i],
       xlim=range(pa1[,"long"]),
       ylim=range(pa1[,"lat"]),
       asp=FALSE,axes=FALSE,box=T,
       legend=FALSE)
}
par(op)
dev.off()
