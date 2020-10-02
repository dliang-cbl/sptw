rm(list=ls())
library(raster)
library(mgcv)

## simulated fishery observations
loc <- read.csv("supp_csp_sim.csv")
## a list of seascape raster stacks
load(file="supp_csp_env.rda")
## each layer is one month
## bathy= bathymetry
## sst = sea surface temp.
## sst2 = quadratic sst
## fpi = frontal probability index
## ssh = sea surface height
## chl = chlorophyll-a index

cat("processing fishery observations.\n")

## load point process codes
source("cellFromXY_nbr.R")
source("ppm_st_glm_6e.R")
source("extract_1.R")
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
ppm1 <- gam(z~bathy+sst+sst2+fpi+ssh,data=pa1,
            weights=wt,
            family=quasi(link="log",variance="mu"))

## summary statistics
summary(ppm1)

## Akaike Information Criterion
qaic(ppm1)
