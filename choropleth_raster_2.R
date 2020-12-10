require(classInt)
## ver 2) allow label to be added to the legend
## reclassify a raster using classInt
## and plot the legend
choropleth.raster <- function(ext,col,leg.title,legend=FALSE,leg.label=NULL,...)
{
  ## ext: raster object
  ## col: color palette
  ## leg.title: title of legend
  ## legend : whether plot legend
  ## ...: additional arguments to classIntervals
  ## value: a legend plot
  ## and return the re-classified raster
  val <- na.omit(values(ext))
  q6 <- classIntervals(var=val,n=length(col),...)
  leg.quant <- findColours(q6,col)
  ## plot legend
  if(legend){
    if(is.null(leg.label)){
      leg.label <- names(attr(leg.quant,"table"))
    }
    oldpar <- par(mar=c(0.1,0.1,0.1,0.1))
    plot(c(0,1),c(0,1),axes=FALSE,type="n",xlab="",ylab="")
    legend("topleft",legend=leg.label,
           fill=attr(leg.quant,"palette"),bty="n",
           title=leg.title)
    par(oldpar)
  }
  ## re-classify raster
  val2 <- (match(leg.quant,attr(leg.quant,"palette")))
  ii2 <-  attr(val,"na.action")
  if(is.null(ii2)){
    values(ext) <- val2
  }
  else{
    r <- values(ext)
    r[-ii2] <- val2
    values(ext) <- r
  }
  if(FALSE){
    plot(r,col=col)
    plot(ext,col=col,legend=FALSE)
  }
  ext
}

test <- function()
{
  source("choropleth_raster_2.R")
  require(RColorBrewer)
  r <- raster(ncols=36, nrows=18)
  r[] <- runif(ncell(r))
  r[100] <- NA
  debug(choropleth.raster)
  pal <- brewer.pal(6,"Blues")
  r2 <- choropleth.raster(r,pal,"test",style="quantile",dataPrecision=3)
  plot(r2,col=pal,legend=FALSE)
}

test2 <- function()
{
  source("choropleth_raster_1.R")
  require(RColorBrewer)
  r <- raster(ncols=36, nrows=18)
  r[] <- runif(ncell(r))
  r[100] <- NA
  debug(choropleth.raster)
  pal <- brewer.pal(6,"Blues")
  r2 <- choropleth.raster(r,pal,"test",legend=TRUE,
                          leg.label=c("0.00-0.17","0.17-0.33","0.34-0.50",
                                         "0.51-0.67","0.68-0.83","0.83-1.00"),
                          style="fixed",fixedBreaks=seq(0,1,len=7),
                          dataPrecision=3)
  plot(r2,col=pal,legend=FALSE,zlim=NULL)
}