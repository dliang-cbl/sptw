#' Helper function to visualize sp object with legends
#'
#' \code{visio} uses leaflet to visualize sp or raster objects in a single command.
#'  More details can be found at \link{https://rstudio.github.io/leaflet/}
#' @param obj a spatial object to visualize
#' @param field the field name to visualize
#' @param leg.col color palette
#' @param choropleth whether showing choropleth or numeric color
#' @param bins for choropleth the binary color dividers
#' @param nquantile.color default division for binary color definitions
#' @param na.color The color to return for NA values.
#' @param pt.cex for point data the size of the point
#' @param pt.col for point data the color of the outline
#' @param pt.weight for point data the weight/width of the outline
#' @param Opacity opacity of filled color
#' @param project whether to project raster to EPSG:3857
#' @param legend whether to include legend
#' @param legend.title title of the legend
#' @param bbox west,south,east,north coordinates to bound,make sure only used once
#' @param Interval the spacing in map units between horizontal and vertical lines
#' @param add add to an existing layer 
#' @param ... additional arguments to graphic layer code
#' @return a leaflet widget
#' @export
visio <- function(
  obj,leg.col="white",field=NULL,
  choropleth=TRUE,bins=NA,nquantile.color=4,na.color="#808080",
  pt.cex=4,pt.col="black",pt.weight=1,Opacity=0.9,
  project=TRUE,
  legend=TRUE,legend.position="bottomright",legend.title="",
  bbox=NULL,Interval=NULL,add=NULL,provider=NULL,
  ...)
{
  ## extract response data fields
  if(class(obj)=="RasterLayer"){
    resp__ <- (as.vector(obj))
  }
  else if(class(obj) %in% c("SpatialPointsDataFrame","SpatialPolygonsDataFrame",
                            "SpatialLinesDataFrame")){
    ## extract response data fields
    if(is.null(field)){
      resp__ <- seq(1,nrow(obj@data))
      #stop("field can not be null for spatial objects.")
    }else{
      resp__ <- obj@data[,field]
    }
  } 
  else{
    stop("can't handle this type of spatial object")
  }
  
  ## define filled color palette
  if((class(obj)!="RasterLayer") && (is.null(field))){
    ## vector without field
    pal <- visio.color.helper(
      x=resp__,leg.col = leg.col[1],choropleth = FALSE,
      bins=bins,nquantile.color = nquantile.color,na.color = na.color)
  }else{
    ## raster or vector with field
    pal <- visio.color.helper(
      x=resp__,leg.col = leg.col,choropleth = choropleth,
      bins=bins,nquantile.color = nquantile.color,na.color = na.color)
  }

  ## define base map
  if(is.null(add)){
    if(class(obj)=="RasterLayer"){
      if(!is.null(provider)){
        base <- leaflet() %>% addProviderTiles(provider = provider)
      }else{
        base <- leaflet() %>% addProviderTiles(provider = providers$OpenStreetMap)  ## add base tiles
      }
      
    }
    else if(class(obj) %in% c("SpatialPointsDataFrame","SpatialPolygonsDataFrame",
                              "SpatialLinesDataFrame")){
      if(!is.null(provider)){
        base <- leaflet(obj) %>% addProviderTiles(provider = provider)  ## add base tiles
      }else{
        base <- leaflet(obj) %>% addProviderTiles(provider = providers$OpenStreetMap)  ## add base tiles
      }
      
    }
    map1 <- base %>% addScaleBar()
    if(!is.null(Interval)){
      map1 <- map1 %>% addGraticule(interval=Interval)
    }
    
  }else{
    if(class(add)[1] == "leaflet" & class(add)[2]=="htmlwidget"){
      map1 <- add
    }else{
      stop("base map is not leaflet.\n")
    }
  }
  
  ## add features
  if(class(obj)=="RasterLayer"){
    map.out <- map1 %>% addRasterImage(
      x=obj, colors=pal, opacity=Opacity,project=project,...
    )
  }
  else if(class(obj) == "SpatialPointsDataFrame"){
    ## spatial point data
    radius__ <- pt.cex  ## size of points
    if(is.na(pt.col)){
      ## outline color transparent
      map.out <- map1%>%addCircleMarkers(
        radius=radius__,
        stroke=FALSE,
        fillColor=~pal(resp__),
        fillOpacity=Opacity,
        data=obj,...
      )
    }else{
      map.out <- map1%>%addCircleMarkers(
        radius=radius__,
        stroke=TRUE,
        color=pt.col,
        weight=pt.weight,
        opacity=1,
        fillColor=~pal(resp__),
        fillOpacity=Opacity,
        data=obj,...
      )
    }
  }
  else if(class(obj) == "SpatialPolygonsDataFrame"){
    if(is.na(pt.col)){
      ## outline color transparent
      map.out <- map1 %>% addPolygons(
        stroke=FALSE, 
        fillColor = ~pal(resp__),
        fillOpacity=Opacity,
		data=obj,...
      )
    }else{
      ## outline color provided
      map.out <- map1 %>% addPolygons(
        stroke=TRUE,
        color=pt.col,
        weight=pt.weight,
        opacity=1,
        fillColor = ~pal(resp__),
        fillOpacity = Opacity,
		data=obj,...
      )
    }
  }
  else if(class(obj) =="SpatialLinesDataFrame"){
    map.out <- map1 %>% addPolylines(
      stroke = TRUE, color = ~pal(resp__),
      weight=pt.weight,opacity = 1,fillOpacity = Opacity,
      data=obj,...
    )
  }

  ## add legend
  if(legend){
    map.out <- map.out %>%
      addLegend(position=legend.position,
                pal=pal,values=resp__,
                opacity=Opacity,
                title=legend.title)
  }
 
  if(!is.null(bbox)){
    map.out <- map.out %>% fitBounds(
      bbox[1],bbox[2],bbox[3],bbox[4]
    ) %>% addRectangles(
      lng1=bbox[1],lat1=bbox[2],lng2=bbox[3],lat2=bbox[4],
      weight=1,fill = FALSE
    )
  }
  map.out
}

## helper function to define color in leaflet
## return a color pallette
visio.color.helper <- function(x,leg.col,choropleth,bins,nquantile.color,na.color){
  if(class(x)=="factor"){
    ## categorical data
    pal <- colorFactor(
      palette = leg.col,
      domain = x,
      na.color=na.color
    )
  }else{
    ## quantitative data
    if(choropleth){
      ## discrete color
      if(any(is.na(bins))){
        ## contain missing values
        if(length(bins)>1){
          ## assume correct
          warning("missing value ignored in bins argument.\n")
          bins <- na.omit(bins)
        }else{
          ## assume using quantiles
          bins <- quantile(x,probs=seq(0,1,length.out = nquantile.color + 1),na.rm=TRUE)
		  bins <- unique(bins)
        }
      }
      pal <- colorBin(
        palette = leg.col,
        domain = x,
        bins = bins,
        na.color = na.color
      )
    }else{
      ## continuous color
      pal <- colorNumeric(
        palette=leg.col,
        domain=x,
        na.color = na.color)
    }
  }
  pal
}

show.knitr <- function(m,file=NULL)
{
  if(is.null(file)){
    file <- tempfile(tmpdir=".",fileext=".png")
  }
  mapshot(m,file=file)
  knitr::include_graphics(substr(file,3,nchar(file)))
}
