climVal <- function(modelled, observed, region = NULL, partial_overlap = TRUE, ...) {
  require("asbio")
  require("MetricsWeighted")
  require("hydroGOF")
  stopifnot(exists("m.stat"))
  stopifnot(c(class(modelled), class(observed)) != "raster")
  stopifnot(inherits(region, c("Extent","sf", "SpatialPolygons", "NULL")))
  ## if region is supplied, mask the data to the correct region
  if (!is.null(region)) {
    rs = crop(stack(modelled, observed), region)
    ## if shapefile, mask to boundary
  }
  if (inherits(region, c("sf", "sp", "SpatialPolygons")) & partial_overlap) {
      region_raster <- rasterize(region, y = raster(res = res(rs),
                                                    ext = extent(rs)),
                                 getCover = TRUE)
      region_raster[region_raster == 0] <- NA
      rs = na.omit(as.data.frame(mask(rs, region_raster), xy = TRUE))
    } else if (inherits(region, c("sf", "sp")) & !partial_overlap) {
      rs = na.omit(as.data.frame(mask(rs, region), xy = TRUE))
    } else {
    rs = na.omit(as.data.frame(rs, xy = TRUE))
  }
  ## calculate weights based on latitudes
  rs$weights = cos(rs$y  * (pi/180))
  ## Calculate pb correlation
  pb = unlist(unname(round(asbio::r.pb(X = rs[, 3], Y = rs[, 4])[1],2)))
  ## M stat
  M = round(m.stat(mod = rs[, 3], obs = rs[, 4])/10,1)
  ## weighted RMSE
  rmse = round(MetricsWeighted::rmse(actual = rs[, 4], predicted = rs[, 3],
                               w = rs$weights),2)
  ## rSD
  rSD = round(sd(rs[, 3]) / sd(rs[, 4]), 2)
  ## md
  md = round(hydroGOF::md(sim = rs[, 3], obs = rs[, 4]) , 2)
  ## pbias
  pbias = round(hydroGOF::pbias(sim = rs[, 3], obs = rs[, 4]), 2)
  return(c("rpb" = pb, "M" = M, "wRMSE" = rmse, "rSD" = rSD, "md" = md,
           "pBias" = pbias))
}
