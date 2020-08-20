# https://pdfs.semanticscholar.org/0df2/9e127c60a59707415053acb1e70f5d0147ac.pdf
##%######################################################%##
#                                                          #
####          Watterson, I. G., Hirst, A. C.,           ####
####         & Rotstayn, L. D. (2013). A skill          ####
####               score based evaluation               ####
####          of simulated Australian climate.          ####
####           Australian Meteorological and            ####
####       Oceanographic Journal, 63(1), 181-190.       ####
#                                                          #
##%######################################################%##

m.stat <- function(mod, obs, ...) {
  # model field = x; observed field = y
  # mse the mean square error between x and y,
  # V and G are spatial variance and domain mean of the respective fields
  require("Metrics")
  obs <- na.omit(obs)
  mod <- na.omit(mod)
  stopifnot(length(obs) == length(mod))
  se = (obs - mod)^2
  mse = mean(se)
  ## variance and mean of modelled
  Vx <- var(mod)
  Gx <- mean(mod)
  ## variance and mean of observed
  Vy <- var(obs)
  Gy <- mean(obs)
  ## calc M-stat
  M <- (2 / pi) * asin(1 - mse / (Vx + Vy + ((Gx - Gy)^2)))
  ## multiply by 1000
  return(round(M*1000, 0))
}


## TESTS ##
# obs <- c(-37.3, -93.3, -45.9, 101.9, -7.9, 84.6, -29.4, 34.2, -129.2, 46.5, -67,
#          -177, -1, -28.8, -81.3, -172.2, 196.6, -3.6, -109.8, 17.6, -60.6, -57.1,
#          119.1, 151.5, 251.1, 212.1, 108.8, -105.4, -3.7, 152.5, -58.5, -254.3,
#          -9.5, -164.4, 271.8, 9, -73.6, 124.6, 13.8, -183.1, -138.7, -140, -199.7,
#          118.3, -20.7, 245.1, 261.5, -49.3, -13.5, 196.8)
# mod <- c(-2.6, -6.5, -3.2, 7.1, -0.5, 5.9, -2, 2.4, -8.9, 3.2, -4.6, -12.3, -0.1,
#          -2, -5.6, -11.9, 13.6, -0.2, -7.6, 1.2, -4.2, -4, 8.2, 10.5, 17.4, 14.7,
#          7.5, -7.3, -0.3, 10.6, -4.1, -17.6, -0.7, -11.4, 18.8, 0.6, -5.1, 8.6, 1,
#          -12.7, -9.6, -9.7, -13.8, 8.2, -1.4, 17, 18.1, -3.4, -0.9, 13.6)
# plot(mod, obs); cor.test(mod, obs) ## perfect correlation, order of mag diff in values
# m.stat(mod, obs) ## 99 = low score due to differences
# m.stat(mod*10, obs) ## 773 = high score, no order of mag difference
# m.stat(obs, obs) ## 1000 = perfect score, modelled = observed
# m.stat(obs, obs/2) ## ~ 594, modelled is half of observed
# clim <- getData("worldclim", res = 10, var = "tmax")
# clim <- resample(crop(clim, extent(112,153,-43,-10)),
#                  y = raster(res = 0.5, ext = extent(112,153,-43,-10)))
# plot(clim[[1]][]/10, clim[[2]][]/10)
# cor.test(clim[[1]][]/10, clim[[2]][]/10)
# m.stat(mod = clim[[2]][]/10, obs = clim[[1]][]/10) ## 859
# futClim <- getData("CMIP5", var="tmax", res=10, rcp = 45, model='AC', year = 50)
# futClim <- resample(crop(futClim, extent(112,153,-43,-10)),
#                     y = raster(res = 0.5, ext = extent(112,153,-43,-10)))
# plot(futClim[[1]][], clim[[1]][])
# spplot(stack(futClim[[1]], clim[[1]]), main = "Future Vs Obs")
# cor.test(futClim[[1]][], clim[[1]][]) ## 0.98
# m.stat(mod = futClim[[1]][], obs = clim[[1]][]) ## 798
