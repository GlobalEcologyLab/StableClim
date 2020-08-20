library(ncdf4)
library(raster)
library(tidyverse)
library(sf)
library(plotrix)
source("04_m_stat.R")
source("04_climate_validation.R")

temp_clim <- raster::stack(
  "cru_ts4.03.1901.2018.tmp.dat_climatology.nc",
  "ensMeanClimatology.nc") %>%
  as.data.frame(., xy = TRUE) %>%
  tbl_df() %>%
  dplyr::rename("CRU" = 3, "Ensemble" = 4) %>%
  na.omit()

prec_clim <- raster::stack(
  "ru_ts4.03.1901.2018.pre.dat_climatology.nc",
  "ensMeanClimatology_pr.nc") %>%
  as.data.frame(., xy = TRUE) %>%
  tbl_df() %>%
  dplyr::rename("CRU" = 3, "Ensemble" = 4) %>%
  na.omit()

taylor.diagram(ref = temp_clim$CRU,
               model = temp_clim$Ensemble, col = "darkgreen",
               normalize = TRUE)
taylor.diagram(ref = prec_clim$CRU,
               model = prec_clim$Ensemble, col = "darkblue",
               normalize = TRUE, add = TRUE)

m.stat(mod = temp_clim$Ensemble, obs = temp_clim$CRU)
m.stat(mod = prec_clim$Ensemble, obs = prec_clim$CRU)

temp_clim %>%
  ggscatterstats(., x = "CRU", y = "Ensemble", type = "r",
                 method = "lm",
                 marginal.type = "densigram",
                 xfill = "#0072B2",
                 yfill = "#009E73",
                 xlab = "Observed temperature (°C/Day)", ylab = "Modelled temperature (°C/Day)")

prec_clim %>%
  ggscatterstats(., x = "CRU", y = "Ensemble", type = "r",
                 method = "lm",
                 marginal.type = "densigram",
                 xfill = "#0072B2",
                 yfill = "#009E73",
                 xlab = "Observed precipitation (mm/month)", ylab = "Modelled precipitation (mm/month)")

# Global PBIAS for precip
pbias(x = prec_clim$Ensemble, y = prec_clim$CRU)

## Temperature
abs_diff_df <- temp_clim %>%
  mutate(abs_diff = (round(CRU,0)-round(Ensemble,0))) %>%
  dplyr::select(1,2,5) %>%
  rasterFromXYZ(crs = crs(raster()))
spplot(abs_diff_df, zlim = c(-max(abs(abs_diff_df[]), na.rm = TRUE),
                             max(abs(abs_diff_df[]), na.rm = TRUE)))

# Weighted RMSE
w <- na.omit(mask(init(abs_diff_df, "y"), abs_diff_df)[])
weights <- cos(w  * (pi/180))
MetricsWeighted::rmse(actual = na.omit(temp_clim$CRU),
                      predicted = na.omit(temp_clim$Ensemble),
                      w = weights)

## rSD
sd(temp_clim$Ensemble)/sd(temp_clim$CRU)
sd(prec_clim$Ensemble)/sd(prec_clim$CRU)

## md
hydroGOF::md(sim = temp_clim$Ensemble, obs = temp_clim$CRU)
hydroGOF::md(sim = prec_clim$Ensemble, obs = prec_clim$CRU)

## Correlation at global and lat-band scales
cru_ts_mean <- raster("cru_ts4.03.1901.2018.tmp.dat_climatology.nc")
hist_ts_mean <- mask(raster("ensMeanClimatology.nc"),
                     cru_ts_mean)
cru_pr_mean <- raster("cru_ts4.03.1901.2018.pre.dat_climatology.nc")
hist_pr_mean <- mask(raster("ensMeanClimatology_pr.nc"), cru_pr_mean)

# Global
climVal(modelled = hist_ts_mean, observed = cru_ts_mean, region = extent(c(-180,180,-90,90)))
climVal(modelled = hist_pr_mean, observed = cru_pr_mean, region = extent(c(-180,180,-90,90)))

# High northern lats
climVal(modelled = hist_ts_mean, observed = cru_ts_mean,
        region = extent(c(-180,180,50,90)),
        partial_overlap = FALSE)
climVal(modelled = hist_pr_mean, observed = cru_pr_mean,
        region = extent(c(-180,180,50,90)),
        partial_overlap = FALSE)

high_north <- crop(stack(cru_ts_mean, hist_ts_mean),
                   y = c(-180,180,50,90))
names(high_north) <- c("cru", "hist")
high_north_pr <- crop(stack(cru_pr_mean, hist_pr_mean),
                   y = c(-180,180,50,90))
names(high_north_pr) <- c("cru", "hist")
high_north <- high_north %>%
  as.data.frame(xy = TRUE) %>%
  tbl_df() %>%
  mutate(Var = "TS") %>%
  bind_rows(., as.data.frame(high_north_pr, xy = TRUE)) %>%
  mutate(Var = if_else(is.na(Var), "PR", "TS")) %>%
  na.omit() %>%
  select(-x)
high_north
high_north$Var <- factor(high_north$Var, levels = c("TS", "PR"))
p_north <- grouped_scatter(high_north, x = cru, y = hist, type = "robust",
                                  method = "lm",
                                  marginal.type = "densigram",
                                  nboot = 100,
                                  point.alpha = 0.7,
                                  point.size = 2,
                                  grouping.var = Var,
                                  axes.range.restrict = FALSE,
                                  xfill = "#0072B2",
                                  yfill = "#009E73",
                                  xlab = "Observed mean",
                                  ylab = "Modelled mean",
                                  caption.size = 3,
                                  title.size = 3)
p_north
# Mid northern lats
climVal(modelled = hist_ts_mean, observed = cru_ts_mean,
        region = extent(c(-180,180,20,50)),
        partial_overlap = FALSE)
climVal(modelled = hist_pr_mean, observed = cru_pr_mean,
        region = extent(c(-180,180,20,50)),
        partial_overlap = FALSE)
mid_north <- crop(stack(cru_ts_mean, hist_ts_mean),
                   y = c(-180,180,20,50))
names(mid_north) <- c("cru", "hist")
mid_north_pr <- crop(stack(cru_pr_mean, hist_pr_mean),
                      y = c(-180,180,20,50))
names(mid_north_pr) <- c("cru", "hist")
mid_north <- mid_north %>%
  as.data.frame(xy = TRUE) %>%
  tbl_df() %>%
  mutate(Var = "TS") %>%
  bind_rows(., as.data.frame(mid_north_pr, xy = TRUE)) %>%
  mutate(Var = if_else(is.na(Var), "PR", "TS")) %>%
  na.omit() %>%
  select(-x)
mid_north
mid_north$Var <- factor(mid_north$Var, levels = c("TS", "PR"))
p_midnorth <- grouped_scatter(mid_north, x = cru, y = hist, type = "robust",
                           method = "lm",
                           marginal.type = "densigram",
                           nboot = 100,
                           point.alpha = 0.7,
                           point.size = 2,
                           grouping.var = Var,
                           axes.range.restrict = FALSE,
                           xfill = "#0072B2",
                           yfill = "#009E73",
                           xlab = "Observed mean",
                           ylab = "Modelled mean",
                           caption.size = 3,
                           title.size = 3)
p_midnorth
# High tropics
climVal(modelled = hist_ts_mean, observed = cru_ts_mean,
        region = extent(c(-180,180,-20,20)),
        partial_overlap = FALSE)
climVal(modelled = hist_pr_mean, observed = cru_pr_mean,
        region = extent(c(-180,180,-20,20)),
        partial_overlap = FALSE)

high_trops <- crop(stack(cru_ts_mean, hist_ts_mean),
                   y = c(-180,180,-20,20))
names(high_trops) <- c("cru", "hist")
high_trops_pr <- crop(stack(cru_pr_mean, hist_pr_mean),
                      y = c(-180,180,-20,20))
names(high_trops_pr) <- c("cru", "hist")
high_trops <- high_trops %>%
  as.data.frame(xy = TRUE) %>%
  tbl_df() %>%
  mutate(Var = "TS") %>%
  bind_rows(., as.data.frame(high_trops_pr, xy = TRUE)) %>%
  mutate(Var = if_else(is.na(Var), "PR", "TS")) %>%
  na.omit() %>%
  select(-x)
high_trops
high_trops$Var <- factor(high_trops$Var, levels = c("TS", "PR"))
p_trops <- grouped_scatter(high_trops, x = cru, y = hist, type = "robust",
                              method = "lm",
                              marginal.type = "densigram",
                              nboot = 100,
                              point.alpha = 0.7,
                              point.size = 2,
                              grouping.var = Var,
                              axes.range.restrict = FALSE,
                              xfill = "#0072B2",
                              yfill = "#009E73",
                              xlab = "Observed mean",
                              ylab = "Modelled mean",
                              caption.size = 3,
                              title.size = 3)
p_trops

# Mid southern lats
climVal(modelled = hist_ts_mean, observed = cru_ts_mean,
        region = extent(c(-180,180,-50,-20)),
        partial_overlap = FALSE)
climVal(modelled = hist_pr_mean, observed = cru_pr_mean,
        region = extent(c(-180,180,-50,-20)),
        partial_overlap = FALSE)

mid_south <- crop(stack(cru_ts_mean, hist_ts_mean),
                   y = c(-180,180,-50,-20))
names(mid_south) <- c("cru", "hist")
mid_south_pr <- crop(stack(cru_pr_mean, hist_pr_mean),
                      y = c(-180,180,-50,-20))
names(mid_south_pr) <- c("cru", "hist")
mid_south <- mid_south %>%
  as.data.frame(xy = TRUE) %>%
  tbl_df() %>%
  mutate(Var = "TS") %>%
  bind_rows(., as.data.frame(mid_south_pr, xy = TRUE)) %>%
  mutate(Var = if_else(is.na(Var), "PR", "TS")) %>%
  na.omit() %>%
  select(-1)
mid_south
mid_south$Var <- factor(mid_south$Var, levels = c("TS", "PR"))
p_midsouth <- grouped_scatter(mid_south, x = cru, y = hist, type = "robust",
                           method = "lm",
                           marginal.type = "densigram",
                           nboot = 100,
                           point.alpha = 0.7,
                           point.size = 2,
                           grouping.var = Var,
                           axes.range.restrict = FALSE,
                           xfill = "#0072B2",
                           yfill = "#009E73",
                           xlab = "Observed mean",
                           ylab = "Modelled mean",
                           caption.size = 3,
                           title.size = 3)
p_midsouth

# Scatter plots
plot_list <- list(p_north, p_midnorth, p_trops, p_midsouth)
names(plot_list) <- c("p_north", "p_midnorth", "p_trops", "p_midsouth")
plot_list

lapply(seq_along(plot_list), function(x) {
  ggsave(filename = paste0(names(plot_list)[x], "_50yrClim.pdf"), plot = plot_list[[x]], device = "pdf",
         path = "",
         width = 8, height = 8, units = "in",
         dpi = 500)
})

## Import the IPCC regions
ipcc_reg <- st_read("WG1AR5_Annex1_DissolveRegions.shp")
ipcc_reg
unique(ipcc_reg$RegionName)

reg <- c("High Latitudes", "Mediterranean and Sahara",
         "North America (East)", "Southern Africa and West Indian Ocean",
         "Australia and New Zealand")

ipcc_reg <- ipcc_reg[ipcc_reg$RegionName %in% reg, ]

for (name in ipcc_reg$RegionName) {
  poly <- ipcc_reg[ipcc_reg$RegionName == name, ]
  print(name)
  print(climVal(modelled = hist_ts_mean, observed = cru_ts_mean,
          region = poly, partial_overlap = TRUE))
}

for (name in ipcc_reg$RegionName) {
  poly <- ipcc_reg[ipcc_reg$RegionName == name, ]
  print(name)
  print(climVal(modelled = hist_pr_mean, observed = cru_pr_mean,
                region = poly, partial_overlap = TRUE))
}

## Import the Wallace Realms
wall_realms <- st_read("WallaceRegions.shp")
wall_realms$Descriptio <- NULL
wall_realms <- st_zm(wall_realms)
unique(wall_realms$HoltRealm)

plot(wall_realms[2])

reg <- c("Neotropical", "Oriental", "Palearctic")

plot(cru_ts_mean, col = bpy.colors(100))
lines(as_Spatial(wall_realms[wall_realms$HoltRealm %in% reg, ]))

wall_realms <- wall_realms[wall_realms$HoltRealm %in% reg, ]

for (name in unique(wall_realms$HoltRealm)) {
  poly <- wall_realms[wall_realms$HoltRealm == name, ]
  poly <- st_union(poly)
  poly <- as_Spatial(poly)
  print(name)
  print(climVal(modelled = hist_ts_mean, observed = cru_ts_mean,
                region = poly, partial_overlap = TRUE))
}

for (name in unique(wall_realms$HoltRealm)) {
  poly <- wall_realms[wall_realms$HoltRealm == name, ]
  poly <- st_union(poly)
  poly <- as_Spatial(poly)
  print(name)
  print(climVal(modelled = hist_pr_mean, observed = cru_pr_mean,
                region = poly, partial_overlap = TRUE))
}
