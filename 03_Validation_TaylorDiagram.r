library(ncdf4)
library(raster)
library(plotrix)

cru_ts <- raster("cru_ts4.03.1901.2018.tmp.dat_climatology.nc")
cru_pr <- raster("cru_ts4.03.1901.2018.pre.dat_climatology.nc")

cmipModels_ts <- lapply(list.files("C:/tmp/cdo_outputs/climatology/",
                         "_ts_", full.names = TRUE), FUN = function(x) {
                           f <- mask(raster(x), mask = cru_ts)
                         })
names(cmipModels_ts) <- paste0("X", 1:length(cmipModels_ts))
cmipModels_pr <- lapply(list.files("C:/tmp/cdo_outputs/climatology/",
                                   "_pr_", full.names = TRUE), FUN = function(x) {
                                     f <- mask(raster(x), mask = cru_pr)
                                   })
names(cmipModels_pr) <- paste0("X", 1:length(cmipModels_pr))

# Greens temp
greens <- paletteer_c("pals::ocean.algae", 25)[-c(1, 21:25)]

# blue prec
blues <- paletteer_c("pals::ocean.ice", 25)[-c(1, 21:25)]

taylor.diagram(ref=rnorm(30,sd=3),
               model = rnorm(30,sd=3)+rnorm(30)/2,
               ref.sd = TRUE,
               main = NULL,
               normalize = TRUE, bg = NA, col = NA, pch = 19)
lapply(names(cmipModels_ts), FUN = function(x) {
  print(x)
  i <- as.integer(gsub("X", "", x))
  f <- cmipModels_ts[[x]]
  taylor.diagram(values(cru_ts), values(f), normalize = TRUE, add = TRUE,
                 col = greens[i], pch = 19, pcex = 1.5, bg = "black")
})
lapply(names(cmipModels_pr), FUN = function(x) {
  print(x)
  i <- as.integer(gsub("X", "", x))
  f <- cmipModels_pr[[x]]
  taylor.diagram(values(cru_pr), values(f), normalize = TRUE, add = TRUE,
                 col = blues[i], pch = 19, pcex = 1.5, bg = "black")
})

temp <- calc(stack(cmipModels_ts, quick = TRUE), mean, na.rm = TRUE)
prec <- calc(stack(cmipModels_pr, quick = TRUE), mean, na.rm = TRUE)

taylor.diagram(values(cru_ts), values(temp), normalize = TRUE, add = TRUE,
               col = greens[1], pch = 17, pcex = 2, bg = "black")
taylor.diagram(values(cru_pr), values(prec), normalize = TRUE, add = TRUE,
               col = blues[1], pch = 17, pcex = 2, bg = "black")

model = cmipModels_ts[1]
f <- raster(model)
f <- mask(f, cru_ts)
model <- values(f)
ref <- values(cru_ts)
taylor.diagram(ref, model, normalize = TRUE)

for (model in cmipModels_ts[-1]) {
  t <- 0
  print(model)
  f <- raster(model)
  f <- mask(f, cru_ts)
  model <- values(f)
  ref <- values(cru_ts)
  oldpar <- taylor.diagram(ref, model, )
  names(df)[3] <- "Value"
  df$Model <- as.factor(model)
  df$Var <- as.factor("Temperature")
  temp_df <- bind_rows(temp_df, df)
  rm(f)
}
f <- raster(ts_files[7], varname = "tmp")
df <- na.omit(as.data.frame(f/10, xy = TRUE))
names(df)[3] <- "Value"
df$Model <- as.factor("CRU")
df$Var <- as.factor("Temperature")
temp_df <- bind_rows(temp_df, df)
