library(raster)
library(tidyverse)
library(ggstatsplot)
library(sf)
library(readxl)

#### ANTARCTICA ####
ant <- st_read("NaturalEarth/110m_cultural/ne_110m_admin_0_countries.shp") %>%
  filter(NAME == "Antarctica")

core_data <- read_xlsx("Shakun_2012_Nature_SI2.xlsx", sheet = 83)
core_data
core_data <- core_data[-1,]
names(core_data) <- c("Temp", "Age", "DROP")
core_data$DROP <- NULL
core_data$Age <- as.integer(core_data$Age)*-1
core_data$Temp <- as.numeric(core_data$Temp)
core_data

# Trace temperature data
trace_temp <- brick("TraCE21ka/paleo_TS_Annual.nc")
trace_temp <- crop(trace_temp, ant)
trace_temp <- mask(trace_temp, ant)
trace_temp <- setZ(trace_temp, as.integer(gsub("X\\.","", names(trace_temp)))*-1, "YearsBP")
core_idx <- which(core_data$Age %in% getZ(trace_temp))
trace_idx <- which(getZ(trace_temp) %in% core_data$Age)

## Subset the trace and core data to the same timesteps
core_data <- core_data[core_idx, ]
trace_temp <- trace_temp[[trace_idx]]
core_data <- core_data[order(core_data$Age), ]

## raster with latitude cell values
w <- init(raster(res = 2.5), "y")
## cosine after transforming to radians
w <- cos(w  * (pi/180))
w <- mask(crop(w, ant), ant)
## multiply values by weights
x <- trace_temp * w

## time series of weighted temperature values
trace_ts <- data.frame(Age = core_data$Age, Temp = cellStats(x, sum) / cellStats(w, sum))

plot_temps <- merge(core_data, trace_ts, by = "Age")
names(plot_temps) <- c("x", "y", "y1")
plot(x = plot_temps$x, y = plot_temps$y,
     ylim = c(range(c(plot_temps$y, plot_temps$y1))),
     type = "l", col = "darkblue")
lines(x = plot_temps$x, y = plot_temps$y1, col = "red")

summary(trace_ts$Temp)
summary(core_data$Temp)

## https://stackoverflow.com/questions/46860333/rolling-regression-on-irregular-time-series
ws <- 150
library(data.table)

snr_core <- setDT(plot_temps)[.(start = x, end = x + ws), on = .(x >= start, x < end),
              c(Slope = as.list(round(coef(lm(y ~ x.x))[2], 5)), Var = round(sd(residuals(lm(y ~ x.x))), 5)), by = .EACHI]
names(snr_core) <- c("Start", "End", "SlopeCore", "VarCore")
snr_core$SNRCore <- round(abs(snr_core$Slope)/snr_core$Var, 5)
snr_core

snr_trace <- setDT(plot_temps)[.(start = x, end = x + ws), on = .(x >= start, x < end),
                              c(Slope = as.list(round(coef(lm(y1 ~ x.x))[2], 5)), Var = round(sd(residuals(lm(y1 ~ x.x))), 5)),
                              by = .EACHI]
names(snr_trace) <- c("Start", "End", "SlopeTrace", "VarTrace")
snr_trace$SNRTrace <- round(abs(snr_trace$Slope)/snr_trace$Var,5)
snr_trace

snr_plot <- merge(snr_core, snr_trace, by = "Start")
snr_plot
snr_plot <- snr_plot[is.finite(rowSums(snr_plot)), ]

summary(snr_plot)

plot(x = snr_plot$Start, y = snr_plot$SNRCore,
     ylim = c(range(c(snr_plot$SNRCore, snr_plot$SNRTrace), na.rm = TRUE, finite = TRUE)),
     type = "l", col = "darkblue")
lines(x = snr_plot$Start, y = snr_plot$SNRTrace, col = "red")

fwrite(snr_plot, "snr_ant_compar.csv")

snr_plot <- gather(snr_plot, key = Layer, value = Value, SNRTrace, SNRCore)
as_tibble(snr_plot); summary(snr_plot)
## Cutting snr_plot$Start into snr_plot$Start_rec
snr_plot$Start_rec <- cut(snr_plot$Start, include.lowest=TRUE,  right=TRUE,
                          breaks=c(-21000, -15000, -11000, -3000, 0))
## Recoding snr_plot$Start_rec into snr_plot$Start_rec_rec
snr_plot$Start_rec <- fct_recode(snr_plot$Start_rec,
               "21-15k BP" = "[-2.1e+04,-1.5e+04]",
               "15-11k BP" = "(-1.5e+04,-1.1e+04]",
               "11-3k BP" = "(-1.1e+04,-3e+03]",
               ">3k BP" = "(-3e+03,0]")
## Reordering snr_plot$Start_rec_rec
snr_plot$Start_rec <- factor(snr_plot$Start_rec, levels=c("21-15k BP", "15-11k BP", "11-3k BP", ">3k BP"))

## Recoding snr_plot$Layer into snr_plot$Layer_rec
snr_plot$Layer_rec <- fct_recode(snr_plot$Layer,
               "TraCE-21" = "SNRTrace",
               "Observed" = "SNRCore")

as_tibble(snr_plot); summary(snr_plot)
names(snr_plot)

snr_plot_df <- snr_plot[, c("Start_rec", "Layer_rec", "Value")]
names(snr_plot_df) <- c("Period", "Layer", "SNR")

fwrite(snr_plot_df, "snr_comparison_antarctica.csv")

ggplot(data = snr_plot_df) +
  facet_wrap(~Period, ncol = 2, scales = "free") +
  geom_violin(aes(x = Layer, y = SNR, group = Layer)) +
  geom_boxplot(aes(x = Layer, y = SNR, group = Layer), fill = NA) +
  geom_point(aes(x = Layer, y = SNR, group = Layer, colour = Layer), alpha = 0.7,
             position = "jitter")

scaleFUN <- function(x) sprintf("%.2f", x)

antarctica_plot <-  ggstatsplot::grouped_ggbetweenstats(
  data = snr_plot_df,
  x = Layer,
  y = SNR,
  grouping.var = Period,
  p.adjust.method = "none",
  conf.level = 0.95,
  var.equal = FALSE,
  mean.plotting = FALSE,
  ggplot.component = list(
    ggplot2::scale_y_continuous(
      labels = scaleFUN
    ),
    ggplot2::scale_colour_manual(
      values = c("#0072B2", "#009E73"),
      aesthetics = c("colour", "fill")
    )),
  type = "p",
  results.subtitle = FALSE,
  xlab = "")
antarctica_plot

# Run PERMANOVA and PERMDISP in PRIMER 6.0.
# See SNR_Comparison_StableClim.pwk inside 03_Validation_PRIMER_SNR_Comparison_StableClim.tar.gz
# P-values from analysis below.

pval <- c(0.916, 0.001, 0.003, 0.368) ## PERMANOVA

pval_permdisp <- c(0.921, 0.014, 0.01, 0.685, 0.039, 0.079, 0.001, 0.003, 0.001,
                   0.91, 0.001, 0.025, 0.001, 0.753, 0.017, 0.001, 0.197, 0.864,
                   0.002, 0.001, 0.203, 0.811, 0.044, 0.093, 0.001, 0.001,
                   0.001, 0.198) ## PERMDISP

p.adjust(pval_permdisp, "BY")
p.adjust(pval, "BY")

#### GREENLAND ####
gland <- st_read("NaturalEarth/110m_cultural/ne_110m_admin_0_countries.shp") %>%
  filter(NAME == "Greenland")

core_data <- read_xlsx("Shakun_2012_Nature_SI2.xlsx", sheet = 4)
core_data
core_data <- core_data[-1,-c(1,2,5)]
names(core_data) <- c("Temp", "Age")
core_data$Age <- as.integer(core_data$Age)*-1
core_data$Temp <- as.numeric(core_data$Temp)
core_data

# Trace temperature data
trace_temp <- brick("TraCE21ka/paleo_TS_Annual.nc")
trace_temp <- crop(trace_temp, gland)
trace_temp <- raster::mask(trace_temp, gland)
trace_temp <- setZ(trace_temp, as.integer(gsub("X\\.","", names(trace_temp)))*-1,
                   "YearsBP")
core_idx <- which(core_data$Age %in% getZ(trace_temp))
trace_idx <- which(getZ(trace_temp) %in% core_data$Age)

## Subset the trace and core data to the same timesteps
core_data <- core_data[core_idx, ]
trace_temp <- trace_temp[[trace_idx]]

core_data <- core_data[order(core_data$Age), ]

## raster with latitude cell values
w <- init(raster(res = 2.5), "y")
## cosine after transforming to radians
w <- cos(w  * (pi/180))
w <- raster::mask(raster::crop(w, gland), gland)
## multiply values by weights
x <- trace_temp * w

## time series of weighted temperature values
trace_ts <- tibble(Age = core_data$Age, Temp = cellStats(x, sum) / cellStats(w, sum))

plot_temps <- as_tibble(merge(core_data, trace_ts, by = "Age"))
names(plot_temps) <- c("x", "y", "y1")
plot(x = plot_temps$x, y = plot_temps$y,
     ylim = c(range(c(plot_temps$y, plot_temps$y1))),
     type = "l", col = "darkblue")
lines(x = plot_temps$x, y = plot_temps$y1, col = "red")

summary(trace_ts$Temp)
summary(core_data$Temp)

## https://stackoverflow.com/questions/46860333/rolling-regression-on-irregular-time-series
ws <- 100
library(data.table)

snr_core <- setDT(plot_temps)[.(start = x, end = x + ws), on = .(x >= start, x < end),
                              c(Slope = as.list(round(coef(lm(y ~ x.x))[2], 5)), Var = round(sd(residuals(lm(y ~ x.x))), 5)),
                              by = .EACHI]
names(snr_core) <- c("Start", "End", "SlopeCore", "VarCore")
snr_core$SNRCore <- round(abs(snr_core$Slope)/snr_core$Var, 5)
snr_core

snr_trace <- setDT(plot_temps)[.(start = x, end = x + ws), on = .(x >= start, x < end),
                               c(Slope = as.list(round(coef(lm(y1 ~ x.x))[2], 5)), Var = round(sd(residuals(lm(y1 ~ x.x))), 5)),
                               by = .EACHI]
names(snr_trace) <- c("Start", "End", "SlopeTrace", "VarTrace")
snr_trace$SNRTrace <- round(abs(snr_trace$Slope)/snr_trace$Var,5)
snr_trace

snr_plot <- merge(snr_core, snr_trace, by = "Start")
snr_plot
snr_plot <- snr_plot[is.finite(rowSums(snr_plot)), ]

summary(snr_plot)

plot(x = snr_plot$Start, y = snr_plot$SNRCore,
     ylim = c(range(c(snr_plot$SNRCore, snr_plot$SNRTrace), na.rm = TRUE, finite = TRUE)),
     type = "l", col = "darkblue")
lines(x = snr_plot$Start, y = snr_plot$SNRTrace, col = "red")

fwrite(snr_plot, "snr_gland_compar.csv")

snr_plot <- gather(snr_plot, key = Layer, value = Value, SNRTrace, SNRCore)
as_tibble(snr_plot); summary(snr_plot)
## Cutting snr_plot$Start into snr_plot$Start_rec
snr_plot$Start_rec <- cut(snr_plot$Start, include.lowest = TRUE,  right = TRUE,
                          breaks = c(-21000, -15000, -11000, -3000, 0))

## Recoding snr_plot$Start_rec
snr_plot$Start_rec <- fct_recode(snr_plot$Start_rec,
               "21-15k BP" = "[-2.1e+04,-1.5e+04]",
               "15-11k BP" = "(-1.5e+04,-1.1e+04]",
               "11-3k BP" = "(-1.1e+04,-3e+03]",
               ">3k BP" = "(-3e+03,0]")

## Reordering snr_plot$Start_rec_rec
snr_plot$Start_rec <- factor(snr_plot$Start_rec, levels=c("21-15k BP", "15-11k BP", "11-3k BP", ">3k BP"))

## Recoding snr_plot$Layer into snr_plot$Layer_rec
snr_plot$Layer_rec <- fct_recode(snr_plot$Layer,
                                 "TraCE-21" = "SNRTrace",
                                 "NGRIP" = "SNRCore")

as_tibble(snr_plot)
names(snr_plot)

snr_plot_df_greenland <- snr_plot[, c("Start_rec", "Layer_rec", "Value")]
names(snr_plot_df_greenland) <- c("Period", "Layer", "SNR")

fwrite(snr_plot_df_greenland, "snr_comparison_greenland.csv")

# Run PERMANOVA and PERMDISP in PRIMER 6.0.
# See SNR_Comparison_StableClim.pwk inside 03_Validation_PRIMER_SNR_Comparison_StableClim.tar.gz
# P-values from analysis below.
pval <- c(0.227, 0.774,0.04,0.107) ## PERMANOVA

pval_permdisp <- c(0.07, 0.044, 0.231, 0.369, 0.127, 0.569, 0.704, 0.001, 0.002,
                   0.002, 0.725, 0.001, 0.213, 0.956, 0.251, 0.001, 0.045, 0.014,
                   0.452, 0.004, 0.181, 0.076, 0.003, 0.501, 0.115, 0.002, 0.321,
                   0.273)  ## PERMDISP

p.adjust(pval_permdisp, "BY")
p.adjust(pval, "BY")

ggplot(data = snr_plot_df_greenland) +
  facet_wrap(~Period, ncol = 2) +
  geom_violin(aes(x = Layer, y = SNR, group = Layer)) +
  geom_boxplot(aes(x = Layer, y = SNR, group = Layer), fill = NA) +
  geom_point(aes(x = Layer, y = SNR, group = Layer, colour = Layer), alpha = 0.7,
             position = "jitter")

greenland_plot <- ggstatsplot::grouped_ggbetweenstats(
  data = snr_plot_df_greenland,
  x = Layer,
  y = SNR,
  grouping.var = Period,
  p.adjust.method = "none",
  conf.level = 0.95,
  var.equal = FALSE,
  mean.plotting = FALSE,
  ggplot.component = list(
    ggplot2::scale_y_continuous(
    limits = c(0, 0.4),
    breaks = seq(from = 0, to = .4, by = 0.1),
    labels = scaleFUN
  ),
  ggplot2::scale_colour_manual(
    values = c("#0072B2", "#009E73"),
    aesthetics = c("colour", "fill")
  )),
  type = "p",
  results.subtitle = FALSE,
  xlab = "")

greenland_plot

ggsave(filename = "antarctica_boxplot_SNR_restrictAxes.pdf", plot = antarctica_plot,
       device = "pdf", path = "",
       width = 200, height = 200, units = "mm", dpi = 300)

ggsave(filename = "greenland_boxplot_SNR_restrictAxes.pdf", plot = greenland_plot,
       device = "pdf", path = "",
       width = 200, height = 200, units = "mm", dpi = 300)


combined_SNR <- combine_plots(antarctica_plot, greenland_plot, ncol = 1)

ggsave(filename = "combined_SNR_restrictAxes.pdf", plot = combined_SNR,
       device = "pdf", path = "",
       width = 200, height = 400, units = "mm", dpi = 300)
