source("renv/activate.R")

suppressPackageStartupMessages(
   {library(data.table)
    library(stringr)
    library(pbapply)
    library(Just.universal)})

data.dir = "/data-coco/TROPOMI"
tropomi.urls.path = file.path(data.dir, "tropomi-urls.txt")

n.workers = 8

pm = function(...) pairmemo(
    directory = file.path(data.dir, "pairmemo"),
    n.frame = 2,
    ...)

download = function(from, to, ...)
    download.update.meta(from, file.path(data.dir, "downloads"), to, ...)

date.first = as.Date("2021-07-01")
  # The first day of the latest version of the TROPOMI
  # nitrogen-dioxide product.
date.last = as.Date("2021-10-23")
  # A somewhat arbitrary date close to when I wrote this code.
dates.all = seq(date.first, date.last, by = 1)

study.bbox = list(
  # Includes part of Canada, the continental US, Mexico, and much of
  # the Carribean (including Puerto Rico).
    lon.min = -125, lon.max = -65,
    lat.min = 14, lat.max = 50)
