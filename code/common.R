if (!exists("renv.activated"))
   {source("renv/activate.R")
    renv.activated = T}

suppressPackageStartupMessages(
   {library(data.table)
    library(stringr)
    library(pbapply)
    library(Just.universal)})

data.dir = "/data-coco/TROPOMI"
tropomi.dir = "/data-coco/TROPOMI/tropomi-repro/extracted"

n.workers = 8

pm = function(...)
   {loadNamespace("jsonlite")
    loadNamespace("digest")
    loadNamespace("fst")
    pairmemo(
        directory = file.path(data.dir, "pairmemo"),
        n.frame = 2,
        ...)}

download = function(from, to, ...)
   {loadNamespace("DBI")
    loadNamespace("RSQLite")
    loadNamespace("digest")
    download.update.meta(from, file.path(data.dir, "downloads"), to, ...)}

date.first = as.Date("2018-04-30")
  # Around the first day of availability of the latest TROPOMI
  # nitrogen-dioxide product.
date.last = as.Date("2021-12-31")
dates.all = seq(date.first, date.last, by = 1)

study.bbox = list(
  # This box is based on the bounding box of a TEMPO shapefile of
  # mysterious origin.
    lon.min = -164, lon.max = -20,
    lat.min = 17, lat.max = 64)
