source("code/common.R")

suppressPackageStartupMessages(
   {library(httr)
    library(ncdf4)})

## * Helpers

mol.m2.to.Pmolecule.cm2 = 6.0221408e4
  # For comparison with Verhoelst et al. (2021).

point.in.quadrilateral = function(
# https://stackoverflow.com/a/2752753
      # Coordinates of the test point
      xp, yp,
      # Coordinates of the quadrilateral vertices
      x1, y1, x2, y2, x3, y3, x4, y4)
    left.of.edge(xp, yp, x1, y1, x2, y2) &
    left.of.edge(xp, yp, x2, y2, x3, y3) &
    left.of.edge(xp, yp, x3, y3, x4, y4) &
    left.of.edge(xp, yp, x4, y4, x1, y1)

left.of.edge = function(
      xp, yp,
      x1, y1, x2, y2)
    (x2 - x1) * (yp - y1) - (xp - x1) * (y2 - y1) > 0

quadrilateral.area = function(
# https://stackoverflow.com/a/1329561
      x1, y1, x2, y2, x3, y3, x4, y4)
    0.5 * abs(
        x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 +
        x3 * y4 - x4 * y3 + x4 * y1 - x1 * y4)

## * Nitrogen dioxide (NO_2) from a satellite

# We use the TROPOMI instrument aboard Sentinel-5P.
# https://documentation.dataspace.copernicus.eu/APIs/OData.html
# User's guide: https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide.pdf

satellite.no2 = function(the.date, filter.by.quality = T)
  # Automatic download is not implemented because methods to subset
  # the data server-side, which we'd prefer, are coming but don't yet
  # exist (as of 17 Apr 2023).
   {keep.quality = .5  # On a scale of [0, 1], where 1 is best.
    vars.to.get = c(
        lon = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds",
        lat = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds",
        quality = "PRODUCT/qa_value",
        no2.mol.m2 = "PRODUCT/nitrogendioxide_tropospheric_column",
        no2.prec.mol.m2 = "PRODUCT/nitrogendioxide_tropospheric_column_precision",
        air.mass.factor = "PRODUCT/air_mass_factor_troposphere",
        cloud.fraction = "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_crb",
        cloud.pressure.Pa = "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_pressure_crb",
        surface.pressure.Pa = "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure")

    files = satellite.file.ids()[
        lubridate::as_date(date.begin) == the.date,
        .(
           time = date.begin,
           path = file.path(
               tropomi.dir,
               sprintf("%d.nc", date.begin)))]
    rbindlist(lapply(seq_len(nrow(files)), function(fi)
       {o = nc_open(files[fi, path])
        on.exit(nc_close(o))
        d = as.data.table(c(
           unlist(rec = F, lapply(vars.to.get[c("lon", "lat")],
               function(vname) lapply(1 : 4, function(i.corner)
                   as.numeric(ncvar_get(o, vname)[i.corner,,])))),
           lapply(vars.to.get[!(vars.to.get %in% vars.to.get[c("lon", "lat")])],
               function(vname) as.numeric(ncvar_get(o, vname)))))
        setnames(d, str_replace(colnames(d),
            "^(lon|lat)([0-9])$", "\\1.c\\2"))
        d = d[with(study.bbox,
            (!filter.by.quality | quality >= keep.quality) &
            (lon.c1 >= lon.min | lon.c2 >= lon.min |
                lon.c3 >= lon.min | lon.c4 >= lon.min) &
            (lon.c1 <= lon.max | lon.c2 <= lon.max |
                lon.c3 <= lon.max | lon.c4 <= lon.max) &
            (lat.c1 >= lat.min | lat.c2 >= lat.min |
                lat.c3 >= lat.min | lat.c4 >= lat.min) &
            (lat.c1 <= lat.max | lat.c2 <= lat.max |
                lat.c3 <= lat.max | lat.c4 <= lat.max))]
        if (nrow(d))
            cbind(time = files[fi, time], d)}))}

pm(fst = T,
satellite.file.ids <- function()
   {product.prefix = "S5P_RPRO_L2__NO2_"
      # A reprocessing of TROPOMI NO_2
      # https://sentinels.copernicus.eu/web/sentinel/-/copernicus-sentinel-5-precursor-full-mission-reprocessed-datasets-further-products-release
    max.results = 1000L

    message("Getting satellite-file IDs")
    d = rbindlist(pblapply(
        unique(lubridate::floor_date(unit = "month",
            seq(date.first, date.last, by = 1))),
        \(date)
           {response = GET(sprintf("%s?$top=%d&$filter=%s",
                "http://catalogue.dataspace.copernicus.eu/odata/v1/Products",
                max.results,
                paste(sep = " and ",
                    with(study.bbox, sprintf("%s((%d %d,%d %d,%d %d,%d %d,%d %d))%s",
                        "OData.CSC.Intersects(area=geography'SRID=4326;POLYGON",
                        lon.min, lat.min, lon.max, lat.min,
                        lon.max, lat.max, lon.min, lat.max,
                        lon.min, lat.min,
                        "')")),
                    sprintf("startswith(Name, '%s')",
                        product.prefix),
                    sprintf("ContentDate/Start ge %sT00:00:00.000Z",
                        max(date.first, date)),
                    sprintf("ContentDate/Start lt %sT00:00:00.000Z",
                        min(date.last + 1, date + months(1))))))
            stop_for_status(response)
            response = jsonlite::fromJSON(simplifyDataFrame = F,
                content(response, type = "text", encoding = "ASCII"))
            assert(!("@odata.nextLink" %in% names(response)))
            rbindlist(lapply(response$value, \(x) data.table(
                date.begin = lubridate::as_datetime(x$ContentDate$Start),
                date.end = lubridate::as_datetime(x$ContentDate$End),
                id = x$Id)))}))

    setkey(d, date.begin)
    d})

earthdata.creds = function()
   {creds = Sys.getenv(names = F,
        c("EARTHDATA_USERNAME", "EARTHDATA_PASSWORD"))
    if (any(creds == ""))
        stop("You need to set the environment variables EARTHDATA_USERNAME and EARTHDATA_PASSWORD. If you don't have an account, you can get one at https://urs.earthdata.nasa.gov/users/new")
    c("--user", paste(collapse = ":", creds),
        "--cookie", "",
        "--location-trusted")}

## * NO_2 from ground stations

# Documentation: https://www.pandonia-global-network.org/wp-content/uploads/2021/10/PGN_DataProducts_Readme_v1-8-4.pdf

pandonia.base.url = "http://data.pandonia-global-network.org"

pandonia.retrieval.codes = c(
    no2.trop = "(nvh|nuh)",
    no2.total = "nvs")

pm(fst = T,
ground.stations.raw <- function()
   {n.beginning.bytes = 1000

    rbindlist(fill = T, unlist(rec = F, pblapply(
        setdiff(pandonia.items(pandonia.base.url),
            c("calibrationfiles", "operationfiles", "robots.txt")),
        function(site)
           unlist(rec = F, lapply(
               pandonia.items(paste0(pandonia.base.url, "/", site)),
               function(stn.id)
                  {p = paste0(pandonia.base.url, "/", site, "/", stn.id)
                   # The data we want comes from level-2 (L2) processing.
                   if (!("L2" %in% pandonia.items(p)))
                       return()
                   lapply(
                       str_subset(pandonia.items(paste0(p, "/", "L2")), sprintf("_r(%s)",
                           paste(collapse = "|", pandonia.retrieval.codes))),
                       function(fname)
                         {# Get the first `n.beginning.bytes` of the file.
                          r = GET(paste0(p, "/L2/", fname),
                              add_headers(Range = paste0("0-",
                                  n.beginning.bytes - 1)))
                          if (status_code(r) == 404)
                             # A broken link. Ignore it.
                              return()
                          stop_for_status(r)
                          # Read and return the file header.
                          x = str_split_fixed(pattern = ": ", n = 2,
                              str_split(pattern = "\n",
                              str_split_fixed(pattern = "\n----------", n = 2,
                                  content(r, "text", encoding = "Windows-1252"))[,1])[[1]])
                          as.data.table(as.list(`names<-`(
                              x[,2], x[,1])))})})))))})

pm(pandonia.items <- function(url)
  # Read the directory at the given web page.
     {r = GET(url)
      stop_for_status(r)
      x = str_match_all(content(r, "text", encoding = "UTF-8"),
          '<a href="([^/"]+)/?">')[[1]][,2]
      x[x != ".."]})

ground.stations = function()
   {d = copy(ground.stations.raw())
    setnames(d, str_replace_all(colnames(d), "[ \\]\\[]+", "."))

    # Filter by location.
    d[, lon := as.numeric(Location.longitude.deg.)]
    d[, lat := as.numeric(Location.latitude.deg.)]
    d = d[with(study.bbox,
        lon.min <= lon & lon <= lon.max &
        lat.min <= lat & lat <= lat.max)]

    # Filter by date.
    d[Data.end.time == "NONE", Data.end.time := NA]
    d[, Data.end.time := lubridate::fast_strptime(lt = F,
        str_replace(Data.end.time, fixed("."), ""),
        "%Y%m%dT%H%M%SZ")]
    d = d[is.na(Data.end.time) |
        lubridate::as_date(Data.end.time) >= date.first]

    # Reduce to one row per station.
    d[, stn := factor(sprintf("%s-%s-%s",
        Short.location.name, Instrument.number, Spectrometer.number))]
    d = d[, by = stn,
        {slice = data.table(
             Short.location.name,
             Instrument.number, Spectrometer.number, Instrument.type,
             lon, lat,
             Location.altitude.m.)
         assert(nrow(unique(slice)) == 1)
         slice = unique(slice)
         fnames = lapply(pandonia.retrieval.codes, function(x)
             {fname = str_subset(File.name, x)
              assert(length(fname) <= 1)
              if (length(fname)) fname else NA_character_})
         c(slice, fnames)}]

    assert(!anyDuplicated(d$stn))
    setkey(d, stn)
    setcolorder(d, "stn")
    d}

pm(fst = T,
ground.obs <- function(no2.kind)
   {keep.qualities = c("assured high quality", "not-assured high quality")
      # Imitating Verhoelst et al. (2021), p. 494
    assert(no2.kind %in% names(pandonia.retrieval.codes))

    rbindlist(pblapply(cl = n.workers, 1 : nrow(ground.stations()), function(stn.i)
       {station = ground.stations()[stn.i]
        if (is.na(station[, get(no2.kind)]))
            return()

        # Download the data file.
        path = with(station, download(
            paste(sep = "/",
                pandonia.base.url,
                Short.location.name,
                sprintf("%s%ss%s",
                    Instrument.type,
                    Instrument.number,
                    Spectrometer.number),
                "L2",
                get(no2.kind)),
            file.path("pandonia", get(no2.kind))))

        # Parse the column descriptions.
        pieces = str_split(readChar(path, file.info(path)$size),
            regex("^----------+\n", multiline = T))[[1]]
        col.descs = str_match_all(
            pieces[2],
            regex("^Column \\d+: (.+)", multiline = T))[[1]][,2]
        max.col = -1 + as.integer(str_match(
            pieces[2],
            regex("^From Column (\\d+):", multiline = T))[,2])
        columns = c(
            time = "date and time",
            quality = "L2 data quality flag for nitrogen dioxide",
            no2 = unname(c(
               no2.trop = "Nitrogen dioxide tropospheric vertical column amount",
               no2.total = "Nitrogen dioxide total vertical column amount")
                   [no2.kind]))
        columns = sapply(columns, function(x)
            str_which(col.descs, fixed(x)))

        # Read the data rows.
        text = pieces[3]
        if (text == "")
            return()
        if (!is.na(max.col))
          # The number of columns can vary. Ultimately we don't want
          # to use any of the extra columns, but they confuse `fread`,
          # so cut them off.
            text = str_replace_all(text,
               regex(multiline = T,
                   sprintf("^(\\S+)(( \\S+){%d}).+", max.col - 1)),
               "\\1\\2")
        d = fread(text = text,
            sep = " ",
            select = unname(columns),
            showProgress = F)
        setnames(d, names(columns))

        # Parse times.
        d[, time := lubridate::fast_strptime(lt = F,
            str_replace(time, "(\\.\\d+)Z$", "Z"),
              # Discard fractional seconds, which can't be read by
              # `fast_strptime` and aren't important for this study.
            "%Y%m%dT%H%M%SZ")]
        assert(!anyNA(d$time))
        # Discard observations outside the study period.
        d = d[
            date.first <= lubridate::as_date(time) &
            lubridate::as_date(time) <= date.last]
        if (!nrow(d))
            return()

        # Get quality codes.
        m = str_match_all(col.descs[columns["quality"]],
            "(\\d+)=([^,]+)")[[1]]
        d[, quality := factor(
            `names<-`(m[,3], m[,2])[as.character(quality)])]
        assert(!anyNA(d$quality))
        # Filter observations by quality.
        d = d[as.character(quality) %in% keep.qualities]
        if (!nrow(d))
            return()

        # Interpret missing-value codes.
        missing.code = as.numeric(str_match(col.descs[columns["no2"]],
            "(\\S+)=retrieval not successful")[,2])
        assert(!is.na(missing.code))
        # Oddly, the code never occurs in the data.
        assert(!any(d$no2 == missing.code))
        # Standardize NO2 units.
        no2.unit = str_match(col.descs[columns["no2"]], "\\[(.+?)\\]")[,2]
        d[, no2 := no2 * (
            if (no2.unit == "moles per square meter")
                1
            else if (no2.unit == "Dobson Units")
                4.4615e-4
                 # Page 17 of https://www.pandonia-global-network.org/wp-content/uploads/2021/10/PGN_DataProducts_Readme_v1-8-4.pdf
            else
                stop())]
        setnames(d, "no2", "no2.mol.m2")

        cbind(
            stn = station[, stn],
            d)}))})

## * Matchup

pm(fst = T,
ground.no2.at.satellite <- function(ground.no2.kind,
        min.dist.minutes = 30,
          # How temporally close a ground observation has to be in
          # order to be matched up to a satellite observation. The
          # default is from Verhoelst et al. (2021), p. 494.
        min.ground.obs = 1L,
        all.times = F)
          # Return all matching observation times, instead of NO_2
          # values.
   {obs = ground.obs(ground.no2.kind)
    stations = ground.stations()[stn %in% obs$stn]
    rbindlist(pblapply(seq_along(dates.all), cl = n.workers, function(date.i)
       {d.satellite = satellite.no2(dates.all[date.i])
        if (!nrow(d.satellite))
            return()
        d.satellite[, i.satellite := .I]
        rbindlist(lapply(stations$stn, function(the.stn)
           {sat = d.satellite[point.in.quadrilateral(
                stations[.(the.stn), lon], stations[.(the.stn), lat],
                lon.c1, lat.c1, lon.c2, lat.c2,
                lon.c3, lat.c3, lon.c4, lat.c4)]
            if (!nrow(sat))
                return()
            rbindlist(lapply(1 : nrow(sat), function(i.sat)
               {os = obs[stn == the.stn]
                # Use the TROPOMI scan start times as the overpass times.
                os[, abs.timediff.s := abs(as.numeric(
                    difftime(time, sat[i.sat, time], units = "secs")))]
                os = os[abs.timediff.s <= 60 * min.dist.minutes]
                if (nrow(os) < min.ground.obs)
                    return()
                if (all.times)
                    data.table(
                        stn = the.stn,
                        time.satellite = sat[i.sat, time],
                        time.ground = os$time)
                else
                    sat[i.sat, .(
                        stn = the.stn,
                        n.ground = nrow(os),
                        no2.ground = os[, mean(no2.mol.m2)],
                        abs.timediff.mean.s = os[, mean(abs.timediff.s)],
                        time.satellite = time,
                        i.satellite,
                        no2.satellite = no2.mol.m2,
                        no2.satellite.prec = no2.prec.mol.m2,
                        satellite.cell.area.deg2 = quadrilateral.area(
                            lon.c1, lat.c1, lon.c2, lat.c2,
                            lon.c3, lat.c3, lon.c4, lat.c4),
                        air.mass.factor,
                        cloud.fraction,
                        cloud.pressure.Pa,
                        surface.pressure.Pa)]}))}))}))})

## * References

# Verhoelst, T., Compernolle, S., Pinardi, G., Lambert, J.-C., Eskes, H. J., Eichmann, K.-U., … Zehner, C. (2021). Ground-based validation of the Copernicus Sentinel-5P TROPOMI NO_2 measurements with the NDACC ZSL-DOAS, MAX-DOAS and Pandonia global networks. Atmospheric Measurement Techniques, 14, 481–510. doi:10.5194/amt-14-481-2021
