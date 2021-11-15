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

## * Nitrogen dioxide (NO_2) from a satellite

# We use the TROPOMI instrument aboard Sentinel-5P.
# https://catalog.data.gov/dataset/sentinel-5p-tropomi-tropospheric-no2-1-orbit-l2-5-5km-x-3-5km-v1-s5p-l2-no2-hir-at-ges-dis
# Example OPeNDAP links: https://tropomi.gesdisc.eosdis.nasa.gov/opendap/S5P_TROPOMI_Level2/S5P_L2__NO2____HiR.1/2020/043/contents.html
# User's guide: https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide

satellite.no2 = function(the.date, filter.by.quality = T)
   {keep.quality = .5  # On a scale of [0, 1], where 1 is best.

    vars.to.get = c(
        lon = "PRODUCT_SUPPORT_DATA_GEOLOCATIONS_longitude_bounds",
        lat = "PRODUCT_SUPPORT_DATA_GEOLOCATIONS_latitude_bounds",
        quality = "PRODUCT_qa_value",
        no2.trop.mol.m2 = "PRODUCT_nitrogendioxide_tropospheric_column",
        no2.strat.mol.m2 = "PRODUCT_SUPPORT_DATA_DETAILED_RESULTS_nitrogendioxide_stratospheric_column",
        no2.total.mol.m2 = "PRODUCT_SUPPORT_DATA_DETAILED_RESULTS_nitrogendioxide_total_column")
    urls = str_subset(
        str_subset(readLines(tropomi.urls.path), "^#", negate = T),
        sprintf("S5P_L2__NO2____HiR.2/%d/%03d/",
            year(the.date), yday(the.date)))
    rbindlist(lapply(urls, function(u)
       {time = str_match(u, "____(\\d+T\\d+)")[,2]
        o = nc_open(download(
            paste0(str_replace(u, "/data//", "/opendap/"), ".nc4?",
                paste(vars.to.get, collapse = ",")),
            sprintf("tropomi_no2/%s.nc", time),
            curl = earthdata.creds()))
        on.exit(nc_close(o))
        d = as.data.table(c(
           unlist(rec = F, lapply(vars.to.get[c("lon", "lat")],
               function(vname) lapply(1 : 4, function(i.corner)
                   as.numeric(ncvar_get(o, vname)[i.corner,,])))),
           lapply(vars.to.get[vars.to.get != vars.to.get[c("lon", "lat")]],
               function(vname) as.numeric(ncvar_get(o, vname)))))
        setnames(d, str_replace(colnames(d),
            "^(lon|lat)([0-9])$", "\\1.c\\2"))
        cbind(
            time = lubridate::fast_strptime(lt = F,
                time, "%Y%m%dT%H%M%S"),
            d[!filter.by.quality | quality >= keep.quality])}))}

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
   {items = function(url)
      # Read the directory at the given web page.
       {r = GET(url)
        stop_for_status(r)
        x = str_match_all(content(r, "text", encoding = "UTF-8"),
            '<a href="([^/"]+)/?">')[[1]][,2]
        x[x != ".."]}

    n.beginning.bytes = 1000

    rbindlist(fill = T, unlist(rec = F, pblapply(
        setdiff(items(pandonia.base.url),
            c("calibrationfiles", "operationfiles", "robots.txt")),
        function(site)
           unlist(rec = F, lapply(
               items(paste0(pandonia.base.url, "/", site)),
               function(stn.id)
                  {p = paste0(pandonia.base.url, "/", site, "/", stn.id)
                   # The data we want comes from level-2 (L2) processing.
                   if (!("L2" %in% items(p)))
                       return()
                   lapply(
                       str_subset(items(paste0(p, "/", "L2")), sprintf("_r(%s)",
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
                str_replace(Short.location.name,
                    "ManhattanNY", "Manhattan"),
                  # In this case, the directory name doesn't quite
                  # match the short location name.
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
            select = unname(columns))
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
        min.dist.hours = .5,
          # How temporally close a ground observation has to be in
          # order to be matched up to a satellite observation. The
          # default is from Verhoelst et al. (2021), p. 494.
        min.ground.obs = NULL)
   {if (is.null(min.ground.obs))
        min.ground.obs = c(
            no2.trop = 3,
            no2.total = 10)[ground.no2.kind]

    obs = ground.obs(ground.no2.kind)
    stations = ground.stations()[stn %in% obs$stn]
    rbindlist(pblapply(seq_along(dates.all), cl = n.workers, function(date.i)
       {d.satellite = satellite.no2(dates.all[date.i])
        d.satellite[, i.satellite := .I]
        rbindlist(lapply(stations$stn, function(the.stn)
           {sat = d.satellite[point.in.quadrilateral(
                stations[.(the.stn), lon], stations[.(the.stn), lat],
                lon.c1, lat.c1, lon.c2, lat.c2,
                lon.c3, lat.c3, lon.c4, lat.c4)]
            if (!nrow(sat))
                return()
            assert(!anyDuplicated(sat$time))
            rbindlist(lapply(1 : nrow(sat), function(i.sat)
               {# Use the TROPOMI scan start times as the overpass times.
                os = obs[
                    stn == the.stn &
                    abs(difftime(time, sat[i.sat, time], units = "hours")) <=
                        min.dist.hours]
                if (nrow(os) < min.ground.obs)
                    return()
                data.table(
                    stn = the.stn,
                    n.ground = nrow(os),
                    time.satellite = sat[i.sat, time],
                    i.satellite = sat[i.sat, i.satellite],
                    no2.trop.satellite = sat[i.sat, no2.trop.mol.m2],
                    no2.strat.satellite = sat[i.sat, no2.strat.mol.m2],
                    no2.total.satellite = sat[i.sat, no2.total.mol.m2],
                    no2.ground = os[i.sat, mean(no2.mol.m2)])}))}))}))})

## * References

# Verhoelst, T., Compernolle, S., Pinardi, G., Lambert, J.-C., Eskes, H. J., Eichmann, K.-U., … Zehner, C. (2021). Ground-based validation of the Copernicus Sentinel-5P TROPOMI NO_2 measurements with the NDACC ZSL-DOAS, MAX-DOAS and Pandonia global networks. Atmospheric Measurement Techniques, 14, 481–510. doi:10.5194/amt-14-481-2021
