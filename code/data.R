source("code/common.R")

suppressPackageStartupMessages(
   {library(httr)
    library(ncdf4)})

## * Nitrogen dioxide (NO_2) from a satellite

# We use the TROPOMI instrument aboard Sentinel-5P.
# https://catalog.data.gov/dataset/sentinel-5p-tropomi-tropospheric-no2-1-orbit-l2-5-5km-x-3-5km-v1-s5p-l2-no2-hir-at-ges-dis
# Example OPeNDAP links: https://tropomi.gesdisc.eosdis.nasa.gov/opendap/S5P_TROPOMI_Level2/S5P_L2__NO2____HiR.1/2020/043/contents.html
# User's guide: https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide

satellite.no2 = function(the.date, filter.by.quality = T)
   {keep.quality = .5  # On a scale of [0, 1], where 1 is best.
    vars.to.get = c(
        lon = "PRODUCT_longitude",
        lat = "PRODUCT_latitude",
        quality = "PRODUCT_qa_value",
        no2.mol.m2 = "PRODUCT_SUPPORT_DATA_DETAILED_RESULTS_nitrogendioxide_total_column")
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
        d = as.data.table(lapply(vars.to.get,
           function(vname) as.numeric(ncvar_get(o, vname))))
        cbind(
            time = lubridate::fast_strptime(lt = F,
                time, "%Y%m%dT%H%M%S"),
            d[!filter.by.quality | quality >= keep.quality])}))}

estimate.satellite.grid.increment = function()
   {n.sample.days = 3
    n.sample.points = 1000

    library(FNN)

    set.seed(15)
    sapply(sample.int(length(dates.all), n.sample.days), function(date.i)
        quantile(drop(get.knn(
            satellite.no2(dates.all[date.i])
                [time == sample(unique(time), 1)]
                [order(
                    (lon - mean(range(lon)))^2 +
                    (lat - mean(range(lat)))^2)]
                [seq_len(n.sample.points),
                   .(lon, lat)],
            k = 1)$nn.dist)))}
satellite.grid.increment.degrees = .051
  # From the result of the above, increased slightly.

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

pm(fst = T,
ground.stations.raw <- function()
   {items = function(url)
      # Read the directory at the given web page.
       {r = GET(url)
        stop_for_status(r)
        x = str_match_all(content(r, "text", encoding = "UTF-8"),
            '<a href="([^/"]+)/?">')[[1]][,2]
        x[x != ".."]}

    retrieval.code = "nvs"  # Total column NO_2
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
                       str_subset(items(paste0(p, "/", "L2")),
                          paste0("_r", retrieval.code)),
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

    d[, stn := factor(sprintf("%s-%s-%s",
        Short.location.name, Instrument.number, Spectrometer.number))]
    stopifnot(!anyDuplicated(d$stn))
    setkey(d, stn)
    setcolorder(d, "stn")
    d}

pm(fst = T,
ground.obs <- function()
   {keep.qualities = c("assured high quality", "not-assured high quality")
      # Imitating Verhoelst et al. (2021), p. 494

    rbindlist(pblapply(cl = n.workers, 1 : nrow(ground.stations()), function(stn.i)
       {# Download the data file.
        path = with(ground.stations()[stn.i], download(
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
                File.name),
            file.path("pandonia", File.name)))

        # Parse the column descriptions.
        pieces = str_split(readChar(path, file.info(path)$size),
            regex("^----------+\n", multiline = T))[[1]]
        col.descs = str_match_all(
            pieces[2],
            regex("^Column \\d+: (.+)", multiline = T))[[1]][,2]
        columns = c(
            time = "date and time",
            quality = "L2 data quality flag for nitrogen dioxide",
            no2 = "Nitrogen dioxide total vertical column amount")
        columns = sapply(columns, function(x)
            str_which(col.descs, fixed(x)))

        # Read the data rows.
        d = fread(text = pieces[3],
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
            stn = ground.stations()[stn.i, stn],
            d)}))})

## * Matchup

pm(fst = T,
ground.no2.at.satellite <- function()
   {min.dist.hours = .5
      # How temporally close a ground observation has to be in order
      # to be matched up to a satellite observation.
    min.ground.obs = 10

    obs = ground.obs()
    stations = ground.stations()[stn %in% obs$stn]
    rbindlist(pblapply(seq_along(dates.all), cl = n.workers, function(date.i)
       {d.satellite = satellite.no2(dates.all[date.i])
        d.satellite[, i.satellite := .I]
        rbindlist(lapply(stations$stn, function(the.stn)
           {sat = d.satellite[
                abs(lon - stations[.(the.stn), lon]) <=
                    satellite.grid.increment.degrees/2 &
                abs(lat - stations[.(the.stn), lat]) <=
                    satellite.grid.increment.degrees/2]
            if (!nrow(sat))
                return()
            if (nrow(sat) > 1)
              # There may be multiple matches in both space and time.
                sat = sat[which.min(
                    (lon - stations[.(the.stn), lon])^2 +
                    (lat - stations[.(the.stn), lat])^2)]
            # Use the TROPOMI scan start time as the overpass time.
            os = obs[
                stn == the.stn &
                abs(difftime(time, sat$time, units = "hours")) <=
                    min.dist.hours]
            if (nrow(os) < min.ground.obs)
                return()
            data.table(
                stn = the.stn,
                n.ground = nrow(os),
                time.satellite = sat[, time],
                i.satellite = sat[, i.satellite],
                no2.satellite = sat[, no2.mol.m2],
                no2.ground = os[, mean(no2.mol.m2)])}))}))})

## * References

# Verhoelst, T., Compernolle, S., Pinardi, G., Lambert, J.-C., Eskes, H. J., Eichmann, K.-U., … Zehner, C. (2021). Ground-based validation of the Copernicus Sentinel-5P TROPOMI NO_2 measurements with the NDACC ZSL-DOAS, MAX-DOAS and Pandonia global networks. Atmospheric Measurement Techniques, 14, 481–510. doi:10.5194/amt-14-481-2021
