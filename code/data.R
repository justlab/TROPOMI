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
# User's guide: https://sentinel.esa.int/documents/247904/2474726/Sentinel-5P-Level-2-Product-User-Manual-Nitrogen-Dioxide.pdf

pm(fst = T,
satellite.values <- \(no2.kind)
   {ids = satellite.file.ids(no2.kind)$id
    paths = satellite.output.path(ids)
    to.get = satellite.file.ids(no2.kind)[!file.exists(paths), id]
    if (length(to.get))
       {message("Files to download: ", scales::comma(length(to.get)))
        dir.create(dirname(paths[1]), showWarnings = F)
        pbwalk(cl = 4L, to.get, \(id)
            satellite.values.for.file(no2.kind, id))}
    message("Reading output files")
    d = rbindlist(pblapply(cl = n.workers, seq_along(ids), \(i)
       {d = qs::qread(paths[i])
        if (is.null(d)) NULL else cbind(d, sat.file.id = factor(ids[i]))}))
    message("Reducing")
    assert(!any(str_detect(levels(d$stn), " ")))
      # We'll use a space to separate station names, so make sure none
      # of them contain spaces themselves.
    d = d[,
        by = .(sat.file.id, sat.ix.x, sat.ix.y),
        .SDcols = setdiff(colnames(d),
            c("sat.file.id", "sat.ix.x", "sat.ix.y", "stn")),
        cbind(head(.SD, 1), stn =
            paste0(sort(as.character(stn)), collapse = " "))]
    message("Cleaning up")
    d[, stn := factor(stn)]
    setkey(d, sat.time, stn)
    setcolorder(d, c("sat.time", "stn", "sat.file.id"))
    d[]})

satellite.output.path = \(id)
    intermediate("tropomi", paste0(id, ".qs"))

pm(fst = T,
satellite.file.ids <- function(ground.no2.kind)
  # Get a list of all satellite files we want via OData.
  # https://documentation.dataspace.copernicus.eu/APIs/OData.html
   {product.prefix = "S5P_RPRO_L2__NO2_"
      # A reprocessing of TROPOMI NO_2
      # https://sentinels.copernicus.eu/web/sentinel/-/copernicus-sentinel-5-precursor-full-mission-reprocessed-datasets-further-products-release
    filename.regex = "_020400_[0-9T]+\\.nc$"
      # This ensures we only get the specified `processor_version`
      # ("020400" in the filename corresponds to `processor_version`
      # 2.4.0).
    max.results = 1000L

    ewkt = st_as_text(
        st_combine(st_as_sf(
            ground.stations(ground.no2.kind)[, .(lon, lat)],
            coords = c(1, 2), crs = crs.lonlat)),
        EWKT = T)

    message("Getting satellite-file IDs")
    d = rbindlist(pblapply(
        unique(lubridate::floor_date(unit = "month",
            seq(date.first, date.last, by = 1))),
        \(date)
           {response = GET(sprintf("%s?$top=%d&$filter=%s",
                "http://catalogue.dataspace.copernicus.eu/odata/v1/Products",
                max.results,
                paste(sep = " and ",
                    sprintf("%s(area=geography'%s')",
                        "OData.CSC.Intersects", ewkt),
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
                filename = x$Name,
                id = x$Id)))}))

    d = d[str_detect(filename, filename.regex)]
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

main.satellite.vars = c("PRODUCT/latitude", "PRODUCT/longitude", "PRODUCT/qa_value", "PRODUCT/nitrogendioxide_tropospheric_column", "PRODUCT/nitrogendioxide_tropospheric_column_precision", "PRODUCT/nitrogendioxide_tropospheric_column_precision_kernel", "PRODUCT/air_mass_factor_troposphere", "PRODUCT/air_mass_factor_total", "PRODUCT/tm5_tropopause_layer_index", "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle", "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_azimuth_angle", "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_zenith_angle", "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_azimuth_angle", "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/geolocation_flags", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/processing_quality_flags", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/number_of_spectral_points_in_retrieval", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/number_of_iterations", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/wavelength_calibration_offset", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/wavelength_calibration_offset_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/wavelength_calibration_stretch", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/wavelength_calibration_stretch_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/wavelength_calibration_chi_square", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_stratospheric_column", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_stratospheric_column_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column_precision_kernel", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_summed_total_column", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_summed_total_column_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/ozone_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/ozone_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/oxygen_oxygen_dimer_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/oxygen_oxygen_dimer_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_liquid_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_liquid_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/ring_coefficient", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/ring_coefficient_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_fraction_crb_nitrogendioxide_window", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_radiance_fraction_nitrogendioxide_window", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/chi_square", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/root_mean_square_error_of_fit", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/degrees_of_freedom", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_stratosphere", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_cloudy", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_clear", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_ghost_column", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_selection_flag", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_fraction_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_fraction_crb_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_pressure_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_pressure_crb_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_height_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_height_crb_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_cloud_albedo_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_scene_albedo", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_scene_albedo_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_apparent_scene_pressure", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_apparent_scene_pressure_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_chi_square", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_continuum_at_reference_wavelength", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_continuum_at_reference_wavelength_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_polynomial_coefficient", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_polynomial_coefficient_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_ring_coefficient", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_ring_coefficient_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_nitrogendioxide_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_nitrogendioxide_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_oxygen_oxygen_dimer_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_oxygen_oxygen_dimer_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_oxygen_oxygen_dimer_slant_column_density_correction_factor", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_ozone_slant_column_density", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_ozone_slant_column_density_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_surface_albedo", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_wavelength_calibration_offset", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_wavelength_calibration_offset_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_wavelength_calibration_stretch", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/O22CLD/o22cld_wavelength_calibration_stretch_precision", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_cloud_fraction_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_cloud_pressure_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_scene_albedo", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_apparent_scene_pressure", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_cloud_albedo_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_surface_albedo", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude_precision", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification", "PRODUCT/SUPPORT_DATA/INPUT_DATA/scaled_small_pixel_variance", "PRODUCT/SUPPORT_DATA/INPUT_DATA/eastward_wind", "PRODUCT/SUPPORT_DATA/INPUT_DATA/northward_wind", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_albedo_nitrogendioxide_window", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_albedo", "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_pressure_crb", "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_crb", "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_albedo_crb", "PRODUCT/SUPPORT_DATA/INPUT_DATA/scene_albedo", "PRODUCT/SUPPORT_DATA/INPUT_DATA/apparent_scene_pressure", "PRODUCT/SUPPORT_DATA/INPUT_DATA/snow_ice_flag", "PRODUCT/SUPPORT_DATA/INPUT_DATA/aerosol_index_354_388")
  # I got these by checking a file for all the variables on the same
  # grid as the nitrogen-dioxode values.
main.windowed.vars = c("PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_cloud_pressure_crb", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_fraction_crb_nitrogendioxide_window", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_radiance_fraction_nitrogendioxide_window", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_stratospheric_column", "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column_precision", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_albedo_nitrogendioxide_window", "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure", "PRODUCT/qa_value")

satellite.values.for.file = \(no2.kind, file.id)
# Download a raw satellite file, save values at station locations, and
# delete the raw file. It's an error to call this function if the
# output file already exists.
   {url.fmt = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products(%s)/$value"
    n.curl.retries = 10L
    zip.path = intermediate("tropomi", paste0(file.id, ".zip"))
    ncdf.path = intermediate("tropomi", paste0(file.id, ".nc"))
    out.path = satellite.output.path(file.id)
    min.quality = .5  # On a scale of [0, 1], where 1 is best.
    n.kernel.layers = 34L
    bounds.vars = c(
        lon = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds",
        lat = "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")
    dv = switch(no2.kind,
        no2.total = "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column",
        stop())
    window.size = 3L
      # Windowed variables will be constructed with a rectangular
      # window that has at most `2*window.size + 1` cells in each
      # dimension.

    if (file.exists(out.path))
        stop(paste0("Already exists: ", out.path))

    if (!file.exists(ncdf.path))
       {if (!file.exists(zip.path))
            # Download the ZIP.
            assert(0 == system2("curl", shQuote(c(
                "--no-progress-meter",
                "--fail", "--anyauth", "--location-trusted",
                "--retry", n.curl.retries,
                "-H", sprintf("Authorization: Bearer %s", copernicus.access.token()),
                sprintf(url.fmt, file.id),
                "--output", zip.path))))
        # Extract the NetCDF from the ZIP.
        assert(0 == system(sprintf("unzip -p %s %s > %s",
            shQuote(zip.path),
            shQuote(str_subset(unzip(zip.path, list = T)$Name, "\\.nc$")),
            shQuote(ncdf.path))))
        assert(0 == unlink(zip.path))}

    # Crack open the NetCDF to retrieve the values of interest.
    o = nc_open(ncdf.path)
    done.with.netcdf = F
    on.exit(
       {nc_close(o)
        if (done.with.netcdf)
            assert(0 == unlink(ncdf.path))})
    done = \(result)
       {qs::qsave(result, out.path)
        done.with.netcdf <<- T
        T}
    mvs = sapply(main.satellite.vars, simplify = F,
        \(vname) ncvar_get(o, vname))
    quality = mvs$"PRODUCT/qa_value"

    # Get longitude and latitude corners.
    d = as.data.table(unlist(rec = F, sapply(simplify = F, bounds.vars,
        function(vname) lapply(1 : 4, function(i.corner)
            as.numeric(ncvar_get(o, vname)[i.corner,,])))))
    setnames(d, str_replace(colnames(d),
        "^(lon|lat)([0-9])$", "\\1.c\\2"))
    # Use them to find the index in the satellite-variable vectors
    # corresponding to each station.
    stations = ground.stations(no2.kind)
    sat.ix = sapply(seq_len(nrow(stations)), \(stn.i)
        d[, which(point.in.quadrilateral(
            stations[stn.i, lon], stations[stn.i, lat],
            lon.c1, lat.c1, lon.c2, lat.c2,
            lon.c3, lat.c3, lon.c4, lat.c4))][1])
    # NA-ify `sat.ix` corresponding to cells that don't meet the
    # quality threshold.
    sat.ix[quality[sat.ix] < min.quality] = NA_real_
    if (all(is.na(sat.ix)))
      # We won't use anything from this file. Save a stub and bail out
      # early.
        return(done(NULL))
    # Reduce `d` to cells with stations and add the column `stn`.
    # (Rows of `d` may be duplicated if a cell has more than one
    # station.)
    d = cbind(
        stations[which(!is.na(sat.ix)), .(stn)],
        d[sat.ix[!is.na(sat.ix)]])
    sat.ix = sat.ix[!is.na(sat.ix)]
    # Find indices of all cells in the window around each selected
    # cell.
    # (It appears that, in the shape `ncvar_get` returns values from
    # these files, the row is the x-coordinate and the column is the
    # y-coordinate, even though that sounds sideways.)
    s = dim(mvs[[1]])
    ary = matrix(seq_len(s[1] * s[2]), nrow = s[1], ncol = s[2])
    sr = row(ary)
    sc = col(ary)
    d[, sat.ix.x := sr[sat.ix]]
    d[, sat.ix.y := sc[sat.ix]]
    sat.window.ix = lapply(sat.ix, \(i) setdiff(y = i, ary[
      # The `setdiff` means our windowed statistics will exclude
      # the center point.
        max(1, sr[i] - window.size) : min(s[1], sr[i] + window.size),
        max(1, sc[i] - window.size) : min(s[2], sc[i] + window.size)]))

    # Now collect more substantive variables.
    # First, time.
    sat.time.str = ncvar_get(o, "PRODUCT/time_utc")[d$sat.ix.y]
    d[, sat.time := lubridate::parse_date_time2(sat.time.str, "YmdHMOS")]
    assert(!anyNA(d$sat.time[!is.na(sat.time.str)]))
      # Assert that every nonmissing string parsed successfully.
    # Now get point values of the variables on the main grid.
    d = cbind(d,
       as.data.table(sapply(main.satellite.vars, simplify = F,
           \(vname) mvs[[vname]][sat.ix])))
    # Now get windowed statistics. Most are computed without quality
    # thresholding, but a few are.
    for (vname in main.windowed.vars)
        d[, paste0(vname, ".wmean") := sapply(sat.window.ix, \(i)
            mean(mvs[[vname]][i], na.rm = T))]
    ok = which(quality >= min.quality)
    d[, paste0(dv, c(".wmean", ".wmad", ".wmiss")) :=
        rbindlist(lapply(sat.window.ix, \(i)
           {v = mvs[[dv]][intersect(ok, i)]
            v = v[!is.na(v)]
            data.frame(
                wmean = mean(v),
                wmad = mean(abs(v - mean(v))),
                wmiss = mean(is.na(mvs[[dv]][i]) | quality[i] < min.quality))}))]
    # Now get each layer of the averaging kernel. We'll represent each
    # layer as its own variable in our result.
    kernel = ncvar_get(o, "PRODUCT/averaging_kernel")
    assert(dim(kernel)[1] == n.kernel.layers)
    for (i.layer in seq_len(n.kernel.layers))
        d[, sprintf("PRODUCT/averaging_kernel.layer%02d", i.layer) :=
            kernel[i.layer,,][sat.ix]]

    # Save the result.
    return(done(d))}

copernicus.access.token = \()
   {url = "https://identity.cloudferro.com/auth/realms/CDSE/protocol/openid-connect/token"

    creds = Sys.getenv(names = F,
        c("COPERNICUS_DATASPACE_EMAIL", "COPERNICUS_DATASPACE_PASSWORD"))
    if (any(creds == ""))
        stop("You need to set the environment variables COPERNICUS_DATASPACE_EMAIL and COPERNICUS_DATASPACE_PASSWORD. If you don't have an account, you can get one at https://dataspace.copernicus.eu")

    r = RETRY("POST", url,
        pause_min = 30, pause_base = 30, pause_cap = 30*60,
        times = n.workers + 1,
        encode = "form", body = list(
            username = creds[1],
            password = creds[2],
            client_id = "cdse-public",
            grant_type = "password"))
    stop_for_status(r)
    jsonlite::fromJSON(content(r, "text", encoding = "ASCII"))$access_token}

## * NO_2 from ground stations

# Documentation: https://www.pandonia-global-network.org/wp-content/uploads/2021/10/PGN_DataProducts_Readme_v1-8-4.pdf

pandonia.base.url = "http://data.pandonia-global-network.org"

pandonia.retrieval.codes = c(
    no2.trop = "(nvh|nuh)",
    no2.total = "nvs")

pm(fst = T,
ground.stations.raw <- function(no2.kind)
   {n.beginning.bytes = 1000

    rbindlist(fill = T, unlist(rec = F, pblapply(
        cl = n.workers,
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
                        str_subset(pandonia.items(paste0(p, "/", "L2")), sprintf("_r%s",
                            pandonia.retrieval.codes[[no2.kind]])),
                        function(fname)
                           {# Get the first `n.beginning.bytes` of the file.
                            r = RETRY("GET", quiet = T,
                                paste0(p, "/L2/", fname),
                                add_headers(Range = paste0("bytes=0-",
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
   {r = RETRY("GET", quiet = T, url)
    stop_for_status(r)
    x = str_match_all(content(r, "text", encoding = "UTF-8"),
        '<a href="([^/"]+)/?">')[[1]][,2]
    x[x != ".."]})

ground.stations.candidate = function(no2.kind)
   {d = copy(ground.stations.raw(no2.kind))
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

    rbindlist(pblapply(cl = n.workers, 1 : nrow(ground.stations.candidate(no2.kind)), function(stn.i)
       {station = ground.stations.candidate(no2.kind)[stn.i]
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
        pieces = str_split(readr::read_file(path),
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
        # Drop negative observations.
        d = d[no2 >= 0]
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

pm(fst = T,
ground.stations <- function(no2.kind)
    ground.stations.candidate(no2.kind)[.(
        sort(unique(ground.obs(no2.kind)$stn)))])

## * Matchup

pm(fst = T,
satellite.ground.matchup <- function(no2.kind,
        min.dist.minutes = 30)
          # How temporally close a ground observation has to be in
          # order to be matched up to a satellite observation. The
          # default is from Verhoelst et al. (2021), p. 494.
   {n.work = 22L

    sat = satellite.values(no2.kind)[, .(stn, sat.time)]
    ground = as.environment(split(ground.obs(no2.kind), by = "stn"))
    rbindlist(pblapply(cl = n.work, 1 : nrow(sat), \(sat.i)
       {sat.obs = sat[sat.i]
        v = unlist(lapply(
            str_split_1(as.character(sat.obs$stn), fixed(" ")),
            \(the.stn) ground[[the.stn]][j = no2.mol.m2, i =
                abs(as.numeric(
                    difftime(time, sat.obs$sat.time, units = "secs")))
                < 60*min.dist.minutes]))
        if (!length(v)) NULL else list(
            sat.i = sat.i,
            n.ground.obs = length(v),
            y.ground = median(v))}))})

## * References

# Verhoelst, T., Compernolle, S., Pinardi, G., Lambert, J.-C., Eskes, H. J., Eichmann, K.-U., … Zehner, C. (2021). Ground-based validation of the Copernicus Sentinel-5P TROPOMI NO_2 measurements with the NDACC ZSL-DOAS, MAX-DOAS and Pandonia global networks. Atmospheric Measurement Techniques, 14, 481–510. doi:10.5194/amt-14-481-2021
