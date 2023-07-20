source("code/data.R")

library(ggplot2)

dv = "y.error"
ivs = c(
    "time.satellite",
    "satellite.x.index",
    "satellite.cell.area.deg2",
    "y.sat.prec",
    "no2.satellite.stratospheric",
    "quality",
    "cloud.pressure.Pa",
    "surface.pressure.Pa",
    "wind.east.m.s",
    "wind.north.m.s",
    "angle.solar.zenith.deg",
    "angle.solar.azimuth.deg",
    "angle.view.zenith.deg",
    "angle.view.azimuth.deg",
    "albedo",
    "cloud.fraction",
    "cloud.radiance.fraction",
    "air.mass.factor.cloudy",
    "air.mass.factor.clear",
    "air.mass.factor.trop")
n.folds = 10L

data.for.modeling = \(no2.kind = "no2.total")
   {no2.unit.factor = 1e6
      # Convert mol / m^2 to μmol / m^2.

    d = satellite.ground.matchup(no2.kind)
    d = cbind(satellite.values(no2.kind)[d$sat.i], d[, .(y.ground)])
    d = d[, .(
        y.ground = no2.unit.factor * y.ground,
        y.sat = no2.unit.factor * switch(no2.kind,
            no2.total = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column`,
            stop()),
        y.sat.prec = no2.unit.factor * switch(no2.kind,
            no2.total = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_total_column_precision`,
            stop()),
        sat.time,
        sat.lon = `PRODUCT/longitude`,
        sat.lat = `PRODUCT/latitude`,
        satellite.x.index = sat.ix.x,
        satellite.cell.area.deg2 = quadrilateral.area(
            lon.c1, lat.c1, lon.c2, lat.c2,
            lon.c3, lat.c3, lon.c4, lat.c4),
        no2.satellite.stratospheric = no2.unit.factor * `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/nitrogendioxide_stratospheric_column`,
        quality = `PRODUCT/qa_value`,
        cloud.pressure.Pa = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/FRESCO/fresco_cloud_pressure_crb`,
        surface.pressure.Pa = `PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure`,
        wind.east.m.s = `PRODUCT/SUPPORT_DATA/INPUT_DATA/eastward_wind`,
        wind.north.m.s = `PRODUCT/SUPPORT_DATA/INPUT_DATA/northward_wind`,
        angle.solar.zenith.deg = `PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle`,
        angle.solar.azimuth.deg = `PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_azimuth_angle`,
        angle.view.zenith.deg = `PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_zenith_angle`,
        angle.view.azimuth.deg =  `PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_azimuth_angle`,
        albedo = `PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_albedo_nitrogendioxide_window`,
        cloud.fraction = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_fraction_crb_nitrogendioxide_window`,
        cloud.radiance.fraction = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/cloud_radiance_fraction_nitrogendioxide_window`,
        air.mass.factor.cloudy = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_cloudy`,
        air.mass.factor.clear = `PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/air_mass_factor_clear`,
        air.mass.factor.trop = `PRODUCT/air_mass_factor_troposphere`)]

    d[, y.error := y.sat - y.ground]

    d[, utc.hour.satellite := as.numeric(difftime(sat.time,
        lubridate::floor_date(sat.time, "day"), "UTC", "hours"))]
    d[, utc.day.satellite := as.numeric(difftime(sat.time,
        lubridate::floor_date(sat.time, "year"), "UTC", "days"))]
    d[, time.satellite := as.numeric(sat.time)]

    # Split `d` into clusters according to the closest station
    # location.
    f = \(x)
        st_as_sf(x, coords = c(1, 2), crs = crs.lonlat)
    d[, cluster := apply(MARGIN = 1, FUN = which.min, st_distance(
        f(d[, .(sat.lon, sat.lat)]),
        f(unique(ground.stations(no2.kind)[, .(lon, lat)]))))]
    # Split the clusters into folds.
    cluster.folds = with.temp.seed(1337L,
        sample(rep(1 : n.folds, len = length(unique(d$cluster)))))
    d[, fold := cluster.folds[match(cluster, unique(cluster))]]

    assert(!anyNA(d))
    d}

correlations = \()
   {cors = pcaPP::cor.fk(
        na.omit(data.for.modeling()[, mget(c(ivs, dv))]))
    cors = melt(
        cbind(as.data.table(cors), v1 = rownames(cors)),
        id.vars = "v1", variable.name = "v2")
    cors = (cors
        [as.character(v1) < as.character(v2)]
        [abs(value) >= 1/4]
        [order(-abs(value))])
    cors[, value := round(value, 3)]
    cors[]}

model.with.ols = \()
   {d = data.for.modeling()

    # Standardize the IVs and the DV, and fill in missing values.
    for (vname in c(ivs, dv))
       {d[, (vname) := (get(vname) - mean(get(vname), na.rm = T))
            / sdn(get(vname), na.rm = T)]
        d[is.na(get(vname)), (vname) := 0]}

    m = lm(reformulate(ivs, dv, intercept = F), data = d)

    list(
        coef = as.data.table(coef(m), keep = T)
            [order(-abs(V2)), .(iv = V1, coef = round(V2, 3))],
        r2 = round(summary(m)$r.squared, 3))}

pm(model.with.xgboost <- \()
   {d = data.for.modeling()

    y.pred = rep(NA_real_, nrow(d))
    d.shap = as.data.table(sapply(ivs, \(x) rep(NA_real_, nrow(d))))

    for (fold.i in seq_len(n.folds))
       {message("Fold ", fold.i)
        fit = fit.xgboost(d[fold != fold.i])
        y.pred[d$fold == fold.i] =
            fit$pred.fun(d[fold == fold.i])
        d.shap[d$fold == fold.i, (ivs) := as.data.table(
            fit$pred.fun(d[fold == fold.i], predcontrib = T)[, ivs])]}

    punl(y.pred, d.shap)})

fit.xgboost = \(d)
   {n.trees = 10L

    xgboost.dart.cvtune(
        d[, mget(c(dv, ivs))],
        dv, ivs,
        n.rounds = n.trees,
        folds = d$fold,
        nthread = n.workers,
        progress = T)}

d.xgb = \()
   {d = data.for.modeling()
    l = model.with.xgboost()
    d[, y.pred := l$y.pred]
    d[, no2.ground.pred := no2.satellite - y.pred]
    d[, paste0("shap.", colnames(l$d.shap)) := l$d.shap]
    d}

summarize.xgboost.results = \()
   {d = data.for.modeling()
    d[, y.ground.pred := y.sat - model.with.xgboost()$y.pred]

    mae = \(x, y) mean(abs(x - y))
    mad = \(x) stats::mad(x, constant = 1)

    d[, .(
        "Cases" = .N,
        "Clusters" = length(unique(cluster)),
        "MAE, raw" = mae(y.ground, y.sat),
        "MAE, corrected" = mae(y.ground, y.ground.pred),
        "Proportion of raw MAE" = mae(y.ground, y.ground.pred) / mae(y.ground, y.sat),
        "Median, ground" = median(y.ground),
        "MAD, ground" = mad(y.ground),
        "MAD, raw" = mad(y.sat),
        "MAD, corrected" = mad(y.ground.pred),
        "Bias, raw" = mean(y.sat - y.ground),
        "Bias, corrected" = mean(y.ground.pred - y.ground),
        "Kendall cor., raw" = pcaPP::cor.fk(y.sat, y.ground),
        "Kendall cor., corrected" = pcaPP::cor.fk(y.ground.pred, y.ground))]}

scatterplot = \(x.var, y.var)
   {d = d.xgb()[
        !is.na(get(x.var)) & !is.na(get(y.var))]
    ggplot(d) +
        geom_point(aes(get(x.var), get(y.var)), alpha = 1/10) +
        xlab(x.var) + ylab(y.var) +
        coord_cartesian(
           xlim = quantile(d[[x.var]], c(.01, .99)),
           ylim = quantile(d[[y.var]], c(.01, .99))) +
        theme_classic()}
