source("code/data.R")

suppressPackageStartupMessages(
   {library(ggplot2)
    library(ParamHelpers)})

dv = "y.error"
all.ivs = c(
    "time.satellite",
    "stn.lon", "stn.lat",
    "stn.dist.m",
    "satellite.x.index",
    "satellite.cell.area.km2",
    "y.sat.se",
    "y.sat.wmean",
    "y.sat.wmad",
    "y.sat.wmiss",
    "no2.satellite.stratospheric",
    "quality",
    "pressure.diff.Pa",
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
    dv.raw = switch(no2.kind,
        no2.total = "nitrogendioxide_total_column",
        stop())

    d = satellite.ground.matchup(no2.kind)
    d = cbind(satellite.values(no2.kind)[d$sat.i], d)
    setnames(d, str_extract(colnames(d), "[^/]+$"))
    assert(!anyDuplicated(colnames(d)))

    d = d[, .(
        y.ground = no2.unit.factor * y.ground,
        n.ground.obs,
        y.sat = no2.unit.factor * get(dv.raw),
        y.sat.se = no2.unit.factor * get(paste0(dv.raw, "_precision")),
          # The `proposed_standard_name` in the user guide is
          # "atmosphere_mole_content_of_nitrogen_dioxide_standard_error".
        y.sat.wmean = no2.unit.factor * get(paste0(dv.raw, ".wmean")),
        y.sat.wmad = no2.unit.factor * get(paste0(dv.raw, ".wmad")),
        y.sat.wmiss = get(paste0(dv.raw, ".wmiss")),
        sat.time,
        sat.lon = longitude,
        sat.lat = latitude,
        satellite.x.index = sat.ix.x,
        satellite.cell.area.km2 = as.numeric(units::set_units(value = "km^2",
            st_area(st_as_sf(crs = crs.lonlat, `[`(
                melt(
                    d[, .(
                        lon.c1, lat.c1, lon.c2, lat.c2,
                        lon.c3, lat.c3, lon.c4, lat.c4,
                        lon.c5 = lon.c1,
                        lat.c5 = lat.c1,
                        polygon = .I)],
                    measure = patterns("lon", "lat")),
                by = "polygon",
                j = .(g =
                    list(st_polygon(list(cbind(value1, value2)))))))))),
        no2.satellite.stratospheric = no2.unit.factor * nitrogendioxide_stratospheric_column,
        quality = qa_value,
        pressure.diff.Pa = fresco_cloud_pressure_crb - surface_pressure,
        wind.east.m.s = eastward_wind,
        wind.north.m.s = northward_wind,
        angle.solar.zenith.deg = solar_zenith_angle,
        angle.solar.azimuth.deg = solar_azimuth_angle,
        angle.view.zenith.deg = viewing_zenith_angle,
        angle.view.azimuth.deg =  viewing_azimuth_angle,
        albedo = surface_albedo_nitrogendioxide_window,
        cloud.fraction = cloud_fraction_crb_nitrogendioxide_window,
        cloud.radiance.fraction = cloud_radiance_fraction_nitrogendioxide_window,
        air.mass.factor.cloudy = air_mass_factor_cloudy,
        air.mass.factor.clear = air_mass_factor_clear,
        air.mass.factor.trop = air_mass_factor_troposphere)]

    assert(all(d$y.sat >= 0) && all(d$y.ground >= 0))
    d[, y.error := y.sat - y.ground]

    d[, utc.hour.satellite := as.numeric(difftime(sat.time,
        lubridate::floor_date(sat.time, "day"), "UTC", "hours"))]
    d[, utc.day.satellite := as.numeric(difftime(sat.time,
        lubridate::floor_date(sat.time, "year"), "UTC", "days"))]
    d[, time.satellite := as.numeric(sat.time)]

    # Split `d` into clusters according to the closest station
    # location.
    stn.clusters = ground.stations(no2.kind)[,
       by = .(lon, lat),
       .(name = factor(paste(collapse = "/", sort(unique(Short.location.name)))))]
    assert(!anyDuplicated(stn.clusters$name))
    levels(stn.clusters$name)[levels(stn.clusters$name) == "AtlantaGA/AtlantaGA-SouthDeKalb"] = "AtlantaGA"
      # Simplify one of the cluster names.
    setkey(stn.clusters, name)
    f = \(x)
        st_as_sf(x, coords = c(1, 2), crs = crs.lonlat)
    dists = st_distance(
        f(d[, .(sat.lon, sat.lat)]),
        f(stn.clusters[, .(lon, lat)]))
    d[, cluster := stn.clusters[apply(dists, 1, which.min), name]]
    d[, stn.dist.m := apply(dists, 1, min)]
    d[, c("stn.lon", "stn.lat") := stn.clusters[.(d$cluster), .(lon, lat)]]

    # Drop the Altzomoni cluster. This station is at a much higher
    # elevation than all the others.
    d = d[cluster != "Altzomoni"]
    d[, cluster := droplevels(cluster)]

    # Split the clusters into folds.
    cluster.folds = with.temp.seed(1337L,
        sample(rep(1 : n.folds, len = length(unique(d$cluster)))))
    d[, fold := cluster.folds[match(cluster, unique(cluster))]]

    assert(!anyNA(d[, -c("y.sat.wmean", "y.sat.wmad")]))
    assert(all(d[,
        is.na(y.sat.wmean) == is.na(y.sat.wmad) &
        is.na(y.sat.wmean) == (y.sat.wmiss == 1)]))
    d}

correlations = \()
   {cors = pcaPP::cor.fk(
        na.omit(data.for.modeling()[, mget(c(all.ivs, dv))]))
    cors = melt(
        cbind(as.data.table(cors), v1 = rownames(cors)),
        id.vars = "v1", variable.name = "v2")
    cors = (cors
        [as.character(v1) < as.character(v2)]
        [abs(value) >= 1/4]
        [order(-abs(value))])
    cors[, value := round(value, 3)]
    cors[]}

pm(model.with.xgboost <- \(ivs = "ALL", cv.unit = "station")
   {if (identical(ivs, "ALL"))
        ivs = all.ivs
    d = data.for.modeling()

    if (cv.unit == "station")
        {}
    else if (cv.unit == "observation")
        d[, fold := with.temp.seed(400L,
            sample(rep(1 : n.folds, len = .N)))]
    else
        stop()

    y.pred = rep(NA_real_, nrow(d))
    d.shap = as.data.table(sapply(ivs, \(x) rep(NA_real_, nrow(d))))
    shap.inter = array(dim = c(nrow(d), length(ivs), length(ivs)))
    dimnames(shap.inter) = list(NULL, ivs, ivs)

    for (fold.i in seq_len(n.folds))
       {message("Fold ", fold.i)
        fit = with.temp.seed(list(ivs, fold.i),
            fit.xgboost(d[fold != fold.i], ivs))
        i = d$fold == fold.i
        y.pred[i] = fit$pred.f(d[i])
        d.shap[i, (ivs) := as.data.table(
            fit$pred.f(d[i], predcontrib = T)[, ivs])]
        shap.inter[i,,] =
            fit$pred.f(d[i], predinteract = T)[, ivs, ivs]}

    punl(y.pred, d.shap, shap.inter)})

fit.xgboost = \(d, ivs, hyperparams = NULL)
   {n.trees = 100L
    n.work = 24L

    if (is.null(hyperparams))
      # Use inner CV to choose hyperparameters.
       {vals = pbsapply(1 : nrow(xgboost.hyperparam.set()), \(p.i)
           {y.pred = rep(NA_real_, nrow(d))
            for (fold.i in sort(unique(d$fold)))
               {fit = fit.xgboost(d[fold != fold.i], ivs,
                    hyperparams = as.list(xgboost.hyperparam.set()[p.i]))
                y.pred[d$fold == fold.i] = fit$pred.f(d[fold == fold.i])}
            mean(abs(y.pred - d[[dv]]))})
        best.i = which.min(vals)
        message(sprintf("Selected hyperparameter vector %d (min %.03f, max %.03f)",
            best.i, min(vals), max(vals)))
        print(xgboost.hyperparam.set()[best.i])
        hyperparams = as.list(xgboost.hyperparam.set()[best.i])}

    model = xgboost::xgboost(
        verbose = 0,
        data = as.matrix(d[, mget(ivs)]),
        label = d[[dv]],
        nrounds = n.trees,
        params = c(
            list(
                nthread = n.work,
                objective = "reg:absoluteerror",
                base_score = median(d[[dv]])),
            hyperparams))

    list(
       model = model,
       pred.f = \(newdata, predcontrib = F, predinteract = F)
          {out = predict(model, as.matrix(newdata[, mget(ivs)]),
               predcontrib = predcontrib, predinteract = predinteract)
           if (!predcontrib && !predinteract)
             # The DV, being `y.sat` minus a nonnegative number, can't
             # meaningfully exceed `y.sat`.
               pmin(newdata$y.sat, out)
           else
               out})}

xgboost.hyperparam.set = \()
   {n.vectors = 100L
    d = makeParamSet(
        makeIntegerParam("max_depth", lower = 3, upper = 9),
        makeDiscreteParam("colsample_bytree", values = seq(.5, 1, by = .1)),
        makeNumericParam("eta", lower = .01, upper = .5),
        makeNumericParam("gamma", lower = -7, upper = 6, trafo = \(x) 2^x),
        makeNumericParam("lambda", lower = -10, upper = 10, trafo = \(x) 2^x),
        makeNumericParam("alpha", lower = -10, upper = 10, trafo = \(x) 2^x))
    d = as.data.table(
       with.temp.seed(as.integer(n.vectors), generateDesign(
           n.vectors, par.set = d, trafo = T,
           fun = lhs::maximinLHS)))
    d[, colsample_bytree := as.numeric(as.character(colsample_bytree))]
    for (dcol in colnames(d))
       {# Round off the selected parameters to a few significant
        # figures.
        if (is.double(d[[dcol]]))
            d[, (dcol) := signif(get(dcol), 2)]}
    d}

pm(fst = T, mem = T,
d.xgb <- \(...)
   {d = data.for.modeling()
    l = model.with.xgboost(...)
    d[, y.error.pred := l$y.pred]
    d[, y.ground.pred := y.sat - y.error.pred]
    d[, paste0("shap.", colnames(l$d.shap)) := l$d.shap]
    d})

summarize.xgboost.results = \(by.expr = NULL, ...)
   {mae = \(x, y) mean(abs(x - y))
    mad = \(x) mae(x, mean(x))

    d.xgb(...)[, keyby = (if (!is.null(by.expr)) eval(by.expr)), .(
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
        "Bias, corrected" = mean(y.ground.pred - y.ground))]}

pretty.cv.summary = \(...)
   {d = summarize.xgboost.results(...)
    for (vname in colnames(d)) d[, (vname) := {v = get(vname);
        if (vname == "Year")
            as.character(v)
        else if (vname == "Month")
            month.abb[v]
        else if (str_detect(tolower(vname), "cor\\.|proportion"))
            sprintf("%.02f", v)
        else if (is.integer(v))
            scales::comma(v)
        else if (is.numeric(v))
            round(v)
        else
            as.character(v)}]
    d[]}

mheatmap = \(x.var, y.var, binwidth)
   {d = d.xgb()[
        !is.na(get(x.var)) & !is.na(get(y.var)),
        .(x = get(x.var), y = get(y.var))]
    ggplot(d) +
        stat_bin2d(aes(x, y),
            binwidth = binwidth) +
        scale_fill_gradient(guide = "none",
            low = "#dddddd", high = "black") +
        xlab(x.var) + ylab(y.var) +
        theme_classic()}

empirical.error.envelope = \(..., target.coverage = 2/3)
  # Find an "error envelope" with an additive and multiplicative
  # term in the fashion of Dark Target validation
  # https://darktarget.gsfc.nasa.gov/validation
   {f = envelope.tester(...)
    optim(c(.05, .1), \(p)
        sum(abs(target.coverage - f(p[1], p[2]))))}

envelope.tester = \(obs, pred, threshold = 200)
   {get.coverage = \(slice, add, mult) slice[, mean(
        pred >= obs - add - mult*obs &
        pred <= obs + add + mult*obs)]
    d = data.table(obs, pred)
    d.lo = d[obs <= threshold]
    d.hi = d[obs > threshold]
    \(add, mult) c(
        lo = get.coverage(d.lo, add, mult),
        hi = get.coverage(d.hi, add, mult))}
