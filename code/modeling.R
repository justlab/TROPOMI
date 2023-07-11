source("code/data.R")

library(ggplot2)

dv = "no2.error"
ivs = c(
    "time.satellite",
    "abs.timediff.mean.s",
    "satellite.x.index",
    "satellite.cell.area.deg2",
    "no2.satellite.prec",
    "no2.satellite.stratospheric",
    setdiff(names(satellite.vars), c(
        "lon", "lat", "no2.mol.m2", "no2.prec.mol.m2",
        "no2.stratospheric.mol.m2")))
n.folds = 10L

data.for.modeling = \()
   {d = copy(ground.no2.at.satellite("no2.total"))

    d[, no2.error := no2.satellite - no2.ground]
    for (vname in colnames(d))
       if (startsWith(vname, "no2."))
           d[, (vname) := get(vname) * mol.m2.to.Pmolecule.cm2]

    d[, utc.hour.satellite := as.numeric(difftime(time.satellite,
        lubridate::floor_date(time.satellite, "day"), "UTC", "hours"))]
    d[, utc.day.satellite := as.numeric(difftime(time.satellite,
        lubridate::floor_date(time.satellite, "year"), "UTC", "days"))]
    d[, time.satellite := as.numeric(time.satellite)]

    d[, c("lon.stn", "lat.stn") :=
        ground.stations()[.(d$stn), .(lon, lat)]]
    stns = unique(d$stn)
    stn.folds = with.temp.seed(1337L,
        sample(rep(1 : n.folds, len = length(stns))))
    d[, fold := stn.folds[match(stn, unique(stn))]]

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
    setnames(d,
       c("stn", "no2.ground", "no2.satellite"),
       c("site", "y.ground", "y.sat"))
    d[, y.ground.pred := y.sat - model.with.xgboost()$y.pred]

    mse = \(x, y) mean((x - y)^2)

    d[, .(
        "Cases" = .N,
        "Sites" = length(unique(site)),
        "RMSE, raw" = sqrt(mse(y.ground, y.sat)),
        "RMSE, corrected" = sqrt(mse(y.ground, y.ground.pred)),
        "Proportion of raw MSE" = mse(y.ground, y.ground.pred) / mse(y.ground, y.sat),
        "Median, ground" = median(y.ground),
        "SD, ground" = sdn(y.ground),
        "SD, raw" = sdn(y.sat),
        "SD, corrected" = sdn(y.ground.pred),
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
