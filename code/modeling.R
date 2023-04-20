source("code/data.R")

relationships = function()
   {d = copy(ground.no2.at.satellite("no2.trop"))
    d[, no2.error := no2.satellite - no2.ground]
    d[, time.satellite := as.numeric(time.satellite)]
    d[, pressure.diff.Pa := cloud.pressure.Pa - surface.pressure.Pa]
    d[, c("lon.stn", "lat.stn") :=
        ground.stations()[.(d$stn), .(lon, lat)]]
    ivs = c(
       "time.satellite",
       "abs.timediff.mean.s",
       "lon.stn",
       "lat.stn",
       "satellite.cell.area.deg2",
       "no2.satellite.prec",
       "air.mass.factor",
       "cloud.fraction",
       "pressure.diff.Pa")
    dv = "no2.error"

    cors = cor(d[, mget(c(ivs, dv))], method = "kendall")
    cors = melt(
        cbind(as.data.table(cors), v1 = rownames(cors)),
        id.vars = "v1", variable.name = "v2")
    cors = (cors
        [as.character(v1) < as.character(v2)]
        [abs(value) >= 1/8]
        [order(-abs(value))])
    cors[, value := round(value, 3)]

    # Standardize the IVs and the DV.
    for (vname in c(ivs, dv))
      d[, (vname) := (get(vname) - mean(get(vname))) / sdn(get(vname))]

    m = lm(reformulate(ivs, dv, intercept = F), data = d)

    list(
       cor = cors,
       coef = as.data.table(coef(m), keep = T)
           [order(-abs(V2)), .(iv = V1, coef = round(V2, 3))],
       r2 = round(summary(m)$r.squared, 3))}
