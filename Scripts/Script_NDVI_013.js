var massifs = ee.FeatureCollection("projects/hale-dynamo-473013-s7/assets/MASSIFS_FORESTIERS_S_013");
Map.centerObject(massifs, 7);
Map.addLayer(massifs, {color: "green"}, "Massifs forestiers");

var ndviCollection = ee.ImageCollection("MODIS/061/MOD13Q1")
  .filterDate("2010-01-01", "2026-01-01")
  .select("NDVI")
  .map(function(img) {
    return img.multiply(0.0001).copyProperties(img, ["system:time_start"]);
  });

var yearsList = ee.List.sequence(2010, 2025);

var annualNDVI = ee.ImageCollection(yearsList.map(function(y) {
  y = ee.Number(y);
  var start = ee.Date.fromYMD(y, 1, 1);
  var end = start.advance(1, "year");
  var img = ndviCollection.filterDate(start, end).mean();
  return img.set("year", y).set("system:time_start", start.millis());
}));

var globalStats = annualNDVI.map(function(img) {
  var stat = img.reduceRegion({
    reducer: ee.Reducer.mean().combine({reducer2: ee.Reducer.minMax(), sharedInputs: true}).combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
    geometry: massifs.geometry(),
    scale: 250,
    maxPixels: 1e9,
    tileScale: 4
  });
  return ee.Feature(null, {
    year: img.get("year"),
    ndvi_mean: stat.get("NDVI_mean"),
    ndvi_min: stat.get("NDVI_min"),
    ndvi_max: stat.get("NDVI_max"),
    ndvi_stdDev: stat.get("NDVI_stdDev")
  });
});
print("Statistiques NDVI globales 2010-2025 (zone d'etude complete) :", globalStats);

var chartGlobal = ui.Chart.feature.byFeature(globalStats, "year", ["ndvi_mean", "ndvi_min", "ndvi_max"])
  .setChartType("LineChart")
  .setOptions({
    title: "Evolution du NDVI (2010-2025) - Zone d'etude complete (25 massifs)",
    hAxis: {title: "Annee", format: "####"},
    vAxis: {title: "NDVI"},
    lineWidth: 2,
    pointSize: 4
  });
print(chartGlobal);

var rowsNested = yearsList.map(function(y) {
  y = ee.Number(y);
  var img = ee.Image(annualNDVI.filter(ee.Filter.eq("year", y)).first());
  var reduced = img.reduceRegions({
    collection: massifs,
    reducer: ee.Reducer.mean(),
    scale: 250,
    tileScale: 4
  }).map(function(f) {
    return f.set("year", y, "ndvi_mean", f.get("mean"));
  });
  return reduced.toList(reduced.size());
});
var perMassifFeatures = ee.FeatureCollection(rowsNested.flatten());
print("Statistiques NDVI par massif 2010-2025 (echantillon) :", perMassifFeatures.limit(10));

var chartByMassif = ui.Chart.feature.groups({
  features: perMassifFeatures,
  xProperty: "year",
  yProperty: "ndvi_mean",
  seriesProperty: "nom_massif"
}).setChartType("LineChart")
  .setOptions({
    title: "Evolution du NDVI moyen par massif (2010-2025)",
    hAxis: {title: "Annee", format: "####"},
    vAxis: {title: "NDVI moyen"},
    lineWidth: 2,
    pointSize: 2,
    width: 900,
    height: 500,
    curveType: "function"
  });
print(chartByMassif);

Export.table.toDrive({
  collection: perMassifFeatures,
  description: "NDVI_par_massif_2010_2025",
  folder: "GEE_exports",
  fileNamePrefix: "NDVI_par_massif_2010_2025",
  fileFormat: "CSV",
  selectors: ["nom_massif", "gid", "year", "ndvi_mean"]
});

Export.table.toDrive({
  collection: globalStats,
  description: "NDVI_global_2010_2025",
  folder: "GEE_exports",
  fileNamePrefix: "NDVI_global_2010_2025",
  fileFormat: "CSV"
});
