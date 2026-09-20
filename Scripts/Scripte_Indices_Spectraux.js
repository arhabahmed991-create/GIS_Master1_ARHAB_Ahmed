var massifs = ee.FeatureCollection("projects/hale-dynamo-473013-s7/assets/MASSIFS_FORESTIERS_S_013");

// MOD09A1 : reflectance de surface MODIS 500m, composites 8 jours.
// Contient Rouge, Proche-Infrarouge, Bleu, Vert, SWIR1, SWIR2 -> permet de calculer
// NDVI, NDMI, NDWI, NBR, SAVI et EVI de facon continue sur toute la periode 2010-2025
// (contrairement a MOD13Q1 qui ne fournit que NDVI/EVI sans les bandes Vert et SWIR1).
// Note: pas de Map.addLayer ici pour eviter de saturer le quota API avec des tuiles
// en plus des calculs de statistiques (qui sont la partie utile de ce script).
var EPS = 1e-6;

var col = ee.ImageCollection("MODIS/061/MOD09A1")
  .filterDate("2010-01-01", "2026-01-01")
  .map(function(img) {
    var scale = 0.0001;
    var red = img.select("sur_refl_b01").multiply(scale);
    var nir = img.select("sur_refl_b02").multiply(scale);
    var blue = img.select("sur_refl_b03").multiply(scale);
    var green = img.select("sur_refl_b04").multiply(scale);
    var swir1 = img.select("sur_refl_b06").multiply(scale);
    var swir2 = img.select("sur_refl_b07").multiply(scale);

    var ndvi = nir.subtract(red).divide(nir.add(red).add(EPS)).clamp(-1, 1).rename("NDVI");
    var ndmi = nir.subtract(swir1).divide(nir.add(swir1).add(EPS)).clamp(-1, 1).rename("NDMI");
    var ndwi = green.subtract(nir).divide(green.add(nir).add(EPS)).clamp(-1, 1).rename("NDWI");
    var nbr = nir.subtract(swir2).divide(nir.add(swir2).add(EPS)).clamp(-1, 1).rename("NBR");
    var savi = nir.subtract(red).divide(nir.add(red).add(0.5).add(EPS)).multiply(1.5).clamp(-1.5, 1.5).rename("SAVI");
    var evi = nir.subtract(red).divide(nir.add(red.multiply(6)).subtract(blue.multiply(7.5)).add(1).add(EPS)).multiply(2.5).clamp(-1.5, 1.5).rename("EVI");

    var indices = ee.Image.cat([ndvi, ndmi, ndwi, nbr, savi, evi]);
    var validMask = red.gte(0).and(nir.gte(0)).and(blue.gte(0)).and(green.gte(0)).and(swir1.gte(0)).and(swir2.gte(0));
    return indices.updateMask(validMask).set("system:time_start", img.get("system:time_start"));
  });

var yearsList = ee.List.sequence(2010, 2025);
var indexNames = ["NDVI", "NDMI", "NDWI", "NBR", "SAVI", "EVI"];

var annualIndices = ee.ImageCollection(yearsList.map(function(y) {
  y = ee.Number(y);
  var start = ee.Date.fromYMD(y, 1, 1);
  var end = start.advance(1, "year");
  var img = col.filterDate(start, end).mean();
  return img.set("year", y).set("system:time_start", start.millis());
})).sort("year");

// ===================== STATISTIQUES PAR MASSIF (tous indices) - calcule en premier =====================
var rowsNested = yearsList.map(function(y) {
  y = ee.Number(y);
  var img = ee.Image(annualIndices.filter(ee.Filter.eq("year", y)).first());
  var reduced = img.reduceRegions({
    collection: massifs,
    reducer: ee.Reducer.mean(),
    scale: 500,
    tileScale: 4
  }).map(function(f) {
    return f.set("year", y);
  });
  return reduced.toList(reduced.size());
});
var perMassifFeatures = ee.FeatureCollection(rowsNested.flatten());

Export.table.toDrive({
  collection: perMassifFeatures,
  description: "Indices_par_massif_2010_2025",
  folder: "GEE_exports",
  fileNamePrefix: "Indices_par_massif_2010_2025",
  fileFormat: "CSV",
  selectors: ["nom_massif", "gid", "year", "NDVI", "NDMI", "NDWI", "NBR", "SAVI", "EVI"]
});

// ===================== STATISTIQUES GLOBALES (toute la zone d'etude) =====================
var globalStats = annualIndices.map(function(img) {
  var stat = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: massifs.geometry(),
    scale: 500,
    maxPixels: 1e9,
    tileScale: 4,
    bestEffort: true
  });
  return ee.Feature(null, {
    year: img.get("year"),
    NDVI: stat.get("NDVI"),
    NDMI: stat.get("NDMI"),
    NDWI: stat.get("NDWI"),
    NBR: stat.get("NBR"),
    SAVI: stat.get("SAVI"),
    EVI: stat.get("EVI")
  });
});

Export.table.toDrive({
  collection: globalStats,
  description: "Indices_global_2010_2025",
  folder: "GEE_exports",
  fileNamePrefix: "Indices_global_2010_2025",
  fileFormat: "CSV"
});

// ===================== Affichage (apres avoir soumis les exports) =====================
print("Statistiques globales des indices 2010-2025 (zone d'etude complete) :", globalStats);

var chartGlobal = ui.Chart.feature.byFeature(globalStats, "year", indexNames)
  .setChartType("LineChart")
  .setOptions({
    title: "Evolution des indices spectraux (2010-2025) - Zone d'etude complete",
    hAxis: {title: "Annee", format: "####"},
    vAxis: {title: "Valeur de l'indice"},
    lineWidth: 2,
    pointSize: 3,
    width: 900,
    height: 500
  });
print(chartGlobal);

print("Indices par massif 2010-2025 (echantillon) :", perMassifFeatures.limit(10));


// ============ Comparaison entre massifs (moyenne 2010-2025) ============
var meanIndicesImage = annualIndices.mean();
var perMassifMean = meanIndicesImage.reduceRegions({
  collection: massifs,
  reducer: ee.Reducer.mean(),
  scale: 500,
  tileScale: 4
});

Export.table.toDrive({
  collection: perMassifMean,
  description: 'Indices_moyenne_par_massif_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'Indices_moyenne_par_massif_2010_2025',
  fileFormat: 'CSV',
  selectors: ['nom_massif','gid','NDVI','NDMI','NDWI','NBR','SAVI','EVI']
});

var chartByMassif = ui.Chart.feature.byFeature(perMassifMean, 'nom_massif', indexNames)
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Comparaison des indices spectraux moyens par massif (2010-2025)',
    hAxis: {title: 'Massif', slantedText: true, slantedTextAngle: 45},
    vAxis: {title: "Valeur moyenne de l'indice"},
    width: 1200,
    height: 600
  });
print('Comparaison par massif (moyenne 2010-2025) :', perMassifMean);
print(chartByMassif);