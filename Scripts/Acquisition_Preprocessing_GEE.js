var massifs = ee.FeatureCollection('projects/hale-dynamo-473013-s7/assets/MASSIFS_FORESTIERS_S_013');
var geom = massifs.geometry();

// Surface de chaque massif (ha), utile pour calculer un pourcentage de surface touchee.
var massifsArea = massifs.map(function(f) {
  return f.set('surface_massif_ha', f.geometry().area(1).divide(10000));
});

// ============ NBR annuel par pixel (2010-2025), identique aux scripts precedents ============
var EPS = 1e-6;
var col = ee.ImageCollection('MODIS/061/MOD09A1')
  .filterDate('2010-01-01', '2026-01-01')
  .map(function(img) {
    var scale = 0.0001;
    var nir = img.select('sur_refl_b02').multiply(scale);
    var swir2 = img.select('sur_refl_b07').multiply(scale);
    var nbr = nir.subtract(swir2).divide(nir.add(swir2).add(EPS)).clamp(-1, 1).rename('NBR');
    var validMask = nir.gte(0).and(swir2.gte(0));
    return nbr.updateMask(validMask).set('system:time_start', img.get('system:time_start'));
  });

var years = ee.List.sequence(2010, 2025).getInfo();

var annualNBR = {};
years.forEach(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end = start.advance(1, 'year');
  annualNBR[y] = col.filterDate(start, end).mean().rename('NBR');
});

// ============ dNBR par pixel pour chaque periode, classification en classes de severite ============
// Seuils USGS classiques (echelle -1..1) :
//   <= 0.10           : pas de changement significatif / non brule
//   0.10 - 0.27       : severite faible
//   0.27 - 0.44       : severite moderee-faible
//   0.44 - 0.66       : severite moderee-elevee
//   > 0.66            : severite elevee
// NB : seuils calibres a l'origine pour Landsat en comparaison pre/post-feu sur un evenement unique.
// Ici, applique a un dNBR "roulant" annee par annee sur MODIS 500m : a interpreter comme un proxy
// de perte de vigueur vegetale (feu possible, mais aussi secheresse/stress hydrique severe).
var SEUIL_PERTURBATION = 0.1;
var SEUILS = [0.1, 0.27, 0.44, 0.66];

var classifiedBands = [];
var valueBands = [];
var periodLabels = [];
var perMassifBurnedArea = null;
var globalBurnedArea = null;

for (var i = 0; i < years.length - 1; i++) {
  var y1 = years[i];
  var y2 = years[i + 1];
  var bandTag = y1 + '_' + y2;
  var label = y1 + '-' + y2;
  periodLabels.push(label);

  var dnbr = annualNBR[y1].subtract(annualNBR[y2]).rename('dNBR_' + bandTag);
  valueBands.push(dnbr);

  var cls = ee.Image(0)
    .where(dnbr.gt(SEUILS[0]), 1)
    .where(dnbr.gt(SEUILS[1]), 2)
    .where(dnbr.gt(SEUILS[2]), 3)
    .where(dnbr.gt(SEUILS[3]), 4)
    .updateMask(dnbr.mask())
    .rename('classe_' + bandTag)
    .toByte();
  classifiedBands.push(cls);

  var burnedMask = dnbr.gt(SEUIL_PERTURBATION);
  var areaImg = ee.Image.pixelArea().divide(10000).updateMask(burnedMask).rename('surface_brulee_ha');

  var reducedMassif = areaImg.reduceRegions({
    collection: massifsArea,
    reducer: ee.Reducer.sum(),
    scale: 500,
    tileScale: 4
  }).map((function(lbl) {
    return function(f) {
      var surf = ee.Number(f.get('sum'));
      var totalHa = ee.Number(f.get('surface_massif_ha'));
      return f.set({
        periode: lbl,
        surface_brulee_ha: surf,
        pct_surface_massif: surf.divide(totalHa).multiply(100)
      });
    };
  })(label));
  perMassifBurnedArea = (perMassifBurnedArea === null) ? reducedMassif : perMassifBurnedArea.merge(reducedMassif);

  var statGlobal = areaImg.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom,
    scale: 500,
    tileScale: 4,
    bestEffort: true
  });
  var fGlobal = ee.Feature(null, {periode: label, year_debut: y1, surface_brulee_ha: statGlobal.get('surface_brulee_ha')});
  globalBurnedArea = (globalBurnedArea === null) ? ee.FeatureCollection([fGlobal]) : globalBurnedArea.merge(ee.FeatureCollection([fGlobal]));
}

var classifiedImage = ee.Image.cat(classifiedBands).clip(geom);
var valueImage = ee.Image.cat(valueBands).clip(geom);

// Carte de severite maximale jamais atteinte sur toute la periode 2010-2025 (par pixel)
var classifiedBandsRenamed = classifiedBands.map(function(img) { return img.rename('classe'); });
var severiteMax = ee.ImageCollection(classifiedBandsRenamed).max().rename('classe_max').clip(geom);
// Frequence de perturbation : nombre de periodes (sur 15) ou le pixel a depasse le seuil
var frequencePerturbation = ee.ImageCollection(classifiedBandsRenamed.map(function(c) { return c.gte(1); }))
  .sum().rename('freq_perturbation').toByte().clip(geom);

// ============ Exports (avant les print/chart pour eviter les soucis de quota) ============
Export.image.toDrive({
  image: classifiedImage,
  description: 'Cartes_classes_dNBR_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'Cartes_classes_dNBR_2010_2025',
  region: geom,
  scale: 500,
  maxPixels: 1e9,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: valueImage,
  description: 'dNBR_valeurs_raster_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'dNBR_valeurs_raster_2010_2025',
  region: geom,
  scale: 500,
  maxPixels: 1e9,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: severiteMax,
  description: 'Carte_severite_max_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'Carte_severite_max_2010_2025',
  region: geom,
  scale: 500,
  maxPixels: 1e9,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: frequencePerturbation,
  description: 'Frequence_perturbation_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'Frequence_perturbation_2010_2025',
  region: geom,
  scale: 500,
  maxPixels: 1e9,
  fileFormat: 'GeoTIFF'
});

Export.table.toDrive({
  collection: perMassifBurnedArea,
  description: 'Surface_brulee_par_massif_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'Surface_brulee_par_massif_2010_2025',
  fileFormat: 'CSV',
  selectors: ['nom_massif', 'gid', 'periode', 'surface_massif_ha', 'surface_brulee_ha', 'pct_surface_massif']
});

Export.table.toDrive({
  collection: globalBurnedArea,
  description: 'Surface_brulee_globale_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'Surface_brulee_globale_2010_2025',
  fileFormat: 'CSV',
  selectors: ['periode', 'year_debut', 'surface_brulee_ha']
});

// ============ Affichage (apres avoir soumis les exports) ============
var palette = ['ffffcc', 'fed976', 'fd8d3c', 'e31a1c', '800026'];
Map.centerObject(geom, 8);
Map.addLayer(severiteMax, {min: 0, max: 4, palette: palette}, 'Severite max 2010-2025');
Map.addLayer(massifs.style({color: '2222ff', fillColor: '00000000', width: 1}), {}, 'Massifs (contours)');

print('Surface brulee globale par periode :', globalBurnedArea);
var chartSurface = ui.Chart.feature.byFeature(globalBurnedArea, 'periode', ['surface_brulee_ha'])
  .setChartType('ColumnChart')
  .setOptions({
    title: 'Surface potentiellement brulee/perturbee (dNBR > 0.1) par periode - zone complete',
    hAxis: {title: 'Periode', slantedText: true, slantedTextAngle: 45},
    vAxis: {title: 'Surface (ha)'},
    width: 1000,
    height: 400
  });
print(chartSurface);

print('Surface brulee par massif (echantillon) :', perMassifBurnedArea.limit(10));
