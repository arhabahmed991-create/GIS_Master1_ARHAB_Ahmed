var massifs = ee.FeatureCollection('projects/hale-dynamo-473013-s7/assets/MASSIFS_FORESTIERS_S_013');

// MOD09A1 : reflectance de surface MODIS 500m, composites 8 jours.
// Calcul du NBR (Normalized Burn Ratio) pour chaque annee 2010-2025,
// puis du dNBR entre chaque paire d'annees consecutives (2010-2011, 2011-2012, ..., 2024-2025).
// dNBR(N -> N+1) = NBR(N) - NBR(N+1)
// dNBR positif = baisse du NBR = perte de vigueur vegetale / signal de feu ou de stress.
// dNBR negatif = hausse du NBR = reprise / croissance de la vegetation.
// NB : ceci est un dNBR "roulant" annee par annee sur toute la serie (proxy de variation
// interannuelle), different du dNBR classique USGS (pre-feu moins post-feu sur un evenement unique).
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

var yearsList = ee.List.sequence(2010, 2025);

var annualNBR = ee.ImageCollection(yearsList.map(function(y) {
  y = ee.Number(y);
  var start = ee.Date.fromYMD(y, 1, 1);
  var end = start.advance(1, 'year');
  var img = col.filterDate(start, end).mean();
  return img.set('year', y).set('system:time_start', start.millis());
})).sort('year');

// ============ NBR global annuel ============
var globalNBR = ee.FeatureCollection(yearsList.map(function(y) {
  y = ee.Number(y);
  var img = ee.Image(annualNBR.filter(ee.Filter.eq('year', y)).first());
  var stat = img.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: massifs.geometry(),
    scale: 500,
    tileScale: 4,
    bestEffort: true
  });
  return ee.Feature(null, {year: y, NBR: stat.get('NBR')});
}));

// ============ NBR par massif annuel ============
var rowsNested = yearsList.map(function(y) {
  y = ee.Number(y);
  var img = ee.Image(annualNBR.filter(ee.Filter.eq('year', y)).first());
  var reduced = img.reduceRegions({
    collection: massifs,
    reducer: ee.Reducer.mean(),
    scale: 500,
    tileScale: 4
  }).map(function(f) {
    return f.set('year', y);
  });
  return reduced.toList(reduced.size());
});
var perMassifNBR = ee.FeatureCollection(rowsNested.flatten());

// ============ dNBR global (paires d'annees consecutives) ============
var yearPairs = ee.List.sequence(2010, 2024);
var globalDNBR = ee.FeatureCollection(yearPairs.map(function(y1) {
  y1 = ee.Number(y1);
  var y2 = y1.add(1);
  var f1 = ee.Feature(globalNBR.filter(ee.Filter.eq('year', y1)).first());
  var f2 = ee.Feature(globalNBR.filter(ee.Filter.eq('year', y2)).first());
  var nbr1 = ee.Number(f1.get('NBR'));
  var nbr2 = ee.Number(f2.get('NBR'));
  var label = ee.String(y1.format('%d')).cat('-').cat(y2.format('%d'));
  return ee.Feature(null, {
    periode: label,
    year_debut: y1,
    NBR_debut: nbr1,
    NBR_fin: nbr2,
    dNBR: nbr1.subtract(nbr2)
  });
}));

// ============ dNBR par massif (paires d'annees consecutives) ============
// Utilisation de dictionnaires (nom_massif -> NBR) par annee plutot qu'un Join,
// et d'une valeur sentinelle pour eviter toute erreur si une valeur venait a manquer
// pour un massif / une annee donnee (pixels masques).
var massifNames = massifs.aggregate_array('nom_massif');
var massifGids = massifs.aggregate_array('gid');
var SENTINEL = -9999;

var yearDictList = yearsList.map(function(y) {
  y = ee.Number(y);
  var fc = perMassifNBR.filter(ee.Filter.eq('year', y));
  var names = fc.aggregate_array('nom_massif');
  var vals = fc.aggregate_array('NBR');
  var n = names.size();
  var m = vals.size();
  return ee.Algorithms.If(n.eq(m), ee.Dictionary.fromLists(names, vals), ee.Dictionary.fromLists(names.slice(0, m), vals));
});

var pairIndices = ee.List.sequence(0, 14);
var rowsDNBRNested = pairIndices.map(function(idx) {
  idx = ee.Number(idx);
  var y1 = ee.Number(yearsList.get(idx));
  var y2 = ee.Number(yearsList.get(idx.add(1)));
  var dict1 = ee.Dictionary(yearDictList.get(idx));
  var dict2 = ee.Dictionary(yearDictList.get(idx.add(1)));
  var label = ee.String(y1.format('%d')).cat('-').cat(y2.format('%d'));
  var rows = massifNames.zip(massifGids).map(function(pair) {
    pair = ee.List(pair);
    var name = pair.get(0);
    var gid = pair.get(1);
    var nbr1 = ee.Number(dict1.get(name, SENTINEL));
    var nbr2 = ee.Number(dict2.get(name, SENTINEL));
    var dnbrRaw = nbr1.subtract(nbr2);
    var valid1 = nbr1.neq(SENTINEL);
    var valid2 = nbr2.neq(SENTINEL);
    var bothValid = valid1.and(valid2);
    return ee.Feature(null, {
      nom_massif: name,
      gid: gid,
      periode: label,
      NBR_debut: ee.Algorithms.If(valid1, nbr1, null),
      NBR_fin: ee.Algorithms.If(valid2, nbr2, null),
      dNBR: ee.Algorithms.If(bothValid, dnbrRaw, null)
    });
  });
  return rows;
});
var perMassifDNBR = ee.FeatureCollection(rowsDNBRNested.flatten());

// ============ Exports (avant les print/chart pour eviter les soucis de quota) ============
Export.table.toDrive({
  collection: perMassifNBR,
  description: 'NBR_annuel_par_massif_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'NBR_annuel_par_massif_2010_2025',
  fileFormat: 'CSV',
  selectors: ['nom_massif', 'gid', 'year', 'NBR']
});

Export.table.toDrive({
  collection: globalNBR,
  description: 'NBR_annuel_global_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'NBR_annuel_global_2010_2025',
  fileFormat: 'CSV',
  selectors: ['year', 'NBR']
});

Export.table.toDrive({
  collection: perMassifDNBR,
  description: 'dNBR_annuel_par_massif_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'dNBR_annuel_par_massif_2010_2025',
  fileFormat: 'CSV',
  selectors: ['nom_massif', 'gid', 'periode', 'NBR_debut', 'NBR_fin', 'dNBR']
});

Export.table.toDrive({
  collection: globalDNBR,
  description: 'dNBR_annuel_global_2010_2025',
  folder: 'GEE_exports',
  fileNamePrefix: 'dNBR_annuel_global_2010_2025',
  fileFormat: 'CSV',
  selectors: ['periode', 'year_debut', 'NBR_debut', 'NBR_fin', 'dNBR']
});

// ============ Affichage (apres avoir soumis les exports) ============
print('NBR global annuel 2010-2025 :', globalNBR);
var chartNBR = ui.Chart.feature.byFeature(globalNBR, 'year', ['NBR'])
  .setChartType('LineChart')
  .setOptions({
    title: 'Evolution du NBR global (2010-2025)',
    hAxis: {title: 'Annee', format: '####'},
    vAxis: {title: 'NBR'},
    lineWidth: 2,
    pointSize: 4,
    width: 900,
    height: 400
  });
print(chartNBR);

print('dNBR global par periode (N -> N+1) :', globalDNBR);
var chartDNBR = ui.Chart.feature.byFeature(globalDNBR, 'periode', ['dNBR'])
  .setChartType('ColumnChart')
  .setOptions({
    title: 'dNBR annuel global (periode N -> N+1)',
    hAxis: {title: 'Periode', slantedText: true, slantedTextAngle: 45},
    vAxis: {title: 'dNBR'},
    width: 1000,
    height: 400
  });
print(chartDNBR);

print('dNBR par massif (echantillon) :', perMassifDNBR.limit(10));
