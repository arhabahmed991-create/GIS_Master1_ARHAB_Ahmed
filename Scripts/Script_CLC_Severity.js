var massifs = ee.FeatureCollection('projects/hale-dynamo-473013-s7/assets/MASSIFS_FORESTIERS_S_013');
var geom = massifs.geometry();

// ============ Cross dNBR (Landsat 30m, 2010-2025) x CLC2018 vegetation classes ============
// "Vegetation" = CLC2018 codes 311-324 (CLC Level-2 groups 3.1 "Forests" and
// 3.2 "Shrub and/or herbaceous vegetation associations"). Excludes 331-335 (sparse/bare
// vegetation, rocks, beaches), agricultural and artificial classes.
var vegCodes = [311, 312, 313, 321, 322, 323, 324];

var clc = ee.Image('COPERNICUS/CORINE/V20/100m/2018').select('landcover');
var vegClass = clc.remap(vegCodes, [1, 2, 3, 4, 5, 6, 7], 0).rename('veg').selfMask().clip(geom);

var EPS = 1e-6;
function maskCloudsShadows(image) {
  var qa = image.select('QA_PIXEL');
  var cloudBitMask = 1 << 3;
  var shadowBitMask = 1 << 4;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(qa.bitwiseAnd(shadowBitMask).eq(0));
  return image.updateMask(mask);
}
function applyScale(image) {
  var optical = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true);
}
function nbrFromBands(image, nirBand, swir2Band) {
  image = applyScale(maskCloudsShadows(image));
  var nir = image.select(nirBand);
  var swir2 = image.select(swir2Band);
  var nbr = nir.subtract(swir2).divide(nir.add(swir2).add(EPS)).clamp(-1, 1).rename('NBR');
  var validMask = nir.gte(0).and(swir2.gte(0));
  return nbr.updateMask(validMask).copyProperties(image, ['system:time_start']);
}
function nbrL57(image) { return nbrFromBands(image, 'SR_B4', 'SR_B7'); }
function nbrL89(image) { return nbrFromBands(image, 'SR_B5', 'SR_B7'); }

var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterBounds(geom).map(nbrL57);
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterBounds(geom).map(nbrL57);
var l8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterBounds(geom).map(nbrL89);
var l9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterBounds(geom).map(nbrL89);
var nbrCol = ee.ImageCollection(l5.merge(l7).merge(l8).merge(l9)).select('NBR');

var years = ee.List.sequence(2010, 2025).getInfo();
var annualNBR = {};
years.forEach(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end = start.advance(1, 'year');
  annualNBR[y] = nbrCol.filterDate(start, end).median().rename('NBR');
});

var SEUILS = [0.10, 0.27, 0.44, 0.66];

for (var i = 0; i < years.length - 1; i++) {
  var y1 = years[i], y2 = years[i + 1];
  var bandTag = y1 + '_' + y2;
  var dnbr = annualNBR[y1].subtract(annualNBR[y2]).rename('dNBR');

  var severity = ee.Image(0)
    .where(dnbr.gt(SEUILS[0]), 1).where(dnbr.gt(SEUILS[1]), 2)
    .where(dnbr.gt(SEUILS[2]), 3).where(dnbr.gt(SEUILS[3]), 4)
    .updateMask(dnbr.mask()).rename('sev');

  // Only vegetation pixels AND severity 1-4 (i.e. actually affected, dNBR > 0.10)
  var typeCode = vegClass.multiply(10).add(severity).rename('type')
    .updateMask(vegClass.mask().and(severity.gte(1)));

  var areaHa = ee.Image.pixelArea().divide(10000).rename('area').updateMask(typeCode.mask());

  var combined = ee.Image.cat([areaHa, typeCode]);
  var reducer = ee.Reducer.sum().group({groupField: 1, groupName: 'type'});
  var grouped = combined.reduceRegion({
    reducer: reducer, geometry: geom, scale: 30, maxPixels: 1e10, tileScale: 8
  });

  var groupsList = ee.List(grouped.get('groups'));
  var features = groupsList.map(function(g) {
    g = ee.Dictionary(g);
    var t = ee.Number(g.get('type'));
    var veg = t.divide(10).floor();
    var sev = t.mod(10);
    return ee.Feature(null, {
      period: bandTag,
      veg_code: veg,
      sev_code: sev,
      area_ha: g.get('sum')
    });
  });
  var fc = ee.FeatureCollection(features);

  Export.table.toDrive({
    collection: fc,
    description: 'CrossTab_' + bandTag,
    folder: 'GEE_exports',
    fileNamePrefix: 'CrossTab_' + bandTag,
    fileFormat: 'CSV'
  });
}

print('Cross-tab exports queued: ' + (years.length - 1) + ' (one per dNBR period, 2010-2011 to 2024-2025)');
