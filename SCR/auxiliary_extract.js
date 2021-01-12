// Author: Zander Venter - zander.venter@nina.no

// This code is part of a workflow to classify land cover over Europe
// For final ELC10 map and link to manuscript see: https://doi.org/10.5281/zenodo.4407051

// Workflow in this script:
  // 1. Import reference locations from cleaned LUCAS dataset
  // 2. Collect auxiliary predictor variables from terrain, climate and nighttime light datasets
  // 3. Extract values for predictors at the reference locations and export as CSV


/*
  // Global Objects ///////////////////////////////////////////////////////////////////////////////
*/

// Import Pflugmacher 30m dataset to define data mask and AOI
  // Dataset from here: https://www.sciencedirect.com/science/article/pii/S0034425718305546
var lc2018 = ee.Image("users/zandersamuel/NINA/Raster/Europe_landcover_2015_RSE_LUCAS");

// Import cleaned LUCAS data points for extraction
  // a combination of Copernicus module from here: https://doi.org/10.5194/essd-2020-178
  // and the theoretical points here: https://www.nature.com/articles/s41597-020-00675-z
var feats = ee.FeatureCollection('projects/nina/ELC10/LUCAS_combined_reference_final');
Map.addLayer(feats, {}, 'features', 0)
print('FEATS length: ', feats.size())

// Define area of interest
var aoi = lc2018.geometry().bounds();
print(aoi, 'aoi')
Map.addLayer(aoi, {}, 'aoi', 0)


/*
  // Auxiliary datasets for extraction ////////////////////////////////////////////////////////////
*/

//// Night time lights
var lightCol = ee.ImageCollection("NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG")
  .filterDate('2018-01-01', '2019-01-01').select('avg_rad');
var light = lightCol.median().rename('light');
Map.addLayer(light, {min:0, max:10},'light', 0)
print(light, 'nighttime light');


///// Terrain
var glob = ee.Image("USGS/GTOPO30").rename('elevation')
var srtm = ee.Image('USGS/SRTMGL1_003').rename('elevation')
var elev = srtm.unmask(glob).unmask(0).rename('elevation');
Map.addLayer(elev, {min:0, max:3000},'elevation', 0)
print(elev, 'elevation');

///// Climate data 
var era5Col = ee.ImageCollection("ECMWF/ERA5/MONTHLY")
  .filterDate('2010-01-01', '2020-01-01')
  .select(['mean_2m_air_temperature', 'total_precipitation'], ['temp', 'precip']);
Map.addLayer(era5Col.select('precip'), {}, 'era5 col', 0)
var era5 = era5Col.mean();
Map.addLayer(era5.select('temp'), {min:250, max:320}, 'era5 temp', 0)
Map.addLayer(era5.select('precip'), {min:0, max:0.5}, 'era5 precip', 0)

var tempSD = era5Col.select('temp').reduce(ee.Reducer.stdDev()).rename('temp_stDev')

var precipSD = era5Col.select('precip').reduce(ee.Reducer.stdDev()).rename('precip_stDev')

var climate = era5
  .addBands(tempSD)
  .addBands(precipSD);
print(climate, 'climate');

/*
  // Extract and export auxiliary data ////////////////////////////////////////////////////////////
*/

///// Stack
var masterStack = elev
  .addBands(climate)
  .addBands(light);
// Multiply by 1000 and convert to integer to reduce file sizes on export
masterStack = masterStack.multiply(1000).round().int()
print(masterStack, 'masterStack');

var table = masterStack.reduceRegions({
  collection: feats,
  reducer: ee.Reducer.mean(),
  scale: 30,
  tileScale: 2
});
table = table.map(function(ft){return ft.setGeometry(null)})
print(table.limit(10))
  
Export.table.toDrive({
  collection: table,
  description: 'auxiliary_data',
  fileFormat: 'CSV'
});
