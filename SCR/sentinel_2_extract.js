// Author: Zander Venter - zander.venter@nina.no

// This code is part of a workflow to classify land cover over Europe
// For final ELC10 map and link to manuscript see: https://doi.org/10.5281/zenodo.4407051

// Workflow in this script:
  // 1. Import reference locations from cleaned LUCAS dataset
  // 2. Collect Sentinel 2 spectral-temporal metrics in an image stack
  // 3. Extract values for S2 metrics at the reference locations and export as CSV


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

// Define data mask
var mask = lc2018.gt(0);

// Import sampling grid to stratify data extraction and export
  // Dataset from here: https://www.eea.europa.eu/data-and-maps/data/eea-reference-grids
var grid = ee.FeatureCollection("projects/nina/MetaComNet/grid_100km_select");
Map.addLayer(grid, {}, 'EAS grid', 0)

// Alternatively....
  // Generate a custom sampling grid to iterate data export over
  // a custom grid allows you to control the number of export pieces
  // if you have a very large reference dataset to extract satellite data over
  // you may want to break the export up into smaller batches.
var samplingGrid = ee.Image.random().multiply(10000000).toInt32()
    .reduceToVectors({
      reducer: ee.Reducer.countEvery(),
      geometry: aoi,
      geometryType: 'bb' ,
      eightConnected: false,
      scale: 250000, // reduce this to get more grid cells to iterate over
      crs: 'EPSG:4326'
    });
samplingGrid = samplingGrid.filterBounds(grid);
Map.addLayer(samplingGrid.style({color:"00FF00", fillColor:"FF000000"}), {}, 'samplingGrid', 0);

// Define Sentinel 2 collection - here I use TOA
var s2Col = 'S2' // or 'S2_SR' for surface reflectance

// Define the percentage of cloud cover below which you want to include
var sceneCloudThreshold = 60;

// Define the pixel cloud mask probability threshold
  // pixels above this cloud probability will be masked
var cloudMaskProbability = 40;

// Define time period - year of 2018 to match LUCAS sampling date
var startDate = '2018-01-01';
var endDate = '2019-01-01';

// Define S2 band common names
var S2_BANDS = ['QA60', 'B1','B2','B3','B4', 'B5', 'B6', 'B7','B8','B11','B12']; // Sentinel bands
var S2_NAMES = ['QA60','cb', 'blue', 'green', 'red', 'R1', 'R2', 'R3','nir','swir1', 'swir2']; // Common names


/*
  // Sentinel processing functions ///////////////////////////////////////////////////////////////////////////
*/

// Function to add spectral indices to Sentinel images
var addIndices = function(image) {
  var ndbi = image.expression(
    '(SWIR - NIR) / (SWIR + NIR)', {
      'SWIR': image.select('swir1'),
      'NIR': image.select('nir'),
    }).rename('ndbi');
  // Add vegetation indices
  var ndvi = image.normalizedDifference(['nir', 'red']).rename('ndvi');
  var nbr = image.normalizedDifference(['nir', 'swir2']).rename("nbr");
  var ndsi = image.normalizedDifference(['green', 'swir1']).rename("ndsi");
  return image.addBands(ndvi).addBands(ndbi).addBands(nbr).addBands(ndsi)
};


//This procedure must be used for proper processing of S2 imagery
function uniqueValues(collection,field){
    var values  =ee.Dictionary(collection.reduceColumns(ee.Reducer.frequencyHistogram(),[field]).get('histogram')).keys();
    return values;
  }
function dailyMosaics(imgs){
  //Simplify date to exclude time of day
  imgs = imgs.map(function(img){
  var d = ee.Date(img.get('system:time_start'));
  var day = d.get('day');
  var m = d.get('month');
  var y = d.get('year');
  var simpleDate = ee.Date.fromYMD(y,m,day);
  return img.set('simpleTime',simpleDate.millis());
  });
  
  //Find the unique days
  var days = uniqueValues(imgs,'simpleTime');
  
  imgs = days.map(function(d){
    d = ee.Number.parse(d);
    d = ee.Date(d);
    var t = imgs.filterDate(d,d.advance(1,'day'));
    var f = ee.Image(t.first());
    t = t.mosaic();
    t = t.set('system:time_start',d.millis());
    t = t.copyProperties(f);
    return t;
    });
    imgs = ee.ImageCollection.fromImages(imgs);
    
    return imgs;
}

var getS2_SR_CLOUD_PROBABILITY = function (aoi, startDate, endDate) { 
  var primary = ee.ImageCollection("COPERNICUS/" + s2Col)
      .filterBounds(aoi)
      .filterDate(startDate, endDate)
      .filterMetadata('CLOUDY_PIXEL_PERCENTAGE', 'less_than', sceneCloudThreshold)
      .select(S2_BANDS, S2_NAMES)
      .map(addIndices);
  var secondary = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
      .filterBounds(aoi)
      .filterDate(startDate, endDate);
  var innerJoined = ee.Join.inner().apply({
    primary: primary,
    secondary: secondary,
    condition: ee.Filter.equals({
      leftField: 'system:index',
      rightField: 'system:index'
    })
  });
  var mergeImageBands = function (joinResult) {
    return ee.Image(joinResult.get('primary'))
          .addBands(joinResult.get('secondary'));
  };
  var newCollection = innerJoined.map(mergeImageBands);
  return ee.ImageCollection(newCollection)
      .map(maskClouds(cloudMaskProbability))
      .sort('system:time_start');
};

var maskClouds = function(cloudProbabilityThreshold){
  return function(_img) {
  var cloudMask = _img.select('probability').lt(cloudProbabilityThreshold);
  return _img.updateMask(cloudMask);
}};


function getSentStack(aoi, startDate, endDate){
  var s2Combo = getS2_SR_CLOUD_PROBABILITY(aoi, startDate, endDate);
  
  var s2Cleaned =s2Combo;
  
  s2Cleaned = dailyMosaics(s2Cleaned); 
  
  var s2Median = s2Cleaned.select(['blue', 'green', 'red', 'R1', 'R2', 'R3','nir','swir1', 'swir2'])
    .reduce('median', 4);
  
  var s2Percs = s2Cleaned.select(['ndvi','ndsi'])
    .reduce(ee.Reducer.percentile([5,25,50,75,95]), 4);
  var stDev = s2Cleaned.select(['nbr'])
    .reduce(ee.Reducer.stdDev(), 4);
    
  var ndviSummer = s2Cleaned.select('ndvi').filter(ee.Filter.calendarRange(6,8, 'month')).median().rename('ndvi_summer');
  var ndviWinter = s2Cleaned.select('ndvi').filter(ee.Filter.calendarRange(12,2, 'month')).median().rename('ndvi_winter');
  var ndviSpring = s2Cleaned.select('ndvi').filter(ee.Filter.calendarRange(9,11, 'month')).median().rename('ndvi_spring');
  var ndviFall = s2Cleaned.select('ndvi').filter(ee.Filter.calendarRange(3,5, 'month')).median().rename('ndvi_fall');
  
  var ndviFocal = s2Cleaned.select('ndvi').reduce('median', 4)
    .reduceNeighborhood(ee.Reducer.stdDev(), ee.Kernel.square(3, 'pixels')).rename('ndvi_texture_sd');
  
  var s2Stack = s2Median
    .addBands(s2Percs)
    .addBands(stDev)
    .addBands(ndviSummer)
    .addBands(ndviWinter)
    .addBands(ndviSpring)
    .addBands(ndviFall)
    .addBands(ndviFocal);
    
  // Multiply by 1000 and convert to integer to reduce file sizes on export
  s2Stack = s2Stack.multiply(1000).round().int();
  //print(s2Stack, 'sentinel stack');

  return s2Stack
}


/*
  // Exporting spectral-temporal metrics ///////////////////////////////////////////////////////////////////////////
*/

// Create a list of grid cells to iterate over
var list= samplingGrid.toList(1000);// here you can use your custom grid, or the "grid" object from L38
var length = list.length().getInfo();
print('length of grid cells for iteration: ', length)

// Only iterating 2 for testing purposes
for (var i = 0; i<2; i++){
  
  // Convert the grid cell into a geometry AOI
  var geo = ee.Feature(list.get(i)).geometry();
  
  // Get the Sentinel 2 stack for this geometry
  var masterStack = getSentStack(geo, startDate, endDate);
  
  // Filter the LUCAS features to this geometry
  var featsSelect = feats.filterBounds(geo)
  
  // Sample the image stack at the LUCAS points
  var table = masterStack.sampleRegions({
    collection: featsSelect,
    scale: 10,
    tileScale: 4
  });
  
  // Submit an export task to your Google Drive
    // NB: you must create a Folder in your Drive named "S2_data"
  Export.table.toDrive({
    collection: table,
    folder: 'S2_data',
    description: s2Col + '_' + String(i),
    fileFormat: 'CSV'
  });
  
}

/*

// After running the script, you will have a whole lot of export tasks ready to run
// Instead of clicking on each one individually, you can:
  // Copy the two functions below
  // Click F12 on your computer
  // In the "Console" tab paste the two functions and hit Enter
  // Then type 'runTaskList()' in the Console and hit Enter
  // Wait a few seconds - roughly 30 seconds wait per 100 pts
  // Then type 'confirmAll()' and hit Enter
  // All the tasks should start running at once.

function runTaskList(){
    var tasklist = document.getElementsByClassName('task local type-EXPORT_FEATURES awaiting-user-config');
    for (var i = 0; i < tasklist.length; i++)
            tasklist[i].getElementsByClassName('run-button')[0].click();
}

function confirmAll() {
    var ok = document.getElementsByClassName('goog-buttonset-default goog-buttonset-action');
    for (var i = 0; i < ok.length; i++)
        ok[i].click();
}

*/