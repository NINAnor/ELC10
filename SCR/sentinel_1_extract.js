// Author: Zander Venter - zander.venter@nina.no

// This code is part of a workflow to classify land cover over Europe
// For final ELC10 map and link to manuscript see: https://doi.org/10.5281/zenodo.4407051

// Workflow in this script:
  // 1. Import reference locations from cleaned LUCAS dataset
  // 2. Collect Sentinel 1 backscatter temporal metrics in an image stack
  // 3. Extract values for S1 metrics at the reference locations and export as CSV


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

// Define whether to perform speckle filtering or not
  // Speckle filter is computationally intensive so here I chose not to do it
  // for efficiency sake. There is also very little increase in classification accuracy
var speckleFilter = false;

// Define whether to perform radiometric slope correction
var slopeCorrection = true;

// Define elevation model to use in radiometric slope correction
var glob = ee.Image("USGS/GTOPO30").rename('elevation')
var srtm = ee.Image('USGS/SRTMGL1_003').rename('elevation')
var elev = srtm.unmask(glob).unmask(0);

// Define time period - year of 2018 to match LUCAS sampling date
var startDate = '2018-01-01';
var endDate = '2019-01-01';


/*
  // Sentinel 1 processing functions ///////////////////////////////////////////////////////////////////////////////
*/

// Mask extreme angels
function maskAngle(img) {
  var ang = img.select(['angle']);
  var first = img.updateMask(ang.gt(30.63993));
  return first.updateMask(ang.lt(44.73993)).select(0).copyProperties(img, ['system:time_start']);
}
// Clip edges of scenes
function maskEdge(img) {
  var mask = img.select(0).unitScale(-25, 5).multiply(255).toByte().connectedComponents(ee.Kernel.rectangle(1,1), 100);
  return img.updateMask(mask.select(0)).copyProperties(img, ['system:time_start']);  
}
//Clip edges
function clipEdge (img) {
          return img.clip(img.geometry().buffer(-10000)).copyProperties(img, ['system:time_start']);
  }
// Convert to natural 
function toNatural(img) {
  img = ee.Image(img)
  return ee.Image(10.0).pow(img.select(0).divide(10.0)).copyProperties(img, ['system:time_start']);
}
// Convert to DB
function toDB(img) {
  return ee.Image(img).log10().multiply(10.0).copyProperties(img, ['system:time_start']);
}

// Sigma Lee filter 
// The RL speckle filter from https://code.earthengine.google.com/2ef38463ebaf5ae133a478f173fd0ab5
// by Guido Lemoine
function RefinedLee(img) {
  // img must be in natural units, i.e. not in dB!
  // Set up 3x3 kernels 
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);

  var mean3 = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  // Use a sample of the 3x3 windows inside a 7x7 windows to determine gradients and directions
  var sample_weights = ee.List([[0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0], [0,1,0,1,0,1,0], [0,0,0,0,0,0,0], [0,1,0,1,0,1,0],[0,0,0,0,0,0,0]]);

  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);

  // Calculate mean and variance for the sampled windows and store as 9 bands
  var sample_mean = mean3.neighborhoodToBands(sample_kernel); 
  var sample_var = variance3.neighborhoodToBands(sample_kernel);

  // Determine the 4 gradients for the sampled windows
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs();
  gradients = gradients.addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs());
  gradients = gradients.addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs());
  gradients = gradients.addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());

  // And find the maximum gradient amongst gradient bands
  var max_gradient = gradients.reduce(ee.Reducer.max());

  // Create a mask for band pixels that are the maximum gradient
  var gradmask = gradients.eq(max_gradient);

  // duplicate gradmask bands: each gradient represents 2 directions
  gradmask = gradmask.addBands(gradmask);

  // Determine the 8 directions
  var directions = sample_mean.select(1).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1);
  directions = directions.addBands(sample_mean.select(6).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2));
  directions = directions.addBands(sample_mean.select(3).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3));
  directions = directions.addBands(sample_mean.select(0).subtract(sample_mean.select(4)).gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4));
  // The next 4 are the not() of the previous 4
  directions = directions.addBands(directions.select(0).not().multiply(5));
  directions = directions.addBands(directions.select(1).not().multiply(6));
  directions = directions.addBands(directions.select(2).not().multiply(7));
  directions = directions.addBands(directions.select(3).not().multiply(8));

  // Mask all values that are not 1-8
  directions = directions.updateMask(gradmask);

  // "collapse" the stack into a singe band image (due to masking, each pixel has just one value (1-8) in it's directional band, and is otherwise masked)
  directions = directions.reduce(ee.Reducer.sum());  

  //var pal = ['ffffff','ff0000','ffff00', '00ff00', '00ffff', '0000ff', 'ff00ff', '000000'];
  //Map.addLayer(directions.reduce(ee.Reducer.sum()), {min:1, max:8, palette: pal}, 'Directions', false);

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));

  // Calculate localNoiseVariance
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  // Set up the 7*7 kernels for directional statistics
  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));

  var diag_weights = ee.List([[1,0,0,0,0,0,0], [1,1,0,0,0,0,0], [1,1,1,0,0,0,0], 
    [1,1,1,1,0,0,0], [1,1,1,1,1,0,0], [1,1,1,1,1,1,0], [1,1,1,1,1,1,1]]);

  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);

  // Create stacks for mean and variance using the original kernels. Mask with relevant direction.
  var dir_mean = img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var = img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));

  dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));

  // and add the bands for rotated kernels
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(img.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var = dir_var.addBands(img.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }

  // "collapse" the stack into a single band image (due to masking, each pixel has just one value in it's directional band, and is otherwise masked)
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var = dir_var.reduce(ee.Reducer.sum());

  // A finally generate the filtered value
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));

  var b = varX.divide(dir_var);

  var result = dir_mean.add(b.multiply(img.subtract(dir_mean)));
  return ee.Image(result.arrayFlatten([['sum']])).copyProperties(img, ['system:time_start']);
}

// Function to wrap all the above into one
var cleanS1 = function(img){
  img = maskAngle(img)
  img = ee.Image(toNatural(img))
  img = ee.Image(RefinedLee(img))
  img = ee.Image(toDB(img))
  return img
}

// Function for radiometric slope correction
// Sentinel-1 image collection
// https://github.com/ESA-PhiLab/radiometric-slope-correction/tree/master/javascript
var slope_correction = function (collection, options){

    // set defaults if undefined options
    options = options || {};
    var model = options.model || 'volume';
    var elevation = options.elevation || ee.Image('USGS/SRTMGL1_003');
    var buffer = options.buffer || 0;

    // we need a 90 degree in radians image for a couple of calculations
    var ninetyRad = ee.Image.constant(90).multiply(Math.PI/180);

    // Volumetric Model Hoekman 1990
    function _volume_model(theta_iRad, alpha_rRad){

      var nominator = (ninetyRad.subtract(theta_iRad).add(alpha_rRad)).tan();
      var denominator = (ninetyRad.subtract(theta_iRad)).tan();
      return nominator.divide(denominator);
    }

    // surface model Ulander et al. 1996
    function _surface_model(theta_iRad, alpha_rRad, alpha_azRad){

      var nominator = (ninetyRad.subtract(theta_iRad)).cos();
      var denominator = alpha_azRad.cos()
        .multiply((ninetyRad.subtract(theta_iRad).add(alpha_rRad)).cos());
      return nominator.divide(denominator);
    }

    // buffer function (thanks Noel)
    function _erode(img, distance) {

      var d = (img.not().unmask(1)
          .fastDistanceTransform(30).sqrt()
          .multiply(ee.Image.pixelArea().sqrt()));

      return img.updateMask(d.gt(distance));
    }

    // calculate masks
    function _masking(alpha_rRad, theta_iRad, proj, buffer){

        // layover, where slope > radar viewing angle
        var layover = alpha_rRad.lt(theta_iRad).rename('layover');

        // shadow
        var shadow = alpha_rRad.gt(ee.Image.constant(-1).multiply(ninetyRad.subtract(theta_iRad))).rename('shadow');

        // combine layover and shadow
        var mask = layover.and(shadow);

        // add buffer to final mask
        if (buffer > 0)
            mask = _erode(mask, buffer);

        return mask.rename('no_data_mask');
   }

    function _correct(image){

        // get image geometry and projection
        var geom = image.geometry();
        var proj = image.select(1).projection();

        // get look direction angle
        var heading = (ee.Terrain.aspect(
            image.select('angle')).reduceRegion(ee.Reducer.mean(), geom, 1000).get('aspect')
            );

        // Sigma0 to Power of input image
        var sigma0Pow = ee.Image.constant(10).pow(image.divide(10.0));

        // Radar geometry
        var theta_iRad = image.select('angle').multiply(Math.PI/180).clip(geom);
        var phi_iRad = ee.Image.constant(heading).multiply(Math.PI/180);

        // Terrain geometry
        var alpha_sRad = ee.Terrain.slope(elevation).select('slope')
            .multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);
        var phi_sRad = ee.Terrain.aspect(elevation).select('aspect')
            .multiply(Math.PI/180).setDefaultProjection(proj).clip(geom);

        // Model geometry

        //reduce to 3 angle
        var phi_rRad = phi_iRad.subtract(phi_sRad);

        // slope steepness in range
        var alpha_rRad = (alpha_sRad.tan().multiply(phi_rRad.cos())).atan();

        // slope steepness in azimuth
        var alpha_azRad = (alpha_sRad.tan().multiply(phi_rRad.sin())).atan();

        // Gamma_nought
        var gamma0 = sigma0Pow .divide(theta_iRad.cos());

               // models
        if (model == 'volume')
          var corrModel = _volume_model(theta_iRad, alpha_rRad);

        if (model == 'surface')
          var corrModel = _surface_model(theta_iRad, alpha_rRad, alpha_azRad);

        if (model == 'direct')
          var corrModel = _direct_model(theta_iRad, alpha_rRad, alpha_azRad);

        // apply model to derive gamma0_flat
        var gamma0_flat = gamma0.divide(corrModel);

        // transform to dB-scale
        var gamma0_flatDB = (ee.Image.constant(10)
            .multiply(gamma0_flat.log10()).select(['VV', 'VH'])
            );

        // get Layover/Shadow mask
        var mask = _masking(alpha_rRad, theta_iRad, proj, buffer);

        // return gamma_flat plus mask
        return gamma0_flatDB.addBands(mask).addBands(image.select('angle')).copyProperties(image);


    }

    // run correction function and return corrected collection
    return collection.map(_correct);

};


// Function to generate an S1 image stack for AOI 
function getSARstack(aoi){
  
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi)
    .filterDate(startDate, endDate)
    .filterMetadata('resolution_meters', 'equals' , 10)
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.eq('instrumentMode', 'IW'));
  
  if (slopeCorrection){
      s1 = slope_correction(
        s1,
        {'model': 'volume', // correction model - volume recommended for land cover classifiction
        'elevation': elev, // elevation model
        'buffer': 0      // buffer in meter
        }
      );
  }
  
  var asc = s1.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
  var desc = s1.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
  
  var asc_vv = asc
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select(['VV', 'angle']);
  var asc_vh = asc
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VH', 'angle']);
  
  var desc_vv = desc
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .select(['VV', 'angle']);
  var desc_vh = desc
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VH', 'angle']);
  
  if (speckleFilter){
    
    asc_vv = asc_vv.map(cleanS1);
    asc_vh = asc_vh.map(cleanS1);
    desc_vv = desc_vv.map(cleanS1);
    desc_vh = desc_vh.map(cleanS1);
    
  }
  
  var stack = asc_vh.select(0).median().rename('asc_vh_median')
    .addBands(asc_vv.select(0).median().rename('asc_vv_median'))
    .addBands(desc_vv.select(0).median().rename('desc_vv_median'))
    .addBands(desc_vh.select(0).median().rename('desc_vh_median'))
    .addBands(desc_vh.select(0).median().divide(desc_vv.select(0).median()).rename('desc_dpol_median'))
    .addBands(asc_vh.select(0).median().divide(asc_vv.select(0).median()).rename('asc_dpol_median'))
    .addBands(asc_vv.select(0).reduce(ee.Reducer.stdDev()).rename('asc_vv_stDev'))
    .addBands(asc_vh.select(0).reduce(ee.Reducer.stdDev()).rename('asc_vh_stDev'))
    .addBands(desc_vv.select(0).reduce(ee.Reducer.stdDev()).rename('desc_vv_stDev'))
    .addBands(desc_vh.select(0).reduce(ee.Reducer.stdDev()).rename('desc_vh_stDev'))
  
  // Multiply by 1000 and convert to integer to reduce file sizes on export
  stack = stack.multiply(1000).round().int()
  
  return stack
}


/*
  // Exporting backscatter temporal metrics ///////////////////////////////////////////////////////////////////////////
*/

// Create a list of grid cells to iterate over
var list= samplingGrid.toList(1000); // here you can use your custom grid, or the "grid" object from L38
var length = list.length().getInfo();
print('length of grid cells for iteration: ', length)

// Only iterating 2 for testing purposes
for (var i = 0; i<2; i++){
  
  // Convert the grid cell into a geometry AOI
  var geo = ee.Feature(list.get(i)).geometry();
  
  // Get the Sentinel 2 stack for this geometry
  var masterStack = getSARstack(geo);
  
  // Filter the LUCAS features to this geometry
  var featsSelect = feats.filterBounds(geo)
  
  // Submit an export task to your Google Drive
    // NB: you must create a Folder in your Drive named "S1_data"
  var table = masterStack.sampleRegions({
      collection: featsSelect,
      scale: 10,
      tileScale: 4
    });
    
  Export.table.toDrive({
    collection: table,
    folder: 'S1_data',
    description: 'S1_' + String(i),
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