// Author: Zander Venter - zander.venter@nina.no

// This code is part of a workflow to classify land cover over Europe
// For final ELC10 map and link to manuscript see: https://doi.org/10.5281/zenodo.4407051

// Workflow in this script:
  // 1. Import land cover predictions output by "RF_predict.js" script
  // 2. Add them to the map


/*
  // Global Objects & functions ///////////////////////////////////////////////////////////////////////////////
*/

// Define image collection directory where your exported LC predictions are located
var exportColDir = 'users/grazingresearch/Europe_ELC10_v04/'

// Import land cover predictions stored in image collection and mosaic
  // Here I am importing the pre-mosaiced image instead of the exported collection of tiles
//var elc10 = ee.ImageCollection('users/zandersamuel/Europe_ELC10_v04').mosaic();
var elc10 = ee.Image('projects/nina/ELC10/ELC10_2018_10m');
print(elc10, 'elcq0')

// Define a dictionary which will be used to make legend and visualize image on map
var dict = {
  "names": [
    "Artificial land", //1
    "Cropland", //2
    "Woodland", //3
    "Shrubland",//4
    "Grassland", //5
    "Bare land", //6
    "Water", //7
    "Wetland",//8
  ],
  "colors": [
    "#CC0303", 
    "#CDB400", 
    "#235123",
    "#B76124", 
    "#92AF1F",
    "#F7E174",
    "#2019A4",
    "#AEC3D6",
  ]};

// Create a panel to hold the legend widget
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});

// Function to generate the legend
function addCategoricalLegend(panel, dict, title) {
  
  // Create and add the legend title.
  var legendTitle = ui.Label({
    value: title,
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
    }
  });
  panel.add(legendTitle);
  
  var loading = ui.Label('Loading legend...', {margin: '2px 0 4px 0'});
  panel.add(loading);
  
  // Creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
    // Create the label that is actually the colored box.
    var colorBox = ui.Label({
      style: {
        backgroundColor: color,
        // Use padding to give the box height and width.
        padding: '8px',
        margin: '0 0 4px 0'
      }
    });
  
    // Create the label filled with the description text.
    var description = ui.Label({
      value: name,
      style: {margin: '0 0 4px 6px'}
    });
  
    return ui.Panel({
      widgets: [colorBox, description],
      layout: ui.Panel.Layout.Flow('horizontal')
    });
  };
  
  // Get the list of palette colors and class names from the image.
  var palette = dict['colors'];
  var names = dict['names'];
  loading.style().set('shown', false);
  
  for (var i = 0; i < names.length; i++) {
    panel.add(makeRow(palette[i], names[i]));
  }
  
  Map.add(panel);
  
}


/*
  // Display map and legend ///////////////////////////////////////////////////////////////////////////////
*/

// Add the legend to the map
addCategoricalLegend(legend, dict, 'land cover');

// Add ELC10 image to the map
Map.addLayer(elc10, {min:1, max:8, palette:dict['colors']}, 'elc10')

// Center the map on image
Map.centerObject(elc10, 5)