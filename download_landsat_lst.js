// Define the area of interest (AOI) as Austin, Texas with a rectangular buffer
var austinCoords = ee.Geometry.Point([-97.7431, 30.2672]);
var bufferSize = 10000; // 10 km

// Create a rectangular buffer around the point
var aoi = austinCoords.buffer(bufferSize).bounds();

Map.addLayer(aoi, {}, 'AOI - Austin');
Map.centerObject(aoi, 11);

// Define start and end date
var startDate = ee.Date('2014-01-01');
var endDate = ee.Date('2023-12-31');

// Applies scaling factors
function applyScaleFactors(image) {
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBands = image.select('ST_B.*').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBands, null, true);
}

// Function to mask clouds and cloud shadows in Landsat 8 imagery
function cloudMask(image) {
  var cloudShadowBitmask = (1 << 3);
  var cloudBitmask = (1 << 5);
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(cloudShadowBitmask).eq(0)
                .and(qa.bitwiseAnd(cloudBitmask).eq(0));
  return image.updateMask(mask);
}

// Get the list of dates with available Landsat 8 images
var imageCollection = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                        .filterBounds(aoi)
                        .filterDate(startDate, endDate);

var dates = imageCollection.aggregate_array('system:time_start')
                            .map(function(date) {
                              return ee.Date(date).format('YYYY-MM-dd');
                            });

dates.evaluate(function(datesList) {
  // Loop through each date and process images
  datesList.forEach(function(date) {
    var dateStart = ee.Date(date);
    var dateEnd = dateStart.advance(1, 'day');

    // Filter images for the specific date
    var dailyImages = imageCollection
                      .filterDate(dateStart, dateEnd)
                      .map(applyScaleFactors)
                      .map(cloudMask)
                      .median()
                      .clip(aoi);

    // Calculate NDVI
    var ndvi = dailyImages.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');

    // Calculate the minimum and maximum NDVI values within the AOI
    var ndviMin = ee.Number(ndvi.reduceRegion({
      reducer: ee.Reducer.min(),
      geometry: aoi,
      scale: 30,
      maxPixels: 1e9
    }).values().get(0));

    var ndviMax = ee.Number(ndvi.reduceRegion({
      reducer: ee.Reducer.max(),
      geometry: aoi,
      scale: 30,
      maxPixels: 1e9
    }).values().get(0));

    // Calculate the Fraction of Vegetation (FV) using the NDVI values
    var fv = ((ndvi.subtract(ndviMin)).divide(ndviMax.subtract(ndviMin)))
              .pow(ee.Number(2))
              .rename('FV');

    // Calculate Land Surface Emissivity (EM) using the Fraction of Vegetation (FV)
    var em = fv.multiply(ee.Number(0.004)).add(ee.Number(0.986)).rename('EM');

    // Select Thermal Band (Band 10) and Rename It
    var thermal = dailyImages.select('ST_B10').rename('thermal');

    // Calculate Land Surface Temperature (LST) using Planck’s Law
    var lst = thermal.expression(
      '(TB / (1 + (0.00115 * (TB / 1.438)) * log(em))) - 273.15', {
        'TB': thermal.select('thermal'),
        'em': em
      }).rename('LST_Austin_' + date);

    // Add the LST Layer to the Map with Custom Visualization Parameters
    Map.addLayer(lst, {
      min: 10.47,
      max: 60.86,
      palette: [
        '040274', '040281', '0502a3', '0502b8', '0502ce', '0502e6',
        '0602ff', '235cb1', '307ef3', '269db1', '30c8e2', '32d3ef',
        '3be285', '3ff38f', '86e26f', '3ae237', 'b5e22e', 'd6e21f',
        'fff705', 'ffd611', 'ffb613', 'ff8b13', 'ff6e08', 'ff500d',
        'ff0000', 'de0101', 'c21301', 'a71001', '911003'
      ]}, 'Land Surface Temperature ' + date);

    // Export the LST image to Google Drive as GeoTIFF
    Export.image.toDrive({
      image: lst,
      description: 'LST_landsat8_Austin_' + date + '_GeoTIFF',
      folder: 'EarthEngineExports',
      scale: 30,
      region: aoi,
      fileFormat: 'GeoTIFF'
    });
  });
});
